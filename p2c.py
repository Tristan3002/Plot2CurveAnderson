import argparse, csv, json, math, random, time, os
import codecs
from pathlib import Path
from datetime import datetime

MIN_OD = 0.01  # guard for very small ODs

from typing import Any, Optional

def safe_float(x: Any) -> Optional[float]:
    """Return float(x) or None for '', None, or non-convertible values.

    This narrows types for static checkers so callers get Optional[float]
    instead of accidentally passing None/"" into float(...).
    """
    try:
        if x is None or x == "":
            return None
        return float(x)
    except (TypeError, ValueError):
        return None

# ============================== I/O HELPERS ==============================

def read_plate_csv(path):
    """
    Accepts tidy CSV from an instrument export OR a live session CSV created by --init-session.
    Required columns: time_s, well, od600, gfp_rfu
    Optional columns (passed through if present): sample, medium, temp_C
    """
    rows = []
    # Use utf-8-sig to tolerate BOM from PowerShell/Excel
    with open(path, newline='', encoding='utf-8-sig') as f:
        r = csv.DictReader(f)
        needed = {'time_s','well','od600','gfp_rfu'}
        if not needed.issubset(set(r.fieldnames or [])):
            raise SystemExit(f"CSV missing columns. Need: {sorted(needed)}; got: {r.fieldnames}")
        for row in r:
            rec = {
                'time_s': float(row['time_s']),
                'well': str(row['well']).strip(),
                'od600': float(row['od600']),
                'gfp_rfu': float(row['gfp_rfu']),
            }
            # optional metadata if already present (e.g., in session.csv)
            # r.fieldnames may be None for some malformed inputs; guard accordingly
            if r.fieldnames and 'sample' in r.fieldnames:
                rec['sample'] = (row.get('sample') or '').strip()
            if r.fieldnames and 'medium' in r.fieldnames:
                rec['medium'] = (row.get('medium') or '').strip()
            if r.fieldnames and 'temp_C' in r.fieldnames:
                tmp = row.get('temp_C')
                try:
                    rec['temp_C'] = float(tmp) if tmp not in (None, "") else ""
                except Exception:
                    rec['temp_C'] = (tmp or "")
            rows.append(rec)
    return rows

def read_layout_json(path):
    """
    Read a JSON file robustly, tolerating an initial UTF-8 BOM (common on Windows).
    """
    data = Path(path).read_bytes()
    if data.startswith(codecs.BOM_UTF8):
        data = data[len(codecs.BOM_UTF8):]
    try:
        return json.loads(data.decode('utf-8'))
    except json.JSONDecodeError as e:
        raise SystemExit(f"Layout JSON parse error in {path}:\n  {e}")


def copy_layout_to(dest_dir: Path, layout_path: Path):
    dest_dir.mkdir(parents=True, exist_ok=True)
    target = dest_dir / "layout.json"
    if Path(layout_path).resolve() != target.resolve():
        data = read_layout_json(layout_path)
        with open(target, 'w', encoding='utf-8') as f:
            json.dump(data, f, indent=2)
    return target

def write_session_manifest(dest_dir: Path, global_medium, global_temp_c, blank_labels):
    man = {
        "created": datetime.now().isoformat(timespec='seconds'),
        "global_medium": global_medium,
        "global_temp_C": global_temp_c,
        "blank_labels": blank_labels,
        "version": 1
    }
    with open(dest_dir / "session_manifest.json", 'w', encoding='utf-8') as f:
        json.dump(man, f, indent=2)

def load_session_manifest(dest_dir: Path):
    p = dest_dir / "session_manifest.json"
    if p.exists():
        with open(p, encoding='utf-8') as f:
            return json.load(f)
    return {}

def well_sort_key(w):
    # Sort wells like A1, A2, ..., B1 ... by row then column
    w = w.upper().strip()
    row = ''.join([c for c in w if c.isalpha()])
    col = ''.join([c for c in w if c.isdigit()])
    row_num = 0
    for ch in row:
        row_num = row_num * 26 + (ord(ch) - ord('A') + 1)
    col_num = int(col) if col else 0
    return (row_num, col_num)

# ============================ METADATA / BLANKS ============================

def assign_metadata(rows, layout, global_medium=None, global_temp_c=None):
    """Attach sample/medium/temp to each row from layout (with global fallbacks).
       Do not overwrite if row already has 'sample' etc. (e.g., session.csv)."""
    for r in rows:
        if r.get('sample'):  # respect pre-filled metadata
            if 'medium' not in r: r['medium'] = global_medium
            if 'temp_C' not in r or r['temp_C'] == "": r['temp_C'] = global_temp_c
            continue
        meta = layout.get(r['well'])
        if isinstance(meta, dict):
            r['sample'] = meta.get('sample') or meta.get('name') or r['well']
            r['medium'] = meta.get('medium', global_medium)
            r['temp_C'] = meta.get('temp_C', global_temp_c)
        else:
            r['sample'] = (meta if isinstance(meta, str) else r['well'])
            r['medium'] = global_medium
            r['temp_C'] = global_temp_c

def detect_blank_wells(rows, blank_labels):
    """Return a set of wells considered blanks based on sample labels (case-insensitive)."""
    labels = {lbl.strip().lower() for lbl in blank_labels if lbl.strip()}
    blanks = set()
    for r in rows:
        name = str(r.get('sample', '')).strip().lower()
        if name in labels or name == 'blank':
            blanks.add(r['well'])
    return blanks

def timepoint_blank_means(rows, blank_wells):
    """Compute per-time_s mean blank OD and GFP across blank wells."""
    by_t = {}
    for r in rows:
        if r['well'] in blank_wells:
            t = r['time_s']
            by_t.setdefault(t, {'n':0,'od_sum':0.0,'gfp_sum':0.0})
            by_t[t]['n'] += 1
            by_t[t]['od_sum'] += r['od600']
            by_t[t]['gfp_sum'] += r['gfp_rfu']
    means = {}
    for t, agg in by_t.items():
        if agg['n'] > 0:
            means[t] = {
                'od': agg['od_sum']/agg['n'],
                'gfp': agg['gfp_sum']/agg['n'],
            }
    return means

def apply_blank_subtraction(rows, blank_means, warnings):
    """Subtract per-timepoint blank OD & GFP; compute adjusted F/OD. Floors at zero."""
    for r in rows:
        bm = blank_means.get(r['time_s'])
        if bm:
            od_adj = r['od600'] - bm['od']
            gfp_adj = r['gfp_rfu'] - bm['gfp']
            if od_adj < 0:
                warnings.append(f"NEG_OD_AFTER_SUB: well {r['well']} t={r['time_s']}s -> {od_adj:.4f} (floored to 0)")
                od_adj = 0.0
            if gfp_adj < 0:
                warnings.append(f"NEG_GFP_AFTER_SUB: well {r['well']} t={r['time_s']}s -> {gfp_adj:.4f} (floored to 0)")
                gfp_adj = 0.0
            r['od600_adj'] = od_adj
            r['gfp_rfu_adj'] = gfp_adj
            r['f_over_od_adj'] = (gfp_adj / od_adj) if od_adj > 0 else None
        else:
            r['od600_adj'] = ""
            r['gfp_rfu_adj'] = ""
            r['f_over_od_adj'] = ""

# ============================== CORE METRICS ==============================

def compute_f_over_od(rows, warnings):
    """Compute raw F/OD and basic warnings (without blanks)."""
    for r in rows:
        if r['od600'] <= MIN_OD:
            warnings.append(f"LOW_OD: well {r['well']} t={r['time_s']}s od600={r['od600']:.4f} (<= {MIN_OD})")
        r['f_over_od'] = (r['gfp_rfu'] / r['od600']) if r['od600'] > 0 else None
        if r['f_over_od'] is None:
            warnings.append(f"DIV_ZERO: well {r['well']} t={r['time_s']}s od600={r['od600']} caused f_over_od=None")

def per_sample_first_last(rows, use_adjusted: bool):
    """Return dicts: first[sample], last[sample] using chosen F/OD field."""
    first, last = {}, {}
    for r in rows:
        s = r['sample']
        if s not in first or r['time_s'] < first[s]['time_s']:
            first[s] = r
        if s not in last or r['time_s'] > last[s]['time_s']:
            last[s] = r
    return first, last

def growth_rates_alpha(first, last, warnings):
    """Compute rough alpha per sample (log-slope of OD across first->last)."""
    alphas = {}
    for s in last.keys():
        r0, r1 = first[s], last[s]
        if r0['od600'] > MIN_OD and r1['od600'] > MIN_OD and r1['time_s'] > r0['time_s']:
            alphas[s] = (math.log(r1['od600']) - math.log(r0['od600'])) / (r1['time_s'] - r0['time_s'])
        else:
            alphas[s] = None
            warnings.append(f"ALPHA_SKIP: sample {s} insufficient OD or time spacing")
    return alphas

def final_per_well(rows, use_adjusted: bool, exclude_wells=set()):
    """For replicate stats: for each (sample, well) pick the last timepoint and return final F/OD list per sample."""
    by_sw = {}
    for r in rows:
        if r['well'] in exclude_wells:
            continue
        key = (r['sample'], r['well'])
        if key not in by_sw or r['time_s'] > by_sw[key]['time_s']:
            fod = r.get('f_over_od_adj') if use_adjusted else r.get('f_over_od')
            by_sw[key] = {'time_s': r['time_s'], 'final_fod': fod}
    per_sample = {}
    for (s, _w), rec in by_sw.items():
        v = safe_float(rec['final_fod'])
        if v is not None:
            per_sample.setdefault(s, []).append(v)
    return per_sample

def bootstrap_ci(values, n_boot=1000, ci=0.95, seed=None):
    """Return (mean, ci_low, ci_high). If <2 values, CI is None."""
    if not values:
        return None, None, None
    rnd = random.Random(seed) if seed is not None else random
    n = len(values)
    mean_val = sum(values)/n
    if n < 2 or n_boot <= 0:
        return mean_val, None, None
    boots = []
    for _ in range(n_boot):
        sample = [values[rnd.randrange(n)] for _ in range(n)]
        boots.append(sum(sample)/n)
    boots.sort()
    alpha = 1.0 - ci
    lo_idx = max(0, int(math.floor(alpha/2 * len(boots))) )
    hi_idx = min(len(boots)-1, int(math.ceil((1 - alpha/2) * len(boots))) - 1)
    return mean_val, boots[lo_idx], boots[hi_idx]

def build_summary(rows, use_adjusted: bool, warnings, n_boot, ci_level, seed, blank_wells):
    first, last = per_sample_first_last(rows, use_adjusted)
    alphas = growth_rates_alpha(first, last, warnings)
    per_sample_finals = final_per_well(rows, use_adjusted, exclude_wells=blank_wells)

    summary = []
    for s in sorted(last.keys()):
        r0, r1 = first[s], last[s]
        fod0 = (r0.get('f_over_od_adj') if use_adjusted else r0.get('f_over_od'))
        fod1 = (r1.get('f_over_od_adj') if use_adjusted else r1.get('f_over_od'))
        f0f = safe_float(fod0)
        f1f = safe_float(fod1)
        delta = None if (f0f is None or f1f is None) else (f1f - f0f)
        fold = None if (f0f is None or f1f is None or f0f == 0) else (f1f / f0f)
        a = alphas.get(s)

        finals = per_sample_finals.get(s, [])
        mean_f, ci_lo, ci_hi = bootstrap_ci(finals, n_boot=n_boot, ci=ci_level, seed=seed)

        med = r0.get('medium'); tmp = r0.get('temp_C')
        temp_str = "" if tmp in (None, "") else (f"{tmp:.1f}" if isinstance(tmp, (int, float)) else str(tmp))

        row = {
            'sample': s,
            'medium': "" if med is None else str(med),
            'temp_C': temp_str,
            'time_s_final': int(r1['time_s']),
            'initial_f_over_od': "" if (f0f is None) else f"{f0f:.3f}",
            'final_f_over_od': "" if (f1f is None) else f"{f1f:.3f}",
            'delta_f_over_od': "" if delta is None else f"{delta:.3f}",
            'fold_change_f_over_od': "" if fold is None else f"{fold:.3f}",
            'alpha_per_s': "" if (a is None) else f"{a:.6f}",
            'replicates_n': "" if not finals else str(len(finals)),
            'final_f_over_od_mean': "" if mean_f is None else f"{mean_f:.3f}",
            'final_f_over_od_ci95_low': "" if ci_lo is None else f"{ci_lo:.3f}",
            'final_f_over_od_ci95_high': "" if ci_hi is None else f"{ci_hi:.3f}",
        }
        summary.append(row)
    return summary

# =============================== REPORTING ================================

def write_csv(path, rows, fieldnames):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, 'w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in fieldnames})

def write_qc_report(path, warnings, used_blank_sub):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, 'w', encoding='utf-8') as f:
        f.write("# P2C mini QC report\n")
        f.write(f"Generated: {datetime.now().isoformat(timespec='seconds')}\n")
        f.write(f"Blank subtraction: {'ON' if used_blank_sub else 'OFF'}\n")
        f.write(f"Warnings: {len(warnings)}\n\n")
        for w in warnings:
            f.write(f"- {w}\n")

# ============================= INTERACTIVE MODES =============================

def init_session(session_dir: Path, layout_path: Path, global_medium, global_temp_c, blank_labels):
    """Create a session directory with session.csv, layout.json copy, and manifest."""
    session_dir.mkdir(parents=True, exist_ok=True)
    out_csv = session_dir / "session.csv"

    # Copy/normalize layout
    layout_target = copy_layout_to(session_dir, layout_path)

    # Manifest & empty session CSV
    write_session_manifest(session_dir, global_medium, global_temp_c, blank_labels)
    header = ["time_s","well","sample","medium","temp_C","od600","gfp_rfu"]
    with open(out_csv, 'w', newline='', encoding='utf-8') as f:
        w = csv.writer(f)
        w.writerow(header)

    print(f"[init] Created session at {session_dir.resolve()}")
    print(f"[init] Empty CSV: {out_csv}")
    print(f"[init] Layout:   {layout_target}")
    print(f"[init] Manifest: {session_dir/'session_manifest.json'}")
    print("[init] Next: append data during the run with --append-session and a timepoint.")

def parse_two_floats(text):
    """Parse '0.05,120' or '0.05 120' -> (0.05,120). Return None if blank."""
    if text is None: return None
    t = text.strip()
    if t == "": return None
    t = t.replace(",", " ")
    parts = [p for p in t.split() if p]
    if len(parts) != 2:
        raise ValueError("Enter two numbers: <od600> <gfp_rfu> (comma or space separated)")
    return float(parts[0]), float(parts[1])

def append_session(session_csv: Path, time_s: float, wells_subset=None):
    """Interactive prompt to append a timepoint of readings for each well in layout.json."""
    session_dir = session_csv.parent
    layout_path = session_dir / "layout.json"
    if not layout_path.exists():
        raise SystemExit(f"layout.json not found next to {session_csv}. Did you run --init-session?")
    layout = read_layout_json(layout_path)
    manifest = load_session_manifest(session_dir)
    glob_med = manifest.get("global_medium")
    glob_temp = manifest.get("global_temp_C")

    # Build well list in sorted order
    wells = sorted(layout.keys(), key=well_sort_key)
    if wells_subset:
        wanted = {w.strip().upper() for w in wells_subset}
        wells = [w for w in wells if w.upper() in wanted]

    # Ensure CSV exists with header
    exists = session_csv.exists()
    if not exists:
        with open(session_csv, 'w', newline='', encoding='utf-8') as f:
            w = csv.writer(f); w.writerow(["time_s","well","sample","medium","temp_C","od600","gfp_rfu"])

    print(f"[append] Session: {session_csv.resolve()}")
    print(f"[append] Timepoint: {time_s} s  (Ctrl+C to stop early)")

    try:
        with open(session_csv, 'a', newline='', encoding='utf-8') as f:
            w = csv.writer(f)
            for well in wells:
                meta = layout.get(well)
                if isinstance(meta, dict):
                    sample = meta.get("sample") or meta.get("name") or well
                    med = meta.get("medium", glob_med)
                    tmp = meta.get("temp_C", glob_temp)
                else:
                    sample = meta if isinstance(meta, str) else well
                    med = glob_med
                    tmp = glob_temp
                prompt = f"{well} ({sample})  enter 'od,gfp' or Enter to skip: "
                inp = input(prompt)
                try:
                    pair = parse_two_floats(inp)
                except ValueError as e:
                    print(f"  ! {e}  (skipped)")
                    continue
                if pair is None:
                    continue
                od, gfp = pair
                w.writerow([time_s, well, sample, med if med is not None else "", tmp if tmp is not None else "", od, gfp])
                print(f"  ✓ recorded {well}: od={od} gfp={gfp}")
    except KeyboardInterrupt:
        print("\n[append] stopped by user")

def run_logger(out_dir: Path, interval_min: float):
    """
    Minimal interactive note logger.
    Run in a separate terminal while the plate reader runs:
      python p2c.py --out data_processed --log-interval-min 10 --logger-only
    Press Ctrl+C to stop. Writes notes_log.csv with timestamp and elapsed minutes.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    log_path = out_dir / "notes_log.csv"
    header_needed = not log_path.exists()
    print(f"[logger] writing notes to {log_path.resolve()}")
    start = time.time()
    try:
        while True:
            time.sleep(max(1, int(interval_min * 60)))
            ts = datetime.now().isoformat(timespec='seconds')
            elapsed = (time.time() - start) / 60.0
            note = input(f"[{ts}] ~{elapsed:.1f} min. Enter note (or Enter to skip): ").strip()
            with open(log_path, 'a', newline='', encoding='utf-8') as f:
                w = csv.writer(f)
                if header_needed:
                    w.writerow(["timestamp_iso", "elapsed_min", "note"])
                    header_needed = False
                w.writerow([ts, f"{elapsed:.2f}", note])
            print("[logger] saved.")
    except KeyboardInterrupt:
        print("\n[logger] stopped.")

# ================================ CLI MAIN ================================

def main():
    p = argparse.ArgumentParser(description="Plot-2-Curve mini MVP (now with sessions)")
    # Modes
    p.add_argument("--init-session", default=None, help="Create a session folder here (or provide a path). Requires --layout.")
    p.add_argument("--append-session", default=None, help="Path to session.csv (or its folder) to append a timepoint interactively.")
    p.add_argument("--time-s", type=float, default=None, help="Time (seconds) for --append-session (>= 0 allowed).")
    p.add_argument("--time-min", type=float, default=None, help="Time (minutes) for --append-session (>= 0 allowed).")
    p.add_argument("--wells", default=None, help="Optional comma-separated subset of wells to prompt (e.g., A1,A2,B1).")

    # Analysis inputs (legacy/simple usage)
    p.add_argument("--plate", help="CSV: time_s,well,od600,gfp_rfu (or a session.csv)")
    p.add_argument("--layout", help="JSON: well -> sample OR {'sample','medium','temp_C'} (not needed if session.csv already has 'sample')")
    p.add_argument("--out", default="data_processed", help="Output folder")

    # Metadata for init/analyze
    p.add_argument("--medium", default=None, help="Global medium (layout can override)")
    p.add_argument("--temp-c", dest="temp_c", type=float, default=None, help="Global temperature °C (layout can override)")
    p.add_argument("--blank-subtract", action="store_true", help="Subtract per-timepoint mean of blank wells (labelled 'blank' or via --blank-labels)")
    p.add_argument("--blank-labels", default="blank", help="Comma-separated sample labels treated as blanks (case-insensitive). Default: 'blank'")

    # Bootstrap CIs
    p.add_argument("--bootstrap-n", type=int, default=1000, help="Bootstrap resamples for final F/OD mean CI (default 1000)")
    p.add_argument("--ci", type=float, default=0.95, help="Confidence level for bootstrap CI (default 0.95)")
    p.add_argument("--seed", type=int, default=None, help="Random seed for bootstrap (optional)")

    # Logger
    p.add_argument("--log-interval-min", type=float, default=None, help="Start an interactive note logger (minutes between prompts)")
    p.add_argument("--logger-only", action="store_true", help="Run only the note logger and exit")

    args = p.parse_args()

    # -------------------- Mode: logger-only --------------------
    if args.logger_only:
        out_dir = Path(args.out)
        run_logger(out_dir, args.log_interval_min or 10.0)
        return

    # -------------------- Mode: init session -------------------
    if args.init_session:
        session_dir = Path(args.init_session)
        if session_dir.suffix.lower() == ".csv":
            session_dir = session_dir.parent
        if not args.layout:
            raise SystemExit("--init-session requires --layout <layout.json>")
        init_session(
            session_dir=session_dir,
            layout_path=Path(args.layout),
            global_medium=args.medium,
            global_temp_c=args.temp_c,
            blank_labels=[s.strip() for s in (args.blank_labels or "").split(",")]
        )
        return

    # -------------------- Mode: append session -----------------
    if args.append_session:
        session_path = Path(args.append_session)
        # Allow passing a folder; then use ./session.csv inside it
        if session_path.is_dir():
            session_csv = session_path / "session.csv"
        else:
            session_csv = session_path
        if not session_csv.exists() and session_csv.suffix.lower() != ".csv":
            raise SystemExit("--append-session must point to a session.csv or its folder.")
        # allow t = 0 (was previously > 0)
        t_s = args.time_s if args.time_s is not None else (
            args.time_min * 60.0 if args.time_min is not None else None
        )
        if t_s is None or t_s < 0:
            raise SystemExit("Provide --time-s or --time-min (>= 0) for the timepoint to append.")
        wells_subset = [w.strip() for w in args.wells.split(",")] if args.wells else None
        append_session(session_csv=session_csv, time_s=t_s, wells_subset=wells_subset)
        return

    # -------------------- Mode: analyze (default) ---------------
    # Minimal backward-compatible analysis mode (same as before, but now tolerates session.csv with sample/medium/temp_C already present)
    if not args.plate:
        raise SystemExit("Provide --plate (CSV) for analysis, or use --init-session / --append-session modes.")

    out_dir = Path(args.out)
    print("[dbg] running:", Path(__file__).resolve())
    print("[dbg] analyzing to out dir:", out_dir.resolve())

    # Load rows and, if needed, attach metadata from layout
    rows = read_plate_csv(args.plate)

    layout = {}
    if args.layout:
        layout = read_layout_json(args.layout)

    assign_metadata(rows, layout, global_medium=args.medium, global_temp_c=args.temp_c)

    warnings = []
    # Raw F/OD
    compute_f_over_od(rows, warnings)

    # Optional blank subtraction
    use_adjusted = False
    blank_wells = set()
    if args.blank_subtract:
        labels = [s.strip() for s in (args.blank_labels or "").split(",")]
        blank_wells = detect_blank_wells(rows, labels)
        if not blank_wells:
            warnings.append("BLANK_SUBTRACT_REQUESTED_BUT_NO_BLANKS_FOUND")
        else:
            bm = timepoint_blank_means(rows, blank_wells)
            if not bm:
                warnings.append("BLANK_SUBTRACT: found blank wells but no matching timepoints")
            apply_blank_subtraction(rows, bm, warnings)
            use_adjusted = True

    # Build summary (includes bootstrap stats across replicate wells)
    summary = build_summary(
        rows,
        use_adjusted=use_adjusted,
        warnings=warnings,
        n_boot=args.bootstrap_n,
        ci_level=args.ci,
        seed=args.seed,
        blank_wells=blank_wells
    )

    # Write outputs
    summary_cols = [
        "sample","medium","temp_C","time_s_final",
        "initial_f_over_od","final_f_over_od","delta_f_over_od","fold_change_f_over_od","alpha_per_s",
        "replicates_n","final_f_over_od_mean","final_f_over_od_ci95_low","final_f_over_od_ci95_high"
    ]
    tidy_cols = ["time_s","well","sample","medium","temp_C","od600","gfp_rfu","f_over_od","od600_adj","gfp_rfu_adj","f_over_od_adj"]

    write_csv(out_dir / "summary.csv", summary, summary_cols)
    write_csv(out_dir / "tidy_timeseries.csv", rows, tidy_cols)
    write_qc_report(out_dir / "qc_report.txt", warnings, used_blank_sub=use_adjusted)

    # quick peek
    print(f"[OK] Wrote {out_dir / 'summary.csv'}")
    with open(out_dir / "summary.csv", encoding='utf-8') as fh:
        for i, line in enumerate(fh):
            print("[summary.csv]", line.strip())
            if i >= 1:
                break
    print(f"[OK] Wrote {out_dir / 'tidy_timeseries.csv'}")
    print(f"[OK] Wrote {out_dir / 'qc_report.txt'}")
    print("[OK] Done.")

if __name__ == "__main__":
    main()

