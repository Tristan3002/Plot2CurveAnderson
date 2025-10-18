import argparse, csv, json, shutil
from pathlib import Path

def read_plate_csv(path):
    rows = []
    with open(path, newline='', encoding='utf-8') as f:
        r = csv.DictReader(f)
        needed = {'time_s','well','od600','gfp_rfu'}
        if not needed.issubset(set(r.fieldnames or [])):
            raise SystemExit(f"CSV missing columns. Need: {sorted(needed)}; got: {r.fieldnames}")
        for row in r:
            rows.append({
                'time_s': float(row['time_s']),
                'well': row['well'].strip(),
                'od600': float(row['od600']),
                'gfp_rfu': float(row['gfp_rfu']),
            })
    return rows

def read_layout_json(path):
    with open(path, encoding='utf-8') as f:
        return json.load(f)  # e.g., {"A1":"J23100","A2":"J23101"}

def compute_summary(rows, layout):
    # tag rows with sample + GFP/OD
    for r in rows:
        r['sample'] = layout.get(r['well'], r['well'])
        r['f_over_od'] = (r['gfp_rfu'] / r['od600']) if r['od600'] > 0 else None

    # first and last timepoint per sample
    first, last = {}, {}
    for r in rows:
        s = r['sample']
        if s not in first or r['time_s'] < first[s]['time_s']:
            first[s] = r
        if s not in last or r['time_s'] > last[s]['time_s']:
            last[s] = r

    summary = []
    for s in sorted(last.keys()):
        f0 = first[s].get('f_over_od')
        f1 = last[s].get('f_over_od')
        delta = (None if (f0 is None or f1 is None) else (f1 - f0))
        fold = (None if (f0 in (None, 0)) else (f1 / f0))
        summary.append({
            'sample': s,
            'time_s_final': int(last[s]['time_s']),
            'initial_f_over_od': "" if f0 is None else f"{f0:.3f}",
            'final_f_over_od': "" if f1 is None else f"{f1:.3f}",
            'delta_f_over_od': "" if delta is None else f"{delta:.3f}",
            'fold_change_f_over_od': "" if fold is None else f"{fold:.3f}",
        })
    return summary, rows

def write_csv(path, rows, fieldnames):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, 'w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in fieldnames})

def main():
    p = argparse.ArgumentParser(description="Plot-2-Curve mini MVP")
    p.add_argument("--plate", required=True, help="CSV: time_s,well,od600,gfp_rfu")
    p.add_argument("--layout", required=True, help="JSON: well -> sample")
    p.add_argument("--out", default="data_processed", help="Output folder")
    p.add_argument("--clean", action="store_true", help="Delete output folder before writing")
    args = p.parse_args()
    print("[dbg] running:", Path(__file__).resolve())
    print("[dbg] writing to out dir:", Path(args.out).resolve())

    out_dir = Path(args.out)
    if args.clean and out_dir.exists():
        shutil.rmtree(out_dir)

    plate_rows = read_plate_csv(args.plate)
    layout = read_layout_json(args.layout)
    summary, all_rows = compute_summary(plate_rows, layout)

    # Write outputs
    summary_cols = ["sample","time_s_final","initial_f_over_od","final_f_over_od","delta_f_over_od","fold_change_f_over_od"]
    write_csv(out_dir / "summary.csv", summary, summary_cols)
    write_csv(out_dir / "tidy_timeseries.csv", all_rows, ["time_s","well","sample","od600","gfp_rfu","f_over_od"])

    # Quick peek so you can see the new columns in the console
    print(f"[OK] Wrote {out_dir / 'summary.csv'}")
    with open(out_dir / "summary.csv", encoding="utf-8") as fh:
        for i, line in enumerate(fh):
            print("[summary.csv]", line.strip())
            if i >= 1:  # header + first row
                break
    print(f"[OK] Wrote {out_dir / 'tidy_timeseries.csv'}")
    print("[OK] Done.")
if __name__ == "__main__":
    main()