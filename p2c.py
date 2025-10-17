import argparse, csv, json
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
    # tag rows with sample name and GFP/OD
    for r in rows:
        r['sample'] = layout.get(r['well'], r['well'])
        r['f_over_od'] = (r['gfp_rfu'] / r['od600']) if r['od600'] > 0 else None

    # keep the last timepoint per sample
    last = {}
    for r in rows:
        s = r['sample']
        if s not in last or r['time_s'] > last[s]['time_s']:
            last[s] = r

    # small summary table
    summary = []
    for s, r in sorted(last.items()):
        summary.append({
            'sample': s,
            'time_s': int(r['time_s']),
            'final_f_over_od': "" if r['f_over_od'] is None else f"{r['f_over_od']:.3f}"
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
    args = p.parse_args()

    plate_rows = read_plate_csv(args.plate)
    layout = read_layout_json(args.layout)
    summary, all_rows = compute_summary(plate_rows, layout)

    out_dir = Path(args.out)
    write_csv(out_dir / "summary.csv", summary, ["sample","time_s","final_f_over_od"])
    write_csv(out_dir / "tidy_timeseries.csv", all_rows,
              ["time_s","well","sample","od600","gfp_rfu","f_over_od"])

    print(f"[OK] Wrote {out_dir / 'summary.csv'}")
    print(f"[OK] Wrote {out_dir / 'tidy_timeseries.csv'}")
    print("[OK] Done.")

if __name__ == "__main__":
    main()
