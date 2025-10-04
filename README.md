# Gene Expression Simulator — README

A lightweight Python simulator for **transcription–translation dynamics** with plotting utilities. It sweeps a promoter library, integrates simple mRNA/protein ODEs, and exports tidy CSVs and publication-ready figures (PNG + a combined PDF). Legends are laid out in a **dedicated right-side panel** so they never collide with axis labels.

---

## Features

* **Deterministic ODE model** of mRNA and protein per promoter.
* **Promoter library** (iGEM J231xx plus “CustomStrong/Weak”) with relative strengths.
* **Translation throughput** accounts for initiation ceiling and ribosome traffic.
* **RBS heuristics** (GC%, melting temp, SD match, spacing → binding probability).
* **Batch outputs**:

  * `promoter_dynamics.csv` — timecourses per promoter.
  * `promoter_summary.csv` — steady states, t50s, etc.
  * `figures.pdf` — all figures aggregated.
  * Individual PNGs for every plot.
* **Plotting**:

  * Stacked per-promoter panels (`plot_stacked_promoter_panels`)
  * **Overlay of N promoters** in one figure (`plot_overlay_promoters`)
  * “All promoters” overlay (`plot_all_promoters_overlay`)
  * Paged gallery for the selected set (`plot_gallery_promoters`)
* **Robust legends**: right-side legend panel (*never overlaps axes/labels*).

---

## Quick start

```bash
# 1) Create & activate a virtual environment (recommended)
python -m venv .venv
source .venv/bin/activate        # Windows: .venv\Scripts\activate

# 2) Install dependencies
pip install numpy pandas matplotlib

# 3) Run
python simulate.py
```

You’ll be prompted for:

* **Promoters to test** (comma-separated, empty = all)
* **Overlay promoters** (optional; comma-separated list of 2+)
* Whether to **plot ALL** promoters together (y/n)

All outputs land in a time-stamped folder like `sim_results_20250130_142233/`.

> **Headless runs (servers/CI):**
> Set `HEADLESS=1` to force a non-interactive backend and avoid GUI windows:
>
> ```bash
> HEADLESS=1 python simulate.py
> ```

---

## Outputs

Inside each `sim_results_<timestamp>/`:

* `promoter_dynamics.csv` — columns:

  * `Time_s`, `Promoter`, `Promoter_relative_strength`, `mRNA_molecules`, `Protein_molecules`, `mRNA_eq_molecules`, `Protein_eq_molecules`, `Steady_state_protein_molecules`, `Translation_rate_k_eff_per_mRNA_per_s`, `RBS_binding_probability`
* `promoter_summary.csv` — one row per promoter:

  * `Rel_strength`, `k_tx_RNAs_per_s`, `m_eq`, `p_eq`, `t50_mRNA_s`, `t50_protein_s`, `mRNA_final`, `Protein_final`, `beta_per_mRNA_per_s`, `P_bind`
* `figures.pdf` — every figure appended
* PNGs for each figure (gallery pages, overlays, etc.)

---

## How the model works (one paragraph)

For each promoter, transcription is `dm/dt = k_tx − δ_m m` and protein is `dp/dt = β m − α p`. The transcription rate is `k_tx = k_tx_baseline × rel_strength`. The effective translation rate `β` accounts for initiation (scaled by RBS binding probability) and ribosome traffic (`min(k_init, v_nt/footprint)`). RBS heuristics (GC%, approximate SD pairing, spacing) give a coarse **binding probability** used to scale initiation.

---

## Usage patterns

### Run everything interactively

Just execute `python simulate.py` and respond to prompts.

### Non-interactive overlay (example in code)

At the “Overlay promoters …” prompt, you can paste:

```
J23100,J23105,J23118,CustomStrong
```

This produces a single figure with **solid** lines for mRNA (left y) and **dashed** lines for protein (right y), color-coded per promoter. The legend is placed in a dedicated right-hand panel.

### Paged gallery

The gallery auto-renders **all selected promoters** in stacks of 4 per page so you can scan the entire set quickly.

---

## Legend system (no overlap, adjustable)

The overlay and “all promoters” plots use:

* `_panel_width_from_cols(ncols)` — maps number of legend columns → panel width (fraction of figure).
* `_add_right_side_legend(fig, handles, labels, ncols=..., title=..., fontsize=..., title_fontsize=..., gutter=0.06)` — shrinks the plotting area and creates a **separate Axes** to host the legend on the right. This guarantees the legend never collides with the right y-axis label/ticks.

**Fine-tuning tips:**

* Increase `gutter` (e.g., `0.08`) to push the legend further to the right.
* Increase `legend_cols` to reduce the legend’s height and make better use of the panel width.
* If you have very long promoter names, expand `_panel_width_from_cols` mapping (e.g., bump values by `+0.02`).

---

## Key parameters you may want to tweak

* `k_tx_baseline` (RNAs/s/cell), `promoter_strengths` dict
* `aa_per_sec`, `ribosome_footprint_nt`, `k_init_max`
* mRNA half-life (`mRNA_t_half`) and doubling time (`doubling_time`)
* Simulation window (`t_max`) and step (`dt`)
* Plot options: log scales, legend columns/title/font

---

## Troubleshooting

* **“No display name and no $DISPLAY”** (Linux/CI): run with `HEADLESS=1`.
* **Legend still feels tight**: raise `gutter` or increase `legend_cols`. For very large lists, try `legend_cols=3` or `4`.
* **Missing promoters**: if you type a promoter not in `promoter_strengths`, the code ignores it (or backfills by simulating if found in the library).

---

## Repo layout (suggested)

```
.
├── simulate.py                 # the script (your current file)
├── README.md                   # this readme
├── requirements.txt            # optional pinning
└── sim_results_YYYYmmdd_HHMMSS/ (generated per run)
```

**requirements.txt (optional pinning):**

```
numpy>=1.23
pandas>=1.5
matplotlib>=3.6
```

---

## Notes & disclaimer

This is a **didactic** kinetic model / plotting tool. It makes simplifying assumptions and is **not** a substitute for experimental characterization or detailed biophysical modeling. No wet-lab procedures are included or implied.

---

## Citation

If this tool helps your work, please cite the repository and include the timestamped release commit.
