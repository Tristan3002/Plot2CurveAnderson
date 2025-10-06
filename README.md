---

# Gene Expression Simulator — README

A lightweight Python simulator for **transcription–translation dynamics** with plotting utilities. It sweeps a promoter library, integrates mRNA → immature protein → **mature** protein ODEs (with sfGFP-style maturation), and exports tidy CSVs plus publication-ready figures (PNGs + a combined PDF). Legends render in a **dedicated right-side panel** so they never collide with axes.

---

## Deliverables (for wet-lab companion)

1. **P2C “Anderson” report (sfGFP)**

   * **`figures.pdf`** → contains the *Anderson k_tx ladder* figure.
   * **`anderson_k_tx_ladder.csv`** → machine-readable **ranked ladder** with `Promoter`, `Rel_strength`, `k_tx_RNAs_per_s`, `Rank`.
   * Use the **top/bottom ranks** as your **shortlist of promoter swaps**.

2. **CSV exports for Learn → Design-2**

   * **`promoter_summary.csv`** supplies per-condition values:

     * **`k_tx_RNAs_per_s`** → *k_tx*
     * **`m_eq`** → *m_eq*
     * **`Protein_mature_eq`** (alias of p_eq) → *protein_eq*
   * Bonus columns: timing metrics (t50/t80), maturation lag, detection threshold time, etc.
   * **`promoter_dynamics.csv`** gives full timecourses (if trajectories are needed).

> Minimal ingest example:
>
> ```python
> import pandas as pd
> s = pd.read_csv("sim_results_*/promoter_summary.csv")
> df_design2 = s[["Promoter","k_tx_RNAs_per_s","m_eq","Protein_mature_eq"]].rename(
>     columns={"k_tx_RNAs_per_s":"k_tx","Protein_mature_eq":"protein_eq"}
> )
> ```

---

## Features

* **Deterministic ODE model** per promoter:

  * mRNA: `dm/dt = k_tx − δ_m m`
  * Protein (immature→mature): initiation-limited translation + first-order maturation
* **Promoter library** (Anderson J23100–J23118 + `CustomStrong`/`CustomWeak`) with relative strengths
* **Translation throughput** with initiation ceiling & ribosome traffic cap
* **RBS heuristics** (GC%, Tm, SD match, spacing → binding probability)
* **Resource coupling (optional)**: global ribosome throttle, β_eff(t) = β / (1 + M_tot/K_R)
* **Custom CDS support**: paste a FASTA/sequence or file path; validated start/stop and frame
* **Batch outputs** (timestamped folder):

  * `promoter_dynamics.csv` — timecourses
  * `promoter_summary.csv` — steady states, t50/t80, maturation lag, detection time, etc.
  * `fluorescence_dashboard.csv` — slim panel for quick comparisons
  * `anderson_k_tx_ladder.csv` — ranked ladder (shortlist source)
  * `figures.pdf` — every figure appended
  * Individual PNGs for each plot

---

## Quick Start (copy–paste)

> **Windows (PowerShell)**

```powershell
# 0) Open VS Code terminal in your project folder
cd "C:\path\to\your\project"      # put your real path in quotes

# 1) Create a virtual environment
py -m venv .venv

# 2) Install dependencies into that venv (activation optional)
.\.venv\Scripts\python -m pip install --upgrade pip setuptools wheel
.\.venv\Scripts\python -m pip install "numpy>=1.26,<3" "pandas>=2.0,<3" "matplotlib>=3.7,<4"

# 3) Run (headless avoids GUI popups)
$env:HEADLESS = "1"
.\.venv\Scripts\python .\Plot2CurveV8.py
```

> **macOS / Linux**

```bash
cd /path/to/your/project
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip setuptools wheel
python -m pip install "numpy>=1.26,<3" "pandas>=2.0,<3" "matplotlib>=3.7,<4"
HEADLESS=1 python Plot2CurveV8.py
```

**Prompts you’ll see:**

* “Enter promoters to test …” → comma-separated, or press Enter for **all**
* “Overlay promoters …” → optional, comma-separated list of **2+**
* “Plot ALL promoters together (y/n)?” → optional “gallery” overlay

> **Tip (Windows):** Don’t type `<your-project-folder>` literally. Angle brackets are **not** syntax in PowerShell—use your real path in quotes.

---

## VS Code setup (recommended)

1. **Interpreter**: Command Palette → **Python: Select Interpreter** → pick
   `.\.venv\Scripts\python.exe` (Windows) or `./.venv/bin/python` (macOS/Linux).

2. *(Optional)* Quiet a benign Pylance warning (compiled libs): **.vscode/settings.json**

   ```json
   {
     "python.defaultInterpreterPath": "${workspaceFolder}/.venv/Scripts/python.exe",
     "python.analysis.diagnosticSeverityOverrides": {
       "reportMissingModuleSource": "information"
     }
   }
   ```

   > macOS/Linux path: `"${workspaceFolder}/.venv/bin/python"`

3. *(Optional)* One-click headless run: **.vscode/launch.json**

   ```json
   {
     "version": "0.2.0",
     "configurations": [
       {
         "name": "Run promoter sim (headless)",
         "type": "python",
         "request": "launch",
         "program": "${file}",
         "console": "integratedTerminal",
         "cwd": "${workspaceFolder}",
         "env": { "HEADLESS": "1" }
       }
     ]
   }
   ```

---

## Environment variables (power users)

* `HEADLESS=1` — force non-interactive backend (no GUI windows)
* `CDS_SEQ` — DNA string / FASTA / file path for **custom CDS**
  *(validation: optional leading ATG allowed; length multiple of 3; ends with TAA/TAG/TGA)*
* `THROTTLE=0` — disable global ribosome resource coupling (default ON)
* `K_R=100.0` — ribosome budget scale (mRNA molecules/cell); smaller → stronger throttling

Examples:

```bash
# macOS/Linux
CDS_SEQ=./my_gene.fasta K_R=60 HEADLESS=1 python Plot2CurveV8.py
# Windows (PowerShell)
$env:CDS_SEQ="C:\path\gene.fasta"; $env:K_R="60"; $env:HEADLESS="1"; .\.venv\Scripts\python .\Plot2CurveV8.py
```

---

## Outputs

All outputs are written to a new folder like `sim_results_YYYYmmdd_HHMMSS/` in your project.

### CSV columns (high-level)

* **`promoter_dynamics.csv`**

  * `Time_s`, `Promoter`, `Promoter_relative_strength`
  * `mRNA_molecules`, `Protein_immature_molecules`, `Protein_molecules`
  * `mRNA_eq_molecules`, `Protein_eq_molecules` (mature), `Steady_state_protein_molecules`
  * `Translation_rate_k_eff_per_mRNA_per_s` (β), `RBS_binding_probability`, `k_tx_RNAs_per_s`
  * (if throttling) `M_tot_all_promoters`, `beta_eff_current`, `throttle_phi`, `K_R`, …

* **`promoter_summary.csv`**

  * **Design-2 keys:** `k_tx_RNAs_per_s` (*k_tx*), `m_eq` (*m_eq*), `Protein_mature_eq` (*protein_eq*)
  * Extras: `Protein_immature_eq`, `t50_mRNA_s`, `t50_protein_s` (robust fallback), `t10_p_m_s`,
    **`t80_p_m_s`**, `t90_p_m_s`, `maturation_lag_s`, `fraction_mature_at_1800s`, `AUC_p_m_0_tmax`,
    detection panel: `t_detect_s_at_thresh`, `detect_thresh_molecules`

* **Other**

  * `anderson_k_tx_ladder.csv` — ranked ladder + `Rank` (shortlist source)
  * `fluorescence_dashboard.csv` — compact view for quick comparisons
  * `figures.pdf` — all figures appended; includes *Anderson k_tx ladder*

---

## Plots

* **Stacked per-promoter panels**: mRNA (left y) + mature protein (right y)
* **Overlay of N promoters**: solid = mRNA, dashed = protein; right-side legend panel
* **All-promoters overlay**: optional, with right-side legend
* **Gallery**: paged stacks (4 per page) so *all* are visible
* **Maturation stacks**: immature vs mature per promoter
* **Maturation lag vs brightness**: scatter with top points annotated
* **Brightness & speed bars**: steady-state brightness and **t80**
* **Detection time vs brightness**: time to reach an **absolute molecule** threshold

---

## How the model works (short)

For each promoter, transcription uses `k_tx = k_tx_baseline × rel_strength`. mRNA decays with half-life `mRNA_t_half`. Protein forms in two stages: initiation-limited translation yields **immature** protein, which matures with rate `k_mat` (sfGFP-like). Growth dilution at rate `α = ln(2)/doubling_time` applies to both protein pools. Translation throughput uses `β = min(k_init_max·P_bind, v_nt/footprint)` where `P_bind` comes from RBS heuristics (SD match & spacing). **Optional** global coupling throttles β by `φ(t) = 1/(1 + M_tot/K_R)` where `M_tot` is the total mRNA across promoters.

---

## Typical tweaks

* **Promoter set**: edit `promoter_strengths` or pass names at the prompt
* **Window**: `t_max` and `dt` (dt auto-stabilizes vs `δ_m`, `α`, `k_mat`)
* **Biology**: `mRNA_t_half`, `doubling_time`, `aa_per_sec`, `ribosome_footprint_nt`, `k_init_max`
* **Legend layout**: `legend_cols`, `legend_fontsize`, `gutter`, `_panel_width_from_cols`

---

## Troubleshooting

* **Pylance “reportMissingModuleSource”** in VS Code
  Select your venv (**Python: Select Interpreter**) and, if desired, add the `diagnosticSeverityOverrides` snippet above.
* **“The '<' operator is reserved for future use.”**
  You pasted `<your-project-folder>` literally. Use your real path in quotes:

  ```powershell
  cd "C:\Users\you\Projects\gene-sim"
  ```
* **No plots on servers / WSL** → set `HEADLESS=1`
* **No CSVs found** → check the printed `sim_results_*` path (outputs go to the working directory)
* **Long runs / no detection time** → increase `t_max` or lower `frac_of_max` in detection plotting

---

## Suggested repo layout

```
.
├── Plot2CurveV8.py                 # main script
├── README.md
├── .vscode/
│   ├── settings.json               # set interpreter + optional diagnostic tweak
│   └── launch.json                 # one-click headless run
└── sim_results_YYYYmmdd_HHMMSS/    # generated per run
```

---

## Citation

If this tool helps your work, please cite the repository and the timestamped release/commit.

---

**Happy simming!**
