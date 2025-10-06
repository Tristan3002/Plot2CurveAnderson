Here’s a tightened, battle-tested README you can drop in your repo. I folded in the exact Windows/PowerShell steps, VS Code setup (including the `reportMissingModuleSource` fix), headless runs, and the key environment variables your script supports.

---

# Gene Expression Simulator — README

A lightweight Python simulator for **transcription–translation dynamics** with plotting utilities. It sweeps a promoter library, integrates mRNA → immature protein → **mature** protein ODEs (with sfGFP-style maturation), and exports tidy CSVs plus publication-ready figures (PNGs + a combined PDF). Legends render in a **dedicated right-side panel** so they never collide with axes.

---

## Features

* **Deterministic ODE model** per promoter:

  * mRNA: `dm/dt = k_tx − δ_m m`
  * Protein (immature→mature): initiation-limited translation and first-order maturation
* **Promoter library** (Anderson J23100–J23118 + `CustomStrong`/`CustomWeak`) with relative strengths
* **Translation throughput** with initiation ceiling & ribosome traffic cap
* **RBS heuristics** (GC%, Tm, SD match, spacing → binding probability)
* **Resource coupling (optional)**: global ribosome throttle, β_eff(t) = β / (1 + M_tot/K_R)
* **Custom CDS support**: paste a FASTA/sequence or file path; validated start/stop and frame
* **Batch outputs** (timestamped folder):

  * `promoter_dynamics.csv` — timecourses
  * `promoter_summary.csv` — steady states, t50/t80, maturation lag, etc.
  * `fluorescence_dashboard.csv` — slim panel for quick comparisons
  * `anderson_k_tx_ladder.csv` — ranked ladder
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

2. *(Optional)* **Quiet a benign Pylance warning** (e.g., compiled libs):
   Create **.vscode/settings.json**:

   ```json
   {
     "python.defaultInterpreterPath": "${workspaceFolder}/.venv/Scripts/python.exe",
     "python.analysis.diagnosticSeverityOverrides": {
       "reportMissingModuleSource": "information"
     }
   }
   ```

   > macOS/Linux path: `"${workspaceFolder}/.venv/bin/python"`

3. *(Optional)* **One-click run, headless**: **.vscode/launch.json**

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
* `CDS_SEQ` — paste a DNA string / FASTA / file path to use a **custom CDS**
  *(validation: optional leading ATG is allowed; length must be multiple of 3; ends with TAA/TAG/TGA)*
* `THROTTLE=0` — disable global ribosome resource coupling (default is ON)
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

All outputs are written to a new folder like `sim_results_YYYYmmdd_HHMMSS/` in your project (or current working directory if you launch from elsewhere).

### CSV columns (high-level)

* **`promoter_dynamics.csv`**

  * `Time_s`, `Promoter`, `Promoter_relative_strength`
  * `mRNA_molecules`, `Protein_immature_molecules`, `Protein_molecules`
  * `mRNA_eq_molecules`, `Protein_eq_molecules` (mature), `Steady_state_protein_molecules`
  * `Translation_rate_k_eff_per_mRNA_per_s` (β), `RBS_binding_probability`, `k_tx_RNAs_per_s`
  * (when throttling is enabled) `M_tot_all_promoters`, `beta_eff_current`, `throttle_phi`, `K_R`, etc.

* **`promoter_summary.csv`**

  * `Rel_strength`, `k_tx_RNAs_per_s`, `m_eq`, `Protein_mature_eq`, `Protein_immature_eq`
  * `t50_mRNA_s`, `t50_protein_s` (numerical fallback if needed), `t10_p_m_s`, **`t80_p_m_s`**, `t90_p_m_s`
  * `maturation_lag_s` (= t50_protein − t50_mRNA), `fraction_mature_at_1800s`, `AUC_p_m_0_tmax`
  * **Detection panel**: `t_detect_s_at_thresh`, `detect_thresh_molecules` (auto-threshold ~10% of library max)

---

## Plots

* **Stacked per-promoter panels**: mRNA (left y) + mature protein (right y)
* **Overlay of N promoters**: solid = mRNA, dashed = protein; right-side legend panel
* **All-promoters overlay**: optional, with right-side legend
* **Gallery**: paged stacks (4 per page) so *all* are visible
* **Maturation stacks**: immature vs mature per promoter
* **Maturation lag vs brightness**: scatter with top points annotated
* **Brightness & speed bars**: steady-state brightness and **t80**
* **Detection time vs brightness**: reaches an **absolute molecule** threshold

---

## How the model works (short)

For each promoter, transcription uses `k_tx = k_tx_baseline × rel_strength`. mRNA decays with half-life `mRNA_t_half`. Protein forms in two stages: initiation-limited translation yields **immature** protein, which matures with rate `k_mat` (sfGFP-like). Growth dilution at rate `α = ln(2)/doubling_time` applies to both protein pools. Translation throughput uses `β = min(k_init_max·P_bind, v_nt/footprint)` where `P_bind` comes from RBS heuristics (SD match & spacing). **Optional** global coupling throttles β by `φ(t) = 1/(1 + M_tot/K_R)` where `M_tot` is the total mRNA across promoters.

---

## Typical tweaks

* **Promoter set**: edit `promoter_strengths` or pass names at the prompt
* **Window**: `t_max` and `dt` (dt auto-keeps stability vs `δ_m`, `α`, `k_mat`)
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
* **No plots on servers / WSL**
  Set `HEADLESS=1`.
* **No CSVs found**
  Look inside the printed `sim_results_*` path after the run. If you launched from another directory, outputs land in that **working directory**.
* **Long runs / no detection time**
  Increase `t_max`, or lower the detection threshold (`frac_of_max`) in `plot_detection_time_vs_brightness`.

---

## Suggested repo layout

```
.
├── Plot2CurveV8.py            # main script
├── README.md
├── .vscode/
│   ├── settings.json          # set interpreter + optional diagnostic tweak
│   └── launch.json            # one-click headless run
└── sim_results_YYYYmmdd_HHMMSS/   # generated per run
```

---

## Citation

If this tool helps your work, please cite the repository and the timestamped release/commit.

---

**Happy simming!** If something’s unclear, open an issue or ping the maintainer with your OS, Python version, and the exact terminal output.
