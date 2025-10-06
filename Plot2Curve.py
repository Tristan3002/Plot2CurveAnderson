# -*- coding: utf-8 -*-
import os
from datetime import datetime

# ---------- Backend selection BEFORE importing pyplot ----------
import matplotlib
# Toggle this to disable GUI windows globally (or set env HEADLESS=1)
SHOW_PLOTS = os.environ.get("HEADLESS", "0") != "1"

if not SHOW_PLOTS:
    matplotlib.use("Agg")  # headless
# If SHOW_PLOTS is True, let Matplotlib pick an interactive backend (e.g., TkAgg)
# You can force one by uncommenting: matplotlib.use("TkAgg")

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.lines import Line2D  # for custom legend handles

# ============== Output folder & PDF collector ==============
RUN_STAMP = datetime.now().strftime("%Y%m%d_%H%M%S")
OUTDIR = os.path.join(os.getcwd(), f"sim_results_{RUN_STAMP}")
os.makedirs(OUTDIR, exist_ok=True)
PDF_PATH = os.path.join(OUTDIR, "figures.pdf")
_pdf = PdfPages(PDF_PATH)

def _save_figure(fig, base_name, outdir=OUTDIR, pdf=_pdf):
    """Save fig as PNG and also append to the PDF."""
    png_path = os.path.join(outdir, f"{base_name}.png")
    fig.savefig(png_path, dpi=300, bbox_inches='tight')
    if pdf is not None:
        pdf.savefig(fig, bbox_inches='tight')
    print(f"Saved: {os.path.relpath(png_path)}")

# ---------- Legend helper: put a dedicated legend panel on the right ----------
def _panel_width_from_cols(ncols: int) -> float:
    """
    Return a figure-fraction width for a right-side legend panel based on cols.
    Tuned to comfortably fit long labels without colliding with the right y-label.
    """
    # Wider panel as columns increase
    return {1: 0.18, 2: 0.26, 3: 0.34, 4: 0.42}.get(int(ncols), 0.42)

def _add_right_side_legend(fig, handles, labels, *, ncols=1, title="Promoters",
                           fontsize=8, title_fontsize=9, gutter=0.06):
    """
    Shrinks the plotting area and adds a separate Axes at the right solely for the legend.
    This guarantees no overlap with axis labels/ticks.
    """
    panel_w = _panel_width_from_cols(ncols)
    # Make space for gutter + legend panel
    fig.subplots_adjust(right=1 - (panel_w + gutter))
    sp = fig.subplotpars  # current subplot geometry after adjustment

    # Create legend Axes spanning the plot height
    # Slight insets on left/right so the frame doesn't touch the canvas edge
    leg_ax = fig.add_axes([1 - panel_w + 0.01, sp.bottom, panel_w - 0.02, sp.top - sp.bottom])
    leg_ax.axis("off")
    leg_ax.legend(handles, labels, loc="center left", frameon=True, ncol=ncols,
                  title=title, fontsize=fontsize, title_fontsize=title_fontsize)
    return leg_ax

# --- Helper Functions (unchanged except noted) ---
def calc_gc_content(seq):
    gc = seq.count('G') + seq.count('C')
    return gc / len(seq) * 100 if len(seq) > 0 else 0

def melting_temp(seq):
    at = seq.count('A') + seq.count('T')
    gc = seq.count('G') + seq.count('C')
    return 2 * at + 4 * gc

def check_sd_best_match(seq, anti_sd="CCUCCU"):
    seq_rna = seq.upper().replace('T', 'U')
    anti_sd = anti_sd.upper()
    pairs = {"A": "U", "U": "A", "G": "C", "C": "G"}
    max_score = 0
    for i in range(len(seq_rna) - len(anti_sd) + 1):
        window = seq_rna[i:i+len(anti_sd)]
        score = sum(1 for a, b in zip(window, anti_sd) if pairs.get(a) == b)
        if score > max_score:
            max_score = score
    return max_score

def check_spacing(rbs_end_pos, start_codon_pos):
    return start_codon_pos - rbs_end_pos - 1

def estimate_binding_probability(sd_score, spacing):
    sd_max = 6
    ideal_spacing = 7
    spacing_score = max(0, 1 - abs(spacing - ideal_spacing) / 5)
    probability = 0.6 * (sd_score / sd_max) + 0.4 * spacing_score
    return min(1.0, max(0.0, round(probability, 3)))

# ---- initiation-limited translation throughput (proteins per mRNA per s)
def compute_k_eff(cds_len_nt, aa_per_sec, footprint_nt, P_bind, k_init_max=1.0):
    v_nt = 3 * aa_per_sec                 # nt/s along CDS
    f_max = v_nt / footprint_nt           # s^-1 per mRNA (traffic limit)
    k_init = k_init_max * P_bind          # s^-1 per mRNA (RBS dependent)
    return min(k_init, f_max)

# ===================== Custom CDS support (NEW) =====================
def _clean_dna(s: str) -> str:
    return "".join(c for c in s.upper() if c in "ACGT")

def _read_seq_from_path(path: str) -> str:
    with open(path, "r", encoding="utf-8") as f:
        lines = f.read().splitlines()
    # Allow FASTA: concatenate non-header lines
    seq = "".join(line for line in lines if not line.startswith(">"))
    return _clean_dna(seq)

def _parse_user_cds_payload(raw: str) -> str:
    """
    Accepts:
      - raw DNA sequence (with or without leading ATG),
      - FASTA text,
      - or a filesystem path to a file containing the sequence.
    Returns the CDS payload *after* ATG (i.e., excludes the start codon),
    and requires a valid stop codon at the end.
    """
    raw = raw.strip()
    if not raw:
        return ""

    # If looks like a path, try to read it
    if os.path.exists(raw):
        seq = _read_seq_from_path(raw)
    else:
        # If FASTA pasted, strip headers; else treat as raw DNA
        if ">" in raw:
            seq = _clean_dna("".join(line for line in raw.splitlines() if not line.startswith(">")))
        else:
            seq = _clean_dna(raw)

    if not seq:
        return ""

    # If user included ATG, strip it off to get the payload-after-start
    if seq.startswith("ATG"):
        seq = seq[3:]

    # Basic integrity checks for payload:
    # - multiple of 3
    # - ends with a stop codon
    if len(seq) % 3 != 0:
        raise AssertionError("Custom CDS length (after removing ATG if present) is not a multiple of 3.")
    if seq[-3:] not in {"TAA", "TAG", "TGA"}:
        raise AssertionError("Custom CDS must end with a valid stop codon (TAA/TAG/TGA).")

    return seq

def _maybe_get_custom_cds(default_payload: str) -> tuple[str, str]:
    """
    Returns (cds_payload, cds_name). If user provides nothing or on error,
    falls back to the provided default_payload (sfGFP).
    """
    # Allow env var override (useful for non-interactive runs)
    env_cds = os.environ.get("CDS_SEQ", "").strip()

    try:
        if env_cds:
            user = env_cds
        else:
            try:
                user = input(
                    'Custom CDS input (paste DNA/FASTA or path). '
                    'Press ENTER to keep sfGFP: '
                ).strip()
            except EOFError:
                user = ""

        if user:
            payload = _parse_user_cds_payload(user)
            if payload:
                print(f"Using custom CDS (length {len(payload)+3} nt including ATG).")
                return payload, "CustomCDS"
    except Exception as e:
        print(f"⚠️ Custom CDS rejected: {e}. Falling back to sfGFP.")

    print("Using default sfGFP CDS.")
    return default_payload, "sfGFP"
# ===================================================================

# --- Sequence and translation setup ---
rbs = "AAAGAGGAGAA"
spacer = "ATACTAG"
start_codon = "ATG"
sfGFP_cds_start = (
"CGTAAAGGCGAAGAGCTGTTC ACTGGTGTCGTCCCTATTCTGGTGGA ACTGGATGGTGATGTCAACGGT"
"CATAAGTTTTCCGTGCGTGGCGAGGGTGAAGGTGACGCAACTAATGGTAAACTGACGCTGAAGTTCATCTGTA"
"CTACTGGTAAACTGCCGGTACCTTGGCCGACTCTGGT AACGACGCTGACTTATGGTGTTCAGTGCTTTGCTC"
"GTTATCCGGACCATATGAAGCAGCATGACTTCTTCAAGTCCGCCATGCCGGAAGGCTATGTGCAGGAACGCAC"
"GAT TTCCTTTAAGGATGACGGCACGTACAAAACGCGTGCGGAAGTGAAATTTGAAGGCGATACCCTGGTAAA"
"CCGCATTGAGCTGAAAGGCATTGACTTTAAAGAAGACGGCAATATCCTGGGCCATAAGCTGGAATACAATTTT"
"AACAGCCACAATGTTTACATCACCGCCGATAAACAAAAAAATGGCATTAAAGCGAATTTTAAAATTCGCCACA"
"ACGTGGAGGATGGCAGCGTGCAGCTGGCTGATCACTACCAGCAAAACACTCCAATCGGTGATGGTCCTGTTCT"
"GCTGCCAGACAATCACTATCTGAGCACGCAAAGCGTTCTGTCTAAAGATCCGAACGAGAAACGCGATCATATG"
"GTTCTGCTGGAGTTCGTAACCGCAGCGGGCATCACGCATGGTATGGATGAACTGTACAAATGATGA"
).replace("\n","").replace(" ", "")

rbs = rbs.upper(); spacer = spacer.upper(); start_codon = start_codon.upper()
sfGFP_cds_start = sfGFP_cds_start.upper()

# === NEW: choose CDS payload (after ATG) ===
cds_payload, cds_name = _maybe_get_custom_cds(sfGFP_cds_start)

# CDS integrity checks (including ATG)
cds_full = start_codon + cds_payload
stop_codons = {"TAA", "TAG", "TGA"}
assert len(cds_full) % 3 == 0, "CDS length is not a multiple of 3."
assert cds_full[-3:] in stop_codons, "CDS does not end with a valid stop codon."

sequence = rbs + spacer + start_codon + cds_payload
gc = calc_gc_content(sequence)
tm = melting_temp(rbs + spacer)
sd_score = check_sd_best_match(rbs)
spacing = check_spacing(len(rbs) - 1, len(rbs) + len(spacer))
binding_prob = estimate_binding_probability(sd_score, spacing)

# Lengths (bp/nt)
mRNA_length_bp = len(spacer + start_codon + cds_payload)
cds_len_nt = len(start_codon + cds_payload)

# Kinetic parameters (ABSOLUTE units)
aa_per_sec = 20                # aa/s per ribosome
ribosome_footprint_nt = 30     # nt occlusion
k_init_max = 1.0               # s^-1 per mRNA (strong RBS ceiling)

# Translation throughput (proteins per mRNA per second)
beta = compute_k_eff(cds_len_nt, aa_per_sec, ribosome_footprint_nt, binding_prob, k_init_max)

# --- Promoter library (expanded) ---
promoter_strengths = {
    "J23100": 1.0, "J23101": 0.7, "J23102": 0.86, "J23103": 0.01, "J23104": 0.72,
    "J23105": 0.24, "J23106": 0.47, "J23107": 0.36, "J23108": 0.51, "J23109": 0.04,
    "J23110": 0.33, "J23111": 0.58, "J23112": 0.00, "J23113": 0.01, "J23114": 0.10,
    "J23115": 0.15, "J23116": 0.16, "J23117": 0.17, "J23118": 0.56,
    "CustomStrong": 1.5, "CustomWeak": 0.1
}

# --- Optional: choose a subset interactively (empty input = all) ---
try:
    user_in = input('Enter promoters to test (comma-separated, e.g., "J23100,J23105"; empty = all): ').strip()
except EOFError:
    user_in = ""  # environments without stdin

if user_in:
    if user_in.lower() in ("all", "*"):
        chosen_promoters = promoter_strengths
        print("Running all promoters.")
    else:
        requested = [s.strip() for s in user_in.split(",") if s.strip()]
        chosen_promoters = {k: promoter_strengths[k] for k in requested if k in promoter_strengths}
        if not chosen_promoters:
            print("No valid promoter names provided; running all.")
            chosen_promoters = promoter_strengths
else:
    chosen_promoters = promoter_strengths
    print("Running all promoters.")

# Degradation/dilution (s^-1)
mRNA_t_half = 6*60
mRNA_decay = np.log(2)/mRNA_t_half
doubling_time = 30*60
alpha = np.log(2)/doubling_time

# ABSOLUTE transcription rate baseline (RNAs/s/cell)
k_tx_baseline = 0.02  # ~1 transcript every 50 s for baseline promoter

t_max = 5000
dt = min(5.0, 0.02 / min(mRNA_decay, alpha))
time = np.arange(0, t_max + dt, dt)

# ---------- simulation function ----------
def simulate_promoters(promoters_dict):
    """Return a dataframe of timecourses for the supplied promoters."""
    recs = []
    for promoter, rel_strength in promoters_dict.items():
        m = np.zeros_like(time, dtype=float)  # mRNA molecules per cell
        p = np.zeros_like(time, dtype=float)  # protein molecules per cell
        k_tx = k_tx_baseline * rel_strength   # RNAs/s/cell

        for i in range(1, len(time)):
            dm_dt = k_tx - mRNA_decay * m[i-1]
            m[i] = m[i-1] + dm_dt * dt
        for i in range(1, len(time)):
            dp_dt = beta * m[i-1] - alpha * p[i-1]
            p[i] = p[i-1] + dp_dt * dt

        m_eq = k_tx / mRNA_decay
        p_eq = beta * m_eq / alpha

        for t, mv, pv in zip(time, m, p):
            recs.append({
                "Time_s": t,
                "Promoter": promoter,
                "Promoter_relative_strength": rel_strength,
                "mRNA_molecules": mv,
                "Protein_molecules": pv,
                "mRNA_eq_molecules": m_eq,
                "Protein_eq_molecules": p_eq,
                "Steady_state_protein_molecules": p_eq,
                "Translation_rate_k_eff_per_mRNA_per_s": beta,
                "RBS_binding_probability": binding_prob,
                "k_tx_RNAs_per_s": k_tx  # NEW: per-timepoint k_tx in dynamics CSV
            })
    return pd.DataFrame(recs)

# --- Run initial simulation for chosen_promoters ---
df = simulate_promoters(chosen_promoters)
uniq_proms = sorted(df["Promoter"].unique())
print(f"Simulated promoters ({len(uniq_proms)}): " + ", ".join(uniq_proms[:10]) + (" …" if len(uniq_proms) > 10 else ""))

# --- Save CSVs (dynamics) ---
dyn_csv = os.path.join(OUTDIR, "promoter_dynamics.csv")
df.to_csv(dyn_csv, index=False)
print(f"CSV saved: {os.path.relpath(dyn_csv)}")

# ============== Summary table per promoter ==============
def _t50_from_trace(t, y, y_target):
    """Linear-interpolated first-passage time to y_target. Returns np.nan if never crosses."""
    y = np.asarray(y); t = np.asarray(t)
    if y_target <= 0 or np.all(y < y_target):
        return np.nan
    idx = np.searchsorted(y, y_target)  # monotone increasing in this model
    if idx == 0:
        return float(t[0])
    if idx >= len(y):
        return np.nan
    t0, t1 = t[idx-1], t[idx]
    y0, y1 = y[idx-1], y[idx]
    if y1 == y0:
        return float(t1)
    return float(t0 + (y_target - y0) * (t1 - t0) / (y1 - y0))

def build_summary_table(df):
    rows = []
    for name, g in df.sort_values("Time_s").groupby("Promoter"):
        g = g.copy()
        rs = float(g["Promoter_relative_strength"].iloc[0])
        m_eq = float(g["mRNA_eq_molecules"].iloc[0])
        p_eq = float(g["Protein_eq_molecules"].iloc[0])
        m_last = float(g["mRNA_molecules"].iloc[-1])
        p_last = float(g["Protein_molecules"].iloc[-1])
        t_arr = g["Time_s"].values
        tm50 = _t50_from_trace(t_arr, g["mRNA_molecules"].values, 0.5*m_eq) if m_eq > 0 else np.nan
        tp50 = _t50_from_trace(t_arr, g["Protein_molecules"].values, 0.5*p_eq) if p_eq > 0 else np.nan
        rows.append({
            "Promoter": name,
            "Rel_strength": rs,
            "k_tx_RNAs_per_s": k_tx_baseline * rs,
            "m_eq": m_eq,
            "p_eq": p_eq,
            "t50_mRNA_s": tm50,
            "t50_protein_s": tp50,
            "mRNA_final": m_last,
            "Protein_final": p_last,
            "beta_per_mRNA_per_s": beta,
            "P_bind": binding_prob
        })
    out = pd.DataFrame(rows).sort_values("Rel_strength", ascending=False)
    return out

summary_df = build_summary_table(df)
sum_csv = os.path.join(OUTDIR, "promoter_summary.csv")
summary_df.to_csv(sum_csv, index=False)
with pd.option_context('display.max_rows', 20, 'display.max_columns', None, 'display.width', 120):
    print("\nQuick summary (top 12 by strength):")
    print(summary_df.head(12).to_string(index=False, float_format=lambda x: f"{x:.3g}"))
print(f"CSV saved: {os.path.relpath(sum_csv)}")

# --- Anderson k_tx ladder (ranked; optional shortlist) ---
# If you want a shortlist (e.g., specific "promoter swaps [44–45]"), populate this list.
SHORTLIST = []  # e.g., ["J23100", "J23105", "J23110"]
ladder = summary_df[["Promoter", "Rel_strength", "k_tx_RNAs_per_s"]].copy()
if SHORTLIST:
    ladder = ladder[ladder["Promoter"].isin(SHORTLIST)]
ladder = ladder.sort_values("k_tx_RNAs_per_s", ascending=False).reset_index(drop=True)
ladder["Rank"] = ladder.index + 1

ladder_csv = os.path.join(OUTDIR, "anderson_k_tx_ladder.csv")
ladder.to_csv(ladder_csv, index=False)
print(f"CSV saved: {os.path.relpath(ladder_csv)}")

# Add a simple ladder plot to the PDF (and a PNG)
fig, ax = plt.subplots(figsize=(8, 0.4 * max(len(ladder), 1) + 2))
ax.barh(ladder["Promoter"], ladder["k_tx_RNAs_per_s"])
ax.invert_yaxis()
ax.set_xlabel("k_tx (RNAs/s/cell)")
ax.set_title(f"Anderson k_tx ladder — {cds_name}")
_save_figure(fig, "anderson_k_tx_ladder", OUTDIR, _pdf)
plt.close(fig)

# ================= Plots (save PNG + PDF, and optionally SHOW) =================
def plot_stacked_promoter_panels(df, promoters_to_plot, *, outdir=OUTDIR, pdf=_pdf,
                                 save_name=None, show=SHOW_PLOTS):
    available = sorted(set(df["Promoter"].unique()))
    keep = [p for p in promoters_to_plot if p in available]
    missing = [p for p in promoters_to_plot if p not in available]
    if missing:
        print(f"⚠️ Not found in current dataframe: {', '.join(missing)}")
        print(f"Available: {', '.join(available)}")
        return
    n = len(keep)
    fig, axes = plt.subplots(nrows=n, ncols=1, figsize=(9, 3.6*n), sharex=True)
    if n == 1:
        axes = [axes]
    for ax, name in zip(axes, keep):
        g = df[df["Promoter"] == name].sort_values("Time_s")
        ax.plot(g["Time_s"], g["mRNA_molecules"], label="mRNA", color="tab:blue", lw=1.8)
        ax.set_ylabel("mRNA (molecules/cell)", color="tab:blue")
        ax.tick_params(axis='y', labelcolor="tab:blue")
        ax2 = ax.twinx()
        ax2.plot(g["Time_s"], g["Protein_molecules"], label="Protein", color="tab:green", lw=1.8)
        ax2.set_ylabel("Protein (molecules/cell)", color="tab:green")
        ax2.tick_params(axis='y', labelcolor="tab:green")
        ax.set_title(f"{name} dynamics")
        h1,l1=ax.get_legend_handles_labels(); h2,l2=ax2.get_legend_handles_labels()
        ax.legend(h1+h2, l1+l2, loc="upper left")
        ax.grid(True, alpha=0.25)
    axes[-1].set_xlabel("Time (s)")
    plt.tight_layout()
    base = save_name or ("stacked_" + "_".join(keep[:2]) + ("" if n <= 2 else "_etc"))
    _save_figure(fig, base, outdir, pdf)
    if show and SHOW_PLOTS:
        plt.show()
    else:
        plt.close(fig)

# ---------- Overlay ANY number of promoters (legend panel to avoid overlap) ----------
def plot_overlay_promoters(df, promoters, *, outdir=OUTDIR, pdf=_pdf,
                           save_name=None, show=SHOW_PLOTS,
                           legend_cols=None, legend_fontsize=8,
                           legend_title="Promoters", logy_mrna=False, logy_protein=False):
    """Overlay N promoters. Color encodes promoter; solid = mRNA (left y), dashed = protein (right y)."""
    if isinstance(promoters, str):
        promoters = [promoters]
    all_avail = set(df["Promoter"].unique())
    keep = [p for p in promoters if p in all_avail]
    if len(keep) < 2:
        print("Need at least two valid promoters to overlay.")
        return

    palette = plt.get_cmap("tab20")(np.linspace(0, 1, max(3, len(keep))))
    color_map = {p: palette[i % len(palette)] for i, p in enumerate(keep)}

    fig, ax = plt.subplots(figsize=(10, 5.6))
    ax2 = ax.twinx()

    for p in keep:
        g = df[df["Promoter"] == p].sort_values("Time_s")
        c = color_map[p]
        ax.plot(g["Time_s"], g["mRNA_molecules"], color=c, lw=1.8, label=p)                # mRNA solid
        ax2.plot(g["Time_s"], g["Protein_molecules"], color=c, lw=1.6, ls="--",
                 label="_nolegend_")                                                       # protein dashed

    ax.set_ylabel("mRNA (molecules/cell)")
    ax2.set_ylabel("Protein (molecules/cell)", labelpad=10)
    ax.set_xlabel("Time (s)")
    if logy_mrna: ax.set_yscale("log")
    if logy_protein: ax2.set_yscale("log")
    ax.grid(True, alpha=0.25)
    ax.set_title(f"Overlay comparison ({len(keep)} promoters)")

    # Style legend (solid vs dashed) inside the axes
    style_handles = [Line2D([0], [0], color="k", lw=2, ls="-"),
                     Line2D([0], [0], color="k", lw=2, ls="--")]
    ax.legend(style_handles, ["mRNA (left y)", "Protein (right y)"], loc="upper left", frameon=True)

    # Promoter legend in dedicated side panel
    if legend_cols is None:
        n = len(keep)
        legend_cols = 1 if n <= 16 else (2 if n <= 32 else 3)
    prom_handles = [Line2D([0], [0], color=color_map[p], lw=2) for p in keep]
    _add_right_side_legend(fig, prom_handles, keep,
                           ncols=legend_cols, title=legend_title,
                           fontsize=legend_fontsize, title_fontsize=legend_fontsize+1)

    base = save_name or ("overlay_" + "_".join(keep[:6]) + ("_etc" if len(keep) > 6 else ""))
    _save_figure(fig, base, outdir, pdf)
    if show and SHOW_PLOTS:
        plt.show()
    else:
        plt.close(fig)

# Back-compat wrapper (two promoters) calls the N-promoter version
def plot_overlay_two_promoters(df, p1, p2, **kwargs):
    return plot_overlay_promoters(df, [p1, p2], **kwargs)

def plot_all_promoters_overlay(df, promoters=None, show_mrna=True, show_protein=True,
                               top_n=None, logy=False,
                               legend_cols=None, legend_fontsize=8, legend_title="Promoters",
                               *, outdir=OUTDIR, pdf=_pdf, save_name="all_promoters_overlay",
                               show=SHOW_PLOTS):
    """All-promoter overlay with a dedicated right-side legend panel."""
    all_proms = sorted(df["Promoter"].unique())
    promoters = all_proms if (promoters is None) else [p for p in promoters if p in all_proms]
    if top_n is not None and top_n < len(promoters):
        end = df.sort_values("Time_s").groupby("Promoter")["Protein_molecules"].last()
        promoters = end.loc[promoters].sort_values(ascending=False).head(top_n).index.tolist()
    if not promoters:
        print("No promoters to plot.")
        return

    palette = plt.get_cmap("tab20")(np.linspace(0, 1, max(3, len(promoters))))
    color_map = {p: palette[i % len(palette)] for i, p in enumerate(promoters)}

    n_rows = 2 if (show_mrna and show_protein) else 1
    fig, axes = plt.subplots(n_rows, 1, figsize=(10, 3.6*n_rows), sharex=True)
    if n_rows == 1:
        axes = [axes]

    if show_mrna:
        ax = axes[0]
        for p in promoters:
            g = df[df["Promoter"] == p].sort_values("Time_s")
            ax.plot(g["Time_s"], g["mRNA_molecules"], lw=1.3, alpha=0.9, label=p, color=color_map[p])
        ax.set_ylabel("mRNA (molecules/cell)")
        if logy: ax.set_yscale("log")
        ax.grid(True, alpha=0.25)
        ax.set_title("All promoters — overlay (mRNA and protein)")

    if show_protein:
        axp = axes[-1]
        for p in promoters:
            g = df[df["Promoter"] == p].sort_values("Time_s")
            axp.plot(g["Time_s"], g["Protein_molecules"], lw=1.3, alpha=0.9, label=p, color=color_map[p])
        axp.set_ylabel("Protein (molecules/cell)", labelpad=10)
        axp.set_xlabel("Time (s)")
        if logy: axp.set_yscale("log")
        axp.grid(True, alpha=0.25)

    # Unique labels/handles from bottom axis
    handles, labels = axes[-1].get_legend_handles_labels()
    uniq, seen = [], set()
    for h, l in zip(handles, labels):
        if l not in seen and l != "_nolegend_":
            uniq.append((h, l)); seen.add(l)
    handles, labels = zip(*uniq) if uniq else ([], [])

    if legend_cols is None:
        n_proms = len(labels)
        legend_cols = 1 if n_proms <= 16 else (2 if n_proms <= 32 else 3)

    _add_right_side_legend(fig, list(handles), list(labels),
                           ncols=legend_cols, title=legend_title,
                           fontsize=legend_fontsize, title_fontsize=legend_fontsize+1)

    _save_figure(fig, save_name, outdir, pdf)
    if show and SHOW_PLOTS:
        plt.show()
    else:
        plt.close(fig)

# ---------- Paged gallery so ALL selected promoters are shown ----------
def plot_gallery_promoters(df, promoters, per_fig=4, *, outdir=OUTDIR, pdf=_pdf, show=SHOW_PLOTS):
    """
    Paged gallery: stacked panels of mRNA (left y) + protein (right y)
    for 'per_fig' promoters per figure. Saves PNGs, adds to PDF, and can show in GUI.
    """
    all_avail = set(df["Promoter"].unique())
    proms = [p for p in promoters if p in all_avail]
    if not proms:
        print("No valid promoters to plot in gallery.")
        return

    for i in range(0, len(proms), per_fig):
        chunk = proms[i:i+per_fig]
        fig, axes = plt.subplots(nrows=len(chunk), ncols=1,
                                 figsize=(9, 3.6*len(chunk)), sharex=True)
        if len(chunk) == 1:
            axes = [axes]

        for ax, name in zip(axes, chunk):
            g = df[df["Promoter"] == name].sort_values("Time_s")
            ax.plot(g["Time_s"], g["mRNA_molecules"], label="mRNA", color="tab:blue", lw=1.8)
            ax.set_ylabel("mRNA (molecules/cell)", color="tab:blue")
            ax.tick_params(axis='y', labelcolor="tab:blue")
            ax2 = ax.twinx()
            ax2.plot(g["Time_s"], g["Protein_molecules"], label="Protein", color="tab:green", lw=1.8)
            ax2.set_ylabel("Protein (molecules/cell)", color="tab:green")
            ax2.tick_params(axis='y', labelcolor="tab:green")
            ax.set_title(f"{name} dynamics")
            h1,l1 = ax.get_legend_handles_labels()
            h2,l2 = ax2.get_legend_handles_labels()
            ax.legend(h1+h2, l1+l2, loc="upper left")
            ax.grid(True, alpha=0.25)

        axes[-1].set_xlabel("Time (s)")
        plt.tight_layout()
        page = i // per_fig + 1
        _save_figure(fig, f"gallery_promoters_page_{page:02d}", outdir, pdf)
        if show and SHOW_PLOTS:
            plt.show()
        else:
            plt.close(fig)

# --------- Auto-generate a quicklook (shows ALL selected promoters, paged) ----------
selected_promoters = list(chosen_promoters.keys())
if selected_promoters:
    plot_gallery_promoters(df, selected_promoters, per_fig=4, show=SHOW_PLOTS)

# Optional interactive compare: overlay N promoters
try:
    ov_in = input('Overlay promoters (comma-separated; 2+ names; empty = skip): ').strip()
except EOFError:
    ov_in = ""
if ov_in:
    parts = [s.strip() for s in ov_in.split(",") if s.strip()]
    current = set(df["Promoter"].unique())
    missing = [p for p in parts if p not in current and p in promoter_strengths]
    if missing:
        df = pd.concat([df, simulate_promoters({k: promoter_strengths[k] for k in missing})],
                       ignore_index=True)
    if len(parts) >= 2:
        save_stub = "overlay_" + "_".join(parts[:6]) + ("_etc" if len(parts) > 6 else "")
        plot_overlay_promoters(df, parts, save_name=save_stub, show=SHOW_PLOTS)
    else:
        print("Need at least two valid promoter names for overlay. Skipping.")

# Optional: ALL-IN-ONE overlay for all promoters
try:
    allplot_in = input('Plot ALL promoters together (y/n, default n)? ').strip().lower()
except EOFError:
    allplot_in = ""
if allplot_in in ("y", "yes"):
    current = set(df["Promoter"].unique())
    missing = [k for k in promoter_strengths if k not in current]
    if missing:
        df = pd.concat([df, simulate_promoters({k: promoter_strengths[k] for k in missing})],
                       ignore_index=True)
    plot_all_promoters_overlay(df, promoters=None, show_mrna=True, show_protein=True,
                               top_n=None, logy=False, save_name="all_promoters_overlay",
                               show=SHOW_PLOTS)

# Finalize PDF and print summary
_pdf.close()
print("\n========================")
print(f"CDS: {cds_name} | length (incl. ATG) = {cds_len_nt} nt | GC% (UTR+CDS) = {gc:.2f}")
print(f"All outputs saved in: {os.path.relpath(OUTDIR)}")
print(f"PDF of all figures:   {os.path.relpath(PDF_PATH)}")
print("========================") 
