# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

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

# CDS integrity checks
cds_full = start_codon + sfGFP_cds_start
stop_codons = {"TAA", "TAG", "TGA"}
assert len(cds_full) % 3 == 0, "CDS length is not a multiple of 3."
assert cds_full[-3:] in stop_codons, "CDS does not end with a valid stop codon."

sequence = rbs + spacer + start_codon + sfGFP_cds_start
gc = calc_gc_content(sequence)
tm = melting_temp(rbs + spacer)
sd_score = check_sd_best_match(rbs)
spacing = check_spacing(len(rbs) - 1, len(rbs) + len(spacer))
binding_prob = estimate_binding_probability(sd_score, spacing)

# Lengths (bp/nt)
mRNA_length_bp = len(spacer + start_codon + sfGFP_cds_start)
cds_len_nt = len(start_codon + sfGFP_cds_start)

# Kinetic parameters (ABSOLUTE units)
aa_per_sec = 20                # aa/s per ribosome
ribosome_footprint_nt = 30     # nt occlusion
k_init_max = 1.0               # s^-1 per mRNA (strong RBS ceiling)

# Translation throughput (proteins per mRNA per second)
beta = compute_k_eff(cds_len_nt, aa_per_sec, ribosome_footprint_nt, binding_prob, k_init_max)

# --- Promoter library (expanded) ---
promoter_strengths = {  # relative multipliers of k_tx_baseline
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

# accept "all"/"*", and print what we'll do
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

t_max = 1000
# Adaptive dt (cap at 5 s); include the final point at 1000 s
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
                "RBS_binding_probability": binding_prob
            })
    return pd.DataFrame(recs)

# --- Run initial simulation for chosen_promoters
df = simulate_promoters(chosen_promoters)

# --- Save CSV ---
df.to_csv("promoter_dynamics.csv", index=False)
print("CSV saved as promoter_dynamics.csv")

# === Stacked comparison util ===
def plot_stacked_promoter_panels(df, promoters_to_plot):
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
    plt.show()

# === Overlay comparison util ===
def plot_overlay_two_promoters(df, p1, p2):
    assert p1 != p2, "Choose two different promoters."
    g1 = df[df["Promoter"] == p1].sort_values("Time_s")
    g2 = df[df["Promoter"] == p2].sort_values("Time_s")
    if g1.empty or g2.empty:
        print("One or both promoters not in dataframe.")
        return
    fig, ax = plt.subplots(figsize=(9, 4.5))
    ax.plot(g1["Time_s"], g1["mRNA_molecules"], color="tab:blue", label=f"mRNA ({p1})", lw=1.8)
    ax.plot(g2["Time_s"], g2["mRNA_molecules"], color="tab:orange", label=f"mRNA ({p2})", lw=1.8, ls="--")
    ax.set_ylabel("mRNA (molecules/cell)")
    ax2 = ax.twinx()
    ax2.plot(g1["Time_s"], g1["Protein_molecules"], color="tab:green", label=f"Protein ({p1})", lw=1.8)
    ax2.plot(g2["Time_s"], g2["Protein_molecules"], color="tab:red", label=f"Protein ({p2})", lw=1.8, ls="--")
    ax2.set_ylabel("Protein (molecules/cell)")
    ax.set_xlabel("Time (s)")
    ax.grid(True, alpha=0.25)
    h1,l1=ax.get_legend_handles_labels(); h2,l2=ax2.get_legend_handles_labels()
    ax.legend(h1+h2, l1+l2, loc="upper left")
    plt.title(f"Overlay comparison: {p1} vs {p2}")
    plt.tight_layout()
    plt.show()

# === NEW: One clean comparison plot (all promoters together) ===
def plot_all_promoters_overlay(df, promoters=None, show_mrna=True, show_protein=True, top_n=None):
    """
    One figure with all promoters together.
    - If top_n is set, keeps the top_n promoters by final protein level (reduces clutter).
    """
    # choose promoters
    all_proms = sorted(df["Promoter"].unique())
    if promoters is None:
        promoters = all_proms
    else:
        promoters = [p for p in promoters if p in all_proms]

    if top_n is not None and top_n < len(promoters):
        # keep top_n by endpoint protein
        end = df.sort_values("Time_s").groupby("Promoter")["Protein_molecules"].last()
        keep = end.loc[promoters].sort_values(ascending=False).head(top_n).index.tolist()
        promoters = keep

    if not promoters:
        print("No promoters to plot.")
        return

    # color map consistent across panels
    palette = plt.cm.tab20(np.linspace(0, 1, max(3, len(promoters))))
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
        ax.grid(True, alpha=0.25)
        ax.set_title("All promoters — overlay (mRNA and protein)")

    if show_protein:
        axp = axes[-1]
        for p in promoters:
            g = df[df["Promoter"] == p].sort_values("Time_s")
            axp.plot(g["Time_s"], g["Protein_molecules"], lw=1.3, alpha=0.9, label=p, color=color_map[p])
        axp.set_ylabel("Protein (molecules/cell)")
        axp.set_xlabel("Time (s)")
        axp.grid(True, alpha=0.25)

    # single shared legend to the right
    handles, labels = axes[-1].get_legend_handles_labels()
    fig.legend(handles, labels, loc="center left", bbox_to_anchor=(1.02, 0.5), borderaxespad=0.)
    fig.tight_layout(rect=[0, 0, 0.84, 1])
    plt.show()

# Quicklook after simulation (first 1–2 promoters)
default_to_plot = list(chosen_promoters.keys())[:2]
if default_to_plot:
    plot_stacked_promoter_panels(df, default_to_plot)

# --- Optional interactive compare: stacked panels (with backfill)
try:
    cmp_in = input('Compare two promoters (comma-separated, e.g., "J23100,J23105"; empty = skip): ').strip()
except EOFError:
    cmp_in = ""
if cmp_in:
    cmp_proms = [s.strip() for s in cmp_in.split(",") if s.strip()]
    current = set(df["Promoter"].unique())
    missing = [p for p in cmp_proms if p not in current and p in promoter_strengths]
    if missing:
        df = pd.concat([df, simulate_promoters({k: promoter_strengths[k] for k in missing})],
                       ignore_index=True)
    if len(cmp_proms) >= 2:
        plot_stacked_promoter_panels(df, cmp_proms[:2])
    else:
        print("Need at least two valid promoter names to compare. Skipping.")

# --- Optional interactive compare: overlay (with backfill)
try:
    ov_in = input('Overlay two promoters on one plot (comma-separated; empty = skip): ').strip()
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
        plot_overlay_two_promoters(df, parts[0], parts[1])
    else:
        print("Need two valid promoter names for overlay. Skipping.")

# --- Optional: ALL-IN-ONE overlay for all promoters ---
try:
    allplot_in = input('Plot ALL promoters together (y/n, default n)? ').strip().lower()
except EOFError:
    allplot_in = ""
if allplot_in in ("y", "yes"):
    # tip: set top_n (e.g., 8) if the figure is too busy
    plot_all_promoters_overlay(df, promoters=list(chosen_promoters.keys()), show_mrna=True, show_protein=True, top_n=None)
