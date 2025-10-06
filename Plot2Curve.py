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
import math  # analytic t50 fallback

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
    return {1: 0.18, 2: 0.26, 3: 0.34, 4: 0.42}.get(int(ncols), 0.42)

def _add_right_side_legend(fig, handles, labels, *, ncols=1, title="Promoters",
                           fontsize=8, title_fontsize=9, gutter=0.06):
    panel_w = _panel_width_from_cols(ncols)
    fig.subplots_adjust(right=1 - (panel_w + gutter))
    sp = fig.subplotpars
    leg_ax = fig.add_axes([1 - panel_w + 0.01, sp.bottom, panel_w - 0.02, sp.top - sp.bottom])
    leg_ax.axis("off")
    leg_ax.legend(handles, labels, loc="center left", frameon=True, ncol=ncols,
                  title=title, fontsize=fontsize, title_fontsize=title_fontsize)
    return leg_ax

# --- Helper Functions ---
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
    v_nt = 3 * aa_per_sec
    f_max = v_nt / footprint_nt
    k_init = k_init_max * P_bind
    return min(k_init, f_max)

# ===================== Custom CDS support =====================
def _clean_dna(s: str) -> str:
    return "".join(c for c in s.upper() if c in "ACGT")

def _read_seq_from_path(path: str) -> str:
    with open(path, "r", encoding="utf-8") as f:
        lines = f.read().splitlines()
    seq = "".join(line for line in lines if not line.startswith(">"))
    return _clean_dna(seq)

def _parse_user_cds_payload(raw: str) -> str:
    raw = raw.strip()
    if not raw:
        return ""
    if os.path.exists(raw):
        seq = _read_seq_from_path(raw)
    else:
        if ">" in raw:
            seq = _clean_dna("".join(line for line in raw.splitlines() if not line.startswith(">")))
        else:
            seq = _clean_dna(raw)
    if not seq:
        return ""
    if seq.startswith("ATG"):
        seq = seq[3:]
    if len(seq) % 3 != 0:
        raise AssertionError("Custom CDS length (after removing ATG if present) is not a multiple of 3.")
    if seq[-3:] not in {"TAA", "TAG", "TGA"}:
        raise AssertionError("Custom CDS must end with a valid stop codon (TAA/TAG/TGA).")
    return seq

def _maybe_get_custom_cds(default_payload: str) -> tuple[str, str]:
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

# Choose CDS payload (after ATG)
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

# Lengths
mRNA_length_bp = len(spacer + start_codon + cds_payload)
cds_len_nt = len(start_codon + cds_payload)

# Kinetic parameters (ABSOLUTE units)
aa_per_sec = 20
ribosome_footprint_nt = 30
k_init_max = 1.0

# Translation throughput (proteins per mRNA per second)
beta = compute_k_eff(cds_len_nt, aa_per_sec, ribosome_footprint_nt, binding_prob, k_init_max)

# === NEW: protein maturation (sfGFP realism) ===
k_mat = np.log(2)/420.0  # 7 min maturation half-time example

# === Resource coupling (ribosome budget; V1 global throttling) ===
# Env overrides: THROTTLE=0 disables; K_R sets budget scale.
ENABLE_RESOURCE_COUPLING = os.environ.get("THROTTLE", "1") != "0"
K_R = float(os.environ.get("K_R", "100.0"))  # budget scale (mRNA molecules/cell); tune 20–200+

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
    user_in = ""

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
# Include k_mat in step-size control for stability with maturation
dt = min(5.0, 0.02 / min(mRNA_decay, alpha, k_mat))
time = np.arange(0, t_max + dt, dt)

# ---------- simulation function (original: no coupling) ----------
def simulate_promoters(promoters_dict):
    """Return a dataframe of timecourses for the supplied promoters."""
    recs = []
    for promoter, rel_strength in promoters_dict.items():
        m  = np.zeros_like(time, dtype=float)  # mRNA molecules per cell
        p_i = np.zeros_like(time, dtype=float) # immature protein
        p_m = np.zeros_like(time, dtype=float) # mature protein (reported fluorescence)
        k_tx = k_tx_baseline * rel_strength    # RNAs/s/cell

        # mRNA dynamics
        for i in range(1, len(time)):
            dm_dt = k_tx - mRNA_decay * m[i-1]
            m[i] = m[i-1] + dm_dt * dt

        # Protein dynamics with maturation (no throttling)
        for i in range(1, len(time)):
            dp_i = beta * m[i-1] - (alpha + k_mat) * p_i[i-1]
            dp_m = k_mat * p_i[i-1] - alpha * p_m[i-1]
            p_i[i] = p_i[i-1] + dp_i * dt
            p_m[i] = p_m[i-1] + dp_m * dt

        # Steady states
        m_eq   = k_tx / mRNA_decay
        p_i_eq = beta * m_eq / (alpha + k_mat)
        p_m_eq = (k_mat / alpha) * p_i_eq  # = beta*m_eq*k_mat / [alpha*(alpha + k_mat)]

        for t, mv, piv, pmv in zip(time, m, p_i, p_m):
            recs.append({
                "Time_s": t,
                "Promoter": promoter,
                "Promoter_relative_strength": rel_strength,
                "mRNA_molecules": mv,
                "Protein_molecules": pmv,               # mature protein reported
                "Protein_immature_molecules": piv,       # NEW
                "mRNA_eq_molecules": m_eq,
                "Protein_eq_molecules": p_m_eq,          # mature steady state
                "Steady_state_protein_molecules": p_m_eq,
                "Translation_rate_k_eff_per_mRNA_per_s": beta,
                "RBS_binding_probability": binding_prob,
                "k_tx_RNAs_per_s": k_tx
            })
    return pd.DataFrame(recs)

# ---------- Coupled simulator (V1 global scalar throttling) ----------
def simulate_promoters_coupled(promoters_dict, K_R=100.0):
    """
    Simulate all promoters simultaneously with a shared ribosome budget:
      beta_eff(t) = beta / (1 + M_tot(t)/K_R),  M_tot = sum_j m_j
    Returns (df, diag) with per-step diagnostics recorded in df.
    """
    names = list(promoters_dict.keys())
    rs = np.array([promoters_dict[n] for n in names], dtype=float)
    n = len(names); T = len(time)

    # States
    m   = np.zeros((n, T), dtype=float)
    p_i = np.zeros((n, T), dtype=float)
    p_m = np.zeros((n, T), dtype=float)

    k_tx = k_tx_baseline * rs  # RNAs/s/cell

    # Global diagnostics arrays
    M_tot_t    = np.zeros(T, dtype=float)
    beta_eff_t = np.zeros(T, dtype=float)
    phi_t      = np.ones(T, dtype=float)  # φ(t) = β_eff/β

    # mRNA dynamics (independent)
    for i in range(1, T):
        dm = k_tx - mRNA_decay * m[:, i-1]
        m[:, i] = m[:, i-1] + dm * dt

    # Protein dynamics with time-varying beta_eff(t)
    for i in range(1, T):
        M_tot = m[:, i-1].sum()
        M_tot_t[i] = M_tot
        phi = 1.0 / (1.0 + M_tot / K_R)
        beta_eff = beta * phi
        phi_t[i] = phi
        beta_eff_t[i] = beta_eff

        dp_i = beta_eff * m[:, i-1] - (alpha + k_mat) * p_i[:, i-1]
        dp_m = k_mat * p_i[:, i-1] - alpha * p_m[:, i-1]
        p_i[:, i] = p_i[:, i-1] + dp_i * dt
        p_m[:, i] = p_m[:, i-1] + dp_m * dt

    # Steady states using coupled beta_eff at M_tot_eq
    m_eq = k_tx / mRNA_decay
    M_tot_eq = float(m_eq.sum())
    phi_eq = 1.0 / (1.0 + M_tot_eq / K_R)
    beta_eff_eq = beta * phi_eq
    p_i_eq = beta_eff_eq * m_eq / (alpha + k_mat)
    p_m_eq = (k_mat / alpha) * p_i_eq

    recs = []
    for j, name in enumerate(names):
        for i_t, t in enumerate(time):
            recs.append({
                "Time_s": t,
                "Promoter": name,
                "Promoter_relative_strength": float(rs[j]),
                "mRNA_molecules": m[j, i_t],
                "Protein_molecules": p_m[j, i_t],
                "Protein_immature_molecules": p_i[j, i_t],
                "mRNA_eq_molecules": float(m_eq[j]),
                "Protein_eq_molecules": float(p_m_eq[j]),
                "Steady_state_protein_molecules": float(p_m_eq[j]),
                "Translation_rate_k_eff_per_mRNA_per_s": beta,  # base β (reference)
                "RBS_binding_probability": binding_prob,
                "k_tx_RNAs_per_s": float(k_tx[j]),
                # --- GLOBAL DIAGNOSTICS (same for all promoters at a timepoint) ---
                "M_tot_all_promoters": float(M_tot_t[i_t]),
                "beta_eff_current": float(beta_eff_t[i_t]),
                "throttle_phi": float(phi_t[i_t]),
                "phi_eq_global": float(phi_eq),
                "K_R": float(K_R),
                "M_tot_eq": float(M_tot_eq)
            })

    diag = {
        "M_tot_eq": M_tot_eq,
        "phi_eq": phi_eq,
        "phi_mean": float(phi_t[1:].mean()),
        "phi_min": float(phi_t[1:].min())
    }
    return pd.DataFrame(recs), diag

# --- Run initial simulation for chosen_promoters ---
if ENABLE_RESOURCE_COUPLING:
    df, throttle_diag = simulate_promoters_coupled(chosen_promoters, K_R=K_R)
else:
    df = simulate_promoters(chosen_promoters)
    throttle_diag = None

uniq_proms = sorted(df["Promoter"].unique())
print(f"Simulated promoters ({len(uniq_proms)}): " + ", ".join(uniq_proms[:10]) + (" …" if len(uniq_proms) > 10 else ""))
print(f"Ribosome throttling: {'ON' if ENABLE_RESOURCE_COUPLING else 'OFF'}  (K_R={K_R:g})")
if throttle_diag is not None:
    print(f"  M_tot_eq≈{throttle_diag['M_tot_eq']:.2f} mRNA; φ_eq≈{throttle_diag['phi_eq']:.3f}; "
          f"⟨φ(t)⟩≈{throttle_diag['phi_mean']:.3f} (min {throttle_diag['phi_min']:.3f})")

# --- Save CSVs (dynamics) ---
dyn_csv = os.path.join(OUTDIR, "promoter_dynamics.csv")
df.to_csv(dyn_csv, index=False)
print(f"CSV saved: {os.path.relpath(dyn_csv)}")

# ============== Summary helpers ==============
def _t50_from_trace(t, y, y_target):
    """Linear-interpolated first-passage time to y_target. Returns np.nan if never crosses."""
    y = np.asarray(y); t = np.asarray(t)
    if y_target <= 0 or np.all(y < y_target):
        return np.nan
    idx = np.searchsorted(y, y_target)
    if idx == 0:
        return float(t[0])
    if idx >= len(y):
        return np.nan
    t0, t1 = t[idx-1], t[idx]
    y0, y1 = y[idx-1], y[idx]
    if y1 == y0:
        return float(t1)
    return float(t0 + (y_target - y0) * (t1 - t0) / (y1 - y0))

# === NEW: value interpolation at a specific time (for fractions, AUC windows, etc.)
def _interp_at_time(t_arr, y_arr, t_query):
    t = np.asarray(t_arr); y = np.asarray(y_arr)
    if t_query <= t[0]: return float(y[0])
    if t_query >= t[-1]: return float(y[-1])
    idx = np.searchsorted(t, t_query)
    t0, t1 = t[idx-1], t[idx]
    y0, y1 = y[idx-1], y[idx]
    if t1 == t0: return float(y1)
    return float(y0 + (y1 - y0) * (t_query - t0) / (t1 - t0))

# ---------- Analytic fallback helpers for protein t50 (no maturation) ----------
def _protein_fraction(t, alpha, gamma):
    if abs(alpha - gamma) < 1e-12:
        return 1.0 - (1.0 + alpha*t) * math.exp(-alpha*t)
    return 1.0 - (alpha*math.exp(-gamma*t) - gamma*math.exp(-alpha*t)) / (alpha - gamma)

def analytic_t50_protein(alpha, gamma):
    lo, hi = 0.0, 20.0 / max(alpha, gamma)
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        if _protein_fraction(mid, alpha, gamma) < 0.5:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)

# ---------- Numeric time-to-fraction WITH maturation ----------
def numeric_tfraction_protein_mature(k_tx, gamma, alpha, beta, k_mat, frac):
    """
    Simulate m, p_i, p_m with a fine step and return t where p_m reaches 'frac' of its steady state.
    Robust for maturation; used as fallback if coarse trace doesn't cross within main t_max.
    """
    m_eq = k_tx / gamma
    p_m_eq = beta * m_eq * k_mat / (alpha * (alpha + k_mat))
    if p_m_eq <= 0 or frac <= 0:
        return np.nan
    dt_fine = min(1.0, 0.02 / min(gamma, alpha, k_mat))
    t_cap = 20.0 / min(gamma, alpha, k_mat)
    t = 0.0
    m = 0.0; p_i = 0.0; p_m = 0.0
    prev_t, prev_pm = 0.0, 0.0
    target = frac * p_m_eq
    while t < t_cap:
        prev_t, prev_pm = t, p_m
        dm = k_tx - gamma * m
        m += dm * dt_fine
        dp_i = beta * m - (alpha + k_mat) * p_i
        p_i += dp_i * dt_fine
        dp_m = k_mat * p_i - alpha * p_m
        p_m += dp_m * dt_fine
        t += dt_fine
        if p_m >= target:
            if p_m == prev_pm:
                return float(t)
            frac_lin = (target - prev_pm) / (p_m - prev_pm)
            return float(prev_t + frac_lin * (t - prev_t))
    return np.nan

# Back-compat wrapper (keeps any external calls intact)
def numeric_t50_protein_mature(k_tx, gamma, alpha, beta, k_mat):
    return numeric_tfraction_protein_mature(k_tx, gamma, alpha, beta, k_mat, 0.5)

# ============== Summary table per promoter ==============
def build_summary_table(df):
    rows = []
    for name, g in df.sort_values("Time_s").groupby("Promoter"):
        g = g.copy()
        rs   = float(g["Promoter_relative_strength"].iloc[0])
        m_eq = float(g["mRNA_eq_molecules"].iloc[0])
        p_eq = float(g["Protein_eq_molecules"].iloc[0])  # mature steady state
        m_last = float(g["mRNA_molecules"].iloc[-1])
        p_last = float(g["Protein_molecules"].iloc[-1])  # mature
        t_arr  = g["Time_s"].values

        # Effective beta at steady state (if coupling present)
        beta_eff_eq_here = float(g["beta_eff_current"].iloc[-1]) if "beta_eff_current" in g.columns and np.isfinite(g["beta_eff_current"].iloc[-1]) else float(g["beta_eff_eq"].iloc[0]) if "beta_eff_eq" in g.columns else beta
        phi_eq_here = float(g["phi_eq_global"].iloc[0]) if "phi_eq_global" in g.columns else (beta_eff_eq_here / beta)

        # t50 from traces
        tm50 = _t50_from_trace(t_arr, g["mRNA_molecules"].values, 0.5*m_eq) if m_eq > 0 else np.nan
        tp50 = _t50_from_trace(t_arr, g["Protein_molecules"].values, 0.5*p_eq) if p_eq > 0 else np.nan

        # Fallbacks if protein trace didn't cross within t_max
        if np.isnan(tp50) and p_eq > 0:
            if k_mat <= 0:
                tp50 = analytic_t50_protein(alpha, mRNA_decay)
            else:
                k_tx_here = k_tx_baseline * rs
                tp50 = numeric_t50_protein_mature(k_tx_here, mRNA_decay, alpha, beta_eff_eq_here, k_mat)

        # === Additional fluorescence-centric metrics ===
        t10_pm = _t50_from_trace(t_arr, g["Protein_molecules"].values, 0.1*p_eq) if p_eq > 0 else np.nan
        t80_pm = _t50_from_trace(t_arr, g["Protein_molecules"].values, 0.8*p_eq) if p_eq > 0 else np.nan
        t90_pm = _t50_from_trace(t_arr, g["Protein_molecules"].values, 0.9*p_eq) if p_eq > 0 else np.nan

        # Fallback for t80 if needed (ensures populated speed metric)
        if np.isnan(t80_pm) and p_eq > 0:
            k_tx_here = k_tx_baseline * rs
            t80_pm = numeric_tfraction_protein_mature(k_tx_here, mRNA_decay, alpha, beta_eff_eq_here, k_mat, 0.8)

        maturation_lag = (tp50 - tm50) if (not np.isnan(tp50) and not np.isnan(tm50)) else np.nan
        frac_mature_1800 = (_interp_at_time(t_arr, g["Protein_molecules"].values, 1800.0)/p_eq) if p_eq > 0 else np.nan
        auc_pm = float(np.trapz(g["Protein_molecules"].values, t_arr))  # molecules*s

        # Immature steady state (for CSV completeness)
        p_i_eq = beta_eff_eq_here * m_eq / (alpha + k_mat)

        rows.append({
            "Promoter": name,
            "Rel_strength": rs,
            "k_tx_RNAs_per_s": k_tx_baseline * rs,
            "m_eq": m_eq,
            "p_eq": p_eq,                               # mature steady state
            "Protein_mature_eq": p_eq,                  # alias for readability
            "Protein_immature_eq": p_i_eq,
            "t50_mRNA_s": tm50,
            "t50_protein_s": tp50,
            "t10_p_m_s": t10_pm,
            "t80_p_m_s": t80_pm,                        # NEW: preferred speed metric
            "t90_p_m_s": t90_pm,
            "maturation_lag_s": maturation_lag,
            "fraction_mature_at_1800s": frac_mature_1800,
            "AUC_p_m_0_tmax": auc_pm,
            "mRNA_final": m_last,
            "Protein_final": p_last,                    # mature at t_max
            "beta_per_mRNA_per_s": beta,
            "beta_eff_eq": beta_eff_eq_here,
            "phi_eq": phi_eq_here,
            "K_R": float(g["K_R"].iloc[0]) if "K_R" in g.columns else np.nan,
            "M_tot_eq_global": float(g["M_tot_eq"].iloc[0]) if "M_tot_eq" in g.columns else np.nan,
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

# === Slim fluorescence dashboard CSV (now includes t80) ===
dash_cols = [
    "Promoter","Rel_strength","k_tx_RNAs_per_s",
    "Protein_mature_eq","Protein_immature_eq",
    "t10_p_m_s","t50_protein_s","t80_p_m_s","t90_p_m_s",
    "maturation_lag_s","fraction_mature_at_1800s","AUC_p_m_0_tmax"
]
fluor_csv = os.path.join(OUTDIR, "fluorescence_dashboard.csv")
summary_df[dash_cols].to_csv(fluor_csv, index=False)
print(f"CSV saved: {os.path.relpath(fluor_csv)}")

# --- Anderson k_tx ladder (ranked; optional shortlist) ---
SHORTLIST = []
ladder = summary_df[["Promoter", "Rel_strength", "k_tx_RNAs_per_s"]].copy()
if SHORTLIST:
    ladder = ladder[ladder["Promoter"].isin(SHORTLIST)]
ladder = ladder.sort_values("k_tx_RNAs_per_s", ascending=False).reset_index(drop=True)
ladder["Rank"] = ladder.index + 1

ladder_csv = os.path.join(OUTDIR, "anderson_k_tx_ladder.csv")
ladder.to_csv(ladder_csv, index=False)
print(f"CSV saved: {os.path.relpath(ladder_csv)}")

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
        ax2.plot(g["Time_s"], g["Protein_molecules"], label="Protein (mature)", color="tab:green", lw=1.8)
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

def plot_overlay_promoters(df, promoters, *, outdir=OUTDIR, pdf=_pdf,
                           save_name=None, show=SHOW_PLOTS,
                           legend_cols=None, legend_fontsize=8,
                           legend_title="Promoters", logy_mrna=False, logy_protein=False):
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
                 label="_nolegend_")                                                       # mature dashed

    ax.set_ylabel("mRNA (molecules/cell)")
    ax2.set_ylabel("Protein (molecules/cell)", labelpad=10)
    ax.set_xlabel("Time (s)")
    if logy_mrna: ax.set_yscale("log")
    if logy_protein: ax2.set_yscale("log")
    ax.grid(True, alpha=0.25)
    ax.set_title(f"Overlay comparison ({len(keep)} promoters)")

    style_handles = [Line2D([0], [0], color="k", lw=2, ls="-"),
                     Line2D([0], [0], color="k", lw=2, ls="--")]
    ax.legend(style_handles, ["mRNA (left y)", "Protein (right y)"], loc="upper left", frameon=True)

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

def plot_all_promoters_overlay(df, promoters=None, show_mrna=True, show_protein=True,
                               top_n=None, logy=False,
                               legend_cols=None, legend_fontsize=8, legend_title="Promoters",
                               *, outdir=OUTDIR, pdf=_pdf, save_name="all_promoters_overlay",
                               show=SHOW_PLOTS):
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

    handles, labels = axes[-1].get_legend_handles_labels()
    uniq, seen = [], set()
    for h, l in zip(handles, labels):
        if l not in seen and l != "_nolegend_":
            uniq.append((h, l)); seen.add(l)
    if uniq:
        handles, labels = zip(*uniq)
    else:
        handles, labels = [], []

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

def plot_gallery_promoters_v2(df, promoters, per_fig=4, *, outdir=OUTDIR, pdf=_pdf, show=SHOW_PLOTS):
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
            ax2.plot(g["Time_s"], g["Protein_molecules"], label="Protein (mature)", color="tab:green", lw=1.8)
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

# === NEW: Fluorescence-first visualizations =====================
def plot_maturation_stacks(df, promoters, per_fig=4, *, outdir=OUTDIR, pdf=_pdf, show=SHOW_PLOTS):
    """
    Stacked panels showing immature (p_i) vs mature (p_m) protein for selected promoters.
    """
    all_avail = set(df["Promoter"].unique())
    proms = [p for p in promoters if p in all_avail]
    if not proms:
        print("No valid promoters for maturation stacks.")
        return

    for i in range(0, len(proms), per_fig):
        chunk = proms[i:i+per_fig]
        fig, axes = plt.subplots(nrows=len(chunk), ncols=1,
                                 figsize=(9, 3.2*len(chunk)), sharex=True)
        if len(chunk) == 1:
            axes = [axes]
        for ax, name in zip(axes, chunk):
            g = df[df["Promoter"] == name].sort_values("Time_s")
            ax.plot(g["Time_s"], g["Protein_immature_molecules"], lw=1.6, label="Protein (immature)")
            ax.plot(g["Time_s"], g["Protein_molecules"], lw=1.8, label="Protein (mature)")
            ax.set_ylabel("Protein (molecules/cell)")
            ax.set_title(f"{name} — maturation")
            ax.grid(True, alpha=0.25)
            ax.legend(loc="upper left")
        axes[-1].set_xlabel("Time (s)")
        plt.tight_layout()
        page = i // per_fig + 1
        _save_figure(fig, f"maturation_stacks_page_{page:02d}", outdir, pdf)
        if show and SHOW_PLOTS:
            plt.show()
        else:
            plt.close(fig)

def plot_maturation_lag_scatter(summary_df, *, outdir=OUTDIR, pdf=_pdf, show=SHOW_PLOTS):
    """
    Scatter of maturation lag (t50_p_m - t50_mRNA) vs steady-state mature brightness.
    """
    sdf = summary_df.copy()
    sdf = sdf[np.isfinite(sdf["maturation_lag_s"]) & np.isfinite(sdf["Protein_mature_eq"])]
    if sdf.empty:
        print("No data for maturation lag scatter.")
        return
    fig, ax = plt.subplots(figsize=(7.5, 5.2))
    ax.scatter(sdf["maturation_lag_s"], sdf["Protein_mature_eq"], s=40, alpha=0.8)
    ax.set_xlabel("Maturation lag (s) = t50_protein - t50_mRNA")
    ax.set_ylabel("Steady-state mature protein (molecules/cell)")
    ax.set_title("Maturation lag vs brightness")
    ax.grid(True, alpha=0.25)
    # Annotate a few most-bright points for readability
    top = sdf.sort_values("Protein_mature_eq", ascending=False).head(8)
    for _, row in top.iterrows():
        ax.annotate(row["Promoter"], (row["maturation_lag_s"], row["Protein_mature_eq"]),
                    xytext=(5,5), textcoords="offset points", fontsize=8)
    _save_figure(fig, "maturation_lag_scatter", outdir, pdf)
    if show and SHOW_PLOTS:
        plt.show()
    else:
        plt.close(fig)

def plot_brightness_and_speed_bars(summary_df, *, outdir=OUTDIR, pdf=_pdf, show=SHOW_PLOTS, top_n=20):
    """
    Side-by-side bars: p_m,eq (brightness) and t80_p_m (speed) for top-N by brightness.
    """
    sdf = summary_df.copy()
    sdf = sdf.sort_values("Protein_mature_eq", ascending=False).head(top_n)
    if sdf.empty:
        print("No data for brightness/speed bars.")
        return
    fig, axes = plt.subplots(1, 2, figsize=(12, 0.4*max(len(sdf),1) + 2.5))
    # Brightness
    axes[0].barh(sdf["Promoter"], sdf["Protein_mature_eq"])
    axes[0].invert_yaxis()
    axes[0].set_xlabel("p_m,eq (molecules/cell)")
    axes[0].set_title("Brightness (mature steady state)")
    # Speed (t80)
    axes[1].barh(sdf["Promoter"], sdf["t80_p_m_s"])
    axes[1].invert_yaxis()
    axes[1].set_xlabel("t80 (s)")
    axes[1].set_title("Time to 80% fluorescence")
    plt.tight_layout()
    _save_figure(fig, "brightness_and_speed_bars", outdir, pdf)
    if show and SHOW_PLOTS:
        plt.show()
    else:
        plt.close(fig)

# === NEW: throttling diagnostics ==============================
def plot_throttling_diagnostics(df, *, outdir=OUTDIR, pdf=_pdf, show=SHOW_PLOTS):
    if "M_tot_all_promoters" not in df.columns or "throttle_phi" not in df.columns:
        return
    first = df["Promoter"].iloc[0]
    g = df[df["Promoter"] == first].sort_values("Time_s")
    fig, ax1 = plt.subplots(figsize=(9, 4.2))
    ax1.plot(g["Time_s"], g["M_tot_all_promoters"], lw=1.8, label="M_tot(t)")
    ax1.set_ylabel("Total mRNA, M_tot(t)")
    ax1.set_xlabel("Time (s)")
    ax2 = ax1.twinx()
    ax2.plot(g["Time_s"], g["throttle_phi"], lw=1.8, ls="--", label="φ(t) = β_eff/β")
    ax2.set_ylabel("Throttle factor φ(t)")
    ax1.grid(True, alpha=0.25)
    ax1.set_title(f"Global ribosome throttling (K_R={K_R:g})")
    _save_figure(fig, "throttling_diagnostics", outdir, pdf)
    if show and SHOW_PLOTS: plt.show()
    else: plt.close(fig)

# === NEW: detection-time vs brightness (amplitude-sensitive) ==================
def _label_offset(i):
    """
    Return a small (dx, dy) in points for label i to reduce overlaps.
    Cycles around a compass with slowly growing radius.
    """
    ring = i // 8 + 1
    step = i % 8
    base = 6 * ring
    compass = [(0,8),(8,0),(0,-8),(-8,0),(6,6),(6,-6),(-6,6),(-6,-6)]
    dx, dy = compass[step]
    return base * (dx/8), base * (dy/8)

def plot_detection_time_vs_brightness(df, summary_df, *, outdir=OUTDIR, pdf=_pdf,
                                      show=SHOW_PLOTS, y_thresh=None, frac_of_max=0.10,
                                      annotate_all=True):
    """
    Scatter of time to reach an absolute mature-protein threshold vs p_m,eq.
    Writes a full CSV (including unreached) and returns the detection DataFrame.
    """
    # Choose a threshold that scales with the library, unless provided
    p_eq_col = "Protein_mature_eq" if "Protein_mature_eq" in summary_df.columns else "p_eq"
    if y_thresh is None:
        p_eq_max = float(summary_df[p_eq_col].max())
        y_thresh = max(1.0, frac_of_max * p_eq_max)  # never below 1 molecule
    y_thresh = float(y_thresh)

    # Helper: first-crossing time for an absolute target
    def _t_cross_abs(t_arr, y_arr, y_target):
        t = np.asarray(t_arr); y = np.asarray(y_arr)
        if np.all(y < y_target): return np.nan
        idx = np.searchsorted(y, y_target)
        if idx == 0: return float(t[0])
        if idx >= len(y): return np.nan
        t0, t1 = t[idx-1], t[idx]
        y0, y1 = y[idx-1], y[idx]
        if y1 == y0: return float(t1)
        return float(t0 + (y_target - y0) * (t1 - t0) / (y1 - y0))

    rows = []
    for name in summary_df["Promoter"]:
        g = df[df["Promoter"] == name].sort_values("Time_s")
        t_detect = _t_cross_abs(g["Time_s"].values, g["Protein_molecules"].values, y_thresh)
        p_eq = float(summary_df.loc[summary_df["Promoter"] == name, p_eq_col].iloc[0])
        rows.append({
            "Promoter": name,
            "t_detect_s": t_detect,
            "reached": np.isfinite(t_detect),
            "threshold_molecules": y_thresh,
            "Protein_mature_eq": p_eq
        })

    det_df = pd.DataFrame(rows)

    # Save full CSV (includes unreached)
    det_csv = os.path.join(OUTDIR, f"detection_time_vs_brightness_thresh_{int(round(y_thresh))}.csv")
    det_df.to_csv(det_csv, index=False)
    print(f"CSV saved: {os.path.relpath(det_csv)}")

    # Plot reached points
    sdf = det_df[det_df["reached"]].copy()
    if sdf.empty:
        print("No points to plot for detection-time vs brightness.")
        return det_df

    fig, ax = plt.subplots(figsize=(7.5, 5.2))
    ax.scatter(sdf["t_detect_s"], sdf["Protein_mature_eq"], s=40, alpha=0.85)
    ax.set_xlabel(f"Time to reach {int(y_thresh)} mature molecules (s)")
    ax.set_ylabel("Steady-state mature protein (molecules/cell)")
    ax.set_title("Detection time vs brightness")
    ax.grid(True, alpha=0.25)

    # Annotate ALL plotted points (offset to reduce collisions)
    if annotate_all:
        for i, row in sdf.reset_index(drop=True).iterrows():
            dx, dy = _label_offset(i)
            ax.annotate(row["Promoter"], (row["t_detect_s"], row["Protein_mature_eq"]),
                        xytext=(dx, dy), textcoords="offset points", fontsize=8,
                        ha="center", va="bottom",
                        arrowprops=dict(arrowstyle='-', lw=0.5, alpha=0.6))

    _save_figure(fig, "detection_time_vs_brightness", outdir, pdf)
    if show and SHOW_PLOTS:
        plt.show()
    else:
        plt.close(fig)

    # Inform if any never reached threshold
    missing = (~det_df["reached"]).sum()
    if missing:
        names = det_df.loc[~det_df["reached"], "Promoter"].tolist()
        print(f"Note: {missing} promoter(s) never reached {int(y_thresh)} molecules within t_max={int(time[-1])} s: {', '.join(names)}.")

    return det_df

# ---------- Paged gallery so ALL selected promoters are shown ----------
def plot_gallery_promoters(df, promoters, per_fig=4, *, outdir=OUTDIR, pdf=_pdf, show=SHOW_PLOTS):
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
            ax2.plot(g["Time_s"], g["Protein_molecules"], label="Protein (mature)", color="tab:green", lw=1.8)
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

# --------- Auto-generate quicklooks ----------
selected_promoters = list(chosen_promoters.keys())
if selected_promoters:
    plot_gallery_promoters(df, selected_promoters, per_fig=4, show=SHOW_PLOTS)

# --- Throttling diagnostics plot (makes impact visible) ---
plot_throttling_diagnostics(df, show=SHOW_PLOTS)

# === NEW: Auto-generate fluorescence-first visuals (top constructs by brightness) ===
try:
    top_bright = summary_df.sort_values("Protein_mature_eq", ascending=False)["Promoter"].head(6).tolist()
    if top_bright:
        plot_maturation_stacks(df, top_bright, per_fig=3, show=SHOW_PLOTS)
        plot_maturation_lag_scatter(summary_df, show=SHOW_PLOTS)
        plot_brightness_and_speed_bars(summary_df, show=SHOW_PLOTS, top_n=min(20, len(summary_df)))
        # Detection-time vs brightness (absolute threshold; set y_thresh for instrument-specific)
        det_df = plot_detection_time_vs_brightness(df, summary_df, show=SHOW_PLOTS, y_thresh=None, frac_of_max=0.10, annotate_all=True)

        # === Merge detection results back into summary & re-save everything ===
        summary_df = summary_df.merge(det_df[["Promoter","t_detect_s","threshold_molecules"]],
                                      on="Promoter", how="left")
        summary_df = summary_df.rename(columns={
            "t_detect_s":"t_detect_s_at_thresh",
            "threshold_molecules":"detect_thresh_molecules"
        })
        # Overwrite summary with new columns
        summary_df.to_csv(sum_csv, index=False)
        print(f"(updated) CSV saved: {os.path.relpath(sum_csv)}")

        # Update dashboard CSV to include detection columns
        dash_cols_extended = dash_cols + ["t_detect_s_at_thresh","detect_thresh_molecules"]
        summary_df[dash_cols_extended].to_csv(fluor_csv, index=False)
        print(f"(updated) CSV saved: {os.path.relpath(fluor_csv)}")

        # Full label table to terminal (all rows)
        with pd.option_context('display.max_rows', None, 'display.max_columns', None, 'display.width', 180):
            label_cols = ["Promoter","Protein_mature_eq","t80_p_m_s","t_detect_s_at_thresh","detect_thresh_molecules"]
            print("\nAll labels & detection times:")
            print(summary_df[label_cols].to_string(index=False, float_format=lambda x: f"{x:.3g}"))
except Exception as e:
    print(f"⚠️ Visualization add-ons skipped due to: {e}")

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
        if ENABLE_RESOURCE_COUPLING:
            # Re-simulate whole set so competition is consistent
            new_set = sorted(current.union(missing))
            df_new, _ = simulate_promoters_coupled({k: promoter_strengths[k] for k in new_set}, K_R=K_R)
            df = df_new
        else:
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
        if ENABLE_RESOURCE_COUPLING:
            new_set = sorted(current.union(missing))
            df_new, _ = simulate_promoters_coupled({k: promoter_strengths[k] for k in new_set}, K_R=K_R)
            df = df_new
        else:
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
