# rwp/experiments/phase_map.py
# ---------------------------------------------
# Geometric Plasticity (preview repo)
# Phase map sweep for ringing vs smooth regimes.
#
# Generates:
#   - results/phase_map/phase_map.csv
#   - results/phase_map/phase_map.png
#
# Method:
#   Uses a lightweight AR(2) surrogate to emulate
#   witness-flux dynamics. For each (alpha, eta),
#   we synthesize an impulse response, compute:
#     • PSD peak gain (dB) excluding DC
#     • Overshoot count (zero-crossing of first difference)
#     • Damping ratio ζ via Hilbert envelope
#   and classify “ringing” if gain ≥ psd_db and
#   overshoots ≥ overs_min.
#
# CLI:
#   python -m rwp.experiments.phase_map \
#     --alphas 0.1,0.2,0.3,0.4,0.5 \
#     --etas 0.02,0.04,0.06,0.08 \
#     --n 512 --psd_db 6 --overs_min 2 --seed 42
#
# Dependencies: numpy, scipy, pandas, matplotlib
# ---------------------------------------------

import argparse
import os
from typing import Tuple, Dict, List

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.signal import lfilter, find_peaks, welch
from scipy.fft import rfftfreq, rfft


def psd_peak(signal: np.ndarray, fs: float = 1.0) -> Dict[str, float]:
    """Return {'peak_freq', 'gain_db'} for non-DC PSD peak.
    Normalizes by std, crops to 99% cumulative energy to reduce tail noise."""
    x = np.asarray(signal, dtype=float)
    x = x - np.mean(x)
    std = np.std(x)
    if std <= 0:
        return {"peak_freq": 0.0, "gain_db": 0.0}
    x /= std

    # Crop to 99% energy
    energy = np.cumsum(x ** 2)
    energy /= energy[-1]
    idx_last = np.argmax(energy >= 0.99) + 1
    if idx_last < 32:
        return {"peak_freq": 0.0, "gain_db": 0.0}
    x = x[:idx_last]

    freqs = rfftfreq(len(x), d=1.0 / fs)
    X = rfft(x)
    Pxx = (np.abs(X) ** 2) / len(x)
    if len(Pxx) <= 2:
        return {"peak_freq": 0.0, "gain_db": 0.0}

    # Exclude DC (index 0)
    Pxx_db = 10.0 * np.log10(Pxx + 1e-12)
    idx_peak = np.argmax(Pxx_db[1:]) + 1
    peak = Pxx_db[idx_peak]
    baseline = np.median(Pxx_db[1:])  # robust to outliers
    return {"peak_freq": float(freqs[idx_peak]), "gain_db": float(peak - baseline)}


def count_overshoots(series: np.ndarray, thresh_std: float = 0.5) -> int:
    """Total significant local extrema (pos+neg) after z-score; threshold in std."""
    x = np.asarray(series, dtype=float)
    mu, sd = np.mean(x), np.std(x)
    if sd == 0:
        return 0
    z = (x - mu) / sd
    peaks_pos, _ = find_peaks(z, height=thresh_std)
    peaks_neg, _ = find_peaks(-z, height=thresh_std)
    return int(len(peaks_pos) + len(peaks_neg))


def estimate_damping_ratio(signal: np.ndarray, peak_freq: float) -> float:
    """Estimate damping ratio ζ using log-envelope slope from analytic signal."""
    if peak_freq <= 0:
        return np.nan
    # Hilbert envelope
    x = np.asarray(signal, dtype=float)
    x = x - np.mean(x)
    # lightweight Hilbert via FFT phase shift
    n = len(x)
    X = np.fft.rfft(x)
    h = np.zeros_like(X, dtype=float)
    if n % 2 == 0:
        h[1:-1] = 2.0
        h[-1] = 1.0
    else:
        h[1:] = 2.0
    analytic = np.fft.irfft(X * (1j * h), n=n) + 1j * np.fft.irfft(X * h, n=n)
    env = np.abs(analytic)
    mask = env > 1e-8
    if np.count_nonzero(mask) < 8:
        return np.nan
    t = np.arange(n, dtype=float)[mask]
    y = np.log(env[mask] + 1e-12)
    # robust linear fit
    A = np.vstack([t, np.ones_like(t)]).T
    slope, _ = np.linalg.lstsq(A, y, rcond=None)[0]
    sigma = -slope
    omega_d = 2 * np.pi * peak_freq
    denom = np.sqrt(sigma ** 2 + omega_d ** 2)
    if denom == 0:
        return np.nan
    zeta = sigma / denom
    if zeta < 0 or zeta > 1:
        return np.nan
    return float(zeta)


def simulate_series(alpha: float, eta: float, n: int, noise: float = 0.0) -> np.ndarray:
    """AR(2) impulse response surrogate parameterized by (alpha, eta).
    Heuristic map: radius r = max(0, 1 - eta), angle θ = α·π (0..π).
    """
    r = max(0.0, 1.0 - eta)
    theta = np.clip(alpha, 0.0, 1.0) * np.pi  # 0..π
    # AR(2): x_t - φ1 x_{t-1} - φ2 x_{t-2} = ε_t
    phi1 = 2.0 * r * np.cos(theta)
    phi2 = -r ** 2
    a = [1.0, -phi1, -phi2]
    b = [1.0]
    imp = np.zeros(n)
    imp[0] = 1.0
    x = lfilter(b, a, imp)
    if noise > 0:
        x = x + noise * np.random.normal(0, 1, size=n)
    return x


def make_phase_map(
    alphas: List[float],
    etas: List[float],
    n: int,
    psd_db_thresh: float,
    overs_min: int,
    seed: int = 42,
) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    rows = []
    for alpha in alphas:
        for eta in etas:
            series = simulate_series(alpha, eta, n, noise=0.02 * rng.standard_normal())
            psd = psd_peak(series, fs=1.0)
            gain_db = psd["gain_db"]
            peak_freq = psd["peak_freq"]
            overs = count_overshoots(series, thresh_std=0.5)
            zeta = estimate_damping_ratio(series, peak_freq)
            is_ringing = int((gain_db >= psd_db_thresh) and (overs >= overs_min))
            rows.append(
                {
                    "alpha": alpha,
                    "eta": eta,
                    "gain_db": gain_db,
                    "peak_freq": peak_freq,
                    "overshoots": overs,
                    "zeta": zeta if np.isfinite(zeta) else "",
                    "is_ringing": is_ringing,
                }
            )
    return pd.DataFrame(rows)


def plot_phase_map(df: pd.DataFrame, out_png: str) -> None:
    # Build grid for imshow
    a_vals = sorted(df["alpha"].unique())
    e_vals = sorted(df["eta"].unique())
    Z = np.zeros((len(a_vals), len(e_vals)))
    for i, a in enumerate(a_vals):
        for j, e in enumerate(e_vals):
            mask = (df["alpha"] == a) & (df["eta"] == e)
            Z[i, j] = df.loc[mask, "is_ringing"].iloc[0] if mask.any() else 0

    plt.figure(figsize=(9, 7))
    im = plt.imshow(
        Z,
        origin="lower",
        aspect="auto",
        extent=[min(e_vals), max(e_vals), min(a_vals), max(a_vals)],
        cmap="RdYlBu_r",
        vmin=0,
        vmax=1,
    )
    plt.xlabel("Learning rate η")
    plt.ylabel("Memory parameter α")
    plt.title("Ringing Boundary Map (1 = ringing, 0 = smooth)")
    cbar = plt.colorbar(im, shrink=0.85)
    cbar.set_label("Ringing (1) vs Smooth (0)")
    # Overlay numeric markers
    for i, a in enumerate(a_vals):
        for j, e in enumerate(e_vals):
            val = int(Z[i, j])
            color = "white" if val > 0.5 else "black"
            plt.text(e, a, f"{val}", ha="center", va="center", fontsize=8, color=color)
    plt.tight_layout()
    plt.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.close()


def parse_list(arg: str, cast=float) -> List[float]:
    return [cast(s.strip()) for s in arg.split(",") if s.strip()]


def main():
    parser = argparse.ArgumentParser(description="Phase map sweep for ringing vs smooth.")
    parser.add_argument(
        "--alphas",
        type=str,
        default="0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9",
        help="Comma-separated α values in [0,1].",
    )
    parser.add_argument(
        "--etas",
        type=str,
        default="0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09",
        help="Comma-separated η values in (0,1).",
    )
    parser.add_argument("--n", type=int, default=512, help="Impulse response length.")
    parser.add_argument("--psd_db", type=float, default=6.0, help="PSD peak threshold (dB).")
    parser.add_argument("--overs_min", type=int, default=2, help="Min overshoots to call ringing.")
    parser.add_argument("--seed", type=int, default=42, help="RNG seed.")
    parser.add_argument(
        "--out_dir", type=str, default="results/phase_map", help="Output directory."
    )
    args = parser.parse_args()

    alphas = parse_list(args.alphas, cast=float)
    etas = parse_list(args.etas, cast=float)
    os.makedirs(args.out_dir, exist_ok=True)

    df = make_phase_map(
        alphas=alphas,
        etas=etas,
        n=args.n,
        psd_db_thresh=args.psd_db,
        overs_min=args.overs_min,
        seed=args.seed,
    )
    csv_path = os.path.join(args.out_dir, "phase_map.csv")
    png_path = os.path.join(args.out_dir, "phase_map.png")
    df.to_csv(csv_path, index=False)
    plot_phase_map(df, png_path)

    # Console summary
    total = len(df)
    ringing = int(df["is_ringing"].sum())
    print(f"Saved CSV -> {csv_path}")
    print(f"Saved PNG -> {png_path}")
    print(f"Ringing detected at {ringing}/{total} grid points.")


if __name__ == "__main__":
    main()
