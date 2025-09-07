# rwp/experiments/hysteresis.py
# -------------------------------------------------------
# Geometric Plasticity (preview repo)
# Hysteresis resonance: Plasticity ON vs OFF
#
# Generates:
#   - results/hysteresis/area_vs_period.csv
#   - results/hysteresis/area_vs_period.png
#   - results/hysteresis/loops_T=<mult>.npz   (Σg_k, R^δ traces for ON/OFF)
#
# Method:
#   1) Estimate τ_geom from a reference response (plasticity OFF).
#   2) Sweep drive periods T_drive = T_mult · (2π τ_geom).
#   3) Simulate a surrogate central-spin/bath model:
#        • Base AR(2) response for the central mode.
#        • N bath fragments with couplings g_k(t).
#        • Plasticity ON: g_k strengthened in proportion to
#          instantaneous correlation with central response (Hebbian-ish).
#   4) Build loops in the plane (Σ g_k,  R_X^δ ) and compute area
#      via shoelace formula. Compare ON vs OFF.
#
# CLI:
#   python -m rwp.experiments.hysteresis \
#     --T 150 --fs 100.0 --N 24 --delta 0.1 \
#     --amplitude 0.12 --seed 42 \
#     --periods 0.3,0.5,0.8,1.0,1.3,2.0,3.0
#
# Dependencies: numpy, scipy, pandas, matplotlib
# -------------------------------------------------------

import argparse
import os
from typing import Tuple, Dict

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.signal import lfilter, correlate, butter, filtfilt, find_peaks


# ------------------------------
# Utilities
# ------------------------------

def autocorr_tau_1e(x: np.ndarray, fs: float = 1.0) -> float:
    """1/e crossing of normalized autocorrelation → τ_geom (seconds)."""
    x = np.asarray(x, float)
    x = x - np.mean(x)
    if np.allclose(x, 0):
        return 1.0
    ac = correlate(x, x, mode="full")
    ac = ac[len(ac)//2:]
    ac = ac / (ac[0] + 1e-12)
    idx = np.argmax(ac < 1.0 / np.e)
    if idx == 0:
        idx = len(ac) - 1
    return idx / fs


def shoelace_area(x: np.ndarray, y: np.ndarray) -> float:
    """Area enclosed by parametric loop (x, y) via shoelace formula.
    If curve is open, this computes polygon area of the polyline closure.
    """
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    if len(x) < 3 or len(y) < 3:
        return 0.0
    return 0.5 * abs(np.dot(x, np.roll(y, -1)) - np.dot(y, np.roll(x, -1)))


def lowpass(x: np.ndarray, fs: float, fc: float) -> np.ndarray:
    """Zero-phase low-pass filter for de-noising traces."""
    if fc <= 0 or fc >= fs/2:
        return x
    b, a = butter(3, fc / (fs / 2.0), btype="low")
    return filtfilt(b, a, x)


# ------------------------------
# Surrogate system
# ------------------------------

def ar2_impulse_response(alpha: float, eta: float, n: int) -> np.ndarray:
    """AR(2) impulse response surrogate parameterized by (alpha, eta)."""
    r = max(0.0, 1.0 - eta)
    theta = np.clip(alpha, 0.0, 1.0) * np.pi  # 0..π
    phi1 = 2.0 * r * np.cos(theta)
    phi2 = -r ** 2
    a = [1.0, -phi1, -phi2]
    b = [1.0]
    imp = np.zeros(n)
    imp[0] = 1.0
    return lfilter(b, a, imp)


def simulate_hysteresis(
    T: float,
    fs: float,
    T_drive: float,
    N: int,
    delta: float,
    amplitude: float,
    plasticity_on: bool,
    seed: int = 42,
    alpha: float = 0.4,
    eta: float = 0.06,
) -> Tuple[np.ndarray, np.ndarray]:
    """Return (Σ g_k(t), R^δ(t)) for a driven bath with/without plasticity.
    - Central response: AR(2) impulse filtered by sinusoidal drive.
    - Bath: N fragments with base couplings g_k(t).
    - Plasticity (ON): Hebbian-like increment ∝ correlation with central response.
    """
    rng = np.random.default_rng(seed)
    t = np.arange(0.0, T, 1.0 / fs)
    n = len(t)
    # Driving signal
    omega = 2 * np.pi / max(T_drive, 1e-6)
    drive = amplitude * np.sin(omega * t)

    # Central AR(2) kernel and response to drive (via convolution)
    kern = ar2_impulse_response(alpha=alpha, eta=eta, n=min(n, 1024))
    x = np.convolve(drive, kern, mode="same")
    x = lowpass(x, fs, fc=min(0.25 * fs, fs / 4.0))

    # Initialize bath couplings
    gk = rng.normal(0.08, 0.02, size=(N,))  # initial coupling strengths
    g_traj = np.zeros((n, N))
    R_traj = np.zeros(n)

    # Online correlation state for each fragment (EMA)
    corr_state = np.zeros(N)
    corr_alpha = 2.0 / (1.0 + fs * 0.2)  # ~0.2 s horizon

    # Fragment carrier signals (slightly detuned sines)
    detunes = rng.normal(0.0, 0.05 * omega, size=N)
    phases = rng.uniform(0, 2 * np.pi, size=N)
    frag = np.sin((omega + detunes)[:, None] * t[None, :] + phases[:, None])

    # Evolution
    for i, ti in enumerate(t):
        # Instantaneous MI-proxy per fragment: |corr(x, frag_k)| (EMA update)
        x_i = x[i]
        fk = frag[:, i]
        corr_state = (1.0 - corr_alpha) * corr_state + corr_alpha * (x_i * fk)

        if plasticity_on:
            # Hebbian-like plasticity with budget projection
            gk += 0.04 * corr_state - 0.01 * gk  # learn - decay
            # Nonnegativity & soft budget
            gk = np.clip(gk, 0.0, None)
            budget = N * 0.15
            scale = min(1.0, budget / (np.sum(gk) + 1e-12))
            gk *= scale

        g_traj[i, :] = gk

        # Redundancy R^δ: count fragments whose running corr magnitude exceeds δ
        # Use normalized corr proxy (bounded by ~1 for sines)
        R_traj[i] = np.sum(np.abs(corr_state) > delta)

    # Outputs
    sum_g = np.sum(g_traj, axis=1)
    R = R_traj
    return sum_g, R


# ------------------------------
# Main experiment
# ------------------------------

def main():
    parser = argparse.ArgumentParser(description="Hysteresis resonance: plasticity ON vs OFF.")
    parser.add_argument("--T", type=float, default=150.0, help="Total simulation time (s).")
    parser.add_argument("--fs", type=float, default=100.0, help="Sampling rate (Hz).")
    parser.add_argument("--N", type=int, default=24, help="Number of bath fragments.")
    parser.add_argument("--delta", type=float, default=0.1, help="Redundancy threshold δ.")
    parser.add_argument("--amplitude", type=float, default=0.12, help="Drive amplitude.")
    parser.add_argument("--seed", type=int, default=42, help="RNG seed.")
    parser.add_argument(
        "--periods",
        type=str,
        default="0.3,0.5,0.8,1.0,1.3,2.0,3.0",
        help="Comma-separated T_mult values.",
    )
    parser.add_argument("--alpha", type=float, default=0.4, help="AR(2) angle parameter α.")
    parser.add_argument("--eta", type=float, default=0.06, help="AR(2) radius parameter via 1-η.")
    parser.add_argument("--out_dir", type=str, default="results/hysteresis", help="Output dir.")
    args = parser.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)

    # 1) Reference run (OFF) to estimate τ_geom
    #    Use a long drive period to avoid entrainment.
    sum_g_ref, _ = simulate_hysteresis(
        T=args.T, fs=args.fs, T_drive=max(args.T * 10.0, 1.0), N=args.N,
        delta=args.delta, amplitude=args.amplitude, plasticity_on=False,
        seed=args.seed, alpha=args.alpha, eta=args.eta
    )
    tau_geom = autocorr_tau_1e(sum_g_ref, fs=args.fs)
    print(f"[hysteresis] Estimated τ_geom ≈ {tau_geom:.3f} s")

    # 2) Sweep T_drive = T_mult · 2π τ_geom
    T_mults = [float(s.strip()) for s in args.periods.split(",") if s.strip()]
    rows = []

    for Tm in T_mults:
        T_drive = Tm * (2.0 * np.pi * tau_geom)
        print(f"[hysteresis] T_mult={Tm:.2f}  ⇒  T_drive={T_drive:.3f} s")

        # ON
        sum_g_on, R_on = simulate_hysteresis(
            T=args.T, fs=args.fs, T_drive=T_drive, N=args.N,
            delta=args.delta, amplitude=args.amplitude, plasticity_on=True,
            seed=args.seed, alpha=args.alpha, eta=args.eta
        )
        area_on = shoelace_area(sum_g_on, R_on)

        # OFF
        sum_g_off, R_off = simulate_hysteresis(
            T=args.T, fs=args.fs, T_drive=T_drive, N=args.N,
            delta=args.delta, amplitude=args.amplitude, plasticity_on=False,
            seed=args.seed, alpha=args.alpha, eta=args.eta
        )
        area_off = shoelace_area(sum_g_off, R_off)

        # Save NPZ for this period
        np.savez(
            os.path.join(args.out_dir, f"loops_T={Tm:.2f}.npz"),
            x_on=sum_g_on, y_on=R_on, x_off=sum_g_off, y_off=R_off,
        )

        rows.append(
            {
                "T_mult": Tm,
                "T_drive": T_drive,
                "area_on": area_on,
                "area_off": area_off,
                "ratio_on_off": area_on / (area_off + 1e-12),
                "tau_geom": tau_geom,
            }
        )

        print(f"  → areas: ON={area_on:.4f}, OFF={area_off:.4f}, ratio={rows[-1]['ratio_on_off']:.2f}")

    # 3) Save CSV
    df = pd.DataFrame(rows).sort_values("T_mult")
    csv_path = os.path.join(args.out_dir, "area_vs_period.csv")
    df.to_csv(csv_path, index=False)
    print(f"[hysteresis] Saved CSV → {csv_path}")

    # 4) Plot
    plt.figure(figsize=(8.2, 6.2))
    plt.plot(df["T_mult"], df["area_on"], marker="o", label="Plasticity ON")
    plt.plot(df["T_mult"], df["area_off"], marker="s", label="Plasticity OFF")
    plt.axvspan(0.8, 1.5, color="gray", alpha=0.20, label="Expected peak window")
    plt.xlabel(r"$T_{\mathrm{mult}}$  (drive period / $2\pi\tau_{\mathrm{geom}}$)")
    plt.ylabel("Hysteresis loop area  (Σg_k vs R$^{\delta}$)")
    plt.title("Hysteresis Resonance: Plasticity ON vs OFF")
    plt.grid(alpha=0.3)
    plt.legend()
    plt.tight_layout()
    png_path = os.path.join(args.out_dir, "area_vs_period.png")
    plt.savefig(png_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"[hysteresis] Saved PNG → {png_path}")

    # Console acceptance hint
    # Check whether ON peak lies near T_mult ≈ 1
    idx_max = int(np.argmax(df["area_on"].values))
    Tm_peak = df["T_mult"].values[idx_max]
    print(f"[hysteresis] ON peak at T_mult ≈ {Tm_peak:.2f} "
          f"(expect ~1.0 and within 0.8–1.5 window).")


if __name__ == "__main__":
    main()
