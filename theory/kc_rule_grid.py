# theory/kc_rule_grid.py
# Geometric Plasticity — Kc comparison grid + heatmaps
# Generates: results/kc_rule/Kc_comparison_grid.csv
#            results/kc_rule/Kc_error_RH_B{...}.png
#            results/kc_rule/Kc_error_DS_B{...}.png
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from math import atan
from scipy.optimize import brentq

# -----------------------------
# Threshold formulas
# -----------------------------

def kc_rh(A, B, Delta):
    """
    Routh–Hurwitz 'boxed' form (Padé(1,1) cubic):
    Kc_RH = [ AB + (2/Δ)(1 + (A+B)Δ/2)( (A+B) + ABΔ/2 ) ] / [ A(1 + Δ/2) ]
    """
    if Delta <= 0:
        # Δ -> 0 limit from RH box
        return A*B / A  # -> B ? but small-Δ asymptotic reduces to AB/A = B; not used in our grids
    num = (A*B) + (2.0/Delta)*(1.0 + (A+B)*Delta/2.0)*((A+B) + (A*B)*Delta/2.0)
    den = A*(1.0 + Delta/2.0)
    return num/den

def kc_ds(A, B, Delta):
    """
    DeepSeek’s denominator form:
    Kc_DS = [ (2/Δ + A + B)( 2(A+B)/Δ + AB ) ] / [ A(2/Δ - B) ]
    (Algebraically close to RH; compare numerically.)
    """
    if Delta <= 0:
        return np.inf  # undefined at Δ=0 in this parameterization
    num = (2.0/Delta + A + B) * (2.0*(A+B)/Delta + A*B)
    den = A*(2.0/Delta - B)
    return num/den

def omega_c_er(A, B, Delta):
    """
    Engineering-rule root: solve arctan(w/A)+arctan(w/B)+w*Δ = 3π/4.
    Monotone LHS in w for Δ>=0, so a single root exists.
    """
    target = 3.0*np.pi/4.0

    def f(w):
        return np.arctan(w/A) + np.arctan(w/B) + w*Delta - target

    # bracket: start near 0, extend to a safe upper bound
    lo = 0.0
    # upper bound grows as 1/Δ; include A,B scales too
    scale = max(A, B, 1.0)
    hi = 100.0*scale + (100.0/Delta if Delta > 0 else 100.0)
    # Ensure sign change
    if f(hi) < 0:
        # increase further if needed
        hi = hi*10.0
    return brentq(f, 1e-9, hi, maxiter=200)

def kc_er(A, B, Delta):
    """
    Engineering-rule Kc from ωc:
    Kc_ER ≈ (1/A) * sqrt( (A^2 + ω^2)(B^2 + ω^2) )
    """
    w = omega_c_er(A, B, Delta)
    return (1.0/A)*np.sqrt((A*A + w*w)*(B*B + w*w))

# -----------------------------
# Grid + CSV + heatmaps
# -----------------------------

def main():
    A_list = [1.0, 2.0, 5.0, 10.0]
    B_list = [0.5, 1.0, 2.0]
    D_list = [0.05, 0.10, 0.20]

    rows = []
    for A in A_list:
        for B in B_list:
            for D in D_list:
                try:
                    rh = kc_rh(A, B, D)
                except Exception:
                    rh = np.nan
                try:
                    ds = kc_ds(A, B, D)
                except Exception:
                    ds = np.nan
                try:
                    er = kc_er(A, B, D)
                except Exception:
                    er = np.nan

                err_rh = abs((rh - er)/er) if np.isfinite(er) and er != 0 else np.nan
                err_ds = abs((ds - er)/er) if np.isfinite(er) and er != 0 else np.nan

                rows.append({
                    "A": A, "B": B, "Δ": D,
                    "Kc_RH": rh, "Kc_DS": ds, "Kc_ER": er,
                    "err_RH": err_rh, "err_DS": err_ds
                })

    df = pd.DataFrame(rows)
    out_dir = "results/kc_rule"
    os.makedirs(out_dir, exist_ok=True)
    csv_path = os.path.join(out_dir, "Kc_comparison_grid.csv")
    df.to_csv(csv_path, index=False)

    # Heatmaps by B (Δ on x, A on y)
    for B in B_list:
        sub = df[df["B"] == B].copy()
        # pivot to 2D arrays
        A_vals = sorted(sub["A"].unique())
        D_vals = sorted(sub["Δ"].unique())
        Z_rh = np.zeros((len(A_vals), len(D_vals)))
        Z_ds = np.zeros((len(A_vals), len(D_vals)))
        for i, A in enumerate(A_vals):
            for j, D in enumerate(D_vals):
                cell = sub[(sub["A"]==A) & (sub["Δ"]==D)]
                Z_rh[i, j] = float(cell["err_RH"]) if not cell.empty else np.nan
                Z_ds[i, j] = float(cell["err_DS"]) if not cell.empty else np.nan

        # RH error heatmap
        plt.figure(figsize=(6,4.5))
        im = plt.imshow(Z_rh, origin="lower", aspect="auto",
                        extent=[min(D_vals), max(D_vals), min(A_vals), max(A_vals)])
        plt.xlabel("Δ")
        plt.ylabel("A")
        plt.title(f"Rel. error (RH vs ER), B={B}")
        cbar = plt.colorbar(im)
        cbar.set_label("abs(Kc_RH - Kc_ER) / Kc_ER")
        plt.tight_layout()
        plt.savefig(os.path.join(out_dir, f"Kc_error_RH_B{B}.png"), dpi=220)
        plt.close()

        # DS error heatmap
        plt.figure(figsize=(6,4.5))
        im = plt.imshow(Z_ds, origin="lower", aspect="auto",
                        extent=[min(D_vals), max(D_vals), min(A_vals), max(A_vals)])
        plt.xlabel("Δ")
        plt.ylabel("A")
        plt.title(f"Rel. error (DS vs ER), B={B}")
        cbar = plt.colorbar(im)
        cbar.set_label("abs(Kc_DS - Kc_ER) / Kc_ER")
        plt.tight_layout()
        plt.savefig(os.path.join(out_dir, f"Kc_error_DS_B{B}.png"), dpi=220)
        plt.close()

    print(f"[ok] wrote {csv_path} and error heatmaps in {out_dir}/")

if __name__ == "__main__":
    main()
