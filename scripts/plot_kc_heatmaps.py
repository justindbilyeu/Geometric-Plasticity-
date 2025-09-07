#!/usr/bin/env python3
"""
plot_kc_heatmaps.py
Generate error heatmaps comparing Kc_RH / Kc_DS vs Engineering Rule.

Inputs:
  data/Kc_comparison_grid.csv

Outputs:
  img/Kc_error_RH.png
  img/Kc_error_DS.png
"""
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

CSV_PATH = os.path.join("data", "Kc_comparison_grid.csv")
OUT_DIR = "img"

def pivot_err(df, err_col, B_value):
    sub = df[df["B"] == B_value].copy()
    # We want Δ on x, A on y; average if duplicates exist
    table = sub.pivot_table(index="A", columns="Δ", values=err_col, aggfunc="mean")
    return table.sort_index(axis=0), table.columns.values, table.index.values

def plot_heatmap(ax, mat, xs, ys, title):
    im = ax.imshow(mat, aspect="auto", origin="lower",
                   extent=[min(xs), max(xs), min(ys), max(ys)])
    ax.set_xlabel("Δ")
    ax.set_ylabel("A")
    ax.set_title(title)
    plt.colorbar(im, ax=ax, label="relative error vs ER")

def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    df = pd.read_csv(CSV_PATH)
    Bs = sorted(df["B"].unique())

    for err_col, fname, nice in [("err_RH", "Kc_error_RH.png", "RH vs ER"),
                                 ("err_DS", "Kc_error_DS.png", "DS vs ER")]:
        fig, axes = plt.subplots(1, len(Bs), figsize=(5*len(Bs), 4), squeeze=False)
        for k, Bv in enumerate(Bs):
            mat, xs, ys = pivot_err(df, err_col, Bv)
            axes[0, k].set_title(f"{nice} (B={Bv})")
            plot_heatmap(axes[0, k], mat.values, mat.columns.values, mat.index.values, axes[0, k].get_title())
        fig.tight_layout()
        out_path = os.path.join(OUT_DIR, fname)
        fig.savefig(out_path, dpi=200, bbox_inches="tight")
        plt.close(fig)
        print(f"Wrote {out_path}")

if __name__ == "__main__":
    main()
