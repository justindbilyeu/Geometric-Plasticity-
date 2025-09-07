---
title: Geometric Plasticity — Preview
---

# Geometric Plasticity — Preview

This page shares **comparison data and reproducible plots** for the **ringing threshold** in Geometric Plasticity (GP), contrasting:
- a **Routh–Hurwitz (RH)** cubic condition with Padé(1,1) delay,
- the **DeepSeek (DS)** closed-form,
- an **Engineering Rule (ER)** solved from the phase condition  
  \[
  \arctan\!\frac{\omega_c}{A} + \arctan\!\frac{\omega_c}{B} + \omega_c \Delta = \frac{3\pi}{4}.
  \]

Jump to: **[Methods](preview_methods.md)** · **[CSV](../data/Kc_comparison_grid.csv)** · **[Plot script](../scripts/plot_kc_heatmaps.py)**

---

## Data & Plots

- **CSV:** [`data/Kc_comparison_grid.csv`](../data/Kc_comparison_grid.csv)  
  Columns: `A,B,Δ,Kc_RH,Kc_DS,Kc_ER,err_RH,err_DS`.

- **Heatmaps (generate locally):**
  1. `python scripts/plot_kc_heatmaps.py`
  2. Images saved to `img/`:
     - `Kc_error_RH.png`
     - `Kc_error_DS.png`

> Observations (from the grid):  
> • RH and DS track ER within ~1–2% across the shown ranges.  
> • Errors shrink as \(A/B\) grows (fast EMA regime).  
> • RH and DS are numerically near-identical ⇒ algebraically equivalent forms.

---

## Reproduce

```bash
# From repo root
python -m pip install matplotlib pandas numpy
python scripts/plot_kc_heatmaps.py
open img/Kc_error_RH.png
