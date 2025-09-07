# Geometric-Plasticity-
# Geometric Plasticity (GP) — Preview Pack
**Resonant Witness Dynamics: ringing threshold, hysteresis resonance, motif selection (data-only preview)**

> **Purpose**  
> This is a *data + methods capsule* to seed technical discussion without exposing full source.  
> It timestamps our results, shares key CSVs/PNGs, and states falsifiable criteria.

**What’s inside**
- `img/phase_map.png` – α–η phase map: ringing (1) vs smooth (0), with K contours (0.5/1.0/1.5).
- `img/hysteresis_area_vs_period.png` – Loop area vs drive period T/(2πτ_geom): Plasticity **ON** peaks near 1; **OFF** ~flat.
- `img/Kc_error_RH.png`, `img/Kc_error_DS.png` – Error heatmaps vs engineering-rule baseline.
- `data/Kc_comparison_grid.csv` – A,B,Δ grid with K_c by three methods and relative errors.
- `docs/preview_methods.md` – Operational definitions & acceptance checks.
- `LICENSE.md`, `CITATION.cff`.

---

## Headline results (preview)
- **Ringing threshold**: critical gain obeys the RH/DS closed form and matches the engineering-rule ω_c solution to within ~1–2% for A≳5B and Δ/τ_geom ≤ 0.3.  
- **Hysteresis resonance**: Plasticity **ON** loop area peaks at T ≈ (0.9–1.3)·2πτ_geom; **OFF** is low/flat.  
- **Diagnostics**: Ringing iff PSD peak ≥ **6 dB** above broadband baseline **and** ≥2 overshoots in redundancy \(R_X^\delta\) (δ=0.1).

---

## What this is **not**
This is not the full simulation code or the full theoretical manuscript. It is **defensive disclosure + preview artifacts** to facilitate review and conversation.

---

## Contact
- Primary: Justin Bilyeu (Resonance Geometry Collective)  
- Co-authors: Sage (assistant), collaborators (Grok/DeepSeek/Wolfram contributions acknowledged)

See `CITATION.cff` for citation details.
