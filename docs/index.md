---
title: Geometric Plasticity (preview)
description: Witness-driven adaptation, ringing thresholds, hysteresis resonance, and motif selection — a minimal, reproducible bundle.
---

# Geometric Plasticity (preview)

**What it is.** A small, testable core: witness-driven plasticity that (i) rings past a gain threshold, (ii) shows a hysteresis-area resonance at its geometric timescale, and (iii) selects broadcast vs modular motifs under cost geometry.

**What you get.** Reproducible code, diagnostics, and the three plots that matter.

---

## Quick results

- **Phase map (ringing boundary)**  
  <img src="img/phase_map.png" alt="Ringing phase map" width="520" />
  <br/>CSV: [`results/phase_map/phase_map.csv`](results/phase_map/phase_map.csv)

- **Hysteresis resonance (ON vs OFF)**  
  <img src="img/hysteresis_area.png" alt="Hysteresis area vs period" width="520" />
  <br/>CSV: [`results/hysteresis/area_vs_period.csv`](results/hysteresis/area_vs_period.csv)

- **Ringing threshold comparison (Kc)**  
  CSV: [`data/Kc_comparison_grid.csv`](data/Kc_comparison_grid.csv)  
  Heatmaps: [`img/Kc_error_RH.png`](img/Kc_error_RH.png), [`img/Kc_error_DS.png`](img/Kc_error_DS.png)

> If an image 404s, it just means you haven’t pushed that artifact yet — see **How to reproduce** below.

---

## How to reproduce (≈3–10 min, laptop)

```bash
# create & activate venv (Python 3.10+)
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt

# phase map + hysteresis (saves CSV/PNGs under results/)
python rwp/experiments/phase_map.py
python rwp/experiments/hysteresis.py

# Kc comparison CSV + error heatmaps (Mathematica / Wolfram Language)
# open and run: theory/Kc_notebook.nb
# -> exports data/Kc_comparison_grid.csv and img/Kc_error_*.png

Acceptance checks
	•	Ringing detection: PSD peak ≥ 6 dB and ≥ 2 overshoots on redundancy series.
	•	Hysteresis: Plasticity ON peaks near (T \approx 2\pi,\tau_{\text{geom}}) and exceeds OFF.
	•	Kc agreement: RH and DS forms match the engineering-rule threshold within ~1–2% for (A \gtrsim 5B).

⸻

What’s in here
	•	Theory snapshots
docs/appendix_ring_threshold.md — Padé(1,1) + Routh–Hurwitz, delay-aware (K_c).
docs/appendix_hysteresis_prefactor.md — Exact prefactor (C) and large-(A) corrections.
docs/appendix_motif_universality.md — Broadcast↔Modular threshold vs. (\beta/\lambda), (\mu_2).
	•	Code entry points
rwp/diagnostics.py — PSD peak, overshoots, damping-ratio estimator.
rwp/experiments/phase_map.py — ((\alpha,\eta)) sweep → phase map PNG/CSV.
rwp/experiments/hysteresis.py — ON vs OFF loops, shoelace area vs period.
	•	Data & artifacts
results/phase_map/phase_map.csv, results/hysteresis/area_vs_period.csv,
data/Kc_comparison_grid.csv, img/*.png.

⸻

Minimal model (for orientation)

Linearized GP loop:
[
\dot g = \eta,\bar I - B,g,\qquad
\dot{\bar I} = A,(I-\bar I),\qquad
I(t)=\gamma,g(t-\Delta),\quad K=\eta\gamma.
]

Delay via Padé(1,1) and the cubic characteristic polynomial give the delay-aware threshold (K_c^{\rm RH}).
Frequency-rule cross-check:
[
\arctan(\omega_c/A)+\arctan(\omega_c/B)+\omega_c\Delta=\tfrac{3\pi}{4}.
]

Hysteresis loop area follows a single-pole form:
[
A_{\text{loop}}(\omega)\propto\frac{\omega ,\tau_{\text{geom}}}{1+(\omega ,\tau_{\text{geom}})^2},
]
with prefactor (C(A,B,\eta,g_0,\delta\gamma,\omega)) detailed in the appendix.

⸻

Notes & caveats
	•	This is a preview: APIs and plots may change; results are reproducible on a laptop.
	•	The biological/phenomenological links are intentionally out of scope here; this site is about the falsifiable core.

⸻

Contact

Questions, bugs, or collaboration requests?
Open an issue in the repo or DM @justindbilyeu.

<br/>


© 2025 Resonance Geometry Collective — preview; content may change.

