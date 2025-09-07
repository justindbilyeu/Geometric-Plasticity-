# Preview Methods (Data-only)

## Model capsule (high-level)
Closed-loop geometric plasticity with:
- EMA witness: \(\bar I_k(t) = (1-\alpha)\bar I_k(t-1) + \alpha I_k(t)\)
- Geometry update: \(g \leftarrow g + \eta[\bar I - \lambda g - \beta L g]\) with budget projection \(\sum c_k g_k^2 \le B\)
- Optional delay Δ (Padé(1,1) in linear analysis)
- Metrics: witness flux \(\Phi_{\text{wit}}(t)=\Delta\sum_k I_k\), redundancy \(R_X^\delta(t)\), \(\tau_{\text{geom}}\) (autocorr 1/e)

## Operational criteria
- **Ringing detection** (either time series qualifies):  
  (i) PSD peak ≥ **6 dB** above baseline (exclude DC), and  
  (ii) ≥ **2** overshoots in \(R_X^{0.1}(t)\).
- **Learning vs. mere ringing** (under periodic drive): per-period loss proxy \(\Delta\mathcal{L}<0\) and redundancy plateaus.
- **Hysteresis loop**: \(\mathcal{A}_{\text{loop}}(T) \propto \frac{\omega\tau}{1+(\omega\tau)^2}\), \(\omega=2\pi/T\). Peak expected near \(T\approx 2\pi\tau_{\text{geom}}\).

## Threshold formulas compared
- **Routh–Hurwitz (RH)** boxed cubic (Padé(1,1)): yields \(K_c^{\text{RH}}(A,B,\Delta)\).  
- **DeepSeek (DS)** variant with \(A(2/\Delta - B)\) denominator: \(K_c^{\text{DS}}\).  
- **Engineering rule (ER)**: solve \(\arctan(\omega_c/A)+\arctan(\omega_c/B)+\omega_c\Delta=\tfrac{3\pi}{4}\), then map to \(K_c^{\text{ER}}\).

## Acceptance checks (preview)
- `data/Kc_comparison_grid.csv`: \(|K_c^{\text{RH}}-K_c^{\text{ER}}|/K_c^{\text{ER}} \le 2\%\) for A≥5B; similar for DS.  
- `img/hysteresis_area_vs_period.png`: ON peak in T/(2πτ) ∈ [0.8, 1.5]; OFF ~flat.

## Provenance
- RH box & asymptotics (our derivation).  
- DS threshold & motif criteria (DeepSeek notes).  
- ER mapping & error heatmaps (Wolfram notebook).

This file is intentionally high-level; implementation details reside in private working repos until manuscript submission.
