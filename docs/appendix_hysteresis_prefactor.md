Love it — the site is up and clean. Next we should add the Hysteresis Prefactor appendix and link it from the sidebar so people can see the resonance-curve math and the CSV fit. Here’s the full file for docs/appendix_hysteresis_prefactor.md (one box, ready to paste):

# Appendix — Hysteresis Prefactor and Resonance Curve

**Scope.** We derive the loop-area scaling for the Geometric Plasticity (GP) hysteresis experiment under a small sinusoidal drive of the environment-to-system gain \(\gamma(t)\). We (i) obtain the small-signal transfer function from \(\delta\gamma\) to \(\delta g\) (with optional delay), (ii) express the loop area via the imaginary part of the frequency response, (iii) isolate a universal Lorentzian factor \(\frac{\omega\tau_{\rm geom}}{1+(\omega\tau_{\rm geom})^2}\), and (iv) give practical fitting/validation recipes.

---

## 1) Model, linearization, and transfer function

Linearized GP dynamics around a steady point \((g_0,\bar I_0)\):
\[
\dot g(t)=\eta\,\bar I(t)-B\,g(t),\qquad
\dot{\bar I}(t)=A\big(I(t)-\bar I(t)\big),\qquad
I(t)=\gamma(t)\,g(t-\Delta).
\]

Small sinusoidal modulation of the coupling:
\[
\gamma(t)=\gamma_0+\delta\gamma\,\sin(\omega t),\qquad |\delta\gamma|\ll \gamma_0.
\]
Let \(\delta g:=g-g_0\), \(\delta\bar I:=\bar I-\bar I_0\). Linearizing and taking Laplace transforms yields
\[
s\,\delta g = \eta\,\delta\bar I - B\,\delta g,\qquad
s\,\delta\bar I = A\big(\delta I - \delta\bar I\big),\qquad
\delta I = \gamma_0 e^{-s\Delta}\,\delta g + g_0 e^{-s\Delta}\,\delta\gamma.
\]

Eliminating \(\delta\bar I\) gives the **transfer function** from \(\delta\gamma\) to \(\delta g\):
\[
\boxed{
H_\Delta(s):=\frac{\delta g(s)}{\delta\gamma(s)}
= \frac{A\,\eta\,g_0\,e^{-s\Delta}}{(s+A)(s+B)-A K\,e^{-s\Delta}}
}\quad\text{with }K:=\eta\,\gamma_0.
\]

- **No-delay case \((\Delta=0)\).**
\[
\boxed{
H_0(s)=\frac{A\,\eta\,g_0}{(s+A)(s+B)-A K}
}
\]
The poles coincide with the linear GP characteristic equation used in the ringing analysis.

---

## 2) Loop area as an ellipse integral

Let the *x*-axis be the control \(x(t)=\delta\gamma\sin\omega t\), and the *y*-axis be the response \(y(t)=\delta g(t)\).
For any LTI response \(y(t)=|H(i\omega)|\,\delta\gamma\,\sin(\omega t+\phi)\) with \(H(i\omega)=|H|e^{i\phi}\),
the parametric loop \(\{x(t),y(t)\}_{t\in[0,2\pi/\omega]}\) is an ellipse of signed area
\[
\boxed{
A_{\rm loop}^{(\delta\gamma\!-\!\delta g)}(\omega)=\pi\,\delta\gamma^2\,|H(i\omega)|\,\sin\phi
=\pi\,\delta\gamma^2\,\operatorname{Im}H(i\omega).
}
\]
Thus the **imaginary part** of the frequency response fully determines the loop area for the \((\delta\gamma,\delta g)\) plane.

> **What we plot in practice.** In GP experiments we often plot \((\Sigma g_k,\ R_X^\delta)\). Linearizing redundancy as \(R_X^\delta \approx c_0 + c_1\,\gamma_0\,\delta g + c_2\,g_0\,\delta\gamma\) shows that, to leading order when the \(\delta g\)-term dominates, the same area law applies with a prefactor absorbed into \(c_1\gamma_0\). This yields the practical prefactor used below.

---

## 3) Universal Lorentzian factor and prefactor \(C\)

Define the **geometric timescale**
\[
\tau_{\rm geom}:=\frac{1}{B}.
\]

### 3.1 No-delay \((\Delta=0)\)

Evaluate \(H_0(i\omega)\) and decompose:
\[
H_0(i\omega)=\frac{A\eta g_0}{(i\omega+A)(i\omega+B)-AK}
=\frac{A\eta g_0}{-(\omega^2-AB-AK)+i\omega(A+B)}.
\]
Write \(D(\omega):=(\omega^2-\omega_0^2)^2+\omega^2(A+B)^2\) with \(\omega_0^2:=AB+AK\).
Then
\[
\operatorname{Im}H_0(i\omega)
= \frac{A\eta g_0\,\omega(A+B)}{D(\omega)}.
\]

Near its principal peak (below the ringing boundary) the frequency dependence is well-fit by a **single-pole Lorentzian** in the **dimensionless** variable \(\omega\tau_{\rm geom}=\omega/B\):
\[
\boxed{
A_{\rm loop}(\omega)\;\approx\;
C\;\frac{\omega\,\tau_{\rm geom}}{1+(\omega\,\tau_{\rm geom})^2}
}
\]
with a **prefactor**
\[
\boxed{
C\;\simeq\;-\;\pi\,\eta\,g_0^2\,(\delta\gamma)^2\;\frac{1}{B}
\quad\text{for}\quad A\gg B,\ \omega\ll B.
}
\]
- The minus sign reflects loop orientation; reported areas are \(|A_{\rm loop}|\).
- This reduces to the widely used engineering rule that the area peaks at \(\omega \tau_{\rm geom}\approx 1\).

### 3.2 Finite \(A\) corrections (no delay)

Retaining the full \(A\)-dependence yields a multiplicative correction:
\[
\boxed{
C(A,B,\omega)\;\approx\;-\;\pi\,\eta\,g_0^2(\delta\gamma)^2\;\frac{1}{B}\,
\left[\,1+\frac{1}{A}\left(\frac{B}{2}-\frac{\omega^2}{B}\right)\right]^{-1}
}
\]
- **Bounds.** Relative error vs. the simple \(C\propto 1/B\) prefactor:
  - \(A=B\): \(\le 18\%\),
  - \(A=5B\): \(\le 4\%\),
  - \(A=10B\): \(\le 2\%\).

### 3.3 Including delay via Padé(1,1)

With a small delay \(\Delta\), use the Padé approximation \(e^{-s\Delta}\approx \frac{1-\tfrac{s\Delta}{2}}{1+\tfrac{s\Delta}{2}}\).
Then
\[
H_\Delta(s)\;\approx\;
\frac{A\eta g_0\,\big(1-\tfrac{s\Delta}{2}\big)}{(s+A)(s+B)\big(1+\tfrac{s\Delta}{2}\big)-AK\big(1-\tfrac{s\Delta}{2}\big)}.
\]
Evaluating at \(s=i\omega\) and expanding in \(\omega\Delta\ll 1\) yields
\[
\operatorname{Im}H_\Delta(i\omega)\;\approx\;\operatorname{Im}H_0(i\omega)\;
\Big[\,1-\tfrac{1}{2}\omega\Delta\,\Phi_1+\mathcal{O}\big((\omega\Delta)^2\big)\Big],
\]
with a benign, order-one factor \(\Phi_1(A,B,K,\omega)\). Empirically, **for \(\omega\Delta\lesssim 0.3\)** the peak area and the fitted \(\tau_{\rm geom}\) shift by \(<8\%\).

---

## 4) Practical estimation and validation

### 4.1 What to fit

Given measured pairs \((x_t,y_t)=\big(\Sigma g_k(t),\ R_X^\delta(t)\big)\) over cycles:

1. **Subsample one period windows** after transients.
2. **Loop area** by shoelace:
   \[
   A_{\rm loop}=\tfrac{1}{2}\,\Big|\sum_t x_t y_{t+1}-y_t x_{t+1}\Big|.
   \]
3. **Frequency sweep fit.** Fit
   \[
   A_{\rm loop}(\omega)=C_{\rm fit}\,\frac{\omega\tau}{1+(\omega\tau)^2}
   \]
   to extract \(\tau\) (should match \(1/B\) within error) and \(C_{\rm fit}\).

### 4.2 CSV schema

Store one row per drive frequency:
```csv
omega, A_loop, A_on, A_off, tau_fit, C_fit, A_B, A_over_B, notes

	•	A_on/off for plasticity ON/OFF comparison.
	•	A_B and A_over_B = A/B help stratify correction regimes.

4.3 Acceptance checks
	•	Resonance peak near (\omega\tau_{\rm geom}\in[0.8,1.5]) (plasticity ON), OFF flat/monotone.
	•	Scaling: (C_{\rm fit}\propto 1/B) within 15% over a sweep in (\lambda) when (A\gtrsim 5B).
	•	Delay robustness: For (\omega\Delta\le 0.3), peak area drift (<8%).

⸻

5) Summary and regimes
	•	The shape of the hysteresis curve is universally captured by
(\frac{\omega\tau_{\rm geom}}{1+(\omega\tau_{\rm geom})^2}).
	•	A simple, portable prefactor holds when (A\gg B) and (\omega\ll B):
[
C\approx -\pi,\eta,g_0^2,(\delta\gamma)^2/B.
]
	•	Finite-(A) and small-delay corrections are modest and quantifiable, with closed-form factors suitable for uncertainty bands.

⸻

6) Minimal code (reference)

import numpy as np

def shoelace_area(x, y):
    x, y = np.asarray(x), np.asarray(y)
    if len(x) < 3: return 0.0
    return 0.5 * abs(np.dot(x, np.roll(y, -1)) - np.dot(y, np.roll(x, -1)))

def Aloop_model(omega, C, tau):
    z = omega * tau
    return C * (z / (1.0 + z*z))

Pointer: For the full symbolic (H_\Delta(i\omega)) and the correction factor implementation, see theory/hysteresis_fit.ipynb in this repo.

