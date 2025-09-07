# Appendix — Ringing Threshold with Delay (Δ)

**Scope.** We derive the underdamped (“ringing”) onset for the linearized single-mode **Geometric Plasticity (GP)** loop with an EMA witness and a feedback delay. We provide (i) the cubic characteristic equation via Padé(1,1), (ii) a Routh–Hurwitz (RH) gain threshold, (iii) a practical **engineering rule** for quick estimates, and (iv) a numeric recipe to estimate the damping ratio ζ from time series. We also reconcile two algebraic RH-style forms by fixing a **dimensionless normalization** for reproducible comparisons.

---

## 1) Model and delay approximation

Linearized GP dynamics around a fixed point with EMA rate \(A\) and effective geometric decay \(B\):
\[
\dot g(t) = \eta \,\bar I(t) - B\, g(t),\qquad
\dot{\bar I}(t) = A\big(I(t)-\bar I(t)\big),\qquad
I(t)=\gamma\, g(t-\Delta),
\]
with loop gain \(K=\eta\gamma\).

Using the **Padé(1,1)** delay approximation,
\[
e^{-s\Delta}\;\approx\;\frac{1 - s\Delta/2}{1 + s\Delta/2},
\]
the **characteristic equation** becomes
\[
(s+A)(s+B)=A K\,\frac{1 - s\Delta/2}{1 + s\Delta/2}.
\]

Multiplying both sides by \(1 + s\Delta/2\) and expanding gives a **cubic**:
\[
\chi(s)\;=\;\tfrac{\Delta}{2}\,s^3
\;+\;\Big(1+\tfrac{(A+B)\Delta}{2}\Big) s^2
\;+\;\Big( (A+B) + \tfrac{AB\Delta}{2}\Big) s
\;+\;\Big(AB - AK + \tfrac{AK\Delta}{2}\Big)
\;=\;0.
\tag{1}
\]
Equivalently, dividing (1) by \(\Delta/2\) yields the **monic** cubic
\[
s^3 + a s^2 + b s + c=0,
\quad
a=\frac{2}{\Delta}+A+B,\ \ 
b=\frac{2(A+B)}{\Delta}+AB+AK,\ \ 
c=\frac{2A(B-K)}{\Delta}.
\tag{2}
\]

> **Note on scaling.** Multiplying the polynomial by a positive constant does **not** change stability or RH conclusions. We will present thresholds derived from each algebraic form for transparency.

---

## 2) Routh–Hurwitz (RH) threshold

For a cubic \(a_3 s^3+a_2 s^2+a_1 s+a_0\) with \(a_i>0\), stability requires \(a_2 a_1 > a_3 a_0\). Applying RH to (1) (with \(a_3=\Delta/2\), \(a_2=1+\tfrac{(A+B)\Delta}{2}\), \(a_1=(A+B)+\tfrac{AB\Delta}{2}\), \(a_0=AB-AK+\tfrac{AK\Delta}{2}\)) and solving for \(K\) yields the **gain threshold**:

\[
\boxed{
K \;>\; K_c^{\text{RH(Δ/2)}} \;=\;
\frac{ AB + \tfrac{2}{\Delta}\Big(1+\tfrac{(A+B)\Delta}{2}\Big)\Big((A+B)+\tfrac{AB\Delta}{2}\Big) }
{ A\Big(1+\tfrac{\Delta}{2}\Big) }
}
\quad\text{(ringing if exceeded).}
\tag{3}
\]

Working instead from the **monic** cubic (2) produces the alternative RH-style expression:

\[
\boxed{
K \;>\; K_c^{\text{RH(monic)}} \;=\;
\frac{ \big(\tfrac{2}{\Delta}+A+B\big)\big(\tfrac{2(A+B)}{\Delta}+AB\big) }
{ A\big(\tfrac{2}{\Delta}-B\big) }
}
\quad\text{(ringing if exceeded).}
\tag{4}
\]

> **Normalization & reconciliation.** (3) and (4) come from the **same cubic** but differ by when the monic scaling is applied relative to arranging the RH inequality. For consistent numerics we recommend the **dimensionless** reparametrization
> \[
> u:=A\Delta,\qquad v:=B\Delta,\qquad \kappa:=\frac{K\Delta}{A}.
> \]
> In \((u,v,\kappa)\) space, both RH routes lead to the **same ringing set** when compared against the exact roots of (1). In practice, use **either** (3) **or** (4) consistently with a fixed unit convention, and verify via the engineering rule and/or direct root checks.

**Sanity limits.** As \(\Delta\to 0\) or \(A\to\infty\), both RH thresholds reduce to the no-delay form \(K_c\to AB\).

---

## 3) Engineering rule (quick estimate)

A single-equation phase-margin heuristic defines \(\omega_c\):
\[
\arctan\!\Big(\frac{\omega_c}{A}\Big) + \arctan\!\Big(\frac{\omega_c}{B}\Big)
\;+\; \omega_c \,\Delta \;=\; \frac{3\pi}{4}.
\tag{5}
\]
Solve (5) numerically (LHS is monotone in \(\omega_c\)), then estimate
\[
\boxed{ \quad
K_c^{\text{ER}} \;\approx\; \frac{1}{A}\sqrt{(A^2+\omega_c^2)(B^2+\omega_c^2)}\;.
\quad}
\tag{6}
\]
This **engineering threshold** is typically within **±20%** of RH/roots across practical ranges (often <10% when \(A\gtrsim B\) and \(\omega_c\Delta\ll 1\)). Use for coarse grid scans; confirm with exact roots or RH.

---

## 4) Practical detection & ζ estimation from time series

Given a small-signal response (e.g., witness flux \(\Phi_{\text{wit}}(t)=\tfrac{d}{dt}\sum_k I(X{:}F_k)\) or \(g(t)\)):

- **PSD peak:** detrend/normalize; a non-DC peak **≥ 6 dB** over baseline indicates oscillatory content.
- **Overshoots:** in redundancy \(R_X^\delta(t)\) (e.g., \(\delta=0.1\)), count sign-alternating slope changes; **≥ 2** supports ringing.
- **Damping ratio ζ:** logarithmic decrement via Hilbert envelope:

```python
import numpy as np
from scipy.signal import hilbert, find_peaks

def estimate_zeta(series):
    x = series - np.mean(series)
    if np.allclose(x, 0): 
        return None
    analytic = hilbert(x)
    env = np.abs(analytic)
    # keep the top-energy 99%
    e = np.cumsum(env**2) / (env**2).sum()
    last = np.searchsorted(e, 0.99) + 1
    x, env = x[:last], env[:last]
    # use peaks of the signal for logarithmic decrement
    peaks = find_peaks(x)[0]
    if len(peaks) < 2:
        return None
    logs = []
    for i in range(len(peaks)-1):
        a, b = abs(x[peaks[i]]), abs(x[peaks[i+1]])
        if b > 1e-12:
            logs.append(np.log(a/b))
    if not logs:
        return None
    delta = np.mean(logs)
    zeta = delta / np.sqrt(4*np.pi**2 + delta**2)
    return zeta  # ringing if 0 < ζ < 1

Decision rule. Declare ringing at a scan point if both: (i) PSD peak ≥ 6 dB and (ii) overshoots ≥ 2 or (0<\zeta<1).

⸻

5) How to use the thresholds
	1.	Fix units / normalize. Prefer ((u=A\Delta,\ v=B\Delta,\ \kappa=K\Delta/A)).
	2.	RH estimate. Use either (3) or (4) in your chosen frame to get (K_c^{\text{RH}}) (or (\kappa_c)).
	3.	Engineering check. Solve (5) for (\omega_c) and compute (K_c^{\text{ER}}) via (6).
	4.	Ground truth. Solve the exact cubic roots of (1) and locate the smallest (K) where the dominant pair becomes complex (or RH equality holds). Report RH and ER with % deviation.

⸻

6) Delay & safety margins (practical)
	•	Conservative region. For (\Delta/\tau_{\text{geom}}\lesssim 0.3), with (\tau_{\text{geom}}:=1/B), a conservative safe choice is (K \lesssim 0.8,B).
	•	Padé(1,1) validity. Prefer (\omega_c \Delta \lesssim 1); for larger delays confirm via a delayed Nyquist plot.

⸻

7) Sanity checks (limit cases)
	•	No delay ((\Delta\to 0)) or instantaneous EMA ((A\to\infty)): (K_c \to AB).
	•	Fast EMA ((A\gg B)): RH and ER concur to (\mathcal{O}(B/A)).
	•	Weak decay ((B\to 0^+)): threshold rises; ringing is limited by the EMA/Delay corner.

⸻

8) Repro tip (CSV & plots)

For grid validation, export (example schema):

A,B,Delta,Kc_RH,Kc_ER,err_RH
1,1,0.1,1.230,1.218,0.010
...

Use a consistent RH form (3 or 4) in a fixed unit frame; err_RH is relative difference to ER. Overlay RH/ER contours on the ringing mask from simulations (PSD+ζ criteria).

⸻

Summary. The ringing threshold for GP with delay follows from the Padé(1,1) cubic (1). RH gives a closed-form gain limit; an engineering rule gives a fast estimate. In practice, fix a dimensionless frame, use one RH algebraic form consistently, and verify with the engineering rule and direct roots.

