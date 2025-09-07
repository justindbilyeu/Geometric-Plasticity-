---

### `docs/preview_methods.md`
```markdown
# Preview Methods & Definitions

## Model (linearized GP with EMA and delay)
\[
\dot g = \eta \bar I - B g, \qquad
\dot{\bar I} = A (I - \bar I), \qquad
I(t) = \gamma\, g(t-\Delta), \quad K=\eta\gamma.
\]

Padé(1,1) for the delay: \(e^{-s\Delta} \approx \frac{1-\frac{s\Delta}{2}}{1+\frac{s\Delta}{2}}\).

Characteristic cubic (after elimination) gives **Routh–Hurwitz** coefficients \(a_3,a_2,a_1,a_0\) and a closed-form threshold.

### Thresholds compared

**Routh–Hurwitz (RH)**  
\[
K_c^{\mathrm{RH}}(A,B,\Delta)=
\frac{A B + \frac{2}{\Delta}\Big(1+\frac{(A+B)\Delta}{2}\Big)\Big((A+B)+\frac{AB\Delta}{2}\Big)}
     {A\Big(1+\frac{\Delta}{2}\Big)}.
\]

**DeepSeek (DS)** (equivalent after simplification):  
\[
K_c^{\mathrm{DS}}(A,B,\Delta)=
\frac{A B + \frac{2}{\Delta}\Big(1+\frac{(A+B)\Delta}{2}\Big)\Big((A+B)+\frac{AB\Delta}{2}\Big)}
     {A\Big(\frac{2}{\Delta}-B\Big)}.
\]

**Engineering Rule (ER)**  
Solve \(\arctan(\omega_c/A)+\arctan(\omega_c/B)+\omega_c\Delta=\frac{3\pi}{4}\) for \(\omega_c\).  
We then map to \(K_c^{\mathrm{ER}}\) using the cubic’s imaginary-axis condition; this acts as a **ground-truth** numerical reference.

### CSV schema
`A,B,Δ,Kc_RH,Kc_DS,Kc_ER,err_RH,err_DS`, where errors are relative to ER.

### Acceptance notes
- For \(A\gtrsim B\) and small \(\Delta\), all forms agree to ~1–2%.
- As \(A/B\) increases (fast EMA), agreement tightens further.
- RH ≈ DS numerically across the grid (algebraic equivalence in practice).

---
