# Resonance Fold Operator (RFO) — Overview

**Purpose.** The RFO formalizes the subspace where pointer information is redundantly recorded across environment fragments. GP’s plasticity drives the dynamics toward this subspace.

**Definition (operator-theoretic).** Let \(X\) be the pointer observable on \(S\) with eigenbasis \(\{|x\rangle\}\). An **RFO isometry** \(V:\mathcal H_S\to \mathcal H_S\otimes \mathcal H_{E^\star}\) satisfies
\(V|x\rangle = |x\rangle \otimes \bigotimes_k |r_k(x)\rangle\),
with (approximately) orthogonal record states \(|r_k(x)\rangle\) across fragments \(F_k\).
The **folded subspace** \(\mathcal C\) is spanned by these branch states, and the **RFO projector** is
\(P_{\rm RFO} := V V^\dagger\) (partial isometry: \(V^\dagger V=\mathbb I_S\), \(VV^\dagger=P_{\mathcal C}\)).

**Fold distance.** For a joint state \(\rho_{SF}\), \(d_{\rm fold}(\rho)=1-\mathrm{Tr}(P_{\rm RFO}\rho)\) quantifies “how unfolded” the state is.

**Data-driven approximation.** Given simulation access to \(\rho_{SF}\) and pointer \(X\):
1) Condition fragments on \(X=x\) to get \(\rho_{F_k|x}\).
2) Let \(R_k(x)\) be the top-eigenspace projector of \(\rho_{F_k|x}\).
3) Define \(\tilde P_{\rm RFO} := \sum_x |x\rangle\langle x| \otimes \bigotimes_k R_k(x)\).
Then \(d_{\rm fold} \approx 1-\mathrm{Tr}(\tilde P_{\rm RFO}\rho_{SF})\).

**Why it matters.**  
- Links redundancy to a concrete subspace/code picture (SBS-like).  
- Supplies a scalar **objective** GP can implicitly minimize.  
- Bridges to QEC: \(\mathcal C\) behaves like a code space for \(X\)-records.

**Predictions.**  
- Below the ringing threshold \(K_c\), \(d_{\rm fold}(t)\) decreases to a plateau.  
- Hysteresis resonance modulates \(d_{\rm fold}\) with a peak near \(T\approx 2\pi\tau_{\rm geom}\).  
- Motif transitions (β/λ sweep) change the structure of \(\tilde P_{\rm RFO}\) (broadcast → modular).
