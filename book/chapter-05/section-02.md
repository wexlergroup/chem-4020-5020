---
jupytext:
  formats: md:myst
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.11.5
kernelspec:
  display_name: Python 3
  language: python
  name: python3
---


# 5.2. Third Law

[Course-wide Conventions & Notation](../notation.md)

## Overview and Learning Objectives

The first and second laws tell us *how* energy and entropy change but leave the entropy scale with an undetermined reference point — the integration constant in $\Delta S = \int \delta q_{\mathrm{rev}}/T$ is never fixed. The third law resolves this by specifying what happens to entropy as $T \to 0$. Using Boltzmann's formula (Section 4.3), we explain why an ideal crystal has $S(0) = 0$, examine how ground-state degeneracy or frozen-in disorder produces a nonzero *residual entropy*, and note two further consequences: heat capacities must vanish as $T \to 0$, and absolute zero itself is unattainable in a finite number of steps.

With the entropy reference point fixed, we can compute **absolute** (third-law) entropies $S^\circ(T)$ for any substance by integrating $C_P/T$ from $0$ to $T$ and adding any phase-transition contributions. These absolute entropies are the basis for tabulated $\Delta_r S^\circ$ and $\Delta_r G^\circ$ values that we will put to use in Section 5.3.

---

Learning objectives:

- State Planck's formulation of the third law and define what is meant by a pure crystalline substance and an ideal crystal.
- Use $S = k_{\mathrm B}\ln\Omega$ to justify $S(0) = 0$ for an ideal crystal with a unique ground state.
- Compute residual entropy $S_{\mathrm{res}} = k_{\mathrm B}\ln\Omega_0$ (or per mole, $R\ln\Omega_0$) for degenerate or disordered ground states.
- Explain why $C_V$ and $C_P$ must vanish as $T \to 0$ and connect this to the quantum behavior of heat capacities (Section 2.6).
- Describe, qualitatively, why absolute zero is unattainable.

## Core Ideas and Derivations

### Third law of thermodynamics

#### Planck's statement

**Planck's statement of the third law:**

> As $T \to 0\,\mathrm{K}$, the entropy of any **pure crystalline substance** tends to a **constant**.

A key special case is the **ideal (perfect) crystal** — one with no defects (no vacancies, impurities, grain boundaries, or stacking faults). For an ideal crystal, the constant is taken to be zero:

```{math}
:label: eq:third-law
S(0\,\mathrm{K}) = 0 \quad \text{(ideal crystal)}.
```

#### Why the ideal crystal has $S(0) = 0$

Recall from Section 4.3 that Boltzmann's formula gives the entropy of an isolated system as

```{math}
S = k_{\mathrm B}\ln\Omega,
```

where $\Omega$ is the number of microstates compatible with the macroscopic constraints. An ideal crystal at $0\,\mathrm{K}$ has a unique ground-state arrangement: every atom sits at its designated lattice site in a single, perfectly ordered configuration. That means $\Omega_0 = 1$, so

```{math}
S(0) = k_{\mathrm B}\ln 1 = 0.
```

The third law, from the statistical-mechanical viewpoint, is the statement that a system with a unique ground state has zero entropy at absolute zero.

#### Real crystals: defects and disorder

Real materials can deviate from this picture. Even at very low temperature, a crystal may contain — or "freeze in" — structural or chemical disorder: vacancies, impurities, grain boundaries, stacking faults, or orientational disorder of molecular units. These imperfections increase the number of distinct microscopic arrangements compatible with the same macroscopic state (i.e., increase $\Omega_0$), leading to a nonzero entropy as $T \to 0$.

---

### Residual entropy

If a system retains **more than one** accessible microscopic arrangement as $T \to 0$, its entropy approaches a nonzero constant called the **residual entropy**:

```{math}
:label: eq:residual-entropy
S_{\mathrm{res}} = \lim_{T\to 0}\,S(T) = k_{\mathrm B}\ln\Omega_0.
```

#### Example: orientational disorder in solid CO

Solid CO at very low temperature provides a classic example. The CO molecule has a small dipole moment, and the two orientations CO and OC are nearly isoenergetic. If each of $N$ molecules can be frozen in either orientation with roughly equal probability, the ground-state multiplicity is

```{math}
\Omega_0 = 2^N,
```

and the residual entropy is

```{math}
S_{\mathrm{res}} = k_{\mathrm B}\ln(2^N) = Nk_{\mathrm B}\ln 2.
```

Per mole, using $R = N_Ak_{\mathrm B}$:

```{math}
S_{\mathrm{res,m}} = R\ln 2 \approx 5.76\ \mathrm{J\,mol^{-1}\,K^{-1}}.
```

```{admonition} Comparison with experiment
:class: dropdown

The experimentally measured residual entropy of CO is approximately $4.6\ \mathrm{J\,mol^{-1}\,K^{-1}}$, somewhat less than the $R\ln 2 = 5.76\ \mathrm{J\,mol^{-1}\,K^{-1}}$ predicted by the fully random two-orientation model. The discrepancy indicates that CO is not *completely* disordered at low temperature — neighboring molecules have weak orientational correlations that partially reduce the number of accessible configurations. The simple $2^N$ estimate captures the right order of magnitude but overcounts somewhat.
```

#### Reconciling residual entropy with the third law

Planck's statement says $S \to$ **a constant** as $T \to 0$, not that the constant must be zero.

- For an **ideal crystal**: $\Omega_0 = 1 \Rightarrow S(0) = 0$.
- For a crystal with **frozen-in disorder**: $\Omega_0 > 1 \Rightarrow S(0) = k_{\mathrm B}\ln\Omega_0 > 0$.

Residual entropy is not a violation of the third law — it reflects that the system is not a perfectly ordered crystal at $0\,\mathrm{K}$. The third law's practical utility lies in establishing a universal reference point for entropy: any ideal crystal has $S(0) = 0$, and deviations from zero are ascribed to measurable disorder.

---

### Consequence: heat capacities vanish as $T \to 0$

The third law has a direct implication for heat capacities. At constant pressure, the entropy at temperature $T$ is related to the entropy at $0\,\mathrm{K}$ by

```{math}
S(T) = S(0) + \int_0^T \frac{C_P(T')}{T'}\,dT'.
```

For this integral to converge (i.e., for $S(T)$ to be finite), $C_P$ must go to zero at least as fast as $T$ does. More precisely:

```{math}
C_P \to 0 \quad \text{as} \quad T \to 0.
```

The same argument applies to $C_V$. This is a nontrivial prediction: it says that the ability of any substance to absorb heat must diminish as absolute zero is approached.

```{admonition} Connection to Chapter 2
:class: tip

In Section 2.6 we saw this behavior explicitly. The Einstein model of an atomic crystal gives a heat capacity

$$
C_V = 3Nk_{\mathrm B}\left(\frac{\Theta_{\mathrm{vib}}}{T}\right)^2 \frac{e^{\Theta_{\mathrm{vib}}/T}}{(e^{\Theta_{\mathrm{vib}}/T}-1)^2},
$$

which goes to zero exponentially as $T \to 0$ because the quantum vibrational levels "freeze out." The classical (equipartition) prediction $C_V = 3Nk_{\mathrm B}$, by contrast, is temperature-independent and violates the third law. The resolution is that equipartition is a high-temperature approximation; quantum mechanics is needed at low $T$.
```

---

### Consequence: unattainability of absolute zero

A closely related statement — sometimes called the **Nernst form** of the third law — is that no finite sequence of thermodynamic operations can bring a system to exactly $T = 0\,\mathrm{K}$.

The intuitive argument is as follows. Cooling typically works by removing entropy (e.g., by adiabatic demagnetization or expansion). As $T \to 0$, the entropy of the system approaches a minimum, and each successive cooling step removes less and less entropy — the "gap" between successive stages shrinks and the process converges but never quite reaches zero. In practice, temperatures below $1\ \mathrm{nK}$ have been achieved in laboratory settings, but $T = 0$ remains a limit, not an endpoint.

---

### Absolute entropies

With $S(0) = 0$ established for ideal crystals, we can compute the **absolute (third-law) entropy** of any substance at temperature $T$ and pressure $P^\circ$:

```{math}
:label: eq:absolute-entropy
S^\circ(T) = \int_0^T \frac{C_P(T')}{T'}\,dT' + \sum_{\text{transitions}} \frac{\Delta H_{\mathrm{trs}}}{T_{\mathrm{trs}}},
```

where the sum runs over any phase transitions (melting, boiling, etc.) between $0$ and $T$, each contributing $\Delta H_{\mathrm{trs}}/T_{\mathrm{trs}}$ to the entropy.

These absolute entropies are tabulated in standard references (NIST WebBook, JANAF tables) as $S^\circ(298.15\ \mathrm{K})$. Combined with standard enthalpies of formation, they provide the reaction entropies and Gibbs energies that we will use in Section 5.3 to analyze the Haber–Bosch process.

---

## Worked Examples

### Example 1: Residual entropy of solid CO

**Problem.** Solid CO has $N$ molecules per crystal, each of which can adopt either of two orientations (CO or OC) as $T \to 0$. Estimate the molar residual entropy assuming complete disorder.

**Solution.**

If all $2^N$ orientational configurations are equally accessible, then $\Omega_0 = 2^N$ and

```{math}
S_{\mathrm{res}} = k_{\mathrm B}\ln(2^N) = Nk_{\mathrm B}\ln 2.
```

Per mole:

```{math}
S_{\mathrm{res,m}} = R\ln 2 = (8.314\ \mathrm{J\,mol^{-1}\,K^{-1}})\ln 2 = 5.76\ \mathrm{J\,mol^{-1}\,K^{-1}}.
```

**Result.** Even as $T \to 0$, a disordered crystal can retain a finite entropy determined by the degeneracy of its ground state. The experimental value ($\approx 4.6\ \mathrm{J\,mol^{-1}\,K^{-1}}$) is somewhat smaller, indicating partial orientational ordering.

---

### Example 2: Residual entropy of ice (Pauling model)

**Problem.** In Linus Pauling's 1935 model of ice, each oxygen atom is tetrahedrally coordinated to four neighbors, with two hydrogen atoms nearby (covalent O–H bonds) and two farther away (hydrogen bonds). Pauling showed that the number of configurations consistent with the "ice rules" (exactly two H atoms close to each O) gives a ground-state multiplicity of approximately $\Omega_0 \approx (3/2)^N$ for $N$ water molecules. Estimate the molar residual entropy.

**Solution.**

With $\Omega_0 = (3/2)^N$:

```{math}
S_{\mathrm{res}} = k_{\mathrm B}\ln\!\left(\frac{3}{2}\right)^N = Nk_{\mathrm B}\ln\!\left(\frac{3}{2}\right).
```

Per mole:

```{math}
S_{\mathrm{res,m}} = R\ln\!\left(\frac{3}{2}\right) = (8.314)\ln(1.5) = 3.37\ \mathrm{J\,mol^{-1}\,K^{-1}}.
```

**Result.** The predicted residual entropy is $3.37\ \mathrm{J\,mol^{-1}\,K^{-1}}$, in remarkable agreement with the experimental value of approximately $3.41\ \mathrm{J\,mol^{-1}\,K^{-1}}$. The key insight is that $\Omega_0$ is not simply $2^N$ (two positions per H atom without constraints) but is reduced by the ice rules, which require exactly two short and two long O–H distances around each oxygen. Pauling's combinatorial analysis captures this constraint, and the agreement with experiment confirms that the proton disorder in ice is well described by the ice rules even at the lowest accessible temperatures.

## Concept Checks

1. Why does Planck's statement allow $S(0)$ to approach a *constant* without requiring that constant to be zero?
2. What physical features distinguish an ideal crystal from a real crystal, and how do these features affect $\Omega_0$?
3. How can a system have residual entropy without violating the third law?
4. Why does the third law require $C_P \to 0$ as $T \to 0$, and how is this consistent with the quantum heat capacities derived in Section 2.6?
5. In the Pauling ice model, why is $\Omega_0 = (3/2)^N$ rather than $2^N$ or $1$?

## Key Takeaways

- The third law states that $S \to$ a constant as $T \to 0$; for ideal crystals, that constant is zero, fixing the entropy scale.
- Boltzmann's formula explains $S(0) = 0$ when the ground state is unique ($\Omega_0 = 1$).
- Residual entropy ($S_{\mathrm{res}} = k_{\mathrm B}\ln\Omega_0 > 0$) arises from ground-state degeneracy or frozen-in disorder, not from a failure of the third law.
- The third law requires heat capacities to vanish at absolute zero — a prediction confirmed by quantum statistical mechanics (Section 2.6) and violated by classical equipartition.
- Absolute (third-law) entropies $S^\circ(T)$ are computed by integrating $C_P/T$ from $0$ to $T$ and provide the standard entropy values used in thermochemistry.
