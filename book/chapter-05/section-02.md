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

See also: [Course-wide Conventions & Notation](../notation.md)

## Overview and Learning Objectives

The first and second laws tell us *how* energy and entropy change but leave the entropy scale with an undetermined reference point — the integration constant in $\Delta S = \int \delta q_{\mathrm{rev}}/T$ is never fixed. The third law resolves this by specifying what happens to entropy as $T \to 0$. Using Boltzmann's formula ([Section 4.3](../chapter-04/section-03.md)), we explain why an ideal crystal has $S(0) = 0$, examine how ground-state degeneracy or frozen-in disorder produces a nonzero *residual entropy*, and note two further consequences: heat capacities must vanish as $T \to 0$, and absolute zero itself is unattainable in a finite number of steps.

With the entropy reference point fixed, we can compute **absolute** (third-law) entropies $S^\circ(T)$ for any substance by integrating $C_P/T$ from $0$ to $T$ and adding any phase-transition contributions. These absolute entropies are the basis for tabulated $\Delta_r S^\circ$ and $\Delta_r G^\circ$ values that we will put to use in [Section 5.3](section-03.md).

---

Learning objectives:

- State Planck's formulation of the third law and define what is meant by a pure crystalline substance and an ideal crystal.
- Use $S = k_{\mathrm B}\ln\Omega$ to justify $S(0) = 0$ for an ideal crystal with a unique ground state.
- Compute residual entropy $S_{\mathrm{res}} = k_{\mathrm B}\ln\Omega_0$ (or per mole, $R\ln\Omega_0$) for degenerate or disordered ground states.
- Explain why $C_V$ and $C_P$ must vanish as $T \to 0$ and connect this to the quantum behavior of heat capacities ([Section 2.6](../chapter-02/section-06.md)).
- Describe, qualitatively, why absolute zero is unattainable.

## Core Ideas and Derivations

### Third Law of Thermodynamics

#### Planck's Statement

**Planck's statement of the third law:**

> As $T \to 0\,\mathrm{K}$, the entropy of any **pure crystalline substance** tends to a **constant**.

A key special case is the **ideal (perfect) crystal** — one with no defects (no vacancies, impurities, grain boundaries, or stacking faults). For an ideal crystal, the constant is taken to be zero:

```{math}
:label: eq:third-law
S(0\,\mathrm{K}) = 0 \quad \text{(ideal crystal)}.
```

#### Why the Ideal Crystal Has $S(0) = 0$

Recall from [Section 4.3](../chapter-04/section-03.md) that Boltzmann's formula gives the entropy of an isolated system as

```{math}
S = k_{\mathrm B}\ln\Omega,
```

where $\Omega$ is the number of microstates compatible with the macroscopic constraints. An ideal crystal at $0\,\mathrm{K}$ has a unique ground-state arrangement: every atom sits at its designated lattice site in a single, perfectly ordered configuration. That means $\Omega_0 = 1$, so

```{math}
S(0) = k_{\mathrm B}\ln 1 = 0.
```

The third law, from the statistical-mechanical viewpoint, is the statement that a system with a unique ground state has zero entropy at absolute zero.

#### Real Crystals: Defects and Disorder

Real materials can deviate from this picture. Even at very low temperature, a crystal may contain — or "freeze in" — structural or chemical disorder: vacancies, impurities, grain boundaries, stacking faults, or orientational disorder of molecular units. These imperfections increase the number of distinct microscopic arrangements compatible with the same macroscopic state (i.e., increase $\Omega_0$), leading to a nonzero entropy as $T \to 0$.

---

### Residual Entropy

If a system retains **more than one** accessible microscopic arrangement as $T \to 0$, its entropy approaches a nonzero constant called the **residual entropy**:

```{math}
:label: eq:residual-entropy
S_{\mathrm{res}} = \lim_{T\to 0}\,S(T) = k_{\mathrm B}\ln\Omega_0.
```

#### Example: Orientational Disorder in Solid CO

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

```{admonition} Comparison with Experiment
:class: dropdown

The experimentally measured residual entropy of CO is approximately $4.6\ \mathrm{J\,mol^{-1}\,K^{-1}}$, somewhat less than the $R\ln 2 = 5.76\ \mathrm{J\,mol^{-1}\,K^{-1}}$ predicted by the fully random two-orientation model. The discrepancy indicates that CO is not *completely* disordered at low temperature — neighboring molecules have weak orientational correlations that partially reduce the number of accessible configurations. The simple $2^N$ estimate captures the right order of magnitude but overcounts somewhat.
```

```{admonition} Residual Entropy and the Third Law
:class: important

Return to Planck's statement at the start of this section and to the definition of the residual entropy, Eq. {eq}`eq:residual-entropy`. Which values of $\Omega_0$ does the statement allow? Is solid CO a counterexample to Planck's statement, to Eq. {eq}`eq:third-law`, or to neither?
```

The third law's practical utility lies in establishing a universal reference point for entropy: any ideal crystal has $S(0) = 0$.

---

### Consequence: Heat Capacities Vanish as $T \to 0$

The third law has a direct implication for heat capacities. At constant pressure, the entropy at temperature $T$ is related to the entropy at $0\,\mathrm{K}$ by

```{math}
S(T) = S(0) + \int_0^T \frac{C_P(T')}{T'}\,dT'.
```

For this integral to converge (i.e., for $S(T)$ to be finite), $C_P$ must vanish as $T \to 0$. If $C_P$ approached a nonzero constant, the integrand $C_P/T'$ would grow as $1/T'$, and the integral would diverge logarithmically at its lower limit. Therefore:

```{math}
C_P \to 0 \quad \text{as} \quad T \to 0.
```

The same argument applies to $C_V$. This is a nontrivial prediction: it says that the ability of any substance to absorb heat must diminish as absolute zero is approached.

````{admonition} Connection to Module 2
:class: important

[Section 2.6](../chapter-02/section-06.md) derived the heat capacity of the Einstein model of an atomic crystal, which in terms of $\Theta_{\mathrm{vib}} = \hbar\omega/k_{\mathrm B}$ reads

```{math}
C_V = 3Nk_{\mathrm B}\left(\frac{\Theta_{\mathrm{vib}}}{T}\right)^2 \frac{e^{\Theta_{\mathrm{vib}}/T}}{(e^{\Theta_{\mathrm{vib}}/T}-1)^2}.
```

Find its limit as $T \to 0$, and compare it with the classical (equipartition) value $C_V = 3Nk_{\mathrm B}$. Which of the two is compatible with the third law?
````

---

### Consequence: Unattainability of Absolute Zero

A closely related statement — sometimes called the **Nernst form** of the third law — is that no finite sequence of thermodynamic operations can bring a system to exactly $T = 0\,\mathrm{K}$.

The intuitive argument is as follows. Cooling to very low temperature typically alternates two steps, as in adiabatic demagnetization: an isothermal step that removes entropy from the system, and an adiabatic step that lowers the temperature at constant entropy. As $T \to 0$, the entropy of the system approaches the same limiting value whatever the magnetic field (or other control variable), so each successive pair of steps removes less entropy and lowers the temperature by less. The sequence converges but never quite reaches zero. In practice, temperatures below $1\ \mathrm{nK}$ have been achieved in laboratory settings, but $T = 0$ remains a limit, not an endpoint.

---

### Absolute Entropies

With $S(0) = 0$ established for ideal crystals, we can compute the **absolute (third-law) entropy** of any substance at temperature $T$ and pressure $P^\circ$:

```{math}
:label: eq:absolute-entropy
S^\circ(T) = \int_0^T \frac{C_P(T')}{T'}\,dT' + \sum_{\text{transitions}} \frac{\Delta H_{\mathrm{trs}}}{T_{\mathrm{trs}}},
```

where the sum runs over any phase transitions (melting, boiling, etc.) between $0$ and $T$, each contributing $\Delta H_{\mathrm{trs}}/T_{\mathrm{trs}}$ to the entropy.

These absolute entropies are tabulated in standard references (NIST WebBook, JANAF tables) as $S^\circ(298.15\ \mathrm{K})$. Combined with standard enthalpies of formation, they provide the reaction entropies and Gibbs energies that we will use in [Section 5.3](section-03.md) to analyze the Haber–Bosch process.

---

## Worked Examples

### Example 1: Entropy of a Crystal at Low Temperature

**Problem.** Calorimetric measurements cannot reach $0\,\mathrm{K}$, so the lowest part of the integral in Eq. {eq}`eq:absolute-entropy` is obtained by extrapolation. At sufficiently low temperature, the heat capacity of a nonmetallic crystal is proportional to $T^3$ (the Debye $T^3$ law): $C_P = aT^3$. A crystal has a measured molar heat capacity of $0.24\ \mathrm{J\,mol^{-1}\,K^{-1}}$ at $10\ \mathrm{K}$, and the $T^3$ law holds below this temperature. Find its molar entropy at $10\ \mathrm{K}$, assuming an ideal crystal.

**Solution.**

For an ideal crystal $S(0) = 0$, and no phase transition lies below $10\ \mathrm{K}$, so

```{math}
S(T) = \int_0^T \frac{C_P(T')}{T'}\,dT' = \int_0^T a\,T'^2\,dT' = \frac{aT^3}{3} = \frac{C_P(T)}{3}.
```

At $T = 10\ \mathrm{K}$:

```{math}
S_{\mathrm{m}}(10\ \mathrm{K}) = \frac{0.24\ \mathrm{J\,mol^{-1}\,K^{-1}}}{3} = 0.080\ \mathrm{J\,mol^{-1}\,K^{-1}}.
```

**Result.** Where the $T^3$ law holds, the entropy is one third of the heat capacity at the same temperature, so a single measured value fixes the contribution of the unmeasured range. The integral converges because $C_P$ vanishes as $T \to 0$.

---

### Example 2: Residual Entropy of Ice (Pauling Model)

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

**Result.** The predicted residual entropy is $3.37\ \mathrm{J\,mol^{-1}\,K^{-1}}$, in remarkable agreement with the experimental value of approximately $3.41\ \mathrm{J\,mol^{-1}\,K^{-1}}$. Without the ice rules, the $2N$ hydrogen atoms of $N$ molecules, each with two possible positions, would have $2^{2N} = 4^N$ arrangements. The ice rules reduce this number, and the agreement with experiment confirms that the proton disorder in ice is well described by the ice rules even at the lowest accessible temperatures.

```{admonition} Pauling's Count
:class: important

Reconstruct Pauling's estimate. Each hydrogen atom lies between two neighboring oxygen atoms and has two possible positions, close to one oxygen atom or close to the other. (a) Around one oxygen atom, the four neighboring hydrogen atoms have $2^4 = 16$ arrangements. How many of them place exactly two hydrogen atoms close to that oxygen atom, as the ice rules require? (b) Treat the $N$ oxygen atoms as independent, and combine the fraction from (a) with the $4^N$ unconstrained arrangements to obtain $\Omega_0$.
```

## Concept Checks

1. Why does Planck's statement allow $S(0)$ to approach a *constant* without requiring that constant to be zero?
2. What physical features distinguish an ideal crystal from a real crystal, and how do these features affect $\Omega_0$?
3. How can a system have residual entropy without violating the third law?
4. Why does the third law require $C_P \to 0$ as $T \to 0$, and how is this consistent with the quantum heat capacities derived in Section 2.6?
5. In the Pauling ice model, why is $\Omega_0 = (3/2)^N$ rather than $4^N$ or $1$?

## Key Takeaways

- The third law states that $S \to$ a constant as $T \to 0$; for ideal crystals, that constant is zero, fixing the entropy scale.
- Boltzmann's formula explains $S(0) = 0$ when the ground state is unique ($\Omega_0 = 1$).
- Residual entropy ($S_{\mathrm{res}} = k_{\mathrm B}\ln\Omega_0 > 0$) arises from ground-state degeneracy or frozen-in disorder.
- The third law requires heat capacities to vanish at absolute zero, a prediction that can be tested against the heat capacities of [Section 2.6](../chapter-02/section-06.md).
- Absolute (third-law) entropies $S^\circ(T)$ are computed by integrating $C_P/T$ from $0$ to $T$ and provide the standard entropy values used in thermochemistry.
