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

The third law fixes the entropy reference point by specifying how entropy behaves as temperature approaches absolute zero. Using Boltzmann’s microscopic definition, the section explains why an ideal crystal has $S(0)=0$, introduces residual entropy when the ground state is degenerate or disordered, and uses Haber–Bosch as a reminder that thermodynamic favorability is not the same as reaction speed.

In this section we introduce the **third law of thermodynamics** and connect it to the statistical-mechanical meaning of entropy.

- State **Planck’s (third-law) statement** and explain what “pure crystalline substance” and “ideal crystal” mean.
- Use **Boltzmann’s equation** to interpret why an **ideal crystal** has $S(0\,\mathrm{K})=0$.
- Define and compute **residual entropy** for a system whose ground state is *degenerate* (multiple microstates at $T\to 0$).
- Use the **Haber–Bosch process** to articulate the difference between **thermodynamics** (driving force/equilibrium) and **kinetics** (rate/activation barrier).

---

Learning objectives:

- State Planck’s formulation of the third law and define what is meant by a pure crystalline substance and an ideal crystal.
- Use $S=k_{\mathrm B}\ln\Omega$ to justify $S(0)=0$ for an ideal crystal with a unique ground state.
- Compute residual entropy $S_{\mathrm{res}}=k_{\mathrm B}\ln\Omega_0$ (or per mole $R\ln\Omega_0$) for degenerate/disordered ground states.
- Contrast thermodynamic driving force ($\Delta G$, equilibrium) with kinetics (activation barriers, rate, catalysts).

## Core Ideas and Derivations

### 5.2.1 Review: Boltzmann’s equation for entropy

A useful microscopic interpretation of entropy is

```{math}
S = k_B \ln \Omega
```

where:

- $k_B$ is Boltzmann’s constant, and
- $\Omega$ is the number of accessible microstates consistent with the macroscopic constraints.

This “counting” perspective is the bridge between statistical mechanics and thermodynamics.

---

### 5.2.2 Third law of thermodynamics

#### Planck’s statement

**Planck’s statement of the third law**:

> As $T \to 0\,\mathrm{K}$, the entropy of any **pure crystalline substance** tends to a **constant**.

A key special case is the **ideal (perfect) crystal**, meaning **no defects** (no vacancies, impurities, grain boundaries, stacking faults, etc.). In that case the constant is taken to be zero:

```{math}
S(0\,\mathrm{K}) = 0 \quad \text{(ideal crystal)}
```

#### Why the “ideal crystal” has $S(0)=0$

Using Boltzmann’s equation:

- An **ideal crystal** at $0\,\mathrm{K}$ is assumed to have a **unique ground-state arrangement**.
- That means $\Omega_0 = 1$, so

```{math}
S(0) = k_B \ln(1) = 0
```

#### Real crystals: defects and disorder

Real materials can deviate from the ideal picture. Even at very low temperature, a crystal may contain (or “freeze in”) structural or chemical disorder such as:

- **vacancies**
- **impurities**
- **grain boundaries**
- **stacking faults**

These imperfections can increase the number of distinct microscopic arrangements compatible with the same macroscopic state, i.e. increase $\Omega_0$, leading to a non-zero constant entropy as $T\to 0$.

---

### 5.2.3 Residual entropy

If a system approaches $T\to 0\,\mathrm{K}$ while retaining **more than one** accessible (or frozen-in) microscopic arrangement, then the entropy approaches a **non-zero constant**. This is often called **residual entropy**.

If the ground state is $\Omega_0$-fold degenerate, then

```{math}
S_{\text{res}} = \lim_{T\to 0} S(T) = k_B \ln \Omega_0
```

#### Example: orientational disorder in solid CO

In the lecture notes, solid CO at very low temperature is used as an example where the entropy at $0\,\mathrm{K}$ can be **non-zero** because molecules can be arranged in more than one orientation in the lattice (e.g., **CO** vs **OC** “head–tail” orientations). If each molecule has two equally likely orientations, then for $N$ molecules:

- microstates: $\Omega_0 = 2^N$
- residual entropy:

```{math}
S_{\text{res}}
= k_B \ln(2^N)
= Nk_B \ln 2
```

Per molecule, this corresponds to:

```{math}
s_{\text{res}} = k_B \ln 2
```

and per mole (optional conversion), since $R = N_A k_B$:

```{math}
S_{\text{res,m}} = R \ln 2
```

#### Reconciling residual entropy with the third law

Planck’s statement says $S\to$ **a constant** as $T\to 0$, not automatically *zero*.

- For an **ideal crystal**, $\Omega_0=1\Rightarrow S(0)=0$.
- For crystals with **frozen-in disorder**, $\Omega_0>1\Rightarrow S(0)=k_B\ln\Omega_0>0$.

So residual entropy is not a violation of the third law—it reflects that the system is not a perfect, uniquely ordered crystal at $0\,\mathrm{K}$.

---

### 5.2.4 Haber–Bosch process: thermodynamics feeds the world

A major industrial example where thermodynamics and kinetics both matter is the **Haber–Bosch process** for ammonia production:

```{math}
\mathrm{N_2(g) + 3H_2(g) \rightarrow 2NH_3(g)}
```

- $\mathrm{N_2}$ can be obtained from **air**.
- $\mathrm{H_2}$ is not “free”; it must be produced (e.g., by splitting water or via **steam methane reforming**, which produces CO and H$_2$).

Historical context from the notes:

- Process developed/invented in **1909**
- Recognized with a **1918 Nobel Prize**

#### Typical operating conditions

Industrial conditions are chosen as a compromise between equilibrium yield and reaction rate:

- **$400\text{–}500\,\text{°C}$**
- **$150\text{–}350\,\mathrm{bar}$**

---

### 5.2.5 Thermodynamics meets its match: kinetics

#### Thermodynamics (driving force)

For the Haber–Bosch reaction, the notes list standard reaction quantities at $298.15\,\mathrm{K}$:

- $\Delta H_r^{\circ}(298.15\,\mathrm{K}) \approx -92\,\mathrm{kJ\,mol^{-1}}$
- $\Delta S_r^{\circ}(298.15\,\mathrm{K}) \approx -198\,\mathrm{J\,K^{-1}\,mol^{-1}}$
- $\Delta G_r^{\circ}(298.15\,\mathrm{K}) \approx -33\,\mathrm{kJ\,mol^{-1}}$

You can see the consistency via:

```{math}
\Delta G = \Delta H - T\Delta S
```

With $\Delta S<0$, the $-T\Delta S$ term becomes **positive** (opposes spontaneity), which is one reason higher temperature can reduce the equilibrium favorability for ammonia synthesis even if it helps the rate.

Also note:

- The reaction goes from **4 moles of gas** to **2 moles of gas**, which helps explain why **high pressure** favors ammonia (Le Châtelier’s principle).

#### Kinetics (rate, activation barrier)

Even if $\Delta G < 0$ (thermodynamically favorable), the reaction may proceed **extremely slowly** if the **activation energy** is large.

The lecture notes include a reaction-coordinate sketch showing:

- a large barrier for the **uncatalyzed** pathway, and
- a much smaller barrier for the **catalyzed** pathway.

Key idea:

- A **catalyst lowers the activation barrier** (increases rate),
- but does **not** change $\Delta G$ for the overall reaction.

---

## Worked Example

### Residual entropy from orientational disorder

Suppose a crystal contains $N$ CO molecules, each of which can be frozen in two orientations (CO vs. OC) even as $T\to 0$.
If all orientations are accessible/frozen in, the ground-state multiplicity is

```{math}
\Omega_0 = 2^N.
```

Boltzmann’s equation gives

```{math}
S_{\mathrm{res}} = k_{\mathrm B}\ln\Omega_0 = k_{\mathrm B}\ln(2^N)=Nk_{\mathrm B}\ln 2.
```

Per mole (using $R=N_Ak_{\mathrm B}$):

```{math}
S_{\mathrm{res,m}} = R\ln 2 = (8.314)\ln 2 = 5.76\ \mathrm{J\,mol^{-1}\,K^{-1}}.
```

**Result.** Even as $T\to 0$, a disordered crystal can retain a finite entropy set by the degeneracy of its ground state.

## Concept Checks

1. Why does Planck’s statement allow $S(0)$ to approach a *constant* without requiring that constant to be zero?
2. What physical features distinguish an ideal crystal from a real crystal?
3. How can a system have residual entropy without violating the third law?
4. In the Haber–Bosch example, why does a catalyst affect the rate but not the equilibrium constant?

## Key Takeaways

- Third law: $S\to$ constant as $T\to 0$; ideal crystals set that constant to 0 by convention.
- Boltzmann counting explains $S(0)=0$ when the ground state is unique ($\Omega_0=1$).
- Residual entropy occurs when $\Omega_0>1$ due to degeneracy or frozen-in disorder.
- Thermodynamics tells what is favored at equilibrium; kinetics tells how fast it happens.
