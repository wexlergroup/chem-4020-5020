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


# 4.3. Microscopic View of Entropy

See also: [Course-wide Conventions & Notation](../notation.md)

## Overview and Learning Objectives

[Section 4.1](section-01.md) and [Section 4.2](section-02.md) developed entropy entirely from macroscopic reasoning: a state function defined by $dS = \delta q_{\mathrm{rev}}/T$, used to derive engine efficiency limits and the direction of spontaneous heat flow. This section shifts to the **microscopic** side, connecting entropy to the statistical-mechanical framework from Module 2. The payoff is substantial: entropy becomes a measure of *how many ways* a system can realize a macrostate, and its connection to the partition function $Q$ provides a direct bridge between microscopic energy levels and macroscopic thermodynamic potentials.

We introduce three formulas for entropy — Boltzmann's $S = k_{\mathrm{B}}\ln\Omega$ (isolated systems), the more general Gibbs form $S = -k_{\mathrm{B}}\sum p_i\ln p_i$ (any ensemble), and the canonical identity $S = U/T + k_{\mathrm{B}}\ln Q$ (the Gibbs form evaluated for the canonical ensemble) — and show that they are mutually consistent. The section culminates in the identification of the Helmholtz free energy $A = -k_{\mathrm{B}}T\ln Q$, which was foreshadowed in [Section 2.3](../chapter-02/section-03.md) and [Section 3.3](../chapter-03/section-03.md).

Learning objectives:

- State the second-law condition $dS \ge 0$ for isolated systems and interpret equilibrium as $dS = 0$.
- Use Boltzmann's formula $S = k_{\mathrm{B}}\ln\Omega$ and explain why the logarithm ensures extensivity.
- Compute entropy from a probability distribution via $S = -k_{\mathrm{B}}\sum_i p_i\ln p_i$.
- Derive the canonical identity $S = U/T + k_{\mathrm{B}}\ln Q$.
- Identify the Helmholtz free energy $A = U - TS = -k_{\mathrm{B}}T\ln Q$ as a consequence of the canonical entropy formula.

## Core Ideas and Derivations

### Entropy Increases Until the System Reaches Equilibrium

The second law for an **isolated system** ([Section 4.2](section-02.md)) requires

```{math}
dS \ge 0,
```

with equality only at equilibrium. The qualitative picture is:

- **Out of equilibrium:** entropy increases spontaneously ($dS > 0$).
- **At equilibrium:** entropy reaches its maximum value subject to the constraints and stops changing ($dS = 0$).

```{code-cell} ipython3
:tags: [hide-input]

import numpy as np
import matplotlib.pyplot as plt

# Schematic "approach to equilibrium" curve
t = np.linspace(0, 10, 400)
S_i = 1.0
S_f = 2.0

# Smooth approach to plateau + slight wiggle to mimic "spontaneous process" steps
S = S_f - (S_f - S_i) * np.exp(-t/2.0) + 0.03*np.sin(5*t)*np.exp(-t/2.0)

fig, ax = plt.subplots(figsize=(4.5, 3))
ax.plot(t, S)
ax.axhline(S_f, linestyle="--")
ax.set_xlabel("time (schematic)")
ax.set_ylabel("entropy (schematic)")
ax.set_title("Entropy increases until equilibrium (isolated system)", fontsize=10, pad=10)
ax.text(0.4, S_i+0.05, r"$S_i$")
ax.text(7.5, S_f+0.03, r"$S_f$")
ax.set_ylim(S_i-0.1, S_f+0.2)
ax.grid(True, alpha=0.3)
plt.show()
plt.close(fig)
```

Entropy rises from $S_i$ to $S_f$ during a spontaneous process and then becomes constant at equilibrium.

A concrete example from [Section 4.2](section-02.md) illustrates this. Two subsystems at different temperatures ($T_A > T_B$) in an insulated box exchange heat, and we showed that

```{math}
dS = dU_A\!\left(\frac{1}{T_A} - \frac{1}{T_B}\right) > 0
```

as long as $T_A \neq T_B$. Energy flows from $A$ to $B$, raising $T_B$ and lowering $T_A$, until $T_A = T_B$. At that point the factor $(1/T_A - 1/T_B)$ vanishes, $dS = 0$, and the system has reached thermal equilibrium at the entropy maximum.

---

### The Clausius Inequality

The entropy definition $dS = \delta q_{\mathrm{rev}}/T$ applies specifically to reversible heat transfer. A more general statement, valid for both reversible and irreversible processes, is the **Clausius inequality**:

```{math}
dS \ge \frac{\delta q}{T}.
```

- **Equality** holds for a **reversible** process: $dS = \delta q_{\mathrm{rev}}/T$.
- **Strict inequality** holds for an **irreversible** process: $dS > \delta q/T$.

Here $T$ is the temperature at which the heat is transferred to the system. For heat exchanged with a reservoir, it is the temperature of the reservoir.

The physical interpretation is that irreversibility *produces* additional entropy beyond the entropy "carried in" by heat flow. For an isolated system ($\delta q = 0$), the Clausius inequality reduces to $dS \ge 0$, recovering the second law.

The Clausius inequality also shows why the Carnot engine ([Section 4.2](section-02.md)) sets the maximum efficiency. Consider any engine, reversible or not, that in each cycle absorbs heat $q_{\mathrm{hot}} > 0$ from the hot reservoir and releases heat to the cold reservoir ($q_{\mathrm{cold}} < 0$). For the Carnot engine these are $q_{AB}$ and $q_{CD}$. The working substance returns to its initial state, so its entropy change over the cycle is zero, and the Clausius inequality, summed over the cycle, gives

```{math}
0 = \Delta S_{\mathrm{cycle}} \ge \frac{q_{\mathrm{hot}}}{T_{\mathrm{hot}}} + \frac{q_{\mathrm{cold}}}{T_{\mathrm{cold}}}.
```

Because $q_{\mathrm{hot}} > 0$, this rearranges to $q_{\mathrm{cold}}/q_{\mathrm{hot}} \le -T_{\mathrm{cold}}/T_{\mathrm{hot}}$. The efficiency, defined as in [Section 4.2](section-02.md), is therefore

```{math}
\eta = \frac{q_{\mathrm{hot}} + q_{\mathrm{cold}}}{q_{\mathrm{hot}}}
= 1 + \frac{q_{\mathrm{cold}}}{q_{\mathrm{hot}}}
\le 1 - \frac{T_{\mathrm{cold}}}{T_{\mathrm{hot}}} = \eta_{\mathrm{Carnot}}.
```

Equality holds only for a reversible engine. An irreversible engine releases more heat to the cold reservoir for the same heat absorbed, so it delivers less work.

```{admonition} Practical takeaway
:class: tip
Entropy is a **state function**, but $\delta q$ is not. To compute $\Delta S$ for any process — reversible or irreversible — choose a *convenient reversible path* between the same initial and final states and evaluate $\int \delta q_{\mathrm{rev}}/T$. This is the strategy we used in [Section 4.1](section-01.md), Example 2 (free expansion).
```

---

### Entropy and the Number of Microstates: Boltzmann's Formula

We now connect entropy to the microscopic picture. For an **isolated system** (microcanonical ensemble), the fundamental postulate of statistical mechanics ([Section 2.1](../chapter-02/section-01.md)) assigns equal probability $p_i = 1/\Omega$ to each of the $\Omega$ accessible microstates. The entropy of such a system is

```{math}
S = k_{\mathrm{B}} \ln \Omega.
```

Here $\Omega$ is the number of microstates compatible with the macroscopic constraints (fixed $U$, $V$, $N$). You can think of $\Omega$ as the degeneracy of the macrostate: fix the total energy and other extensive quantities, then count how many distinct microscopic configurations share those values.

```{admonition} Scope of Boltzmann's formula
:class: note
Boltzmann's $S = k_{\mathrm{B}}\ln\Omega$ applies to isolated systems where all accessible microstates are equally probable (the microcanonical ensemble). For systems in thermal contact with a reservoir — where microstates are *not* equally probable — we need the more general Gibbs formula introduced below.
```

#### Why the Logarithm?

A key property of thermodynamic entropy is that it is **extensive**: for two independent subsystems $A$ and $B$,

```{math}
S_{\mathrm{total}} = S_A + S_B.
```

But the number of microstates is *multiplicative* for independent systems:

```{math}
\Omega_{\mathrm{total}} = \Omega_A\,\Omega_B.
```

The logarithm converts the product into a sum:

```{math}
S_{\mathrm{total}}
= k_{\mathrm{B}} \ln(\Omega_A\,\Omega_B)
= k_{\mathrm{B}} \ln \Omega_A + k_{\mathrm{B}} \ln \Omega_B
= S_A + S_B.
```

Without the logarithm, entropy would not be additive for independent subsystems.

---

### Entropy in Terms of Probabilities: The Gibbs Entropy

The Gibbs (or Gibbs–Shannon) formula generalizes Boltzmann's result from "counting equally likely microstates" to working with an arbitrary probability distribution $\{p_i\}$ over microstates $i = 1, \dots, M$:

```{math}
S = -k_{\mathrm{B}} \sum_{i=1}^{M} p_i\,\ln p_i.
```

A microstate with $p_i = 0$ contributes nothing to the sum, because $p\ln p \to 0$ as $p \to 0$. For two independent subsystems $A$ and $B$, the probability that $A$ is in microstate $i$ while $B$ is in microstate $j$ is the product $p_i^{(A)}\,p_j^{(B)}$.

This formula has two important features:

1. **It reduces to Boltzmann's formula in the microcanonical case.** If all $\Omega$ accessible microstates are equally likely ($p_i = 1/\Omega$), then

   ```{math}
   S
   = -k_{\mathrm{B}} \sum_{i=1}^{\Omega} \frac{1}{\Omega}\ln\!\left(\frac{1}{\Omega}\right)
   = -k_{\mathrm{B}}\ln\!\left(\frac{1}{\Omega}\right)
   = k_{\mathrm{B}}\ln\Omega.
   ```

2. **It applies when microstates are not equally likely** — for instance, when the system is in contact with a heat bath (canonical ensemble). In this case, $p_i = e^{-\beta E_i}/Q$, and different microstates have different probabilities.

---

### Canonical Ensemble: Entropy in Terms of the Partition Function

We now evaluate the Gibbs entropy for the canonical ensemble. Recall from [Section 2.2](../chapter-02/section-02.md) that a closed system (fixed $N$, $V$) in thermal contact with a reservoir at temperature $T$ has microstate probabilities

```{math}
p_i = \frac{e^{-\beta E_i}}{Q},
\qquad
\beta \equiv \frac{1}{k_{\mathrm{B}} T},
```

where the canonical partition function is

```{math}
Q \equiv \sum_{i=1}^{M} e^{-\beta E_i}.
```

#### Derivation

Starting from the Gibbs entropy,

```{math}
S = -k_{\mathrm{B}} \sum_i p_i\ln p_i,
```

insert $\ln p_i = -\beta E_i - \ln Q$:

```{math}
\begin{aligned}
S
&= -k_{\mathrm{B}} \sum_i p_i\,(-\beta E_i - \ln Q) \\
&= k_{\mathrm{B}}\beta \sum_i p_i E_i \;+\; k_{\mathrm{B}}(\ln Q)\sum_i p_i.
\end{aligned}
```

Using $\sum_i p_i = 1$ and the definition of internal energy from [Section 2.3](../chapter-02/section-03.md),

```{math}
U \equiv \langle E \rangle = \sum_i p_i E_i,
```

we obtain

```{math}
S = k_{\mathrm{B}}\beta U + k_{\mathrm{B}}\ln Q.
```

Substituting $\beta = 1/(k_{\mathrm{B}} T)$:

```{math}
\boxed{
S = \frac{U}{T} + k_{\mathrm{B}}\ln Q.
}
```

This is a central result: it expresses the macroscopic state function $S$ directly in terms of the partition function $Q$ and the internal energy $U$, both of which we know how to compute from microscopic energy levels.

---

### Connection to the Helmholtz Free Energy

Rearranging the boxed result gives

```{math}
U - TS = -k_{\mathrm{B}}T\ln Q.
```

The left-hand side is the **Helmholtz free energy**, $A \equiv U - TS$ — a state function whose natural variables are $T$ and $V$. We therefore have

```{math}
\boxed{
A = -k_{\mathrm{B}}T\ln Q.
}
```

This result was foreshadowed in [Section 2.3](../chapter-02/section-03.md) (where it was stated without derivation) and in [Section 3.3](../chapter-03/section-03.md) (where we noted that $A = U - TS$ would be "convenient at constant $T$ and $V$"). It is arguably the single most important equation in canonical statistical mechanics: **if you can compute $Q$, you can compute $A$, and from $A$ you can derive all other thermodynamic quantities** — $S$, $U$, $P$, $C_V$, and chemical potentials — by taking appropriate derivatives. [Section 5.1](../chapter-05/section-01.md) develops this machinery.

```{admonition} Defining $A$ by the same strategy as $H$
:class: note
In [Section 3.3](../chapter-03/section-03.md), we defined enthalpy as $H = U + PV$ to simplify the First Law at constant pressure. Defining $A = U - TS$ is the same kind of move: a Legendre transform ([Section 5.1](../chapter-05/section-01.md)) that replaces $S$ (hard to control experimentally) with $T$ (easy to control) as an independent variable. Just as $q_P = \Delta H$ at constant pressure, we will see that $A$ plays a similarly natural role at constant temperature and volume.
```

---

### Microscopic Interpretation of Heat (Connection to Section 3.2)

The results above connect naturally to the microscopic interpretation of the First Law developed in [Section 3.2](../chapter-03/section-02.md). Differentiating $U = \sum_i p_i E_i$ gives

```{math}
dU = \sum_i p_i\, dE_i + \sum_i E_i\, dp_i.
```

For a reversible process in a closed system with only $PV$ work, comparing with $dU = \delta q_{\mathrm{rev}} + \delta w_{\mathrm{rev}} = \delta q_{\mathrm{rev}} - P\,dV$, we identified ([Section 3.2](../chapter-03/section-02.md)):

- **Heat**: $\delta q_{\mathrm{rev}} = \sum_i E_i\, dp_i$ — energy exchange via redistributing probabilities among fixed energy levels.
- **Work**: $\delta w_{\mathrm{rev}} = \sum_i p_i\, dE_i$ — energy exchange via shifting the energy levels themselves (e.g., changing volume changes the particle-in-a-box levels).

The Gibbs entropy formula leads to the same identification. Its differential is

```{math}
dS = -k_{\mathrm{B}}\sum_i \left(\ln p_i + 1\right) dp_i = -k_{\mathrm{B}}\sum_i \ln p_i\, dp_i,
```

because $\sum_i dp_i = 0$. Along a reversible path the system stays in equilibrium, so the probabilities remain canonical, with $\ln p_i = -E_i/(k_{\mathrm{B}}T) - \ln Q$. Substituting, and using $\sum_i dp_i = 0$ once more,

```{math}
dS = \frac{1}{T}\sum_i E_i\, dp_i = \frac{\delta q_{\mathrm{rev}}}{T}.
```

This is the definition of entropy from [Section 4.1](section-01.md), now obtained from the microscopic side. In a reversible process the entropy changes only through the heat term. Reversible work alone ($\sum_i E_i\,dp_i = 0$) leaves the entropy unchanged, in agreement with the macroscopic result that $dS = 0$ on a reversible adiabat.

The restriction to reversible processes matters. In the free expansion of Example 2 in [Section 4.1](section-01.md), no heat flows and no work is done, yet the distribution over microstates changes and the entropy of the gas increases.

---

## Worked Examples

### Example 1: Entropy of a Two-State Distribution

**Problem.** A system has two microstates with probabilities $p$ and $1-p$. Compute the Gibbs entropy and find where it is maximized.

**Solution.** The Gibbs entropy is

```{math}
S(p) = -k_{\mathrm{B}}\left[p\ln p + (1-p)\ln(1-p)\right].
```

The derivative $dS/dp = -k_{\mathrm{B}}\left[\ln p - \ln(1-p)\right]$ vanishes at $p = 1/2$, where

```{math}
S = -k_{\mathrm{B}}\left[\tfrac{1}{2}\ln\tfrac{1}{2} + \tfrac{1}{2}\ln\tfrac{1}{2}\right]
= -k_{\mathrm{B}}\ln\tfrac{1}{2}
= k_{\mathrm{B}}\ln 2.
```

To confirm this is a maximum, note that $S(p)$ is concave (its second derivative $d^2S/dp^2 = -k_{\mathrm{B}}/[p(1-p)] < 0$ for $0 < p < 1$) and is zero at the endpoints $p = 0$ and $p = 1$.

**Result.** The maximum entropy $k_{\mathrm{B}}\ln 2$ occurs at $p = 1/2$, where the distribution is most spread out. This is consistent with Boltzmann's formula: two equally likely microstates give $\Omega = 2$ and $S = k_{\mathrm{B}}\ln 2$.

---

### Example 2: Entropy of the Two-State System from Section 2.2

**Problem.** In [Section 2.2](../chapter-02/section-02.md), we analyzed a two-state system with energies $E_1 = 0$ and $E_2 = \varepsilon = 0.010\ \mathrm{eV}$ and found $p_1 = 0.596$ and $p_2 = 0.404$ at $T = 300\ \mathrm{K}$. Compute the entropy two ways: (a) from the Gibbs formula, and (b) from the canonical identity $S = U/T + k_{\mathrm{B}}\ln Q$.

**Solution.**

**(a) Gibbs formula.**

```{math}
S = -k_{\mathrm{B}}\left[0.596\ln(0.596) + 0.404\ln(0.404)\right].
```

Evaluating: $0.596\ln(0.596) = -0.3084$ and $0.404\ln(0.404) = -0.3662$, so

```{math}
S = -k_{\mathrm{B}}(-0.3084 - 0.3662) = 0.675\,k_{\mathrm{B}} = 9.32 \times 10^{-24}\ \mathrm{J/K}.
```

**(b) Canonical identity.**

From [Section 2.2](../chapter-02/section-02.md): $Q = 1.679$, and $U = \varepsilon\,p_2 = (0.010)(0.404) = 4.04 \times 10^{-3}\ \mathrm{eV} = 6.47 \times 10^{-22}\ \mathrm{J}$.

```{math}
S = \frac{U}{T} + k_{\mathrm{B}}\ln Q = \frac{6.47 \times 10^{-22}}{300} + (1.381 \times 10^{-23})\ln(1.679).
```

```{math}
S = 2.16 \times 10^{-24} + 7.16 \times 10^{-24} = 9.32 \times 10^{-24}\ \mathrm{J/K}.
```

**Result.** Both methods give the same answer, as they must — the canonical identity is derived from the Gibbs formula. Notice that the entropy ($0.675\,k_{\mathrm{B}}$) is less than the maximum ($k_{\mathrm{B}}\ln 2 = 0.693\,k_{\mathrm{B}}$) because the distribution is not perfectly uniform: the ground state is slightly more populated than the excited state at 300 K. At higher temperatures, the two probabilities would become more equal, and $S$ would approach $k_{\mathrm{B}}\ln 2$ from below.

---

## Concept Checks

1. How does "entropy increases until equilibrium" translate into a statement about probability distributions over microstates?
2. Why does $S = k_{\mathrm{B}}\ln\Omega$ require the microcanonical assumption $p_i = 1/\Omega$? What goes wrong if you try to apply it to a canonical ensemble?
3. In the Gibbs entropy formula, why does $S = 0$ when $p_i = 1$ for one state and $p_j = 0$ for all others?
4. How can entropy be extensive (additive for independent subsystems) while probability distributions are normalized to 1?
5. Starting from $A = -k_{\mathrm{B}}T\ln Q$, how would you obtain $S$ and $U$ as functions of $T$? (*Hint:* think about which partial derivatives of $A$ give $S$ and $U$.)

## Key Takeaways

- Entropy quantifies how many microstates — or how much probability weight — a macrostate contains.
- Boltzmann's $S = k_{\mathrm{B}}\ln\Omega$ applies to isolated systems (microcanonical ensemble, equal probabilities). The logarithm ensures extensivity.
- The Gibbs formula $S = -k_{\mathrm{B}}\sum p_i\ln p_i$ generalizes Boltzmann's result to any probability distribution and reduces to it in the microcanonical limit.
- In the canonical ensemble, $S = U/T + k_{\mathrm{B}}\ln Q$ ties entropy directly to the partition function.
- Rearranging gives the Helmholtz free energy $A = U - TS = -k_{\mathrm{B}}T\ln Q$ — the master potential from which all canonical thermodynamic quantities can be derived.
- For a reversible process, the Gibbs formula gives $dS = \tfrac{1}{T}\sum_i E_i\,dp_i = \delta q_{\mathrm{rev}}/T$, the definition of [Section 4.1](section-01.md): the entropy changes through the heat term of the microscopic First Law ([Section 3.2](../chapter-03/section-02.md)).
