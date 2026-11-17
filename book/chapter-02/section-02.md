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
numbering:
  enumerator: 2.2.%s
---


# 2.2. Canonical Ensemble

See also: [Course-wide Conventions & Notation](../notation.md)

## Overview and Learning Objectives

The canonical ensemble describes a closed system in thermal contact with a large heat bath at fixed temperature. Its central result is the Boltzmann assignment $p_i \propto e^{-\beta E_i}$, normalized by the partition function. This section derives and interprets that probability distribution and introduces a two-state system as a concrete model.

In a **closed system** (energy exchange allowed, matter exchange forbidden) in thermal contact with a heat bath, the appropriate statistical description is the **canonical ensemble** (panel (b) of [](#fig:ensembles)).

Learning objectives:

- Describe the physical meaning of a heat bath and why its temperature is approximately constant during energy exchange.
- Derive the canonical probability $p_i = e^{-\beta E_i}/Q$ and define $\beta$ and $Q$.
- Use the partition function $Q=\sum_i e^{-\beta E_i}$ to normalize probabilities and compute ratios $p_j/p_i$.
- Apply the canonical distribution to a two-state system and interpret the low- and high-temperature limits.

## Core Ideas and Derivations

### Probability of a Microstate in the Canonical Ensemble

Consider an ensemble of $\mathcal{A}$ identical closed systems in thermal contact with a heat bath at temperature $T$. A heat bath is an environment so large that it can absorb or release energy with negligible change in temperature.

````{admonition} How Can an Environment Absorb or Release Energy Without Changing Its Temperature?
:class: dropdown
The key point is that the temperature of a collection of particles is determined by the average kinetic energy per particle. According to the equipartition theorem, for particles with three translational degrees of freedom, the total kinetic energy is

```{math}
E_{\text{kin}} = \frac{3}{2} N\, k_\text{B}\, T.
```

Thus, for a fixed amount of energy transfer, the corresponding change in temperature is

```{math}
\Delta T = \frac{2}{3}\frac{\Delta E_\text{kin}}{N\, k_\text{B}}.
```

Because the environment has many more particles than the system ($N_\text{env} \gg N_\text{sys}$), the same energy transfer $\Delta E$ produces a much smaller $\Delta T$ in the environment. For example, if $N_\text{env} = 1000N_\text{sys}$ (e.g., $N_\text{sys} = 6.022\times10^{23}$ and $N_\text{env} = 6.022\times10^{26}$ particles), then the same $\Delta E$ yields a temperature change that is one thousand times smaller for the environment.

This illustrates why an environment (a heat bath) can absorb or release energy nearly isothermally: its very large number of particles (and hence high heat capacity) makes the fractional change in average kinetic energy per particle negligible, even while energy is exchanged.

Thus, even though the total energy of the universe (system plus environment) remains constant, the environment's temperature barely changes because any energy loss or gain is diluted among a huge number of particles.
````

#### Population ratios depend only on energy differences

Intuition tells us that a system is more likely to be found in lower-energy microstates.

```{admonition} Building This Intuition
:class: dropdown
Consider a raindrop that can exist at two different elevations: one high up at Lake Itasca in northern Minnesota and one low down in the Mississippi River in St. Louis. At the higher elevation, the raindrop has more gravitational potential energy, making that state less favorable. Consequently, the raindrop is more likely to be found at the lower elevation, where its gravitational potential energy is reduced. In essence, systems tend to prefer lower-energy states, which is why, on a macroscopic scale, water naturally flows downhill.
```

Let $a_1$ and $a_2$ denote the numbers of systems in microstates 1 and 2. We assume the population ratio depends only on the energy difference:

```{math}
\frac{a_2}{a_1} = f \left( E_1, E_2 \right) = f \left( E_1 - E_2 \right),
```

where $E_1$ and $E_2$ are the corresponding microstate energies.

```{admonition} Why Does the Ratio Depend Only on the Energy Difference?
:class: dropdown
The energy of a system is defined relative to an arbitrary (often convenient) reference level. For example,

* The kinetic energy of a moving particle is defined relative to a stationary particle.
* The gravitational potential energy of a raindrop is defined relative to the surface of the Earth.

Therefore, only energy differences—not the absolute zero of energy—enter into probability ratios.
```

#### Finding an Acceptable Form for $f$

Because ratios multiply, we can write

```{math}
\begin{align*}
\frac{a_3}{a_1} &= \frac{a_2}{a_1}\,\frac{a_3}{a_2} \\
f \left( E_1 - E_3 \right) &= f \left( E_1 - E_2 \right)\, f \left( E_2 - E_3 \right).
\end{align*}
```

If $f$ is well behaved (e.g., continuous, measurable), it must be of the form

```{math}
f \left( E_m - E_n \right) = e^{\beta \left( E_m - E_n \right)},
```

where $\beta$ is a constant that this argument does not determine. The dropdown "Alternative derivation roadmap" below shows that $\beta = 1/(k_\text{B} T)$, where $T$ is the temperature of the heat bath. Because $k_\text{B} T$ is an energy, $\beta$ has units of inverse energy. Since $T > 0$, $\beta$ is positive, in agreement with the intuition that lower-energy microstates are more probable.

````{admonition} Checking the Form of $f$
:class: dropdown
Using this form for $f$, the functional equation is satisfied:

```{math}
\begin{align*}
f(E_1-E_2)\,f(E_2-E_3)
&= e^{\beta(E_1-E_2)}\,e^{\beta(E_2-E_3)} \\
&= e^{\beta(E_1-E_3)} \\
&= f(E_1-E_3).
\end{align*}
```
````

### Converting $f$ to a Probability

The same argument applies to any pair of microstates $m$ and $n$, so

```{math}
\frac{a_n}{a_m} = f \left( E_m - E_n \right) = \frac{e^{\beta E_m}}{e^{\beta E_n}}.
```

Separating the indices $m$ and $n$,

```{math}
a_m e^{\beta E_m} = a_n e^{\beta E_n} = \mathcal{C},
```

where $\mathcal{C}$ is a constant. Because the equality holds for every pair of microstates, the product $a_m e^{\beta E_m}$ has the same value for all of them. Therefore, the number of systems in a microstate with energy $E_m$ is

```{math}
a_m = \mathcal{C} e^{-\beta E_m}.
```

The constant $\mathcal{C}$ is determined by the normalization condition

```{math}
\sum_n a_n = \mathcal{A} = \mathcal{C} \sum_n e^{-\beta E_n}.
```

Solving for $\mathcal{C}$ and substituting into the expression for $a_m$ gives

```{math}
:label: canonical-ensemble-probability
p_m = \frac{a_m}{\mathcal{A}} = \frac{e^{-\beta E_m}}{\sum_n e^{-\beta E_n}} = \frac{e^{-\beta E_m}}{Q},
```

where $p_m$ is the probability of finding the system in microstate $m$ and

```{math}
:label: canonical-partition-function
Q = \sum_n e^{-\beta E_n}
```

is the **partition function**. The sum runs over all microstates of the system.

The argument above does not show how the heat bath enters the result. The alternative derivation below starts from the isolated combination of the system and the reservoir, and it identifies $\beta$ with the temperature of the heat bath.

````{admonition} Alternative derivation roadmap: where does $p_i \propto e^{-\beta E_i}$ come from?
:class: dropdown

We will use the canonical-ensemble probability
```{math}
p_i \propto e^{-\beta E_i}
\quad\text{with}\quad
\beta = \frac{1}{k_\text{B} T}
```
throughout this module. Before we use it heavily, here are two quick (complementary) ways to see
**why the exponential “Boltzmann factor” appears**.

---

**Setup (common to both routes)**

- **System** $S$: a *closed* system (fixed $N,V$) with discrete microstates $i=1,2,\dots$ and energies $E_i$.
- **Reservoir / heat bath** $R$: very large, can exchange energy with $S$, and stays at (approximately) constant temperature $T$.
- **Universe** $S+R$ is **isolated**, so the total energy is fixed:
  ```{math}
  E_{\text{tot}} = E_i + E_R.
  ```
- **Weak coupling**: the interaction energy is negligible, so energies add.

The question is: *What is the probability $p_i$ that the system is in microstate $i$?*

---

**Route A: “Counting” argument (microcanonical $\rightarrow$ canonical)**

**Idea:** the system is more likely to be in a microstate $i$ if the reservoir has *many* compatible microstates
when the system has energy $E_i$.

1. **Fundamental postulate for an isolated system (the universe):**  
   all accessible microstates of $S+R$ at fixed $E_{\text{tot}}$ are equally likely.

2. **If $S$ is in microstate $i$, then $R$ must have energy**
   ```{math}
   E_R = E_{\text{tot}} - E_i.
   ```
   The number of compatible universe microstates is therefore the number of reservoir microstates at that energy:
   ```{math}
   \Omega_R(E_{\text{tot}}-E_i).
   ```
   So
   ```{math}
   p_i \propto \Omega_R(E_{\text{tot}}-E_i).
   ```

3. **Convert multiplicity to entropy** using Boltzmann’s definition:
   ```{math}
   S_R(E) = k_\text{B} \ln \Omega_R(E).
   ```
   Taking logs,
   ```{math}
   \ln p_i = \text{const} + \ln \Omega_R(E_{\text{tot}}-E_i)
   = \text{const} + \frac{S_R(E_{\text{tot}}-E_i)}{k_\text{B}}.
   ```

4. **Reservoir is huge $\Rightarrow$ expand $S_R$ in a Taylor series** around $E_{\text{tot}}$.  
   Because $E_i$ is “small” compared to the reservoir’s energy scale,
   ```{math}
   S_R(E_{\text{tot}}-E_i)
   \approx S_R(E_{\text{tot}})
   - E_i\left(\frac{\partial S_R}{\partial E}\right)_{E_{\text{tot}}}
   + \frac{E_i^2}{2}\left(\frac{\partial^2 S_R}{\partial E^2}\right)_{E_{\text{tot}}}
   + \cdots
   ```

5. **Use the thermodynamic definition of temperature (for the reservoir):**
   ```{math}
   \left(\frac{\partial S}{\partial E}\right)_{V,N} = \frac{1}{T}.
   ```
   So the *linear* term becomes $-E_i/T$.

   The *quadratic* term is typically negligible because the reservoir’s heat capacity is enormous:
   ```{math}
   \left(\frac{\partial^2 S}{\partial E^2}\right)_{V,N}
   =
   -\frac{1}{T^2 C_V},
   ```
   so the correction looks like $-E_i^2/(2T^2 C_V)$, which is tiny when $C_V$ is very large.

6. **Keep only the dominant (linear) term**:
   ```{math}
   S_R(E_{\text{tot}}-E_i) \approx S_R(E_{\text{tot}}) - \frac{E_i}{T}.
   ```

7. **Exponentiate to get the Boltzmann factor**:
   ```{math}
   p_i \propto \exp\!\left(\frac{S_R(E_{\text{tot}})}{k_\text{B}}\right)\,
          \exp\!\left(-\frac{E_i}{k_\text{B} T}\right).
   ```
   The first exponential is just a constant (independent of $i$), so the *physics* is:
   ```{math}
   \boxed{p_i \propto e^{-E_i/(k_\text{B} T)} \equiv e^{-\beta E_i}}
   \quad\text{with}\quad
   \beta \equiv \frac{1}{k_\text{B} T}.
   ```

8. **Normalize** to turn “$\propto$” into “$=$”. Define the **partition function**
   ```{math}
   Q \equiv \sum_i e^{-\beta E_i},
   ```
   so
   ```{math}
   \boxed{p_i = \frac{e^{-\beta E_i}}{Q}}.
   ```

---

**Route B: “Maximum entropy” argument (Gibbs entropy + constraints)**

**Idea:** in equilibrium with a heat bath, the system adopts the probability distribution that
maximizes entropy subject to what is fixed.

1. Start from the Gibbs/Shannon entropy functional for a probability distribution $\{p_i\}$:
   ```{math}
   S[\{p_i\}] = -k_\text{B} \sum_i p_i \ln p_i.
   ```

2. Impose the two canonical constraints:
   - normalization: $\sum_i p_i = 1$
   - fixed mean energy: $\sum_i p_i E_i = U$

3. Maximize $S$ with Lagrange multipliers $\alpha$ and $\lambda$:
   ```{math}
   \delta\Big[-k_\text{B}\sum_i p_i\ln p_i
   -\alpha\Big(\sum_i p_i-1\Big)
   -\lambda\Big(\sum_i p_iE_i-U\Big)\Big]=0.
   ```

4. The maximizer satisfies
   ```{math}
   p_i = \text{const}\times e^{-(\lambda/k_\text{B})E_i}.
   ```
   Rename $\beta \equiv \lambda/k_\text{B}$:
   ```{math}
   p_i \propto e^{-\beta E_i}.
   ```

5. Normalize exactly as before:
   ```{math}
   p_i=\frac{e^{-\beta E_i}}{Q},
   \qquad
   Q=\sum_i e^{-\beta E_i}.
   ```

6. Identify $\beta$ with temperature by using the thermodynamic identity
   ```{math}
   \left(\frac{\partial S}{\partial U}\right)_{V,N}=\frac{1}{T}.
   ```
   For the canonical distribution one can show $dS/dU = k_\text{B}\beta$, so
   ```{math}
   \boxed{\beta=\frac{1}{k_\text{B} T}}.
   ```

---

**Quick checks / intuition (why this form makes sense)**

- **Energy differences matter, not absolute zero**:
  ```{math}
  \frac{p_j}{p_i}=e^{-\beta(E_j-E_i)}.
  ```
  Shifting all energies by a constant leaves all probability ratios unchanged.

- **Low $T$ (large $\beta$)**: weights concentrate on the lowest-energy microstates.  
  **High $T$ (small $\beta$)**: many microstates become appreciably populated.

- **Degeneracy is automatic** if we sum over microstates:  
  if many microstates share the same energy, that energy level is more probable because it appears many times in the sum.

---

**Bottom line:** In the canonical ensemble, the exponential appears because the reservoir’s multiplicity grows
exponentially with its entropy, and the entropy changes approximately linearly with energy over the small energy
exchanges relevant to a large heat bath.
````

## Two-State System

Consider a system with two microstates (two “levels”): state 1 with energy $E_1$ and state 2 with energy $E_2$. We label the states so that $E_2 > E_1$, which makes state 1 the lower-energy state.

```{list-table} Chemical Contexts Where a Two-State Approximation Might be Appropriate
:header-rows: 1
:name: two-level-systems

* - Chemical Context
  - State 1
  - State 2
* - Electronic transitions in atoms or molecules
  - Ground state
  - Excited state
* - Donor–acceptor electron transfer
  - Reduced state
  - Oxidized state
* - Molecular isomerization
  - Reactant
  - Product
* - Defects in solids
  - Defect-free
  - Defective
* - Protein folding
  - Folded
  - Unfolded
```

### Partition Function for a Two-State System

The partition function $Q_\text{two-state}$ for a two-state system is

```{math}
\begin{align*}
Q_\text{two-state} &= e^{-\beta E_1} + e^{-\beta E_2} \\
&= e^{-\beta E_1} \left( 1 + e^{-\beta \Delta E} \right),
\end{align*}
```

where $\Delta E = E_2 - E_1$ is the energy difference between the two states.

### Probability of Finding the System in State 1

The probability of finding the system in state 1 is

```{math}
p_1 = \frac{e^{-\beta E_1}}{Q_\text{two-state}} = \frac{1}{1 + e^{-\beta \Delta E}} = \frac{1}{1 + e^{-\frac{\Delta E}{k_\text{B} T}}}.
```

### Probability of Finding the System in State 2

The probability of finding the system in state 2 is

```{math}
p_2 = 1 - p_1.
```

### Partition Function as the Effective Number of Thermally Accessible Microstates

When energies are measured from the lowest state ($E_1 = 0$), every Boltzmann factor $e^{-\beta E_m}$ is at most 1, and the factor for the lowest state equals 1. The partition function therefore lies between 1 and the total number of microstates, and it can be read as the effective number of microstates that are thermally accessible at temperature $T$. [](#fig:two-state-probabilities) shows $Q_\text{two-state}$ and the two probabilities as functions of temperature for $E_1 = 0$ and $\Delta E = 0.01\ \mathrm{eV}$. The dashed lines in panel (a) mark the low- and high-temperature limits of $Q_\text{two-state}$.

```{code-cell} ipython3
:label: fig:two-state-probabilities
:caption: Two-state system with $E_1 = 0$ and $\Delta E = 0.01\ \mathrm{eV}$ as a function of temperature. (a) The partition function $Q_\text{two-state}$ rises from 1, where only state 1 is thermally accessible, toward 2, where both states are. (b) Probabilities $p_1$ and $p_2$ of finding the system in state 1 and in state 2.
:tags: [hide-input]

import sys

import matplotlib.pyplot as plt
import numpy as np
from scipy.constants import eV, k

sys.path.append("..")            # the book folder, which holds figstyle
import figstyle as fs

k_B = k / eV                     # Boltzmann constant in eV/K
E1, E2 = 0.0, 0.01               # energies of the two states in eV

T = np.linspace(1, 1000, 1000)   # K
beta = 1 / (k_B * T)             # 1/eV
Q = np.exp(-beta * E1) + np.exp(-beta * E2)
p1 = np.exp(-beta * E1) / Q
p2 = np.exp(-beta * E2) / Q


def draw(c):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(fs.WIDTH, 3.0), constrained_layout=True)

    for limit, label, va, dy in [(1, "one accessible microstate", "top", -0.035),
                                 (2, "two accessible microstates", "bottom", 0.035)]:
        ax1.axhline(limit, color=c.muted, lw=1.2, ls=(0, (4, 3)))
        ax1.text(1000, limit + dy, label, color=c.muted, ha="right", va=va)
    ax1.plot(T, Q, color=c.ink)
    ax1.set(xlabel="Temperature $T$ (K)", ylabel=r"$Q_{\mathrm{two\text{-}state}}$", xlim=(0, 1000), ylim=(0.8, 2.2),
            yticks=[1.0, 1.5, 2.0])
    fs.panel(ax1, "a", "Partition function")

    ax2.plot(T, p1, color=c.blue)
    ax2.plot(T, p2, color=c.vermilion)
    ax2.text(330, 0.69, "state 1, $p_1$", color=c.blue, va="center")
    ax2.text(330, 0.31, "state 2, $p_2$", color=c.vermilion, va="center")
    ax2.set(xlabel="Temperature $T$ (K)", ylabel="Probability", xlim=(0, 1000), ylim=(0, 1), yticks=[0, 0.5, 1])
    fs.panel(ax2, "b", "Probabilities")
    return fig


fs.show(draw, alt="Two panels against temperature from 0 to 1000 K for a two-state system with an energy gap of 0.01 "
                  "eV. Panel a: the partition function starts at 1 at low temperature, rises steeply up to about "
                  "200 K, and then approaches 2 slowly, reaching about 1.9 at 1000 K. Panel b: the probability of "
                  "state 1 starts at 1 and falls, the probability of state 2 starts at 0 and rises, and the two "
                  "curves approach each other at high temperature.")
```

## Computational Studio: Two-State System

In [](#fig:two-state-studio), you can change the energy difference between the two states and follow the probabilities, the heat capacity, and the partition function as the temperature rises.

```{iframe} https://wexlergroup.github.io/chem-4020-5020/computational-studios/two-level-system.html
:label: fig:two-state-studio
:width: 100%
:title: Two-state system
:placeholder: ../_static/two-level-system/two-level-system.png

Two-state system with $E_1 = 0$ and $E_2 = \Delta E$ as a function of temperature. Panel (a) shows the probabilities $p_1$ and $p_2$ of the two states, panel (b) the heat capacity of one particle, $C_V/k_\text{B}$, and panel (c) the partition function $Q$. Energies and $k_\text{B}T$ are in units of an arbitrary energy $\varepsilon$.
```

You can also open the studio in a new tab: [Two-State System Studio](https://wexlergroup.github.io/chem-4020-5020/computational-studios/two-level-system.html).

## Worked Example

### Two-state system probabilities

A single particle has two states with energies $E_1=0$ and $E_2=\varepsilon$, where $\varepsilon=0.010\ \mathrm{eV}$. Find $p_1$ and $p_2$ at $T=300\ \mathrm{K}$.

Use $p_i = e^{-\beta E_i}/Q$, $Q=1+e^{-\beta\varepsilon}$, and $k_\text{B}=8.617\times 10^{-5}\ \mathrm{eV/K}$.

1. **Compute $\beta\varepsilon$**

   ```{math}
   \beta\varepsilon=\frac{\varepsilon}{k_\text{B}T}
   =\frac{0.010\ \mathrm{eV}}{(8.617\times10^{-5}\ \mathrm{eV/K})(300\ \mathrm{K})}=0.387.
   ```

2. **Compute the Boltzmann factor**

   ```{math}
   e^{-\beta\varepsilon}=e^{-0.387}=0.679.
   ```

3. **Normalize**

   ```{math}
   Q=1+0.679=1.679,
   \qquad
   p_1=\frac{1}{Q}=0.596,
   \qquad
   p_2=\frac{0.679}{Q}=0.404.
   ```

**Result.** At $300\ \mathrm{K}$, about $40\%$ of systems occupy the higher-energy state (state 2) for this energy gap.

## Concept Checks

1. Why must the argument of the exponential in $e^{-\beta E}$ be dimensionless?
2. What happens to $p_2/p_1$ as $T\to 0$? As $T\to\infty$?
3. Which part of the canonical derivation relies on the reservoir being much larger than the system?
4. How does shifting all energies by a constant affect $p_i$, and why?

## Key Takeaways

- In the canonical ensemble, microstate probabilities follow $p_i=e^{-\beta E_i}/Q$ with $\beta=1/(k_\text{B}T)$.
- The partition function $Q$ normalizes probabilities and encodes thermodynamic information.
- Probability ratios depend only on **energy differences**, not on the absolute zero of energy.
- Two-state systems vividly illustrate how temperature controls population of higher-energy states.
