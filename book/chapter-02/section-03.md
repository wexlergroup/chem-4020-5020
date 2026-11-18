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
  enumerator: 2.3.%s
---


# 2.3. Ensemble Averages

See also: [Course-wide Conventions & Notation](../notation.md)

## Overview and Learning Objectives

Once microstate probabilities are known, thermodynamic properties follow as ensemble averages. This section derives canonical expressions for internal energy, heat capacity, and pressure in terms of derivatives of $\ln Q$. A recurring theme is that fluctuations are not a nuisance—they encode response functions such as $C_V$. A **response function** measures how strongly a property of the system changes when a controlled variable changes; $C_V$ is the response of $U$ to a change in $T$.

[Section 2.1](section-01.md) emphasized that *macroscopic properties are expectation values of microscopic properties*, with the expectation taken over the microstate probabilities. [Section 2.2](section-02.md) derived the canonical probability distribution for a closed system’s microstates. Here, we connect these results by computing **ensemble averages** that determine macroscopic properties.

Learning objectives:

- Compute internal energy both as an ensemble average, $U=\langle E\rangle$, and as a derivative, $U=-(\partial\ln Q/\partial\beta)_{N,V}$.
- Relate heat capacity to energy fluctuations via $\sigma_E^2 = k_{\mathrm B}T^2 C_V$.
- Compute pressure from the partition function using $P=k_{\mathrm B}T(\partial\ln Q/\partial V)_{N,T}$ (when applicable).
- Apply these formulas to simple discrete models (e.g., a two-state system).

## Core Ideas and Derivations

### Thermodynamic Equilibrium

Recall from [Section 1.1](../chapter-01/section-01.md) that **thermodynamic equilibrium** is a state of simultaneous mechanical, thermal, and chemical equilibrium. Below are working definitions of each type of equilibrium:

```{glossary}
Thermal contact
: A state in which two systems can exchange energy.

Chemical contact
: A state in which two systems can exchange matter.

Mechanical equilibrium
: A state in which there are no unbalanced forces. Across a movable boundary, the pressure of the system equals that of its surroundings.

Thermal equilibrium
: A state in which there is no net exchange of energy between systems in thermal contact.

Chemical equilibrium
: A state in which there is no net exchange of matter between systems in chemical contact.
```

### Internal Energy

The **internal energy** $U$ of a macroscopic system at thermodynamic equilibrium is defined as the ensemble average $\langle E \rangle$ of the microstate energy $E$:

```{math}
U \;=\; \langle E \rangle \;=\; \sum_{i=1}^M E_i\,p_i,
```

where $E_i$ is the energy of microstate $i$, $p_i$ is the probability of that microstate, and $M$ is the total number of microstates.

For the **canonical ensemble**, using the Boltzmann factor and the partition function $Q$, we have

```{math}
:label: eq:ensemble-average-internal-energy-sum
U \;=\; \frac{1}{Q} \sum_{i=1}^M E_i\,e^{-\beta E_i},
```

where $\beta = 1/(k_{\mathrm B} T)$ and

```{math}
Q \;=\; \sum_{i=1}^M e^{-\beta E_i}.
```

```{admonition} Why can we factor out $Q$?
:class: dropdown
At fixed $N$, $V$, and $T$ (equivalently fixed $\beta$), $Q$ does not depend on the microstate index $i$. Therefore, it is constant with respect to the summation over $i$ and can be pulled outside the sum.
```

#### Partial Derivative of the Partition Function

Taking the partial derivative of $Q$ with respect to $\beta$ gives

```{math}
\left(\frac{\partial Q}{\partial \beta}\right)_{N,V}
=\; \sum_{i=1}^M \frac{\partial}{\partial \beta}\left(e^{-\beta E_i}\right)
=\; -\sum_{i=1}^M E_i\,e^{-\beta E_i}.
```

Substituting this result into the expression for $U$ yields

```{math}
U 
=\; -\frac{1}{Q}\,\left(\frac{\partial Q}{\partial \beta}\right)_{N,V}
=\; -\left(\frac{\partial \ln Q}{\partial \beta}\right)_{N,V}.
```

````{admonition} Where does the natural log come in?
:class: dropdown
Using the chain rule,

```{math}
\frac{\partial \ln Q}{\partial \beta}
=\; \frac{1}{Q} \left(\frac{\partial Q}{\partial \beta}\right).
```
````

#### Example: Two-State System

For a two-state system with energies $E_1$ and $E_2$, the partition function is

```{math}
Q_{\text{two-state}} = e^{-\beta E_1} + e^{-\beta E_2}.
```

Then the internal energy is

```{math}
U_{\text{two-state}}
=\; -\frac{\partial}{\partial \beta}\,\ln\left(e^{-\beta E_1} + e^{-\beta E_2}\right)
=\; \frac{E_1\, e^{-\beta E_1} + E_2\, e^{-\beta E_2}}{e^{-\beta E_1} + e^{-\beta E_2}},
```

which is the ensemble average $\langle E\rangle$.

````{admonition} Consistency Check
:class: dropdown
Using $Q_{\text{two-state}} = e^{-\beta E_1} + e^{-\beta E_2}$,

```{math}
U_{\text{two-state}}
=\; \frac{E_1 e^{-\beta E_1} + E_2 e^{-\beta E_2}}{Q_{\text{two-state}}}
=\; E_1\,p_1 + E_2\,p_2
=\; \langle E\rangle.
```
````

[](#fig:two-state-energy) shows $U_{\text{two-state}}$, measured from $E_1$, as a function of temperature for an energy gap of $0.01\ \mathrm{eV}$, together with the heat capacity at constant volume, which the next subsection defines.

```{code-cell} ipython3
:label: fig:two-state-energy
:caption: Two-state system with an energy gap $E_2 - E_1 = 0.01\ \mathrm{eV}$ ($10\ \mathrm{meV}$) as a function of temperature. (a) Internal energy measured from $E_1$. The dashed line is its high-temperature limit, $(E_2 - E_1)/2$. (b) Heat capacity at constant volume in units of $k_{\mathrm{B}}$.
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

# Internal energy U = -d(ln Q)/d(beta), which is the average of the two energies weighted by their probabilities
Q = np.exp(-beta * E1) + np.exp(-beta * E2)
U = (E1 * np.exp(-beta * E1) + E2 * np.exp(-beta * E2)) / Q

# Heat capacity in units of k_B: C_V / k_B = x^2 e^(-x) / (1 + e^(-x))^2, with x = (E2 - E1) / (k_B T)
x = (E2 - E1) * beta
C_V = x**2 * np.exp(-x) / (1 + np.exp(-x))**2


def draw(c):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(fs.WIDTH, 3.0), constrained_layout=True)

    half_gap = 1e3 * (E2 - E1) / 2           # meV
    ax1.axhline(half_gap, color=c.muted, lw=1.2, ls=(0, (4, 3)))
    ax1.text(1000, half_gap + 0.15, "$(E_2 - E_1)/2$", color=c.muted, ha="right", va="bottom")
    ax1.plot(T, 1e3 * (U - E1), color=c.ink)
    ax1.set(xlabel="Temperature $T$ (K)", ylabel=r"$U_{\mathrm{two\text{-}state}} - E_1$ (meV)", xlim=(0, 1000),
            ylim=(0, 6))
    fs.panel(ax1, "a", "Internal energy")

    ax2.plot(T, C_V, color=c.ink)
    ax2.set(xlabel="Temperature $T$ (K)", ylabel=r"$C_V / k_{\mathrm{B}}$", xlim=(0, 1000), ylim=(0, 0.5))
    fs.panel(ax2, "b", "Heat capacity")
    return fig


fs.show(draw, alt="Two panels against temperature from 0 to 1000 K for a two-state system with an energy gap of 0.01 "
                  "eV. Panel a: the internal energy measured from the lower energy starts at 0, rises steeply up to "
                  "about 200 K, and then approaches a dashed line at half the gap, 5 meV, reaching about 4.7 meV at "
                  "1000 K. Panel b: the heat capacity divided by the Boltzmann constant starts at 0, rises to a "
                  "maximum of about 0.44 near 50 K, and then decreases toward 0 at high temperature.")
```

### Heat Capacity at Constant Volume

From [Section 1.1](../chapter-01/section-01.md), **heat** is energy transferred due to a temperature difference. The **heat capacity** at constant volume, $C_V$, measures how much heat is required to change a system’s temperature at fixed $N$ and $V$:

```{math}
C_V 
\;=\; \left(\frac{\partial U}{\partial T}\right)_{N,V}.
```

Panel (b) of [](#fig:two-state-energy) shows $C_V$ for the two-state system.

#### Heat Capacity and Energy Fluctuations

In the canonical ensemble, the energy of the system is not fixed, because the system exchanges energy with the heat bath ([Section 2.2](section-02.md)); the internal energy $U = \langle E\rangle$ is its average. In statistics, the variance $\sigma_X^2$ of a random variable $X$ is $\sigma_X^2 = \langle (X - \langle X\rangle)^2\rangle$. In statistical mechanics, $\sigma_E^2$ quantifies fluctuations in the total energy:

```{math}
\sigma_E^2 
=\; \langle E^2\rangle - \langle E\rangle^2,
```

where $\langle E^2\rangle = \sum_{i=1}^M E_i^2\,p_i$. One can show that

```{math}
\sigma_E^2 
=\; \left(\frac{\partial^2 \ln Q}{\partial \beta^2}\right)_{N,V}
=\; k_{\mathrm B}\,T^2\,C_V,
```

which implies

```{math}
C_V 
=\; \frac{\sigma_E^2}{k_{\mathrm B}\,T^2}.
```

```{admonition} What does the relationship mean physically?
:class: dropdown
Consider two systems in contact with heat baths at the same temperature $T$, one with a larger $C_V$ than the other. Use the relationship above to compare their energy fluctuations $\sigma_E^2$.
```

````{admonition} Complete Derivation of the $C_V$ and $\sigma_E^2$ Relationship
:class: dropdown
**1. Show that $\left( \partial^2 \ln Q / \partial \beta^2 \right)_{N, V} = \sigma_E^2$**

```{math}
\begin{align*}
\left( \frac{\partial^2 \ln Q}{\partial \beta^2} \right)_{N, V}
&= \left( \frac{\partial}{\partial \beta} \frac{\partial \ln Q}{\partial \beta} \right)_{N, V} \\
&= -\left( \frac{\partial \langle E \rangle}{\partial \beta} \right)_{N, V} \\
&= -\left( \frac{\partial}{\partial \beta} \sum_{i = 1}^M E_i p_i \right)_{N, V} \\
&= -\left( \frac{\partial}{\partial \beta} \frac{1}{Q} \sum_{i = 1}^M E_i e^{-\beta E_i} \right)_{N, V} \\
&= -\left[ \frac{1}{Q} \sum_{i = 1}^M \left( -E_i^2 e^{-\beta E_i} \right)
+ \left( -\frac{1}{Q^2} \right) \left( \frac{\partial Q}{\partial \beta}\right)_{N, V} \sum_{i = 1}^M E_i e^{-\beta E_i} \right] \\
&= \frac{1}{Q} \sum_{i = 1}^M E_i^2 e^{-\beta E_i}
+ \left[ \frac{1}{Q} \left( \frac{\partial Q}{\partial \beta} \right)_{N, V} \right] \frac{1}{Q} \sum_{i = 1}^M E_i e^{-\beta E_i} \\
&= \langle E^2 \rangle + \left[ -\langle E \rangle \right] \langle E \rangle \\
&= \langle E^2 \rangle - \langle E \rangle^2 \\
&= \sigma_E^2
\end{align*}
```

**2. Show that $\partial T / \partial \beta = -k_{\mathrm B} T^2$**

```{math}
\begin{align*}
\beta &= \frac{1}{k_{\mathrm B} T} \\
T &= \frac{1}{k_{\mathrm B} \beta} \\
\frac{\partial T}{\partial \beta} &= -\frac{1}{k_{\mathrm B} \beta^2} \\
&= -\frac{1}{k_{\mathrm B} \left( \frac{1}{k_{\mathrm B} T} \right)^2} \\
&= -k_{\mathrm B} T^2
\end{align*}
```

**3. Show that $\left( \partial^2 \ln Q / \partial \beta^2 \right)_{N, V} = k_{\mathrm B} T^2 C_V$**

```{math}
\begin{align*}
\left( \frac{\partial^2 \ln Q}{\partial \beta^2} \right)_{N, V}
&= -\left( \frac{\partial \langle E \rangle}{\partial \beta} \right)_{N, V} \\
&= -\left( \frac{\partial \langle E \rangle}{\partial T} \frac{\partial T}{\partial \beta} \right)_{N, V} \\
&= k_{\mathrm B} T^2 \left( \frac{\partial \langle E \rangle}{\partial T} \right)_{N, V} \\
&= k_{\mathrm B} T^2 C_V
\end{align*}
```

**4. Combine Steps 1 and 3 and solve for $C_V$**

```{math}
\begin{align*}
\sigma_E^2 &= k_{\mathrm B} T^2 C_V \\
C_V &= \frac{\sigma_E^2}{k_{\mathrm B} T^2}
\end{align*}
```
````

#### Fluctuations and System Size

How large are these fluctuations? For independent particles, both $U$, measured from the lowest energy of the system, and $C_V$ are proportional to the number of particles $N$. The standard deviation $\sigma_E = \sqrt{k_{\mathrm B} T^2 C_V}$ is then proportional to $\sqrt{N}$, so the width of the energy distribution relative to its mean falls as

```{math}
:label: eq:relative-energy-fluctuation
\frac{\sigma_E}{U} \;\propto\; \frac{1}{\sqrt{N}}.
```

[](#fig:energy-fluctuations) shows this for $N$ independent two-state particles with $E_1 = 0$ and $E_2 = \varepsilon$. Each particle is in state 2 with probability $p_2$, so the number of particles in state 2 follows the binomial distribution, and

```{math}
U = N \varepsilon\, p_2,
\qquad
\sigma_E^2 = N \varepsilon^2\, p_1 p_2,
\qquad
\frac{\sigma_E}{U} = \sqrt{\frac{p_1}{N p_2}}.
```

The energy of 10 particles is spread over a wide range around $U$. For one mole of particles, $\sigma_E/U$ is about $10^{-12}$. The energy of a macroscopic system in contact with a heat bath never stops fluctuating, but the fluctuations are far too small to measure, and this is why thermodynamics can speak of *the* internal energy of the system.

```{code-cell} ipython3
:label: fig:energy-fluctuations
:caption: Energy of $N$ independent two-state particles with $E_1 = 0$ and $E_2 = \varepsilon = 0.01\ \mathrm{eV}$ at $300\ \mathrm{K}$, the conditions of the worked example below. (a) Distribution of the energy per particle, $E/N$, for three values of $N$. Each distribution has an area of 1, and the dashed line marks $U/N = 4.04\ \mathrm{meV}$. (b) Relative width $\sigma_E/U$ of the distribution against $N$ on logarithmic axes. The filled points are the three systems of panel (a), and the open point is one mole of particles.
:tags: [hide-input]

from scipy.constants import N_A
from scipy.stats import binom

epsilon = 0.01                   # energy of state 2 in eV; state 1 has the energy 0
T_room = 300                     # K
x_room = epsilon / (k_B * T_room)
p2 = np.exp(-x_room) / (1 + np.exp(-x_room))     # probability of state 2
p1 = 1 - p2

sizes = [10, 100, 1000]          # numbers of particles in panel (a)


def relative_width(N):
    """sigma_E / U for N independent two-state particles."""
    return np.sqrt(p1 / (N * p2))


def draw(c):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(fs.WIDTH, 3.1), constrained_layout=True,
                                   gridspec_kw={"width_ratios": [1.25, 1]})
    colors = [c.blue, c.vermilion, c.green]

    # (a) The energy is E = n * epsilon, where n, the number of particles in state 2, is binomial.
    for N, color in zip(sizes, colors):
        n = np.arange(N + 1)
        step = 1e3 * epsilon / N                               # spacing of the values of E/N in meV
        density = binom.pmf(n, N, p2) / step                   # probability per meV of E/N
        edges = (np.arange(N + 2) - 0.5) * step
        ax1.stairs(density, edges, color=color, lw=1.6, fill=False)
        ax1.stairs(density, edges, color=c.tint(color, 0.8), fill=True, zorder=0)
    ax1.axvline(1e3 * epsilon * p2, color=c.muted, lw=1.2, ls=(0, (4, 3)), zorder=1)
    ax1.text(1e3 * epsilon * p2 + 0.15, 2.72, "$U/N$", color=c.muted, ha="left", va="top")
    ax1.text(7.3, 0.36, "$N = 10$", color=colors[0], ha="left", va="center")
    ax1.text(5.35, 0.95, "$N = 100$", color=colors[1], ha="left", va="center")
    ax1.text(4.75, 2.2, "$N = 1000$", color=colors[2], ha="left", va="center")
    ax1.set(xlabel="Energy per particle $E/N$ (meV)", ylabel="Probability per meV of $E/N$", xlim=(0, 10),
            ylim=(0, 2.8))
    fs.panel(ax1, "a", "Distribution of the energy")

    # (b) sigma_E / U = sqrt(p1 / (N p2)), a straight line of slope -1/2 on logarithmic axes
    N_line = np.logspace(0, 24, 200)
    ax2.plot(N_line, relative_width(N_line), color=c.ink, lw=1.6)
    for N, color in zip(sizes, colors):
        ax2.plot(N, relative_width(N), "o", color=color, ms=6.5, zorder=3)
    ax2.plot(N_A, relative_width(N_A), "o", mfc=c.ground, mec=c.ink, mew=1.6, ms=6.5, zorder=3)
    ax2.annotate("1 mol", (N_A, relative_width(N_A)), xytext=(-9, -4), textcoords="offset points", ha="right",
                 va="top", color=c.ink)
    ax2.text(1e9, 6e-4, r"$\sigma_E/U \propto 1/\sqrt{N}$", color=c.ink, ha="left", va="bottom")
    ax2.set(xscale="log", yscale="log", xlabel="Number of particles $N$", ylabel=r"Relative width $\sigma_E/U$",
            xlim=(1, 1e25), ylim=(1e-13, 10))
    ax2.set_xticks([1, 1e6, 1e12, 1e18, 1e24])
    ax2.set_yticks([1, 1e-3, 1e-6, 1e-9, 1e-12])
    ax2.minorticks_off()
    fs.panel(ax2, "b", "Relative width")
    return fig


fs.show(draw, alt="Two panels. Panel a: three distributions of the energy per particle between 0 and 10 meV, all "
                  "centered on a dashed line at 4.04 meV. For 10 particles the distribution is a low staircase that "
                  "covers almost the whole range. For 100 particles it is a narrower peak, and for 1000 particles a "
                  "tall, sharp peak. Panel b: on logarithmic axes, the relative width of the distribution falls "
                  "along a straight line from about 1.2 at one particle to about 0.38, 0.12, and 0.038 at 10, 100, "
                  "and 1000 particles, and to about 1.6 times ten to the minus twelve at one mole.")
```

### Pressure

The pressure follows from two results that are derived in later sections:

1. $A = -k_{\mathrm B}\,T\ln Q$, where $A$ is the Helmholtz free energy ([Section 4.3](../chapter-04/section-03.md)).
2. $P = -\left(\frac{\partial A}{\partial V}\right)_{N,T}$, where $P$ is pressure ([Section 5.1](../chapter-05/section-01.md)).

Combining these expressions gives

```{math}
P
=\; k_{\mathrm B}\,T \left(\frac{\partial \ln Q}{\partial V}\right)_{N,T}.
```

[Section 2.5](section-05.md) applies this formula to $N$ particles in a cube and recovers the ideal gas law.

## Worked Example

### Two-state internal energy and heat capacity

Consider a two-state system with $E_1=0$ and $E_2=\varepsilon$. Let $x\equiv \beta\varepsilon$.

1. **Partition function**

   ```{math}
   Q=1+e^{-x}.
   ```

2. **Internal energy**

   ```{math}
   U=\frac{0\cdot 1+\varepsilon e^{-x}}{1+e^{-x}}
   =\varepsilon\frac{e^{-x}}{1+e^{-x}}.
   ```

3. **Heat capacity at constant volume**

   Differentiate $U$ with respect to $T$ (holding $N$ and $V$ fixed). Using $x=\varepsilon/(k_{\mathrm B}T)$ and $dx/dT=-x/T$,

   ```{math}
   C_V=\left(\frac{\partial U}{\partial T}\right)_{N,V}
   = k_{\mathrm B}\,x^2\,\frac{e^{-x}}{(1+e^{-x})^2}.
   ```

4. **Numerical illustration**

   For $\varepsilon=0.010\ \mathrm{eV}$ at $T=300\ \mathrm{K}$, $x=0.3868$ and $e^{-x}=0.6792$, so

   ```{math}
   U = 0.010\ \mathrm{eV}\times \frac{0.6792}{1.6792}=4.04\times 10^{-3}\ \mathrm{eV},
   \qquad
   \frac{C_V}{k_{\mathrm B}}=(0.3868)^2\,\frac{0.6792}{1.6792^2}=3.60\times 10^{-2}.
   ```

**Result.** $U$ interpolates between 0 and $\varepsilon/2$, while $C_V$ peaks when thermal energy is comparable to $\varepsilon$.

## Concept Checks

1. Why is $U$ obtained from a derivative of $\ln Q$ rather than from $Q$ directly?
2. What does it mean physically when $C_V$ is large? When it is small?
3. In the two-state model, why does $C_V\to 0$ both as $T\to 0$ and as $T\to\infty$?
4. How do fluctuations in energy relate to the system’s ability to exchange heat with a reservoir?

## Key Takeaways

- Thermodynamic properties are ensemble averages computed from microstate probabilities.
- In the canonical ensemble, $U=-(\partial\ln Q/\partial\beta)_{N,V}$.
- Response functions connect to fluctuations: $\sigma_E^2=k_{\mathrm B}T^2C_V$.
- Simple discrete models already capture key qualitative behaviors of $U(T)$ and $C_V(T)$.
