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

# 2.2. Canonical Ensemble

## Overview

In a **closed system**—one that exchanges energy but not matter with its surroundings—the appropriate statistical description is the **canonical ensemble**.

```{code-cell} ipython3
:tags: [hide-input]

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from myst_nb import glue

# Helper function to plot a system
def plot_system(ax, title, annotations, boundary_color='b'):
    box = mpatches.FancyBboxPatch((0, 0), 1, 1, boxstyle='roundtooth', ec=boundary_color, fc='w')
    ax.add_patch(box)
    ax.set_title(title, fontsize=14)
    ax.text(0.5, 0.5, 'System', ha='center', va='center', fontsize=12)
    ax.text(0.5, -0.65, 'Surroundings', ha='center', va='center', fontsize=12)
    ax.text(0.5, 1.3, 'Boundary', ha='center', va='bottom', fontsize=12, color=boundary_color)
    for annotation in annotations:
        if "arrowprops" in annotation:  # Arrow annotations
            ax.annotate('', **annotation)
        else:  # Text annotations
            ax.text(**annotation)
    ax.set_xlim(-1, 2)
    ax.set_ylim(-1, 2)
    ax.set_aspect('equal')
    ax.axis('off')

# Define annotations for each system
annotations = [
        dict(xy=(-0.6, 0.15), xytext=(0.15, 0.15), arrowprops=dict(arrowstyle='<->', color='r')),
        dict(x=-1, y=0.3, s='Energy', ha='left', va='bottom', fontsize=12, color='r'),
]

fig, ax = plt.subplots(1, 1, figsize=(4, 4))
plot_system(ax, "", annotations)

plt.show()
plt.close(fig)
```

A closed system, exchanging energy but not matter with its surroundings.

## Probability of a Microstate in the Canonical Ensemble

Consider an ensemble of $\mathcal{A}$ closed systems exchanging energy with a heat bath at temperature $T$. A heat bath is an environment that can absorb or release energy without changing its temperature because it is much larger than the system.

````{admonition} How Can an Environment Absorb or Release Energy Without Changing Its Temperature?
:class: dropdown
The key point is that the temperature of a collection of particles is determined by the average kinetic energy per particle. According to the equipartition theorem, for particles with three translational degrees of freedom the total kinetic energy is

```{math}
E_{\text{kin}} = \frac{3}{2} N\, k_\text{B}\, T.
```

Thus, for a fixed amount of energy transfer, the corresponding change in temperature is given by

```{math}
\Delta T = \frac{2}{3}\frac{\Delta E_\text{kin}}{N\, k_\text{B}}.
```

Because the environment has many more particles than the system ($N_\text{env} \gg N_\text{sys}$), the same amount of energy $\Delta E$ will result in a much smaller $\Delta T$ for the environment. In our example, even though the system (with $N_\text{sys} = 6.022\times10^{23}$ particles) absorbs energy and its temperature increases from 300 K to 400 K, the environment (with $N_\text{env} = 6.022\times10^{24}$ particles, i.e. 10 times as many) loses the same amount of energy per mole but its temperature only drops from 400 K to 390 K.

This illustrates that an environment (or heat bath) can absorb or release energy almost isothermally because its very large number of particles (and hence its high heat capacity) means that the fractional change in average kinetic energy per particle is negligible—even though energy is being exchanged.

Thus, even though the total kinetic energy of the universe (system plus environment) remains constant, the environment's temperature barely changes because any energy loss or gain is diluted among a huge number of particles.
````

### Ratio of Numbers of Systems in Two Microstates is a Function of Their Relative Energies

Intuition tells us that a system is more likely to be found in microstates with lower energy.

```{admonition} Building This Intuition
:class: dropdown
Consider a raindrop that can exist at two different elevations: one high up at Lake Itasca in northern Minnesota and one low down in the Mississippi River in St. Louis. At the higher elevation, the raindrop has more gravitational potential energy, making that state less favorable. Consequently, the raindrop is more likely to be found at the lower elevation, where its gravitational potential energy is reduced. In essence, systems tend to "prefer" lower energy states, which is why, on a macroscopic scale, water naturally flows downhill.
```

Mathematically, this intuition asserts that the ratio of the numbers $a_1$ and $a_2$ of systems in two microstates 1 and 2 is given by

```{math}
\frac{a_2}{a_1} = f \left( E_1, E_2 \right) = f \left( E_1 - E_2 \right),
```

where $E_1$ and $E_2$ are the energies of the two microstates. The function $f$ depends only on the difference in energy between the two states.

```{admonition} Why Does the Ratio Depend Only on the Energy Difference?
:class: dropdown
The energy of a system is defined relative to an arbitrary (but often convenient) reference level. For example,

* The kinetic energy of a moving particle is defined relative to a stationary particle.

* The gravitational potential energy of a raindrop is defined relative to the surface of the Earth.

Therefore, the absolute energy of a system is not essential in determining the likelihood of finding the system in a particular microstate, only the difference in energy between two microstates.
```

### Finding an Acceptable Form for $f$

Since $\{ a_1, a_2, a_3, \ldots \}$ is a set of numbers, we can write

```{math}
\begin{align*}
\frac{a_3}{a_1} &= \frac{a_2}{a_1} \times \frac{a_3}{a_2} \\
f \left( E_1 - E_3 \right) &= f \left( E_1 - E_2 \right) \times f \left( E_2 - E_3 \right).
\end{align*}
```

If $f$ is "well-behaved" (i.e., continuous, measurable, etc.), it must be of the form

```{math}
f \left( E_m - E_n \right) = e^{\beta \left( E_m - E_n \right)},
```

where $\beta$ is an undetermined constant, which we will determine later to be $\frac{1}{k_\text{B} T}$.

```{admonition} Checking the Form of $f$
:class: dropdown
Using this form for $f$, we can verify that

```{math}
\begin{align*}
f \left( E_1 - E_3 \right) &= f \left( E_1 - E_2 \right) \times f \left( E_2 - E_3 \right) \\
e^{\beta \left( E_1 - E_3 \right)} &= e^{\beta \left( E_1 - E_2 \right)} \times e^{\beta \left( E_2 - E_3 \right)} \\
{\color{blue} \cancel{e^{\beta E_1}}} \times {\color{magenta} \cancel{e^{-\beta E_3}}} &= {\color{blue} \cancel{e^{\beta E_1}}} \times {\color{red} \cancel{e^{-\beta E_2}}} \times {\color{red} \cancel{e^{-\beta E_2}}} {\color{magenta} \cancel{e^{-\beta E_3}}} \\
1 &= 1.
\end{align*}
```

### Converting $f$ to a Probability

Separating the variables with indices $m$ and $n$,

```{math}
a_m e^{\beta E_m} = a_n e^{\beta E_n} = \mathcal{C},
```

where $\mathcal{C}$ is a constant. Therefore, the number of systems in a microstate with energy $E_m$ is

```{math}
a_m = \mathcal{C} e^{-\beta E_m}.
```

The constant $\mathcal{C}$ is determined by the normalization condition

```{math}
\sum_m a_m = \mathcal{A} = \mathcal{C} \sum_m e^{-\beta E_m},
```

Solving the second equation for $\mathcal{C}$ and substituting it into the equation for $a_m$ gives

```{math}
:label: canonical-ensemble-probability
p_m = \frac{a_m}{\mathcal{A}} = \frac{e^{-\beta E_m}}{\sum_m e^{-\beta E_m}} = \frac{e^{-\beta E_m}}{Q},
```

where $p_m$ is the probability of finding the system in microstate $m$ and $Q$ is called the **partition function**, which is the essential quantity in statistical mechanics.

````{admonition} Alternative derivation roadmap: where does $p_i \propto e^{-\beta E_i}$ come from?
:class: dropdown

We will use the canonical-ensemble probability
```{math}
p_i \propto e^{-\beta E_i}
\quad\text{with}\quad
\beta = \frac{1}{k_B T}
```
throughout this module. Before we use it heavily, here are two quick (complementary) ways to see
**why the exponential “Boltzmann factor” appears**.

---

### Setup (common to both routes)

- **System** $S$: a *closed* system (fixed $N,V$) with discrete microstates $i=1,2,\dots$ and energies $E_i$.
- **Reservoir / heat bath** $R$: very large, can exchange energy with $S$, and stays at (approximately) constant temperature $T$.
- **Universe** $S+R$ is **isolated**, so the total energy is fixed:
  ```{math}
  E_{\text{tot}} = E_i + E_R.
  ```
- **Weak coupling**: the interaction energy is negligible, so energies add.

The question is: *What is the probability $p_i$ that the system is in microstate $i$?*

---

## Route A: “Counting” argument (microcanonical $\rightarrow$ canonical)

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
   S_R(E) = k_B \ln \Omega_R(E).
   ```
   Taking logs,
   ```{math}
   \ln p_i = \text{const} + \ln \Omega_R(E_{\text{tot}}-E_i)
   = \text{const} + \frac{S_R(E_{\text{tot}}-E_i)}{k_B}.
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
   p_i \propto \exp\!\left(\frac{S_R(E_{\text{tot}})}{k_B}\right)\,
          \exp\!\left(-\frac{E_i}{k_B T}\right).
   ```
   The first exponential is just a constant (independent of $i$), so the *physics* is:
   ```{math}
   \boxed{p_i \propto e^{-E_i/(k_B T)} \equiv e^{-\beta E_i}}
   \quad\text{with}\quad
   \beta \equiv \frac{1}{k_B T}.
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

## Route B: “Maximum entropy” argument (Gibbs entropy + constraints)

**Idea:** in equilibrium with a heat bath, the system adopts the probability distribution that
maximizes entropy subject to what is fixed.

1. Start from the Gibbs/Shannon entropy functional for a probability distribution $\{p_i\}$:
   ```{math}
   S[\{p_i\}] = -k_B \sum_i p_i \ln p_i.
   ```

2. Impose the two canonical constraints:
   - normalization: $\sum_i p_i = 1$
   - fixed mean energy: $\sum_i p_i E_i = U$

3. Maximize $S$ with Lagrange multipliers $\alpha$ and $\lambda$:
   ```{math}
   \delta\Big[-k_B\sum_i p_i\ln p_i
   -\alpha\Big(\sum_i p_i-1\Big)
   -\lambda\Big(\sum_i p_iE_i-U\Big)\Big]=0.
   ```

4. The maximizer satisfies
   ```{math}
   p_i = \text{const}\times e^{-(\lambda/k_B)E_i}.
   ```
   Rename $\beta \equiv \lambda/k_B$:
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
   For the canonical distribution one can show $dS/dU = k_B\beta$, so
   ```{math}
   \boxed{\beta=\frac{1}{k_B T}}.
   ```

---

## Quick checks / intuition (why this form makes sense)

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

**Bottom line:** In the canonical ensemble, the exponential appears because *the reservoir’s multiplicity grows
exponentially with its entropy*, and the entropy changes approximately linearly with energy over the small energy
exchanges relevant to a large heat bath.
````

## Two-State System

Consider a system with two states: state 1 with energy $E_1$ and state 2 with energy $E_2$.

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
  - Unfolded
  - Folded
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
p_2 = 1 - p_1
```

### Partition Function as the Effective Number of Thermally Accessible Microstates

```{code-cell} ipython3
:tags: [hide-input]

import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import k, eV
from labellines import labelLines
from myst_nb import glue
from matplotlib.patches import Rectangle

k_B = k / eV  # Boltzmann constant in eV/K

# Define the partition function for a two-state system
def partition_function_two_state(E1, E2, T):
    beta = 1 / (k_B * T)
    return np.exp(-beta * E1) + np.exp(-beta * E2)

# Calculate the partition function for a two-state system
E1 = 0
E2 = 0.01  # Energy difference between the two states in eV
T_values = np.linspace(1, 1000, 1000)
Q_values = [partition_function_two_state(E1, E2, T) for T in T_values]

# Calculate the probabilities of finding the system in each state for a two-state system
p1_values = [np.exp(-1 / (k_B * T) * E1) / Q for T, Q in zip(T_values, Q_values)]
p2_values = [np.exp(-1 / (k_B * T) * E2) / Q for T, Q in zip(T_values, Q_values)]

# Plot the partition function and the probabilities of finding the system in each state
fig, axs = plt.subplots(1, 2, figsize=(8, 4))

# Plot the partition function on axs[0]
axs[0].plot(T_values, Q_values, 'k-')
axs[0].set_xlabel('Temperature (K)')
axs[0].set_ylabel('$Q_{\\text{two-state}}$')
axs[0].grid(True)
axs[0].annotate(
    '$\\rightarrow 1$ accessible\nmicrostate', xy=(40, Q_values[0] + 0.01), xytext=(300, 1.1),
    arrowprops=dict(arrowstyle='->', color='b'),
    bbox=dict(boxstyle='round,pad=0.3', fc='w', ec='b'),
    ha='center', va='center', color='b'
)
axs[0].annotate(
    '$\\rightarrow 2$ accessible\nmicrostates', xy=(T_values[-1], Q_values[-1]), xytext=(750, 1.7),
    arrowprops=dict(arrowstyle='->', color='m'),
    bbox=dict(boxstyle='round,pad=0.3', fc='w', ec='m'),
    ha='center', va='center', color='m'
)
axs[0].set_xlim(0, 1000)
axs[0].set_ylim(1, 2)

# Plot the probabilities of finding the system in each state on axs[1]
p1_line, = axs[1].plot(T_values, p1_values, 'b-', label='State 1')
p2_line, = axs[1].plot(T_values, p2_values, 'r-', label='State 2')
labelLines([p1_line, p2_line], zorder=2.5)
axs[1].set_xlabel('Temperature (K)')
axs[1].set_ylabel('Probability')
axs[1].grid(True)
axs[1].set_ylim(0, 1)  # ensure y-axis spans from 0 to 1

# Add tall outlined rectangles around the probabilities at low and high temperatures.
# For T -> 0: highlight T from 1 to 50 K.
# For T = 1000: highlight T from 950 to 1000 K.
rect_low = Rectangle((1, 0), 50 - 1, 1, edgecolor='b', facecolor='b', linestyle='-', alpha=0.2)
rect_high = Rectangle((950, 0), 1000 - 950, 1, edgecolor='m', facecolor='m', linestyle='-', alpha=0.2)
axs[1].add_patch(rect_low)
axs[1].add_patch(rect_high)
axs[1].annotate(
    'Only state 1 is accessed', xy=(50, 0.95), xytext=(100, 0.95),
    arrowprops=dict(arrowstyle='-', color='b'),
    ha='left', va='center', color='b'
)
axs[1].annotate(
    'Both states are accessed', xy=(950, 0.05), xytext=(900, 0.05),
    arrowprops=dict(arrowstyle='-', color='m'),
    ha='right', va='center', color='m'
)

plt.tight_layout()
plt.show()
plt.close(fig)
```

Partition function of a two-state system and the probabilities of finding it in each state as a function of temperature. The energy difference between the two states is 0.01 eV.

## Computational Studio: Two-Level System

Explore how the canonical probabilities, partition function, and heat capacity respond to changes in the energy gap and temperature range.

<div style="width: 100%; border: 1px solid #cbd5e1; border-radius: 12px; overflow: hidden; box-shadow: 0 4px 12px rgba(0,0,0,0.08);">
  <iframe
    src="https://chem-4020-5020.vercel.app/"
    title="Two-Level System Studio"
    style="width: 100%; height: 850px; border: 0;"
    loading="lazy"
  ></iframe>
</div>

If the embed does not load, you can open the studio in a new tab:
<a href="https://chem-4020-5020.vercel.app/" target="_blank" rel="noopener">
  Two-Level System Studio
</a>.
