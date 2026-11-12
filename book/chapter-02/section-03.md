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

```{code-cell} ipython3
:tags: [hide-input]

import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import k, eV
from labellines import labelLines

kB = k / eV  # Boltzmann constant in eV/K

# Define the partition function for a two-state system
def partition_function_two_state(E1, E2, T):
    beta = 1 / (kB * T)
    return np.exp(-beta * E1) + np.exp(-beta * E2)

# Calculate the partition function for a two-state system
E1 = 0
E2 = 0.01  # Energy of state 2 in eV (equal to the energy gap, since E1 = 0)
T_values = np.linspace(1, 1000, 1000)
Q_values = [partition_function_two_state(E1, E2, T) for T in T_values]

# Calculate the internal energy for a two-state system
beta_values = 1 / (kB * T_values)
ln_Q_values = np.log(Q_values)
U_values = -np.gradient(ln_Q_values, beta_values)

# Calculate the heat capacity in units of kB from the two-state result
# C_V / kB = x^2 e^(-x) / (1 + e^(-x))^2, with x = (E2 - E1) / (kB T)
x_values = (E2 - E1) * beta_values
C_V_values = x_values**2 * np.exp(-x_values) / (1 + np.exp(-x_values))**2

fig, axs = plt.subplots(1, 2, figsize=(8, 4))

# Plot the internal energy, relative to E1, as a function of temperature
ax = axs[0]
ax.plot(T_values, U_values - E1, 'k-')
ax.set_xlabel('Temperature (K)')
ax.set_ylabel('$U_{\\text{two-state}} - E_1$ (eV)')
ax.grid(True)
ax.annotate(
    '$U_{\\text{two-state}} - E_1 \\rightarrow 0$', xy=(40, U_values[0] - E1 + 0.0001), xytext=(300, 0.001),
    arrowprops=dict(arrowstyle='->', color='b'),
    bbox=dict(boxstyle='round,pad=0.3', fc='w', ec='b'),
    ha='center', va='center', color='b'
)
x_values_high_T = np.linspace(0, 1000, 1001)
y_values_high_T = ((E2 - E1) / 2) * np.ones_like(x_values_high_T)
line = ax.plot(x_values_high_T, y_values_high_T, 'm--', label='$U_{\\text{two-state}} - E_1 \\rightarrow (E_2 - E_1) / 2$')
labelLines(line, zorder=2.5)
ax.set_xlim(0, 1000)
ax.set_ylim(0, 0.006)

# Plot the heat capacity at constant volume as a function of temperature
ax = axs[1]
ax.plot(T_values, C_V_values, 'k-')
ax.set_xlabel('Temperature (K)')
ax.set_ylabel('$C_V / k_{\\mathrm{B}}$')
ax.grid(True)
ax.set_xlim(0, 1000)
ax.set_ylim(0, 0.5)

plt.tight_layout()
plt.show()
plt.close(fig)
```

Internal energy (left) and heat capacity at constant volume (right; defined in the next subsection) of a two-state system as functions of temperature, for an energy gap of $0.01\ \mathrm{eV}$.

### Heat Capacity at Constant Volume

From [Section 1.1](../chapter-01/section-01.md), **heat** is energy transferred due to a temperature difference. The **heat capacity** at constant volume, $C_V$, measures how much heat is required to change a system’s temperature at fixed $N$ and $V$:

```{math}
C_V 
\;=\; \left(\frac{\partial U}{\partial T}\right)_{N,V}.
```

The right panel of the figure above shows $C_V$ for the two-state system.

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
