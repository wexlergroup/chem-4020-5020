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

# 4.2. Carnot Cycle

## Overview

This section connects **entropy** to the performance limits of heat engines using the **Carnot cycle**.
The handwritten notes for Class 28 emphasize four linked ideas:

1. **Entropy and reversible heat**  
   The differential definition of entropy is (for a reversible path)

   ```{math}
   dS = \frac{\delta q_{\mathrm{rev}}}{T}.
   ```

   In words: entropy change tracks how much *reversible* heat flows at a given temperature.  
   (This is the starting point written at the top of the notes.)  

2. **The Carnot cycle as a “best possible” engine**  
   A Carnot engine is an idealized, fully reversible cycle operating between two thermal reservoirs:
   a hot reservoir at temperature $T_{\mathrm{hot}}$ and a cold reservoir at $T_{\mathrm{cold}}$.
   The cycle consists of **two reversible isotherms** and **two reversible adiabats** (see the sketch on page 1).

3. **Maximum efficiency for converting heat to work**  
   For any engine operating between the same two reservoir temperatures, the Carnot engine has the
   maximum possible efficiency. The notes derive

   ```{math}
   \eta_{\mathrm{Carnot}} = 1 - \frac{T_{\mathrm{cold}}}{T_{\mathrm{hot}}}.
   ```

4. **Second law = direction of heat flow (entropy increase)**  
   The notes use a two-body “A + B in an insulated box” argument to show that spontaneous heat flow
   must go from higher temperature to lower temperature, because the **total entropy of an isolated system**
   must not decrease.

Throughout, we will use the standard Carnot labeling $A \to B \to C \to D \to A$ shown in the P–V sketch:

- $A\to B$: isothermal expansion at $T_{\mathrm{hot}}$  
- $B\to C$: adiabatic expansion (insulated)  
- $C\to D$: isothermal compression at $T_{\mathrm{cold}}$  
- $D\to A$: adiabatic compression (insulated)  

---

## The Cycle

The figure below plots the Carnot cycle on a **pressure–volume (P–V) diagram** for an ideal gas.
The key geometric features match the lecture sketch:

- The **isotherms** (red/blue in the plot) follow the ideal-gas relation $P = nRT/V$ at fixed temperature.
- The **adiabats** (dashed/dotted in the plot) are *steeper* than the isotherms and follow $PV^\gamma = \text{constant}$
  for a reversible adiabatic process of an ideal gas.

```{code-cell} ipython3
:tags: [hide-input]

import numpy as np
import matplotlib.cm as cm
import matplotlib.pyplot as plt
from scipy.constants import R as R_J_per_K_mol
from labellines import labelLines
from myst_nb import glue

# Constants
R = R_J_per_K_mol / 1000  # J/K/mol to kJ/K/mol (keeps PV values numerically small for plotting)

# Isothermal Process for Ideal Gas
def P_isothermal_ideal_gas(V, T, n):
    return n * R * T / V

# Adiabatic Process for Ideal Gas
def P_adiabatic_ideal_gas(V, V1, P1, gamma=5/3):
    return P1 * (V1 / V) ** gamma

# Intersection of two pressure curves
def intersection(P1, P2, V):
    return V[np.argmin(np.abs(P1 - P2))]

# Plotting the Carnot Cycle
fig, axs = plt.subplots(1, 1, figsize=(4, 4))

# Calculate the pressure for each segment of the cycle
V = np.linspace(0.1, 5, 1000)  # L
T_hot = 298.15  # K
T_cold = 273.15  # K
n = 1  # mol
V_A = 0.5  # L
V_C = 1.0  # L
P_AB = P_isothermal_ideal_gas(V, T_hot, n)
P_DA = P_adiabatic_ideal_gas(V, V_A, P_isothermal_ideal_gas(V_A, T_hot, n))
P_CD = P_isothermal_ideal_gas(V, T_cold, n)
P_BC = P_adiabatic_ideal_gas(V, V_C, P_isothermal_ideal_gas(V_C, T_cold, n))

# Plot the segments
isotherm_1 = axs.plot(V, P_AB, "r-", label=r'A$\rightarrow$B isothermal expansion ($T_{\text{hot}}$)')
adiabat_1 = axs.plot(V, P_BC, "k:", label=r'B$\rightarrow$C adiabatic expansion', zorder=-10)
isotherm_2 = axs.plot(V, P_CD, "b-", label=r'C$\rightarrow$D isothermal compression ($T_{\text{cold}}$)')
adiabat_2 = axs.plot(V, P_DA, "k--", label=r'D$\rightarrow$A adiabatic compression', zorder=-10)

# Legend
axs.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0., title="Steps")

# Mark the points
axs.plot(V_A, P_isothermal_ideal_gas(V_A, T_hot, n), "ks", markerfacecolor='none')
axs.annotate("A", (V_A, P_isothermal_ideal_gas(V_A, T_hot, n)), textcoords="offset points", xytext=(10, 10), ha='center', va='center')
V_B = intersection(P_AB, P_BC, V)
P_B = P_isothermal_ideal_gas(V_B, T_hot, n)
axs.plot(V_B, P_B, "ks", markerfacecolor='none')
axs.annotate("B", (V_B, P_B), textcoords="offset points", xytext=(10, 10), ha='center', va='center')
V_C = intersection(P_BC, P_CD, V)
P_C = P_isothermal_ideal_gas(V_C, T_cold, n)
axs.plot(V_C, P_C, "ks", markerfacecolor='none')
axs.annotate("C", (V_C, P_C), textcoords="offset points", xytext=(-10, -10), ha='center', va='center')
V_D = intersection(P_CD, P_DA, V)
P_D = P_isothermal_ideal_gas(V_D, T_cold, n)
axs.plot(V_D, P_D, "ks", markerfacecolor='none')
axs.annotate("D", (V_D, P_D), textcoords="offset points", xytext=(-10, -10), ha='center', va='center')

axs.set_xlim(0.4, 1.1)
axs.set_ylim(2, 6)

axs.set_xticks([V_A, V_B, V_C, V_D])
axs.set_xticklabels([r"$V_A$", r"$V_B$", r"$V_C$", r"$V_D$"])

P_A = P_isothermal_ideal_gas(V_A, T_hot, n)
axs.set_yticks([P_A, P_B, P_C, P_D])
axs.set_yticklabels([r"$P_A$", r"$P_B$", r"$P_C$", r"$P_D$"])

# Add heat-flow annotations (matches the lecture sketch: heat in on the hot isotherm, heat out on the cold isotherm)
x = (V_B + V_D) / 2
y = (P_B + P_D) / 2
axs.annotate("system\nabsorbs\nheat", (x, y - 0.1), textcoords="offset points", xytext=(50, 50), ha='center', va='center', fontsize=10, color='r', arrowprops=dict(facecolor='r', shrink=0.05, edgecolor='w', width=2, headwidth=8, headlength=8))
axs.annotate("system\nreleases\nheat", (x, y - 0.1), textcoords="offset points", xytext=(-50, -50), ha='center', va='center', fontsize=10, color='b', arrowprops=dict(facecolor='b', shrink=0.05, edgecolor='w', width=2, headwidth=8, headlength=8))

axs.grid()

plt.show()
plt.close(fig)
```

The Carnot cycle for an ideal gas. The isothermal expansion and compression curves follow the ideal-gas equation
(Boyle's-law form at fixed $T$), while the adiabatic expansion and compression curves follow the adiabatic equation
(which produces a steeper curve than the isotherms).

### Step-by-step thermodynamic description (matching the lecture sketch)

| Step | Process type | Thermal contact | Heat $q$ | Entropy change $\Delta S$ | Qualitative work |
| --- | --- | --- | --- | --- | --- |
| $A\to B$ | Reversible isothermal *expansion* at $T_{\mathrm{hot}}$ | Hot reservoir | $q_{AB} > 0$ | $\Delta S_{AB} = q_{AB}/T_{\mathrm{hot}}$ | System does work on surroundings |
| $B\to C$ | Reversible adiabatic *expansion* | Insulated | $q_{BC}=0$ | $\Delta S_{BC}=0$ | System does work; temperature drops |
| $C\to D$ | Reversible isothermal *compression* at $T_{\mathrm{cold}}$ | Cold reservoir | $q_{CD} < 0$ | $\Delta S_{CD} = q_{CD}/T_{\mathrm{cold}}$ | Work is done on the system |
| $D\to A$ | Reversible adiabatic *compression* | Insulated | $q_{DA}=0$ | $\Delta S_{DA}=0$ | Work on system; temperature rises back to $T_{\mathrm{hot}}$ |

Two quick consequences (both appear explicitly in the notes):

- Because adiabatic steps have $q=0$, the reversible definition $dS = \delta q_{\mathrm{rev}}/T$ gives
  $\Delta S = 0$ on the adiabats (**entropy is constant** on a reversible adiabat).
- Over one full cycle, the system returns to its initial state, so

  ```{math}
  \Delta S_{\text{cycle}}
  = \Delta S_{AB} + \Delta S_{BC} + \Delta S_{CD} + \Delta S_{DA}
  = \Delta S_{AB} + \Delta S_{CD}
  = 0,
  ```

  which is exactly the relation written in blue in the notes.

---

## How to Operate a Carnot Cycle/Engine

The notes give a piston-and-weights picture for how to realize each step (page 2).

**Core idea:** you alternate between (i) placing the gas in contact with a thermal reservoir so it can exchange heat *reversibly* at a fixed temperature (isotherms), and (ii) insulating the gas so it cannot exchange heat (adiabats).

A clean way to narrate the sequence is:

1. **$A\to B$: isothermal expansion at $T_{\mathrm{hot}}$**  
   Put the cylinder in contact with a hot thermal source. Remove weights slowly so the piston rises
   quasi-statically. The system **absorbs heat** $q_{AB}$ from the hot reservoir to keep $T$ constant,
   while doing work on the weights.

2. **$B\to C$: adiabatic expansion (insulated)**  
   Remove the thermal contact and insulate the cylinder. Continue removing weights slowly so the piston rises.
   No heat flows ($q=0$), so the work done by the gas comes from its internal energy, and the temperature drops
   from $T_{\mathrm{hot}}$ to $T_{\mathrm{cold}}$.

3. **$C\to D$: isothermal compression at $T_{\mathrm{cold}}$**  
   Put the cylinder in contact with a cold thermal source. Add weights slowly so the piston lowers.
   The surroundings do work on the gas, and the gas **releases heat** $q_{CD}$ to the cold reservoir to keep $T$ constant.

4. **$D\to A$: adiabatic compression (insulated)**  
   Insulate the cylinder again and keep adding weights slowly to compress the gas. With no heat exchange,
   the work done on the gas increases its internal energy and raises the temperature back to $T_{\mathrm{hot}}$,
   returning to state $A$. (This “work $\to$ heat” idea is drawn in the notes.)

The schematic below is a simplified visualization of those operating steps.

```{code-cell} ipython3
:tags: [hide-input]

import numpy as np
import matplotlib.cm as cm
import matplotlib.pyplot as plt
from scipy.constants import R as R_J_per_K_mol
from labellines import labelLines
from myst_nb import glue

fig, axs = plt.subplot_mosaic([[0, 1, 2, 3]], figsize=(8, 4), constrained_layout=True, sharex=True, sharey=True)

axs[0].set_title('(A)')
axs[1].set_title('(B)')
axs[2].set_title('(C)')
axs[3].set_title('(D)')

piston_height = [1.5, 3.5, 4.5, 2.5]
step = ["A", "B", "C", "D"]
for i in range(4):
    # Cylinder
    axs[i].plot([0, 0], [5, 0], 'k-')
    cylinder = axs[i].plot([0, 3], [0, 0], 'k-',label='cylinder')
    axs[i].plot([3, 3], [0, 5], 'k-')
    labelLines(cylinder, xvals=[1.5], fontsize=10, color='k')

    # Thermal source / insulation block
    axs[i].plot([0, 0], [-1, -2], 'k:')
    thermal_source = axs[i].plot([0, 3], [-2, -2], 'k:',label='thermal source')
    axs[i].plot([3, 3], [-2, -1], 'k:')
    axs[i].plot([3, 0], [-1, -1], 'k:')
    labelLines(thermal_source, xvals=[1.5], fontsize=10, color='k')
    if i == 0:
        axs[i].text(1.5, -1.5, 'insulation', fontsize=10, ha='center', va='center')
    if i == 1:
        axs[i].text(1.5, -1.5, r'$T_{\text{hot}}$', fontsize=10, ha='center', va='center')
    if i == 2:
        axs[i].text(1.5, -1.5, 'insulation', fontsize=10, ha='center', va='center')
    if i == 3:
        axs[i].text(1.5, -1.5, r'$T_{\text{cold}}$', fontsize=10, ha='center', va='center')

    # Piston
    axs[i].plot([1.5, 1.5], [piston_height[i], 6], 'k-', lw=2)
    piston = axs[i].plot([0, 3], [piston_height[i], piston_height[i]], 'k-', lw=2, label='piston')
    labelLines(piston, xvals=[1.5], fontsize=10, color='k')

    # P, V labels
    label = '$V_{\\text{' + step[i] + '}}, P_{\\text{' + step[i] + '}}$'
    axs[i].text(1.5, 0.75, label, fontsize=10, ha='center', va='center')

    # Heat transfer (only on isothermal steps)
    if i == 1:
        axs[i].annotate("heat", (0.5, 0.75), textcoords="offset points", xytext=(0, -60), ha='center', va='center', fontsize=10, color='r', arrowprops=dict(facecolor='r', shrink=0.05, edgecolor='w', width=2, headwidth=8, headlength=8))
    if i == 3:
        axs[i].annotate("heat", (0.5, -1.5), textcoords="offset points", xytext=(0, 60), ha='center', va='center', fontsize=10, color='r', arrowprops=dict(facecolor='r', shrink=0.05, edgecolor='w', width=2, headwidth=8, headlength=8))
    
    # Work (qualitative; work occurs on every leg, but the direction differs)
    if i == 1:
        axs[i].annotate("work done by system", (1.5, piston_height[i] - 0.1), textcoords="offset points", xytext=(0, -30), ha='center', va='center', fontsize=10, color='b', arrowprops=dict(facecolor='b', shrink=0.05, edgecolor='w', width=2, headwidth=8, headlength=8))
    if i == 3:
        axs[i].annotate("work done on system", (1.5, piston_height[i] - 1.5), textcoords="offset points", xytext=(0, 30), ha='center', va='center', fontsize=10, color='b', arrowprops=dict(facecolor='b', shrink=0.05, edgecolor='w', width=2, headwidth=8, headlength=8))

    # Axis settings
    axs[i].set_ylim(-3, 6)

plt.show()
```

---

## Maximum Efficiency of Converting Heat to Work (Carnot efficiency)

The notes define the heat quantities on the two isothermal legs:

- $q_{AB}$: heat **absorbed by the system** from the hot reservoir to maintain $T_{\mathrm{hot}}$  
- $q_{CD}$: heat **released by the system** to the cold reservoir to maintain $T_{\mathrm{cold}}$

For one cycle, energy conservation gives

```{math}
W_{\text{net}} = q_{AB} + q_{CD}.
```

Here $q_{AB} > 0$ and $q_{CD} < 0$, so $W_{\text{net}}$ is the net work output (area enclosed on the P–V diagram).

The thermal efficiency is defined as

```{math}
\eta
= \frac{\text{useful work out}}{\text{heat absorbed}}
= \frac{W_{\text{net}}}{q_{AB}}
= \frac{q_{AB} + q_{CD}}{q_{AB}}.
```

Now use the entropy relation on the isotherms (this is exactly the blue statement in the notes):

```{math}
\Delta S_{AB} = \frac{q_{AB}}{T_{\mathrm{hot}}},\qquad
\Delta S_{CD} = \frac{q_{CD}}{T_{\mathrm{cold}}}.
```

Because the full reversible cycle returns the system to its initial state,

```{math}
\Delta S_{\text{cycle}} = \Delta S_{AB} + \Delta S_{CD} = 0
\quad\Rightarrow\quad
\Delta S_{CD} = -\Delta S_{AB}.
```

Substitute $q_{AB} = T_{\mathrm{hot}}\,\Delta S_{AB}$ and $q_{CD} = T_{\mathrm{cold}}\,\Delta S_{CD}$ into the efficiency:

```{math}
\eta
= \frac{T_{\mathrm{hot}}\Delta S_{AB} + T_{\mathrm{cold}}\Delta S_{CD}}{T_{\mathrm{hot}}\Delta S_{AB}}
= \frac{T_{\mathrm{hot}}\Delta S_{AB} - T_{\mathrm{cold}}\Delta S_{AB}}{T_{\mathrm{hot}}\Delta S_{AB}}
= 1 - \frac{T_{\mathrm{cold}}}{T_{\mathrm{hot}}}.
```

**Interpretation:** the maximum possible fraction of hot-reservoir heat that can be converted to work depends
*only* on the two reservoir temperatures. Lowering $T_{\mathrm{cold}}$ or raising $T_{\mathrm{hot}}$ increases the
maximum possible efficiency.

---

## The Second Law: Entropy and the Direction of Heat Flow

The notes (page 3) use a simple composite system:

- Two subsystems $A$ and $B$ are inside an **insulating boundary** (no heat exchange with the surroundings).
- Each subsystem has fixed volume (so $dV_A = dV_B = 0$).
- $A$ and $B$ *can exchange heat with each other* through an internal wall.

Because the composite $A+B$ is isolated,

```{math}
dU_A + dU_B = 0
\quad\Rightarrow\quad
dU_B = -dU_A.
```

Entropy is extensive, so the total entropy is

```{math}
S = S_A + S_B
\quad\Rightarrow\quad
dS = dS_A + dS_B.
```

At fixed volume, the fundamental relation written in the notes,
$dU = T\,dS - P\,dV$, reduces to $dU = T\,dS$. Therefore

```{math}
dS_A = \frac{dU_A}{T_A},\qquad dS_B = \frac{dU_B}{T_B}.
```

So the total entropy change is

```{math}
dS
= \frac{dU_A}{T_A} + \frac{dU_B}{T_B}
= \frac{dU_A}{T_A} - \frac{dU_A}{T_B}
= dU_A\left(\frac{1}{T_A} - \frac{1}{T_B}\right).
```

**Second law statement for an isolated system:** spontaneous change requires $dS \ge 0$, with equality for a reversible exchange.

Now suppose $T_A > T_B$ ($A$ is hotter than $B$). Then
$\tfrac{1}{T_A} - \tfrac{1}{T_B} < 0$.  
To satisfy $dS \ge 0$, we must have $dU_A \le 0$, meaning subsystem $A$ loses internal energy while $B$ gains it.

That is exactly the familiar conclusion:

> **Heat flows spontaneously from hot to cold.**

Attempting to make heat flow from cold to hot without any other change would make $dS < 0$ for the isolated
composite system, which violates the second law.

---

## Summary (what to remember)

- $dS = \delta q_{\mathrm{rev}}/T$ is the entropy bookkeeping rule highlighted in the notes.
- The Carnot cycle consists of two reversible isotherms (heat exchange at fixed $T$) and two reversible adiabats (no heat exchange).
- For a reversible Carnot engine between $T_{\mathrm{hot}}$ and $T_{\mathrm{cold}}$,

  ```{math}
  \eta_{\mathrm{Carnot}} = 1 - \frac{T_{\mathrm{cold}}}{T_{\mathrm{hot}}}.
  ```

- The second law can be phrased as: for an isolated system, $\Delta S \ge 0$, which enforces the direction of spontaneous heat flow.
