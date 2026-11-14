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

See also: [Course-wide Conventions & Notation](../notation.md)

## Overview and Learning Objectives

The Carnot cycle is a fully reversible heat-engine cycle that sets the maximum efficiency any engine can achieve between two reservoir temperatures. Using the entropy definition $dS=\delta q_{\mathrm{rev}}/T$, the cycle relates heat flows to temperature and establishes $\eta_{\mathrm{Carnot}}=1-T_{\mathrm{cold}}/T_{\mathrm{hot}}$. The same reasoning yields a second-law statement about the direction of spontaneous heat flow.

This section connects **entropy** to the performance limits of heat engines using the **Carnot cycle**. It develops four linked ideas:

1. **Entropy and reversible heat.**
   From [Section 4.1](section-01.md), the differential definition of entropy for a reversible path is

   ```{math}
   dS = \frac{\delta q_{\mathrm{rev}}}{T}.
   ```

   Entropy change tracks how much *reversible* heat flows at a given temperature.

2. **The Carnot cycle as a "best possible" engine.**
   A Carnot engine is an idealized, fully reversible cycle operating between two thermal reservoirs:
   a hot reservoir at temperature $T_{\mathrm{hot}}$ and a cold reservoir at $T_{\mathrm{cold}}$.
   The cycle consists of **two reversible isotherms** and **two reversible adiabats**.

3. **Maximum efficiency for converting heat to work.**
   For any engine operating between the same two reservoir temperatures, the Carnot engine has the
   maximum possible efficiency:

   ```{math}
   \eta_{\mathrm{Carnot}} = 1 - \frac{T_{\mathrm{cold}}}{T_{\mathrm{hot}}}.
   ```

4. **Second law as the direction of heat flow.**
   A two-body argument shows that spontaneous heat flow must go from higher temperature to lower temperature,
   because the **total entropy of an isolated system** must not decrease.

Throughout, we use the standard Carnot labeling $A \to B \to C \to D \to A$:

- $A\to B$: isothermal expansion at $T_{\mathrm{hot}}$
- $B\to C$: adiabatic expansion (insulated)
- $C\to D$: isothermal compression at $T_{\mathrm{cold}}$
- $D\to A$: adiabatic compression (insulated)

---

Learning objectives:

- Describe the four reversible steps of the Carnot cycle (two isotherms and two adiabats) on a $P$–$V$ diagram.
- Use entropy changes on reversible isotherms and adiabats to show $\Delta S_{\mathrm{cycle}}=0$.
- Derive the Carnot efficiency $\eta_{\mathrm{Carnot}}=1-T_{\mathrm{cold}}/T_{\mathrm{hot}}$.
- Explain why spontaneous heat flow from hot to cold is required by non-decreasing total entropy.

## Core Ideas and Derivations

### The Cycle

The figure below plots the Carnot cycle on a **pressure–volume ($P$–$V$) diagram** for one mole of a monatomic ideal gas between reservoirs at $500\ \mathrm{K}$ and $300\ \mathrm{K}$, the temperatures of the Worked Example below.
The key geometric features are:

- The **isotherms** (red at $T_{\mathrm{hot}}$, blue at $T_{\mathrm{cold}}$) follow the ideal-gas relation $P = nRT/V$ at fixed temperature. The thin lines continue them beyond the cycle.
- The **adiabats** (dotted and dashed) are *steeper* than the isotherms and follow $PV^\gamma = \text{constant}$
  for a reversible adiabatic process of an ideal gas ([Section 3.2](../chapter-03/section-02.md)).

```{code-cell} ipython3
:tags: [hide-input]

import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import R as R_J

# One mole of a monatomic ideal gas between the reservoirs of the Worked Example
R = R_J / 100        # L bar K^-1 mol^-1
n = 1.0              # mol
gamma = 5 / 3
T_hot, T_cold = 500.0, 300.0   # K
V_A, V_B = 1.0, 2.6            # L

# Reversible adiabats: T V^(gamma - 1) = constant
ratio = (T_hot / T_cold) ** (1 / (gamma - 1))
V_C, V_D = V_B * ratio, V_A * ratio

def P_isotherm(V, T):
    return n * R * T / V

def P_adiabat(V, V0, T0):
    return P_isotherm(V0, T0) * (V0 / V) ** gamma

steps = [
    ("A", "B", V_A, V_B, lambda V: P_isotherm(V, T_hot), "r-",
     r"A$\rightarrow$B isothermal expansion ($T_{\mathrm{hot}}$)"),
    ("B", "C", V_B, V_C, lambda V: P_adiabat(V, V_B, T_hot), "k:",
     r"B$\rightarrow$C adiabatic expansion"),
    ("C", "D", V_C, V_D, lambda V: P_isotherm(V, T_cold), "b-",
     r"C$\rightarrow$D isothermal compression ($T_{\mathrm{cold}}$)"),
    ("D", "A", V_D, V_A, lambda V: P_adiabat(V, V_A, T_hot), "k--",
     r"D$\rightarrow$A adiabatic compression"),
]

fig, ax = plt.subplots(figsize=(6, 4.2))

# Faint isotherms for reference
V_all = np.linspace(0.8, 6.2, 400)
ax.plot(V_all, P_isotherm(V_all, T_hot), "r-", lw=0.6, alpha=0.35)
ax.plot(V_all, P_isotherm(V_all, T_cold), "b-", lw=0.6, alpha=0.35)

# The four steps, each with an arrow that shows its direction
for start, end, V1, V2, P_of_V, style, label in steps:
    V = np.linspace(V1, V2, 200)
    ax.plot(V, P_of_V(V), style, lw=2.2, label=label)
    Vm = np.sqrt(V1 * V2)
    dV = 0.02 * (V2 - V1)
    ax.annotate("", xy=(Vm + dV, P_of_V(Vm + dV)), xytext=(Vm - dV, P_of_V(Vm - dV)),
                arrowprops=dict(arrowstyle="-|>", color=style[0], lw=0, mutation_scale=16))

# Corner points
corners = {"A": (V_A, P_isotherm(V_A, T_hot)), "B": (V_B, P_isotherm(V_B, T_hot)),
           "C": (V_C, P_isotherm(V_C, T_cold)), "D": (V_D, P_isotherm(V_D, T_cold))}
offsets = {"A": (9, 6), "B": (8, 8), "C": (8, 8), "D": (-8, -9)}
for name, (Vc, Pc) in corners.items():
    ax.plot(Vc, Pc, "ks", markerfacecolor="w", zorder=5)
    ax.annotate(name, (Vc, Pc), textcoords="offset points", xytext=offsets[name], ha="center", va="center")

# Heat enters on the hot isotherm and leaves on the cold isotherm
V_in = 1.7
ax.annotate("", xy=(V_in, P_isotherm(V_in, T_hot) + 0.6), xytext=(V_in + 0.7, P_isotherm(V_in, T_hot) + 9),
            arrowprops=dict(arrowstyle="-|>", color="r", lw=1.5, mutation_scale=14))
ax.text(V_in + 0.75, P_isotherm(V_in, T_hot) + 9.3, "heat in\n" + r"($q_{AB} > 0$)", color="r", ha="center", va="bottom")
V_out = 3.8
ax.annotate("", xy=(V_out - 0.6, 1.8), xytext=(V_out, P_isotherm(V_out, T_cold) - 0.5),
            arrowprops=dict(arrowstyle="-|>", color="b", lw=1.5, mutation_scale=14))
ax.text(V_out - 0.7, 1.4, "heat out " + r"($q_{CD} < 0$)", color="b", ha="right", va="center")

ax.set_xlim(0.6, 6.1)
ax.set_ylim(0, 46)
ax.set_xticks([corners[c][0] for c in "ADBC"])
ax.set_xticklabels([r"$V_A$", r"$V_D$", r"$V_B$", r"$V_C$"])
ax.set_yticks([corners[c][1] for c in "ABDC"])
ax.set_yticklabels([r"$P_A$", r"$P_B$", r"$P_D$", r"$P_C$"])
ax.set_xlabel(r"Volume $V$")
ax.set_ylabel(r"Pressure $P$")
ax.grid(alpha=0.4)
ax.legend(loc="upper right", title="Steps", fontsize=8.5)

plt.show()
plt.close(fig)
```

The Carnot cycle for one mole of a monatomic ideal gas between $T_{\mathrm{hot}} = 500\ \mathrm{K}$ and $T_{\mathrm{cold}} = 300\ \mathrm{K}$. The arrows on the curves give the direction of each step. Heat enters the system on the hot isotherm and leaves it on the cold isotherm. The adiabats, which are steeper than the isotherms, exchange no heat.

#### Step-by-Step Thermodynamic Description

| Step | Process type | Thermal contact | Heat $q$ | Entropy change $\Delta S$ | Qualitative work |
| --- | --- | --- | --- | --- | --- |
| $A\to B$ | Reversible isothermal *expansion* at $T_{\mathrm{hot}}$ | Hot reservoir | $q_{AB} > 0$ | $\Delta S_{AB} = q_{AB}/T_{\mathrm{hot}}$ | System does work on surroundings ($w_{AB} < 0$) |
| $B\to C$ | Reversible adiabatic *expansion* | Insulated | $q_{BC}=0$ | $\Delta S_{BC}=0$ | System does work; temperature drops ($w_{BC} < 0$) |
| $C\to D$ | Reversible isothermal *compression* at $T_{\mathrm{cold}}$ | Cold reservoir | $q_{CD} < 0$ | $\Delta S_{CD} = q_{CD}/T_{\mathrm{cold}}$ | Surroundings do work on system ($w_{CD} > 0$) |
| $D\to A$ | Reversible adiabatic *compression* | Insulated | $q_{DA}=0$ | $\Delta S_{DA}=0$ | Surroundings do work on system; $T$ rises ($w_{DA} > 0$) |

Two quick consequences:

- Because adiabatic steps have $q=0$, the reversible definition $dS = \delta q_{\mathrm{rev}}/T$ gives
  $\Delta S = 0$ on the adiabats (**entropy is constant** on a reversible adiabat — such a process is called *isentropic*). Compare the free expansion in Example 2 of [Section 4.1](section-01.md), which is adiabatic but not reversible.
- Over one full cycle, the system returns to its initial state, so

  ```{math}
  \Delta S_{\mathrm{cycle}}
  = \Delta S_{AB} + \Delta S_{BC} + \Delta S_{CD} + \Delta S_{DA}
  = \Delta S_{AB} + \Delta S_{CD}
  = 0,
  ```

  which is exactly the relation we should expect: $S$ is a state function, so its change over any closed cycle must vanish.

---

### How to Operate a Carnot Engine

The Carnot cycle can be realized conceptually with a gas-filled cylinder, a piston loaded with weights, and two thermal reservoirs. The core idea: you alternate between (i) placing the gas in contact with a thermal reservoir so it can exchange heat *reversibly* at a fixed temperature (isotherms), and (ii) insulating the gas so it cannot exchange heat (adiabats).

1. **$A\to B$: isothermal expansion at $T_{\mathrm{hot}}$**
   Put the cylinder in contact with the hot reservoir. Remove weights slowly so the piston rises
   quasi-statically. The system **absorbs heat** $q_{AB}$ from the hot reservoir to keep $T$ constant,
   while doing work on the weights.

2. **$B\to C$: adiabatic expansion (insulated)**
   Remove the thermal contact and insulate the cylinder. Continue removing weights slowly so the piston rises.
   No heat flows ($q=0$), so the work done by the gas comes from its internal energy, and the temperature drops
   from $T_{\mathrm{hot}}$ to $T_{\mathrm{cold}}$.

3. **$C\to D$: isothermal compression at $T_{\mathrm{cold}}$**
   Put the cylinder in contact with the cold reservoir. Add weights slowly so the piston lowers.
   The surroundings do work on the gas, and the gas **releases heat** $q_{CD}$ to the cold reservoir to keep $T$ constant.

4. **$D\to A$: adiabatic compression (insulated)**
   Insulate the cylinder again and keep adding weights slowly to compress the gas. With no heat exchange,
   the work done on the gas increases its internal energy and raises the temperature back to $T_{\mathrm{hot}}$,
   returning to state $A$.

The schematic below illustrates each step, showing the thermal contact, piston position, and energy flows at the end of each process.

```{code-cell} ipython3
:tags: [hide-input]

import matplotlib.pyplot as plt

fig, axs = plt.subplots(1, 4, figsize=(8, 4), constrained_layout=True, sharex=True, sharey=True)

# Each panel shows the END state of the corresponding step.
titles = [r"Step 1: A$\to$B", r"Step 2: B$\to$C", r"Step 3: C$\to$D", r"Step 4: D$\to$A"]
end_state = ["B", "C", "D", "A"]
piston_height = [3.5, 4.5, 2.5, 1.5]   # schematic, not to scale
contact = ["hot reservoir\n" + r"$T_{\mathrm{hot}}$", "insulation",
           "cold reservoir\n" + r"$T_{\mathrm{cold}}$", "insulation"]
contact_color = ["mistyrose", "0.88", "lightcyan", "0.88"]

for i, ax in enumerate(axs):
    h = piston_height[i]
    ax.set_title(titles[i], fontsize=10)

    # Cylinder walls and base
    ax.plot([0, 0, 3, 3], [5.6, 0, 0, 5.6], "k-", lw=1.5)

    # Piston and rod
    ax.plot([0, 3], [h, h], "k-", lw=3)
    ax.plot([1.5, 1.5], [h, 6.2], "k-", lw=2)
    ax.text(0.7, h + 0.3, "piston", fontsize=9, ha="center", va="bottom")

    # State of the gas at the end of the step
    ax.text(2.0, 0.6, rf"$V_{end_state[i]},\ P_{end_state[i]}$", fontsize=10, ha="center", va="center")

    # Thermal contact during the step
    ax.add_patch(plt.Rectangle((0, -2.4), 3, 2.0, facecolor=contact_color[i], edgecolor="k",
                               linestyle="-" if i % 2 == 0 else ":"))
    ax.text(1.5, -1.4, contact[i], fontsize=9, ha="center", va="center")

    arrow = dict(arrowstyle="-|>", lw=1.8, mutation_scale=14)
    if i == 0:   # heat enters the gas; the gas does work
        ax.annotate("", xy=(0.6, 1.0), xytext=(0.6, -0.6), arrowprops=dict(color="r", **arrow))
        ax.text(0.6, 1.1, "heat in", fontsize=9, color="r", ha="center", va="bottom")
        ax.annotate("", xy=(2.3, h + 1.5), xytext=(2.3, h + 0.2), arrowprops=dict(color="b", **arrow))
        ax.text(2.3, h + 1.6, "work done\nby gas", fontsize=9, color="b", ha="center", va="bottom")
    if i == 2:   # heat leaves the gas; work is done on the gas
        ax.annotate("", xy=(0.6, -0.6), xytext=(0.6, 1.0), arrowprops=dict(color="r", **arrow))
        ax.text(0.6, 1.1, "heat out", fontsize=9, color="r", ha="center", va="bottom")
        ax.annotate("", xy=(2.3, h + 0.2), xytext=(2.3, h + 1.5), arrowprops=dict(color="b", **arrow))
        ax.text(2.3, h + 1.6, "work done\non gas", fontsize=9, color="b", ha="center", va="bottom")

    ax.set_xlim(-0.3, 3.3)
    ax.set_ylim(-2.6, 6.3)
    ax.axis("off")

plt.show()
plt.close(fig)
```

Schematic of the four Carnot-cycle steps. Each panel shows the thermal contact during that step and the piston position at its endpoint (not to scale). Heat exchange occurs only during the isothermal steps (Steps 1 and 3); the adiabatic steps (Steps 2 and 4) are insulated. Work is exchanged in all four steps, and the arrows show it for Steps 1 and 3.

---

### Maximum Efficiency of Converting Heat to Work (Carnot Efficiency)

We define the heat quantities on the two isothermal legs from the system's perspective, following the course sign convention ($q > 0$ when heat is absorbed by the system):

- $q_{AB} > 0$: heat **absorbed by the system** from the hot reservoir during step $A \to B$
- $q_{CD} < 0$: heat **released by the system** to the cold reservoir during step $C \to D$

For one complete cycle, $\Delta U_{\mathrm{cycle}} = 0$ (state function returns to its initial value), so the First Law gives

```{math}
0 = (q_{AB} + q_{CD}) + w_{\mathrm{cycle}},
```

where $w_{\mathrm{cycle}}$ is the total work over the cycle in the course convention ($w > 0$ when work is done *on* the system). Since the engine does net work *on the surroundings*, $w_{\mathrm{cycle}} < 0$, and the net work *output* is

```{math}
W_{\mathrm{out}} \equiv -w_{\mathrm{cycle}} = q_{AB} + q_{CD}.
```

The thermal efficiency is defined as the fraction of absorbed heat converted to useful work output:

```{math}
\eta
= \frac{W_{\mathrm{out}}}{q_{AB}}
= \frac{q_{AB} + q_{CD}}{q_{AB}}.
```

Now use the entropy definition on the isotherms:

```{math}
\Delta S_{AB} = \frac{q_{AB}}{T_{\mathrm{hot}}},\qquad
\Delta S_{CD} = \frac{q_{CD}}{T_{\mathrm{cold}}}.
```

Because the full reversible cycle returns the system to its initial state,

```{math}
\Delta S_{\mathrm{cycle}} = \Delta S_{AB} + \Delta S_{CD} = 0
\quad\Rightarrow\quad
\Delta S_{CD} = -\Delta S_{AB}.
```

Substitute $q_{AB} = T_{\mathrm{hot}}\,\Delta S_{AB}$ and $q_{CD} = T_{\mathrm{cold}}\,\Delta S_{CD} = -T_{\mathrm{cold}}\,\Delta S_{AB}$ into the efficiency:

```{math}
\eta
= \frac{T_{\mathrm{hot}}\Delta S_{AB} - T_{\mathrm{cold}}\Delta S_{AB}}{T_{\mathrm{hot}}\Delta S_{AB}}
= 1 - \frac{T_{\mathrm{cold}}}{T_{\mathrm{hot}}}.
```

**Interpretation:** the maximum possible fraction of hot-reservoir heat that can be converted to work depends
*only* on the two reservoir temperatures.

The derivation above gives the efficiency of a fully reversible cycle. [Section 4.3](section-03.md) uses the Clausius inequality to show that no engine operating between the same two reservoirs can have a higher efficiency.

---

### The Second Law: Entropy and the Direction of Heat Flow

We can use the fundamental relation from [Section 4.1](section-01.md) to derive the direction of spontaneous heat flow. Consider a simple composite system:

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

At fixed volume, the fundamental relation $dU = T\,dS - P\,dV$ reduces to $dU = T\,dS$. The relation applies to each subsystem even though the exchange of heat between them is irreversible, because it contains only state functions ([Section 4.1](section-01.md)). Therefore

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

**Second-law statement for an isolated system:** spontaneous change requires $dS \ge 0$, with equality for a reversible (equilibrium) exchange.

Now suppose $T_A > T_B$ ($A$ is hotter than $B$). Then
$\tfrac{1}{T_A} - \tfrac{1}{T_B} < 0$.
To satisfy $dS \ge 0$, we must have $dU_A \le 0$, meaning subsystem $A$ loses internal energy while $B$ gains it.

That is exactly the familiar conclusion:

> **Heat flows spontaneously from hot to cold.**

Attempting to make heat flow from cold to hot without any other change would make $dS < 0$ for the isolated
composite system, which violates the second law.

```{admonition} When Does the Flow of Heat Stop?
:class: important
The expression for $dS$ contains the factor $(1/T_A - 1/T_B)$. Follow this factor, and with it $dS$, as heat flows and the two temperatures approach each other. [Section 4.3](section-03.md) returns to this point when it discusses the approach to equilibrium.
```

---

## Worked Example

### Carnot Efficiency and Reversible Entropy Bookkeeping

An ideal Carnot engine operates between $T_{\mathrm{hot}}=500\ \mathrm{K}$ and $T_{\mathrm{cold}}=300\ \mathrm{K}$. Suppose the system absorbs $q_{AB}=1000\ \mathrm{J}$ from the hot reservoir during the isothermal expansion step.

1. **Maximum efficiency**

   ```{math}
   \eta_{\mathrm{Carnot}} = 1-\frac{T_{\mathrm{cold}}}{T_{\mathrm{hot}}}
   =1-\frac{300}{500}=0.40.
   ```

   So $40\%$ of the absorbed heat is converted to net work output:

   ```{math}
   W_{\mathrm{out}} = \eta\,q_{AB} = 0.40 \times 1000 = 400\ \mathrm{J}.
   ```

   In the course sign convention ($w > 0$ = work done *on* the system), the engine does work on the surroundings, so

   ```{math}
   w_{\mathrm{cycle}} = -W_{\mathrm{out}} = -400\ \mathrm{J}.
   ```

2. **Heat rejected**

   Energy conservation over the cycle ($\Delta U = 0$) gives $q_{AB} + q_{CD} + w_{\mathrm{cycle}} = 0$:

   ```{math}
   q_{CD} = -q_{AB} - w_{\mathrm{cycle}} = -1000 - (-400) = -600\ \mathrm{J}.
   ```

   The negative sign confirms that the system releases 600 J to the cold reservoir, as expected.

3. **Entropy changes of the reservoirs**

   Each reservoir exchanges heat reversibly at its own fixed temperature. From the *reservoir's* perspective, the heat it receives is opposite in sign to what the system exchanges (the system's gain is the reservoir's loss):

   ```{math}
   \Delta S_{\mathrm{hot}} = \frac{-q_{AB}}{T_{\mathrm{hot}}}= \frac{-1000}{500}=-2.0\ \mathrm{J/K},
   ```

   ```{math}
   \Delta S_{\mathrm{cold}} = \frac{-q_{CD}}{T_{\mathrm{cold}}}= \frac{-(-600)}{300}=+2.0\ \mathrm{J/K}.
   ```

   The hot reservoir loses entropy (it gave up heat); the cold reservoir gains entropy (it absorbed heat).

**Result.** $\Delta S_{\mathrm{hot}}+\Delta S_{\mathrm{cold}}=0$, and the engine itself returns to its initial state ($\Delta S_{\mathrm{engine}} = 0$), so $\Delta S_{\mathrm{total}} = 0$. This is the hallmark of a fully reversible process. Any real (irreversible) engine would produce $\Delta S_{\mathrm{total}} > 0$, meaning less work output for the same heat input.

## Concept Checks

1. Why are adiabatic steps isentropic only in the *reversible* case?
2. What feature of the Carnot cycle makes it "best possible" compared to real engines?
3. If $T_{\mathrm{cold}}$ is lowered while $T_{\mathrm{hot}}$ is fixed, how does $\eta_{\mathrm{Carnot}}$ change and why?
4. In the two-body heat-flow argument, what happens to $dS$ as $T_A \to T_B$? What does this tell you about equilibrium?
5. A real engine operating between $T_{\mathrm{hot}} = 500\ \mathrm{K}$ and $T_{\mathrm{cold}} = 300\ \mathrm{K}$ absorbs 1000 J of heat and produces 350 J of work output. Compute $\Delta S_{\mathrm{total}}$ for the reservoirs and verify that $\Delta S_{\mathrm{total}} > 0$.

## Key Takeaways

- A Carnot engine is fully reversible. No engine operating between the same two reservoirs has a higher efficiency ([Section 4.3](section-03.md)).
- Reversible isotherms exchange heat at fixed $T$; reversible adiabats exchange no heat and keep $S$ constant (isentropic).
- The Carnot efficiency depends only on reservoir temperatures: $\eta_{\mathrm{Carnot}}=1-T_{\mathrm{cold}}/T_{\mathrm{hot}}$.
- The entropy bookkeeping rule $dS = \delta q_{\mathrm{rev}}/T$ applied to a cycle yields the efficiency formula directly.
- Second-law reasoning ($dS \geq 0$ for an isolated system) implies that spontaneous heat flow goes from hot to cold.
