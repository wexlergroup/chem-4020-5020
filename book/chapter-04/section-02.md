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
  enumerator: 4.2.%s
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

[](#fig:carnot-cycle) plots the Carnot cycle on a **pressure–volume ($P$–$V$) diagram** for one mole of a monatomic ideal gas between reservoirs at $500\ \mathrm{K}$ and $300\ \mathrm{K}$, the temperatures of the Worked Example below.
The key geometric features are:

- The **isotherms** (orange at $T_{\mathrm{hot}}$, blue at $T_{\mathrm{cold}}$) follow the ideal-gas relation $P = nRT/V$ at fixed temperature. The thin lines continue them beyond the cycle.
- The **adiabats** (black) are *steeper* than the isotherms and follow $PV^\gamma = \text{constant}$
  for a reversible adiabatic process of an ideal gas ([Section 3.2](../chapter-03/section-02.md)).

```{code-cell} ipython3
:label: fig:carnot-cycle
:caption: The Carnot cycle for one mole of a monatomic ideal gas between $T_{\mathrm{hot}} = 500\ \mathrm{K}$ and $T_{\mathrm{cold}} = 300\ \mathrm{K}$. The arrows on the curves give the direction of each step. Heat enters the system on the hot isotherm and leaves it on the cold isotherm. The adiabats, which are steeper than the isotherms, exchange no heat. The thin lines continue the two isotherms beyond the cycle.
:tags: [hide-input]

import sys

import matplotlib.pyplot as plt
import numpy as np
from scipy.constants import R as R_J

sys.path.append("..")            # the book folder, which holds figstyle
import figstyle as fs

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


def draw_cycle(c):
    fig, ax = plt.subplots(figsize=(6.4, 4.1), constrained_layout=True)
    steps = [(V_A, V_B, lambda V: P_isotherm(V, T_hot), c.vermilion),
             (V_B, V_C, lambda V: P_adiabat(V, V_B, T_hot), c.ink),
             (V_C, V_D, lambda V: P_isotherm(V, T_cold), c.blue),
             (V_D, V_A, lambda V: P_adiabat(V, V_A, T_hot), c.ink)]

    # The two isotherms continue beyond the cycle as thin lines
    V_all = np.linspace(0.84, 6.2, 400)
    ax.plot(V_all, P_isotherm(V_all, T_hot), color=c.tint(c.vermilion, 0.55), lw=0.9)
    ax.plot(V_all, P_isotherm(V_all, T_cold), color=c.tint(c.blue, 0.55), lw=0.9)

    # The four steps, each with an arrow that shows its direction
    for V1, V2, P_of_V, color in steps:
        V = np.linspace(V1, V2, 200)
        ax.plot(V, P_of_V(V), color=color, lw=fs.LINE)
        Vm, dV = np.sqrt(V1 * V2), 0.02 * (V2 - V1)
        ax.annotate("", xy=(Vm + dV, P_of_V(Vm + dV)), xytext=(Vm - dV, P_of_V(Vm - dV)),
                    arrowprops=dict(arrowstyle="-|>", color=color, lw=0, mutation_scale=17))

    # The four states
    corners = {"A": (V_A, P_isotherm(V_A, T_hot)), "B": (V_B, P_isotherm(V_B, T_hot)),
               "C": (V_C, P_isotherm(V_C, T_cold)), "D": (V_D, P_isotherm(V_D, T_cold))}
    offsets = {"A": (10, 2), "B": (6, 9), "C": (2, 11), "D": (-8, -9)}
    for name, (Vc, Pc) in corners.items():
        ax.plot(Vc, Pc, "o", ms=6.5, mfc=c.ground, mec=c.ink, mew=1.4, zorder=5)
        ax.annotate(name, (Vc, Pc), textcoords="offset points", xytext=offsets[name], ha="center", va="center",
                    fontweight="bold")

    # Each step is named beside its curve
    ax.text(1.72, 29.5, "isothermal expansion at $T_{\\mathrm{hot}}$\nheat in, $q_{AB} > 0$", color=c.vermilion,
            ha="left", va="center")
    ax.text(3.35, 15.6, "adiabatic expansion, $q_{BC} = 0$", ha="left", va="center")
    ax.text(3.6, 2.4, "isothermal compression at $T_{\\mathrm{cold}}$\nheat out, $q_{CD} < 0$", color=c.blue,
            ha="center", va="center")
    ax.text(1.46, 11.8, "adiabatic\ncompression\n$q_{DA} = 0$", ha="right", va="center")

    ax.set(xlim=(0.6, 6.1), ylim=(0, 46), xlabel="Volume $V$", ylabel="Pressure $P$")
    ax.set_xticks([corners[k][0] for k in "ADBC"], ["$V_A$", "$V_D$", "$V_B$", "$V_C$"])
    ax.set_yticks([corners[k][1] for k in "ABDC"], ["$P_A$", "$P_B$", "$P_D$", "$P_C$"])
    return fig


fs.show(draw_cycle, alt="Pressure against volume for a Carnot cycle with the states A, B, C, and D. From A to B the gas "
                        "expands along the hot isotherm and takes in heat. From B to C it expands along a steeper "
                        "adiabat to the largest volume. From C to D it is compressed along the cold isotherm and gives "
                        "off heat. From D to A it is compressed along a second adiabat back to the highest pressure.")
```

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

[](#fig:carnot-engine) illustrates each step, showing the thermal contact, the position of the piston, and the flows of energy.

```{code-cell} ipython3
:label: fig:carnot-engine
:caption: Schematic of the four steps of the Carnot cycle (not to scale). Each panel shows the thermal contact during the step and the piston at the end of the step. The dashed line is the position of the piston at the start of the step. The number of weights on the piston stands for the pressure of the gas. Heat is exchanged only in the isothermal steps, (a) and (c), and the adiabatic steps, (b) and (d), are insulated. Work is exchanged in all four steps.
:tags: [hide-input]

import sys

import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

sys.path.append("..")            # the book folder, which holds figstyle
import figstyle as fs

HEIGHT = {"A": 1.3, "D": 2.5, "B": 3.4, "C": 4.7}          # schematic height of the piston in each state
WEIGHTS = {"A": 4, "B": 3, "D": 2, "C": 1}                 # more weights on the piston at a higher pressure


def draw_step(ax, c, start, end, contact, heat):
    """One step of the cycle: the cylinder with its piston at the end of the step, on a reservoir or on insulation."""
    h0, h1 = HEIGHT[start], HEIGHT[end]
    expansion = h1 > h0

    # Thermal contact below the cylinder
    face = {"hot": c.tint(c.vermilion, 0.78), "cold": c.tint(c.blue, 0.78), "insulation": c.fill}[contact]
    ax.add_patch(Rectangle((0, -2.2), 3, 1.8, fc=face, ec=c.muted if contact == "insulation" else "none",
                           lw=1.0, ls=(0, (4, 3))))
    label = {"hot": "hot reservoir\n$T_{\\mathrm{hot}}$", "cold": "cold reservoir\n$T_{\\mathrm{cold}}$",
             "insulation": "insulation"}[contact]
    ax.text(1.5, -1.3, label, ha="center", va="center")

    # Cylinder, piston at the start (dashed) and at the end of the step, and the weights on the piston
    ax.plot([0, 0, 3, 3], [5.9, 0, 0, 5.9], color=c.ink, lw=1.4, solid_joinstyle="round")
    ax.plot([0.08, 2.92], [h0, h0], color=c.gray, lw=1.2, ls=(0, (4, 3)))
    ax.plot([0.05, 2.95], [h1, h1], color=c.ink, lw=3.2, solid_capstyle="butt")
    for k in range(WEIGHTS[end]):
        ax.add_patch(Rectangle((0.2 + 0.34 * k, h1 + 0.09), 0.26, 0.3, fc=c.gray, ec="none"))

    # Work: the piston moves from the dashed line to the solid line
    gap = 0.14 if expansion else -0.14
    fs.arrow(ax, (2.45, h0 + gap), (2.45, h1 - 1.6 * gap), c.ink, lw=1.6)
    ax.text(2.2, max(h0, h1) + 0.25, "work done\nby the gas" if expansion else "work done\non the gas",
            ha="center", va="bottom")

    # Heat crosses the base of the cylinder in the isothermal steps
    if heat == "in":
        fs.arrow(ax, (0.4, -0.3), (0.4, 0.8), c.vermilion, lw=1.6)
    if heat == "out":
        fs.arrow(ax, (0.4, 0.8), (0.4, -0.3), c.vermilion, lw=1.6)
    if heat:
        ax.text(0.58, 0.45, "heat " + heat, ha="left", va="center", color=c.vermilion)

    # State of the gas at the end of the step
    ax.text(2.3, 0.45, f"$V_{end}$, $P_{end}$", ha="center", va="center")

    ax.set(xlim=(-0.25, 3.25), ylim=(-2.3, 6.3))
    ax.axis("off")


def draw_engine(c):
    fig, axes = plt.subplots(1, 4, figsize=(fs.WIDTH, 3.3), constrained_layout=True)
    draw_step(axes[0], c, "A", "B", "hot", "in")
    draw_step(axes[1], c, "B", "C", "insulation", None)
    draw_step(axes[2], c, "C", "D", "cold", "out")
    draw_step(axes[3], c, "D", "A", "insulation", None)
    for ax, letter, title in zip(axes, "abcd", ["A → B", "B → C", "C → D", "D → A"]):
        fs.panel(ax, letter, title)
    return fig


fs.show(draw_engine, alt="Four sketches of a gas in a cylinder under a piston that carries weights. In panel a, from A to B, the cylinder "
                         "stands on the hot reservoir, heat enters the gas, and the piston rises as the gas does work. "
                         "In panel b, from B to C, the cylinder stands on insulation, and the piston rises further. In "
                         "panel c, from C to D, the cylinder stands on the cold reservoir, heat leaves the gas, and "
                         "the piston is pushed down as work is done on the gas. In panel d, from D to A, the cylinder "
                         "stands on insulation, and the piston is pushed down to its lowest position.")
```

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

[](#fig:carnot-energy-flow) shows the three energies of this definition for one cycle.

```{code-cell} ipython3
:label: fig:carnot-energy-flow
:caption: Energy flow through a Carnot engine in one cycle. The widths of the three streams are drawn to scale for the engine of the Worked Example below, which runs between $T_{\mathrm{hot}} = 500\ \mathrm{K}$ and $T_{\mathrm{cold}} = 300\ \mathrm{K}$. Of every $1000\ \mathrm{J}$ of heat that the engine absorbs from the hot reservoir, $400\ \mathrm{J}$ leaves as net work and $600\ \mathrm{J}$ leaves as heat to the cold reservoir. Because $\Delta U_{\mathrm{cycle}} = 0$, the width of the stream that enters equals the sum of the widths of the two streams that leave.
:tags: [hide-input]

import sys

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Polygon, Rectangle

sys.path.append("..")            # the book folder, which holds figstyle
import figstyle as fs

T_HOT, T_COLD = 500.0, 300.0                 # K, the reservoirs of the Worked Example
Q_IN = 1000.0                                # J, heat absorbed from the hot reservoir in one cycle
W_OUT = (1 - T_COLD / T_HOT) * Q_IN          # J, net work output
Q_OUT = Q_IN - W_OUT                         # J, heat released to the cold reservoir

WIDTH_PER_JOULE = 2.0 / Q_IN                 # drawn width of a stream per joule
X_LEFT = 1.4                                 # left edge of the stream of heat
Y_SPLIT = 3.6                                # height at which the stream of work turns aside
R_INNER = 0.3                                # inner radius of that turn
HEAD = 0.45                                  # length of the pointed end of a stream


def draw_flow(c):
    fig, ax = plt.subplots(figsize=(5.2, 4.5), constrained_layout=True)
    w_in, w_work, w_out = (WIDTH_PER_JOULE * e for e in (Q_IN, W_OUT, Q_OUT))
    heat, work = c.tint(c.vermilion, 0.3), c.tint(c.ink, 0.45)

    # Reservoirs and engine
    ax.add_patch(Rectangle((0, 4.8), 4.4, 0.7, fc=c.tint(c.vermilion, 0.78), ec="none"))
    ax.add_patch(Rectangle((0, 0), 4.4, 0.7, fc=c.tint(c.blue, 0.78), ec="none"))
    ax.add_patch(FancyBboxPatch((0.12, 1.72), 4.16, 2.16, boxstyle="round,pad=0.12", fc=c.fill, ec=c.ink, lw=1.4))
    ax.text(2.2, 5.15, "hot reservoir, $T_{\\mathrm{hot}}$", ha="center", va="center")
    ax.text(2.2, 0.35, "cold reservoir, $T_{\\mathrm{cold}}$", ha="center", va="center")
    ax.text(0.7, 2.8, "engine", ha="center", va="center")

    # Heat: from the hot reservoir down to the split, and on to the cold reservoir
    ax.add_patch(Rectangle((X_LEFT, Y_SPLIT), w_in, 4.8 - Y_SPLIT, fc=heat, ec="none"))
    ax.add_patch(Polygon([(X_LEFT, Y_SPLIT), (X_LEFT + w_out, Y_SPLIT), (X_LEFT + w_out, 0.7 + HEAD),
                          (X_LEFT + w_out / 2, 0.7), (X_LEFT, 0.7 + HEAD)], fc=heat, ec="none"))

    # Work: turns aside in a quarter circle and leaves the engine to the right
    x_c, r_outer = X_LEFT + w_in + R_INNER, R_INNER + w_work
    angle = np.linspace(np.pi, 1.5 * np.pi, 40)
    outer = np.column_stack([x_c + r_outer * np.cos(angle), Y_SPLIT + r_outer * np.sin(angle)])
    inner = np.column_stack([x_c + R_INNER * np.cos(angle), Y_SPLIT + R_INNER * np.sin(angle)])
    y_low, y_high, x_tip = Y_SPLIT - r_outer, Y_SPLIT - R_INNER, 6.2
    end = [(x_tip - HEAD, y_low), (x_tip, (y_low + y_high) / 2), (x_tip - HEAD, y_high)]
    ax.add_patch(Polygon(np.vstack([outer, end, inner[::-1]]), fc=work, ec="none"))

    # Labels beside the streams
    ax.text(X_LEFT + w_in + 0.15, 4.4, "heat in\n$q_{AB}$", ha="left", va="center", color=c.vermilion)
    ax.text(X_LEFT + w_out + 0.15, 1.2, "heat out\n$|q_{CD}|$", ha="left", va="center", color=c.vermilion)
    ax.text(5.3, y_high + 0.12, "net work out\n$W_{\\mathrm{out}}$", ha="center", va="bottom")

    # Energies of the Worked Example, written on the streams
    ax.text(X_LEFT + w_in / 2, 4.4, f"{Q_IN:.0f} J", ha="center", va="center", color=fs.on(heat))
    ax.text(X_LEFT + w_out / 2, 1.45, f"{Q_OUT:.0f} J", ha="center", va="center", color=fs.on(heat))
    ax.text(5.0, (y_low + y_high) / 2, f"{W_OUT:.0f} J", ha="center", va="center", color=fs.on(work))

    ax.set(xlim=(-0.1, 6.3), ylim=(-0.05, 5.55))
    ax.set_aspect("equal")
    ax.axis("off")
    return fig


fs.show(draw_flow, alt="Flow diagram of a Carnot engine. A wide stream of heat, q AB, runs from the hot reservoir at the top "
                       "into the engine. Inside the engine the stream divides. Six tenths of its width continue down "
                       "to the cold reservoir as heat out, and four tenths turn to the right and leave the engine as "
                       "net work out.")
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

   The negative sign confirms that the system releases 600 J to the cold reservoir, as expected. [](#fig:carnot-energy-flow) shows these three energies to scale.

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
