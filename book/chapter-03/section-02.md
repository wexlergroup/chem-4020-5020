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
  enumerator: 3.2.%s
---


# 3.2. Applications of the First Law

See also: [Course-wide Conventions & Notation](../notation.md)

Run the code: [Notebook for Section 3.2](https://wexlergroup.github.io/chem-4020-5020/lite/notebooks/index.html?path=chapter-03/section-02.ipynb), in which you can change and rerun the Python cells of this section in your browser.

## Overview and Learning Objectives

If you heat a gas in a sealed rigid container versus a piston open to the atmosphere, the temperature changes will differ—even for the same amount of heat. Why? The answer lies in how constraints (constant $V$ vs. constant $P$) redirect energy between heat and work. This section develops a systematic workflow for applying the First Law under different process constraints, then illustrates it with common ideal-gas cases and a microscopic interpretation.

We begin by clarifying the definitions of quasi-static, reversible, and irreversible processes. We then outline a step-by-step method for using the First Law under different constraints, derive results for isochoric, isothermal, and adiabatic ideal-gas processes, and end with a microscopic interpretation that connects probability distributions of microstates with heat and work. Processes at constant pressure are taken up in [Section 3.3](section-03.md).

---

Learning objectives:

- Differentiate quasi-static, reversible, and irreversible processes and explain why reversibility is an idealization.
- Explain why reversible expansion produces the maximum work for a given volume change (and reversible compression requires the minimum work).
- Rewrite the First Law in terms of chosen independent variables and apply constraints (isothermal, isochoric, isobaric, adiabatic).
- Integrate $PV$ work and heat for common ideal-gas processes and relate results to state-variable changes.
- Interpret microscopic changes in $p_i$ and $E_i$ as heat and work contributions to $dU$.

## Core Ideas and Derivations

### Thermodynamic Processes

```{glossary}
Quasi-static process
: A process carried out *infinitesimally slowly* so that the system remains *nearly in equilibrium* at all times. Each intermediate state is well-defined thermodynamically.

Reversible process
: An *idealized* process that is (1) *quasi-static* and (2) free of *dissipative effects*, meaning it can be reversed with *no net change* to the system or surroundings.

Irreversible process
: Any process that violates one or more of the conditions for reversibility (e.g., too rapid, frictional, or dissipative).

Dissipative effects
: Processes that convert ordered energy into disordered (thermal) energy in a way that cannot be undone spontaneously. Common examples include friction, viscous drag, turbulence, electrical resistance, and heat flow across a finite temperature difference.
```

```{admonition} Examples: Reversible vs. irreversible
:class: tip

**Approximately reversible:** Slowly compressing a gas in a well-lubricated, frictionless piston, one infinitesimal step at a time. At each step, the system is in mechanical equilibrium with the surroundings ($P \approx P_{\mathrm{ext}}$), and the process can be reversed by an infinitesimal decrease in external pressure.

**Irreversible:** Puncturing a membrane so that a gas undergoes free expansion into vacuum (as in [Section 3.1](section-01.md)). The gas is *never* in equilibrium during the expansion—there is no well-defined pressure for the system as a whole—and no infinitesimal change to the surroundings will reverse the process.

**Key practical consequence:** For a given initial and final state, a reversible expansion does the *maximum* work on the surroundings ($|w_{\mathrm{rev}}| \geq |w_{\mathrm{irrev}}|$ for expansion). This is because, at every infinitesimal step, the system pushes against the highest possible opposing pressure ($P_{\mathrm{ext}} = P$), maximizing $\int P_{\mathrm{ext}}\,dV$.
```

---

### How to Apply the First Law of Thermodynamics

The First Law states ([Section 3.1](section-01.md), Eq. {eq}`first-law-differential`):

```{math}
dU = \delta q + \delta w,
```

where $dU$ is the change in internal energy, $\delta q$ is the heat added to the system, and $\delta w$ is the work done on the system. To make this law useful for specific processes, follow the five-step workflow below.

**Workflow overview:** (1) Choose independent variables → (2) Rewrite the First Law → (3) Apply constraints → (4) Define the system/EOS → (5) Integrate. Click each step for details.

```{admonition} Step 1. Choose Two Independent Variables
:class: dropdown

For a single-component, single-phase system, specifying two independent variables (plus the amount of substance) completely determines the thermodynamic state.

- **Convenience**: Which variables are easiest to control with your experimental setup?
- **Necessity**: Which variables can you reliably measure or keep constant?
- **Curiosity**: Which variables are most relevant to the phenomenon you wish to explore?

Common choices include $(V, T)$ for gases in rigid containers and $(P, T)$ for processes open to the atmosphere.
```

````{admonition} Step 2. Rewrite the First Law in Terms of Your Chosen Variables
:class: dropdown

Start from $\delta q = dU - \delta w$. For a reversible process with only $PV$ work, $\delta w = -P\,dV$ ([Section 3.1](section-01.md), Eq. {eq}`pv-work`), so

```{math}
\delta q = dU + P\,dV.
```

Now expand $dU$ as a total differential in your chosen variables. For example, with $V$ and $T$:

```{math}
dU = \left(\frac{\partial U}{\partial T}\right)_V dT 
     + \left(\frac{\partial U}{\partial V}\right)_T dV,
```

which gives

```{math}
\delta q = \left(\frac{\partial U}{\partial T}\right)_V dT 
     + \left[\left(\frac{\partial U}{\partial V}\right)_T + P\right] dV.
```

Note that this expression for $\delta q$ has *three* contributions: the response of $U$ to temperature, the response of $U$ to volume, *and* the $P\,dV$ work term. This is not just the total differential of $U$—the work term matters.
````

```{admonition} Step 3. Apply Constraints
:class: dropdown

- **Constant state variables**:
  - Isobaric ($dP=0$), Isochoric ($dV=0$), Isothermal ($dT=0$).
- **Adiabatic**:
  - If the boundary is thermally insulating, $\delta q = 0$.
  - "Adiabatic" means no heat transfer, but the system may still perform or receive work.

Applying a constraint eliminates one or more terms from your Step 2 expression, often making the resulting differential equation separable and integrable.
```

```{admonition} Step 4. Define the System (Equation of State)
:class: dropdown

- Specify an equation of state if known (e.g., $PV=nRT$ for an ideal gas, or the van der Waals EOS from [Section 1.4](../chapter-01/section-04.md)) and an expression for $U$ (e.g., $U = \tfrac{3}{2}Nk_{\mathrm{B}}T$ for a monatomic ideal gas, from [Section 2.5](../chapter-02/section-05.md)).
- Use them to evaluate $P$ and the partial derivatives of $U$ that appear in your expression for $\delta q$.
- This step is where the *physics* of the system enters: different substances yield different results for the same constraint.
```

```{admonition} Step 5. Integrate the First Law
:class: dropdown

- Integrate the differential relationships to link total heat ($q$) or work ($w$) to finite changes in state variables.
- For reversible processes, the integration path in state space is well-defined.
- For irreversible processes, compute $\Delta U$ from state variables (it's path-independent) and use the First Law to relate $q$ and $w$.
```

```{admonition} Why It Matters
:class: tip

Different process constraints (isobaric, isochoric, isothermal, adiabatic) dictate how heat, work, and changes in state variables interrelate. Mastering these constraints reveals how energy flows and how measurable quantities (e.g., temperature, volume, pressure) change under specific conditions.
```

[](#fig:process-paths) shows the four constraints as paths of an ideal gas on a $P$–$V$ diagram. All four paths start from the same state.

% Try this: Change `V_2` from `20.0` to `10.0`, so that every expansion stops at half the final volume. Which of the three paths keeps the largest fraction of its work, and why?

```{code-cell} ipython3
:label: fig:process-paths
:caption: Four reversible paths of one mole of a monatomic ideal gas from the same initial state, $(V_1, T_1) = (5\,\mathrm{L}, 300\,\mathrm{K})$: isochoric heating, and isobaric, isothermal, and adiabatic expansion to $V_2 = 20\,\mathrm{L}$. The area between a path and the volume axis is $|w|$ for that path. The areas overlap: the area under the isobaric path includes the blue and orange areas, and the area under the isothermal path includes the orange area. The isochoric path has no area under it, so $w = 0$.
:tags: [hide-input]

import sys

import matplotlib.pyplot as plt
import numpy as np

sys.path.append("..")            # the book folder, which holds figstyle
import figstyle as fs

R_GAS = 0.083145                 # L bar K^-1 mol^-1
GAMMA = 5 / 3                    # monatomic ideal gas
V_1, V_2, T_1 = 5.0, 20.0, 300.0             # L, L, K (one mole)
P_1 = R_GAS * T_1 / V_1                      # bar
P_HEATED = 1.3 * P_1                         # bar, the end of the isochoric heating (390 K)


def draw_paths(c):
    fig, ax = plt.subplots(figsize=(6.4, 4.0), constrained_layout=True)
    V = np.linspace(V_1, V_2, 400)
    P_isobaric, P_isothermal, P_adiabatic = np.full_like(V, P_1), P_1 * V_1 / V, P_1 * (V_1 / V) ** GAMMA

    # The area under a path is |w|. Each band has the color of the path that bounds it from above.
    ax.fill_between(V, P_isothermal, P_isobaric, color=c.tint(c.green, 0.8), lw=0)
    ax.fill_between(V, P_adiabatic, P_isothermal, color=c.tint(c.blue, 0.75), lw=0)
    ax.fill_between(V, 0, P_adiabatic, color=c.tint(c.vermilion, 0.75), lw=0)

    # Paths, each named at its final state with the condition that defines it
    for P, color, name in [(P_isobaric, c.green, "isobaric\n$dP = 0$"), (P_isothermal, c.blue, "isothermal\n$dT = 0$"),
                           (P_adiabatic, c.vermilion, "adiabatic\n$\\delta q = 0$")]:
        ax.plot(V, P, color=color)
        ax.plot(V[-1], P[-1], "o", color=color, ms=6.5, zorder=5, clip_on=False)
        ax.text(V[-1] + 0.5, P[-1], name, color=color, va="center", clip_on=False)
    ax.plot([V_1, V_1], [P_1, P_HEATED], color=c.purple)
    ax.plot(V_1, P_HEATED, "o", color=c.purple, ms=6.5, zorder=5)
    ax.text(V_1 + 0.45, P_HEATED, "isochoric, $dV = 0$", color=c.purple, va="center")

    # Initial state
    ax.plot(V_1, P_1, "o", color=c.ink, ms=6.5, zorder=6)
    ax.text(V_1 - 0.35, P_1, "initial\nstate", ha="right", va="center")

    ax.set(xlabel="Volume $V$ (L)", ylabel="Pressure $P$ (bar)", xlim=(2.6, 20), ylim=(0, 7.0))
    ax.set_xticks([5, 10, 15, 20])
    return fig


fs.show(draw_paths, alt="Pressure against volume with four paths that start from one state at 5 L and about 5 bar. The "
                        "isochoric path rises vertically at 5 L to about 6.5 bar. The isobaric path runs horizontally to "
                        "20 L. The isothermal path falls to about 1.2 bar at 20 L. The adiabatic path falls more steeply "
                        "to about 0.5 bar at 20 L. The areas under the isobaric, isothermal, and adiabatic paths are "
                        "shaded, and the isobaric path has the largest area.")
```

---

### Applying the Workflow: $V$ and $T$ as Independent Variables

For many systems—especially gases—choosing $V$ and $T$ as independent variables simplifies calculations. Let us walk through the five-step workflow for a monatomic ideal gas. Steps 3 and 4 can be taken in either order. Here the ideal gas is specified first (Step 4), because this simplifies the expression to which the constraints are then applied (Step 3).

#### Step 1–2: First Law in Terms of $V$ and $T$

Starting from the First Law for a reversible process with only $PV$ work,

```{math}
\delta q = dU + P\,dV,
```

we rewrite $dU$ as a total differential in $V$ and $T$:

```{math}
dU = \left(\frac{\partial U}{\partial T}\right)_V dT 
     \;+\; \left(\frac{\partial U}{\partial V}\right)_T dV.
```

Substituting:

```{math}
:label: first-law-VT
\delta q = 
  \underbrace{\left(\frac{\partial U}{\partial T}\right)_V}_{C_V} dT 
  \;+\; \left[\left(\frac{\partial U}{\partial V}\right)_T + P\right] dV.
```

#### Step 4: Ideal-Gas Simplification

For an ideal gas, $U$ depends only on $T$ ([Section 2.5](../chapter-02/section-05.md)), so $\left(\frac{\partial U}{\partial V}\right)_T = 0$. For a monatomic ideal gas specifically,

```{math}
U = \frac{3}{2} N k_{\mathrm{B}} T,
\quad
P = \frac{N k_{\mathrm{B}} T}{V},
\quad
C_V = \left(\frac{\partial U}{\partial T}\right)_V = \frac{3}{2} N k_{\mathrm{B}}.
```

Equation {eq}`first-law-VT` then simplifies to

```{math}
:label: first-law-VT-ideal
\delta q = C_V \, dT + P \, dV
\qquad\text{(ideal gas, reversible process, $PV$-only work).}
```

#### Step 3: Process Constraints

````{admonition} **Prediction exercise**
:class: important

Before we derive the results, predict which quantities vanish for each process. Fill in the blanks:

| Process | $dV = 0$? | $dT = 0$? | $\delta q = 0$? | $\Delta U = ?$ |
|---------|:---------:|:---------:|:---------------:|:--------------:|
| Isochoric heating | | | | |
| Isothermal expansion | | | | |
| Adiabatic expansion | | | | |

```{admonition} Check your predictions
:class: dropdown

| Process | $dV = 0$? | $dT = 0$? | $\delta q = 0$? | $\Delta U$ |
|---------|:---------:|:---------:|:---------------:|:----------:|
| Isochoric heating | **Yes** | No | No | $q$ (= $C_V \Delta T$) |
| Isothermal expansion | No | **Yes** | No | **0** (ideal gas) |
| Adiabatic expansion | No | No | **Yes** | $w$ |

Key insight: each constraint eliminates a different term from Equation {eq}`first-law-VT-ideal`, producing a different relationship among $q$, $w$, and state-variable changes.
```
````

Applying these constraints to Equation {eq}`first-law-VT-ideal`:

```{list-table} Processes for an Ideal Gas ($V$ and $T$ as Independents)
:header-rows: 1
:name: thermodynamic-process-constraints

* - Constraint
  - Condition
  - Resulting $\delta q$
* - Isochoric
  - $dV = 0$
  - $\delta q = C_V\, dT$
* - Isothermal
  - $dT = 0$
  - $\delta q = P\, dV$
* - Adiabatic
  - $\delta q = 0$
  - $C_V\, dT = -P\, dV$
```

#### Step 5: Integrations Under Specific Constraints

1. **Isochoric ($dV = 0$)**

   With no volume change, there is no $PV$ work ($w = 0$). All energy enters as heat:

   ```{math}
   \delta q = C_V\, dT 
   \quad \Longrightarrow \quad
   q = \int_{T_1}^{T_2} C_V \, dT 
     = \frac{3}{2} N k_{\mathrm{B}} \,\Delta T.
   ```

   Since $w = 0$, we also have $\Delta U = q = C_V \Delta T$.

2. **Isothermal ($dT = 0$)**

   For an ideal gas, $\Delta U = 0$ (because $U$ depends only on $T$). The system absorbs heat to offset the work it does:

   ```{math}
   q = \int_{V_1}^{V_2} \frac{N k_{\mathrm{B}} T}{V} \, dV
     = N k_{\mathrm{B}} T \ln\!\bigl(\tfrac{V_2}{V_1}\bigr).
   ```

   And $w = -q = -N k_{\mathrm{B}} T \ln(V_2/V_1)$.

3. **Adiabatic ($\delta q = 0$)**

   Setting $\delta q = 0$ in Equation {eq}`first-law-VT-ideal`:

   ```{math}
   C_V\, dT = -P\, dV = -\frac{N k_{\mathrm{B}} T}{V}\, dV.
   ```

   Separating variables:

   ```{math}
   \frac{C_V}{T}\, dT = -\frac{N k_{\mathrm{B}}}{V}\, dV
   \quad \Longrightarrow \quad
   \int_{T_1}^{T_2} \frac{C_V}{T}\, dT 
   = - \int_{V_1}^{V_2} \frac{N k_{\mathrm{B}}}{V}\, dV.
   ```

   Integrating both sides:

   ```{math}
   \frac{3}{2} N k_{\mathrm{B}} \ln\frac{T_2}{T_1} = -N k_{\mathrm{B}} \ln\frac{V_2}{V_1}.
   ```

   Dividing through by $N k_{\mathrm{B}}$:

   ```{math}
   \ln\frac{T_2}{T_1} = -\frac{2}{3}\ln\frac{V_2}{V_1} = \ln\!\left(\frac{V_1}{V_2}\right)^{2/3}.
   ```

   Exponentiating:

   ```{math}
   :label: adiabatic-TV-relation
   \boxed{T_1 V_1^{2/3} = T_2 V_2^{2/3}}
   \qquad\text{(monatomic ideal gas).}
   ```

   ````{admonition} Where does the exponent $2/3$ come from?
   :class: note

   The exponent arises from the ratio $N k_{\mathrm{B}} / C_V = N k_{\mathrm{B}} / (\frac{3}{2} N k_{\mathrm{B}}) = 2/3$. For a general ideal gas, this ratio is $\gamma - 1$, where $\gamma = C_P / C_V$ is the **heat-capacity ratio** and $C_P$ is the heat capacity at constant pressure ([Section 3.3](section-03.md)). The general adiabatic relation is

   ```{math}
   T V^{\gamma - 1} = \text{const},
   ```

   which, combined with $PV = nRT$, also gives $PV^\gamma = \text{const}$. We will use these relations when we compare isothermal and adiabatic processes below.
   ````

---

### Comparing Isothermal and Adiabatic Expansions

How does the path affect the work? [](#fig:isothermal-adiabatic) compares a reversible isothermal expansion with a reversible adiabatic expansion from the same initial state. The area under each curve is $|w|$ for that process.

% Try this: Change `gamma` from `5/3` to `7/5`, the value of a diatomic gas near room temperature. How does the work of the adiabatic expansion change, and why does the work of the isothermal expansion stay the same?

```{code-cell} ipython3
:label: fig:isothermal-adiabatic
:caption: Reversible isothermal and adiabatic expansions of one mole of a monatomic ideal gas from $(V_1, T_1) = (5\,\mathrm{L}, 300\,\mathrm{K})$ to $V_2 = 20\,\mathrm{L}$. The area under each curve is $|w|$ for that path, and the labels give $w$. The orange area lies under both curves.
:tags: [hide-input]

import sys
from fractions import Fraction

import matplotlib.pyplot as plt
import numpy as np

sys.path.append("..")            # the book folder, which holds figstyle
import figstyle as fs

# Parameters
n = 1.0        # mol
R = 8.314      # J/(mol K)
T1 = 300       # K
V1 = 5.0e-3    # m^3 (5 L)
V2 = 20.0e-3   # m^3 (20 L)
gamma = 5/3    # monatomic ideal gas

P1 = n * R * T1 / V1  # initial pressure

# Volume array
V = np.linspace(V1, V2, 500)
V_L = V * 1e3  # convert to liters for plotting

# Isothermal: P = nRT/V
P_iso = n * R * T1 / V

# Adiabatic: P V^gamma = P1 V1^gamma  =>  P = P1 (V1/V)^gamma
P_adi = P1 * (V1 / V)**gamma

# Convert pressures to bar for plotting
P_iso_bar = P_iso / 1e5
P_adi_bar = P_adi / 1e5
P1_bar = P1 / 1e5

# Work along each path
w_iso = -n * R * T1 * np.log(V2 / V1)  # J
w_adi = -(P1 * V1**gamma / (1 - gamma)) * (V2**(1 - gamma) - V1**(1 - gamma))  # J


def draw(c):
    fig, ax = plt.subplots(figsize=(6.4, 3.7), constrained_layout=True)

    # Areas under the curves
    ax.fill_between(V_L, P_adi_bar, P_iso_bar, color=c.tint(c.blue, 0.75), lw=0)
    ax.fill_between(V_L, 0, P_adi_bar, color=c.tint(c.vermilion, 0.75), lw=0)

    # Curves, each labeled at its final state
    for P_bar, color, name, w in [(P_iso_bar, c.blue, f"isothermal ($T = {T1}$ K)", w_iso),
                                  (P_adi_bar, c.vermilion,
                                   rf"adiabatic ($\gamma = {Fraction(gamma).limit_denominator(9)}$)", w_adi)]:
        ax.plot(V_L, P_bar, color=color)
        ax.plot(V_L[-1], P_bar[-1], "o", color=color, ms=6.5, zorder=5, clip_on=False)
        ax.text(V_L[-1] + 0.5, P_bar[-1], f"{name}\n$w = {w / 1e3:.2f}$ kJ".replace("-", "−"), color=color,
                va="center", clip_on=False)

    # Initial state
    ax.plot(V1 * 1e3, P1_bar, "o", color=c.ink, ms=6.5, zorder=5)
    ax.text(V1 * 1e3 + 0.4, P1_bar, "initial state", va="center")

    ax.set(xlabel="Volume $V$ (L)", ylabel="Pressure $P$ (bar)", xlim=(4, 20), ylim=(0, 5.5))
    return fig


fs.show(draw, alt="Pressure against volume for two expansions of an ideal gas from 5 L and about 5 bar to 20 L. The "
                  "isothermal curve at 300 K ends at about 1.2 bar, and the work is minus 3.46 kJ. The adiabatic "
                  "curve falls more steeply and ends at about 0.5 bar, and the work is minus 2.26 kJ. The area under "
                  "the adiabatic curve is smaller than the area under the isothermal curve.")
```

```{admonition} Comparing the two paths
:class: important

Which expansion does more work on the surroundings? Use the First Law to explain why the two curves separate as the gas expands, and say what happens to the temperature along each path.
```

---

### Microscopic Interpretation of the First Law

As developed in [Section 2.1](../chapter-02/section-01.md), macroscopic thermodynamic properties are ensemble averages over microstates. This perspective gives a powerful physical picture of heat and work at the molecular level.

From a statistical mechanics perspective, the internal energy is

```{math}
U = \sum_{i=1}^M p_i E_i,
```

where $p_i$ is the probability of occupying the $i$-th microstate with energy $E_i$, and $M$ is the total number of accessible microstates. The total differential can be written as

```{math}
:label: dU-microscopic
dU 
= \underbrace{\sum_{i=1}^M E_i\, dp_i}_{\text{changing which states are occupied}} 
  + \underbrace{\sum_{i=1}^M p_i \, dE_i}_{\text{changing the energies of the states}}.
```

For a closed system (constant $N$) with only $PV$ work, the microstate energies depend on the volume, so $dE_i = \left(\frac{\partial E_i}{\partial V}\right)_N dV$. Noting that

```{math}
\sum_{i=1}^M p_i \left(\frac{\partial E_i}{\partial V}\right)_N 
  = \Bigl\langle \bigl(\tfrac{\partial E}{\partial V}\bigr)_N \Bigr\rangle,
```

we identify pressure as

```{math}
P = - \Bigl\langle \bigl(\tfrac{\partial E}{\partial V}\bigr)_N \Bigr\rangle.
```

Hence,

```{math}
dU 
= \sum_{i=1}^M E_i\, dp_i 
  \;-\; P\, dV.
```

Comparing term by term with $dU = \delta q + \delta w = \delta q - P\,dV$, which holds for a reversible process, we identify:

```{math}
\delta q = \sum_{i=1}^M E_i\, dp_i
\qquad\text{and}\qquad
\delta w = -P\,dV = \sum_{i=1}^M p_i\, dE_i.
```

```{admonition} Microscopic meaning of heat and work
:class: note

- **Heat** ($\delta q = \sum E_i\,dp_i$): Energy exchange via *changing the probability distribution* over microstates. The energy levels stay the same, but the system redistributes among them. This is what happens when you bring a cold system into contact with a hot reservoir.

- **Work** ($\delta w = \sum p_i\,dE_i$): Energy exchange via *shifting the energy levels themselves*. The distribution stays the same, but every level moves. This is what happens when you compress a gas—the particle-in-a-box energy levels all increase as the box shrinks.
```

[](#fig:heat-and-work-levels) illustrates these two mechanisms. Panel (a) shows the initial state, panel (b) the result of adding heat (the probabilities change, and the energy levels stay), and panel (c) the result of doing work by compression (the energy levels shift, and the probabilities stay).

% Try this: Keep the temperature fixed during the compression: replace `boltzmann(E_compressed, 2 * kT_initial)` by `boltzmann(E_compressed, kT_initial)`. In which direction do the probabilities of panel (c) shift, and what does that say about the heat in an isothermal compression?

```{code-cell} ipython3
:label: fig:heat-and-work-levels
:caption: Microstate energies (horizontal lines) and occupation probabilities (bars, with the value of $p_i$ beside each bar). (a) Initial state. (b) After heating: the energy levels are unchanged, and the distribution shifts toward higher-energy states ($\sum E_i\,dp_i > 0$). The dashed outlines are the initial probabilities. (c) After compression: the probabilities are unchanged, and every energy level shifts upward ($\sum p_i\,dE_i > 0$). The dashed lines are the initial levels, and each arrow joins a level before and after the compression. The probabilities in each panel are Boltzmann probabilities ([Section 2.2](../chapter-02/section-02.md)): heating doubles $T$ at fixed levels, and the compression doubles the level spacing and $T$ together, which leaves every $p_i$ unchanged.
:tags: [hide-input]

from matplotlib.patches import Rectangle


def boltzmann(E, kT):
    """Boltzmann probabilities for levels E at thermal energy kT."""
    weights = np.exp(-E / kT)
    return weights / weights.sum()


# Energies and kT in units of the initial level spacing
E_initial = np.array([0.5, 1.5, 2.5, 3.5, 4.5])
kT_initial = 1.5
p_initial = boltzmann(E_initial, kT_initial)
p_heated = boltzmann(E_initial, 2 * kT_initial)                  # heating: twice the temperature, same levels
E_compressed = 2 * E_initial                                     # compression: twice the level spacing
p_compressed = boltzmann(E_compressed, 2 * kT_initial)           # and twice the temperature: same probabilities


def draw_levels(c):
    gutter = 0.2                                                 # room for the arrows of panel (c)
    fig, axes = plt.subplots(1, 3, figsize=(fs.WIDTH, 4.0), constrained_layout=True, sharey=True,
                             gridspec_kw=dict(width_ratios=[0.8, 0.8, 0.8 + gutter]))
    bar = 0.44                                                   # height of a probability bar
    small = fs.BASE - 1.5

    def levels_and_bars(ax, E, p, level_color, bar_color, text_color, label_from=None):
        """Levels E with a bar of length p on each, and the value of p beside the bar."""
        ends = p if label_from is None else np.maximum(p, label_from)
        for i in range(len(E)):
            ax.plot([0, 0.8], [E[i], E[i]], color=level_color, lw=1.4, solid_capstyle="butt")
            ax.add_patch(Rectangle((0, E[i]), p[i], bar, fc=bar_color, ec="none"))
            ax.text(ends[i] + 0.02, E[i] + bar / 2, f"{p[i]:.2f}", color=text_color, va="center", fontsize=small)

    # (a) Initial state
    levels_and_bars(axes[0], E_initial, p_initial, c.gray, c.tint(c.gray, 0.45), c.ink)
    axes[0].text(0.78, 8.25, "Boltzmann distribution", ha="right", va="center")

    # (b) After heating: same levels, new probabilities (the dashed outlines are the initial ones)
    levels_and_bars(axes[1], E_initial, p_heated, c.gray, c.tint(c.vermilion, 0.3), c.vermilion,
                    label_from=p_initial)
    for E, p in zip(E_initial, p_initial):
        axes[1].add_patch(Rectangle((0, E), p, bar, fc="none", ec=c.ink, lw=1.0, ls=(0, (3, 2))))
    axes[1].text(0.78, 8.25, r"$\delta q = \sum_i E_i\,dp_i > 0$", ha="right", va="center")

    # (c) After compression: new levels, same probabilities (the dashed lines are the initial levels)
    for E in E_initial:
        axes[2].plot([0, 0.8], [E, E], color=c.ink, lw=1.0, ls=(0, (3, 2)))
    levels_and_bars(axes[2], E_compressed, p_compressed, c.vermilion, c.tint(c.gray, 0.45), c.ink)
    for E_old, E_new, x in zip(E_initial, E_compressed, [-0.035, -0.035, -0.075, -0.115, -0.155]):
        fs.arrow(axes[2], (x, E_old), (x, E_new), c.vermilion, lw=1.2, scale=9)
    axes[2].text(0.78, 8.25, r"$\delta w = \sum_i p_i\,dE_i > 0$", ha="right", va="center")

    for ax, letter, title in zip(axes, "abc", ["Initial state", "After heating", "After compression"]):
        left = -gutter if ax is axes[2] else 0
        ax.plot([0, 0], [0, 10], color=c.muted, lw=0.9)          # zero of the probability
        ax.set(xlim=(left, 0.8), ylim=(0, 10), xticks=[0, 0.25, 0.5, 0.75], yticks=[], xlabel="Probability $p_i$")
        ax.set_xticklabels(["0", "0.25", "0.5", "0.75"])
        ax.spines["left"].set_visible(False)
        ax.spines["bottom"].set_bounds(0, 0.8)
        fs.panel(ax, letter, title)
    axes[0].set_ylabel("Energy")
    return fig


fs.show(draw_levels,
        alt="Three panels, each with five energy levels drawn as horizontal lines and a bar on each level whose "
            "length is the probability of that level. Panel a, the initial state: the levels are equally spaced, "
            "and the probabilities fall from 0.50 on the lowest level to 0.04 on the highest. Panel b, after "
            "heating: the levels are at the same energies, and the probabilities have changed to 0.35, 0.25, 0.18, "
            "0.13, and 0.09, so the higher levels are more populated. Panel c, after compression: every level has "
            "moved up to twice its energy, and the probabilities are the same as in the initial state.")
```

In [](#fig:heat-and-work-interactive), a particle in a one-dimensional box ([Section 2.5](../chapter-02/section-05.md)) takes the place of the schematic levels. Choose a process, move its slider, and follow the heat, the work, and the change in internal energy.

```{iframe} https://wexlergroup.github.io/chem-4020-5020/computational-studios/heat-and-work.html
:label: fig:heat-and-work-interactive
:width: 100%
:title: Heat and work on the energy levels of a particle in a box
:placeholder: ../_static/heat-and-work/heat-and-work.png

Heat and work for a particle in a one-dimensional box. Panel (a) shows the energy levels $E_i$ with a bar for each probability $p_i$, and panel (b) shows the heat $q$, the work $w$, and the change in internal energy $\Delta U$ since the starting state. Energies are in units of $\varepsilon = h^2/(8mL_0^2)$, the energy of the lowest level at the starting length $L_0$.
```

Two things to try:

1. Heat the particle at fixed length, and then compress the box without heat. For each step, say which of $q$ and $w$ changes and why.
2. Go from the starting state to $L = 0.8\,L_0$ and $k_{\mathrm{B}}T = 10\,\varepsilon$ along two routes. On the first route, heat to $6.4\,\varepsilon$ and then compress. On the second route, compress and then heat. Compare $q$, $w$, and $\Delta U$ for the two routes.

````{admonition} **Think–pair–share**
:class: important

In the microscopic picture, heating at *constant volume* changes $\{p_i\}$ while $\{E_i\}$ stay fixed (panel (b) of [](#fig:heat-and-work-levels)). What would the figure look like if you heated at *constant pressure* instead? Would only the probabilities change, only the energies, or both? Why?

```{admonition} Discussion
:class: dropdown

**Both.** At constant pressure, adding heat increases $T$, which shifts the probability distribution toward higher-energy states (just as at constant $V$). But at constant $P$, the gas also *expands*, which changes the volume and therefore shifts the energy levels (just as in compression, but in reverse—the levels *decrease*). A process at constant pressure therefore changes both $\{p_i\}$ and $\{E_i\}$.

[Section 3.3](section-03.md) takes up heating at constant pressure.
```
````

## Worked Example

### Isothermal Reversible Expansion of an Ideal Gas

One mole of an ideal gas expands reversibly and isothermally at $T=300\ \mathrm{K}$ from $V_1=5.0\ \mathrm{L}$ to $V_2=20.0\ \mathrm{L}$. Compute $w$, $q$, and $\Delta U$.

```{admonition} Notation
:class: note

We switch to per-mole notation ($n$, $R$) for this macroscopic calculation, following the course convention (see [Notation](../notation.md): use $k_{\mathrm{B}}$ for per-particle quantities, $R$ for per-mole quantities).
```

**Assumptions.** Ideal gas, reversible process, isothermal; for an ideal gas $U=U(T)$.

1. **Work**

   For a reversible process, $P_{\mathrm{ext}} = P = nRT/V$ at every step:

   ```{math}
   w = -\int_{V_1}^{V_2} P\,dV = -\int_{V_1}^{V_2}\frac{nRT}{V}\,dV
   = -nRT\ln\left(\frac{V_2}{V_1}\right).
   ```

   With $n=1\,\mathrm{mol}$, $R=8.314\ \mathrm{J\,mol^{-1}\,K^{-1}}$, $T=300\ \mathrm{K}$, and $V_2/V_1=4$:

   ```{math}
   w = -(1)(8.314)(300)\ln 4 = -3.46\times10^{3}\ \mathrm{J}=-3.46\ \mathrm{kJ}.
   ```

2. **Internal energy**

   ```{math}
   \Delta U = 0 \quad (\text{isothermal ideal gas: } U = U(T) \text{ and } \Delta T = 0).
   ```

3. **Heat from the First Law**

   ```{math}
   \Delta U = q + w \;\Rightarrow\; 0 = q + w \;\Rightarrow\; q = -w = +3.46\ \mathrm{kJ}.
   ```

**Result.** $w=-3.46\ \mathrm{kJ}$, $q=+3.46\ \mathrm{kJ}$, and $\Delta U=0$. The gas does work on the surroundings, and an equal amount of heat flows in from the thermal reservoir to maintain the temperature.

````{admonition} **Computational check**
:class: important

Verify this result using Python. Compute $w = -nRT\ln(V_2/V_1)$ numerically and confirm it matches the hand calculation. Then plot $P$ vs. $V$ for this isothermal expansion and shade the area under the curve. Does the shaded area match $|w|$?

% Try this: Change the number of points of the numerical integration from `10000` to `5`, and run this cell and the next one. By how many percent does the numerical work differ from the analytical one, and on which side of the curve do the straight segments of the figure lie?

```{code-cell} ipython3
:tags: [hide-input]

n = 1.0        # mol
R = 8.314      # J/(mol K)
T = 300        # K
V1 = 5.0e-3    # m^3 (5 L)
V2 = 20.0e-3   # m^3 (20 L)

# Analytical result
w = -n * R * T * np.log(V2 / V1)
print(f"Work (analytical): w = {w:.1f} J = {w/1e3:.3f} kJ")

# Numerical integration (trapezoidal) for verification
V = np.linspace(V1, V2, 10000)
P = n * R * T / V
w_numerical = -np.trapezoid(P, V)
print(f"Work (numerical):  w = {w_numerical:.1f} J = {w_numerical/1e3:.3f} kJ")
```

```{code-cell} ipython3
:label: fig:isothermal-work
:caption: Reversible isothermal expansion of one mole of an ideal gas at $T = 300\ \mathrm{K}$ from $V_1 = 5.0\ \mathrm{L}$ to $V_2 = 20.0\ \mathrm{L}$. The shaded area under the curve is $|w|$.
:tags: [hide-input]

def draw_isotherm(c):
    fig, ax = plt.subplots(figsize=(5.4, 3.3), constrained_layout=True)
    V_L = V * 1e3      # L
    P_bar = P / 1e5    # bar
    ax.fill_between(V_L, 0, P_bar, color=c.tint(c.blue, 0.75), lw=0)
    ax.plot(V_L, P_bar, color=c.blue)
    ax.plot([V_L[0], V_L[-1]], [P_bar[0], P_bar[-1]], "o", color=c.blue, ms=6.5, zorder=5, clip_on=False)
    ax.text(9.5, 0.75, f"$|w| = {abs(w) / 1e3:.2f}$ kJ", ha="center", va="center")
    ax.text(11.5, 2.75, "$T = 300$ K", color=c.blue, va="center")
    ax.set(xlabel="Volume $V$ (L)", ylabel="Pressure $P$ (bar)", xlim=(4, 21), ylim=(0, 5.5))
    return fig


fs.show(draw_isotherm, alt="Pressure against volume for the isothermal expansion of an ideal gas at 300 K. The curve "
                           "falls from about 5 bar at 5 L to about 1.2 bar at 20 L. The area under the curve between "
                           "the two volumes is shaded, and it equals the magnitude of the work, 3.46 kJ.")
```
````

## Concept Checks

1. Why is reversibility useful even when real processes are irreversible?
2. What distinguishes an isothermal process from an adiabatic process in terms of energy transfer mechanisms?
3. For an ideal gas, why does $\left(\partial U/\partial V\right)_T=0$?
4. In the microscopic expression $dU=\sum_i E_i\,dp_i+\sum_i p_i\,dE_i$, which term corresponds to heat and which to work (under the stated conditions)?
5. An isothermal expansion and an adiabatic expansion start from the same initial state and end at the same final volume. Which produces more work? Which produces a larger temperature change?

## Key Takeaways

- A **reversible** expansion does the maximum work on the surroundings, and a reversible compression requires the minimum work.
- Choosing independent variables and constraints turns the First Law into solvable differential relations.
- Ideal-gas isothermal expansion has $\Delta U=0$ and $q=-w$.
- Adiabatic constraints set $q=0$ and force $T$ to change when work is done; the relation $TV^{\gamma-1} = \text{const}$ governs the temperature–volume coupling.
- On a $P$–$V$ diagram, the magnitude of the reversible work equals the area under the curve, so the work depends on the path that is followed.
- Microscopically, heat changes the probability distribution over fixed energy levels, while work shifts the energy levels themselves.
