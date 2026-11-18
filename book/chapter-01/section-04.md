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
  enumerator: 1.4.%s
---

# 1.4. Real Gases

See also: [Course-wide Conventions & Notation](../notation.md)

Run the code: [Notebook for Section 1.4](https://wexlergroup.github.io/chem-4020-5020/lite/notebooks/index.html?path=chapter-01/section-04.ipynb), in which you can change and rerun the Python cells of this section in your browser.

## Overview and Learning Objectives

Real gases deviate from ideal behavior when finite molecular size and intermolecular forces become important. This section introduces the compressibility factor as a quantitative measure of nonideality and develops the van der Waals equation of state (EOS), including its critical constants and the principle of corresponding states.

At sufficiently high pressures (high densities) or low temperatures, we can no longer treat gases as collections of noninteracting point particles. High pressure/density brings molecules closer together, increasing collision frequency and making both excluded volume and intermolecular interactions more consequential. Lower temperature reduces kinetic energy, so these interactions more strongly perturb $PVT$ behavior. This section develops practical ways to quantify deviations from ideality and introduces the van der Waals EOS.

Learning objectives:

- Define the compressibility factor $Z$ and interpret $Z<1$ vs. $Z>1$ in terms of net attractions/repulsions.
- Write the van der Waals EOS in particle and molar forms and interpret the parameters $a$ and $b$.
- Derive (or use) the van der Waals critical constants $V_{\mathrm m,c}=3b$, $P_c=a/(27b^2)$, $T_c=8a/(27bR)$.
- Non-dimensionalize the van der Waals EOS to obtain the corresponding-states form in reduced variables.

## Core Ideas and Derivations

### Deviations from Ideal Behavior

#### Compressibility Factor

```{code-cell} ipython3
:label: fig:compressibility-water
:caption: Compressibility factor $Z$ of water at 776 K, above its critical temperature of 647 K, as a function of pressure, with the ideal-gas line ($Z=1$) for comparison. Where $Z<1$ (blue), attractive forces dominate, most strongly at the minimum of $Z$. Where $Z>1$ (orange), repulsive forces dominate. The crossing of $Z=1$ marks a “coincidentally ideal” state.
:tags: [hide-input]

import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.constants import R

sys.path.append("..")            # the book folder, which holds figstyle
import figstyle as fs

# Measured isotherm of water at 776.4 K (T_r = 1.2, above the critical temperature), NIST Chemistry WebBook
water = pd.read_csv("../_static/section-04/isothermal-properties-for-water-Tr=1.2.csv")
P_water = water["Pressure (bar)"].to_numpy()
Z_water = P_water * 1e5 * water["Volume (m3/mol)"].to_numpy() / (R * water["Temperature (K)"].to_numpy())
T_water = water["Temperature (K)"].iloc[0]                 # temperature of the isotherm (K), for the label


def draw(c):
    fig, ax = plt.subplots(figsize=(6, 3.9), constrained_layout=True)
    ax.fill_between(P_water, Z_water, 1, where=Z_water <= 1, color=c.tint(c.blue, 0.78), lw=0, interpolate=True)
    ax.fill_between(P_water, Z_water, 1, where=Z_water >= 1, color=c.tint(c.vermilion, 0.78), lw=0, interpolate=True)
    ax.axhline(1, color=c.gray, ls="--", lw=1.4)
    ax.plot(P_water, Z_water, color=c.ink)
    ax.text(3950, 0.97, "ideal gas, $Z = 1$", color=c.gray, ha="right", va="top")
    ax.text(1000, 0.47, f"water at {T_water:.0f} K", color=c.ink, ha="left", va="top")
    P_cross = np.interp(1, Z_water[Z_water.argmin():], P_water[Z_water.argmin():])
    ax.plot(P_cross, 1, "o", mfc=c.ground, mec=c.ink, mew=1.5, ms=6, zorder=5)
    ax.annotate("coincidentally ideal", xy=(P_cross, 1), xytext=(2000, 1.20), ha="center", va="center",
                arrowprops=dict(arrowstyle="-", color=c.muted, lw=0.9, shrinkB=5))
    ax.text(1050, 0.86, "attraction dominates\n$Z < 1$", color=c.blue, ha="center", va="center")
    ax.text(3500, 1.36, "repulsion dominates\n$Z > 1$", color=c.vermilion, ha="right", va="center")
    ax.set(xlabel="Pressure $P$ (bar)", ylabel="Compressibility factor $Z$", xlim=(0, 4000), ylim=(0.3, 1.5))
    return fig


fs.show(draw, alt="Compressibility factor of water at 776 K against pressure from 0 to 4000 bar. The curve starts at 1, "
                  "falls to a minimum near 0.5, rises through 1, and continues to rise. The area between the curve "
                  "and the line Z equals 1 is labeled attraction dominates where Z is below 1 and repulsion dominates "
                  "where Z is above 1.")
```

One convenient measure of nonideality is the compressibility factor $Z$:

```{math}
:label: compressibility-factor-eq
Z = \frac{P V}{N k_\text{B} T} = \frac{P V}{n R T}.
```

```{list-table} Behavior of Water at Different Pressures
:header-rows: 1
:name: water-compressibility-analysis

* - Pressure
  - $Z$
  - Volume
  - Analysis
  - Water
* - High
  - $> 1$
  - $V_\text{real} > V_\text{ideal}$
  - Repulsive interactions
  - Pauli (steric) repulsion
* - Medium
  - $< 1$
  - $V_\text{real} < V_\text{ideal}$
  - Attractive interactions
  - H-bonding and van der Waals interactions
* - Low
  - $\to 1$
  - $V_\text{real} \to V_\text{ideal}$
  - Ideal-gas behavior
  - —
```

```{tip}
Show that $V_\text{real} > V_\text{ideal}$ for $Z > 1$ and $V_\text{real} < V_\text{ideal}$ for $Z < 1$.
```

```{tip}
What might be the origin of attractive interactions in a gas of carbon dioxide or methane molecules?
```

### Van der Waals Fluid

#### Van der Waals Equation of State

The van der Waals equation of state is among the simplest real-gas models (and was central to the work recognized by the [Nobel Prize in Physics 1910](https://www.nobelprize.org/prizes/physics/1910/summary/)):

````{list-table} Forms of the van der Waals Equation of State
:header-rows: 1
:name: van-der-waals-forms

* - 
  - Extensive Volume
  - Intensive Volume
* - Per Particle
  - ```{math}
    \left( P + \frac{a_\text{p} N^2}{V^2} \right) \left( V - b_\text{p} N \right) = N k_\text{B} T
    ```
  - ```{math}
    \left( P + \frac{a_\text{p}}{V_\text{p}^2} \right) \left( V_\text{p} - b_\text{p} \right) = k_\text{B} T
    ```
* - Per Mole
  - ```{math}
    \left( P + \frac{a_\text{m} n^2}{V^2} \right) \left( V - b_\text{m} n \right) = n R T
    ```
  - ```{math}
    :label: van-der-waals-eos
    \left( P + \frac{a_\text{m}}{V_\text{m}^2} \right) \left( V_\text{m} - b_\text{m} \right) = R T,
    ```
````

where $a$ and $b$ are substance-dependent parameters capturing two departures from ideality:

- **$a$**: Net attractive interactions between particles.
- **$b$**: Volume associated with finite particle size (excluded volume).

[](#fig:pair-potential) relates the two parameters to the potential energy of a pair of molecules.

```{code-cell} ipython3
:label: fig:pair-potential
:caption: Schematic potential energy of two molecules as a function of the distance between them. The steep rise at short distance keeps the molecules apart, and the parameter $b$ accounts for the volume that they exclude. The shallow well at larger distance is the net attraction between the molecules, which the parameter $a$ accounts for.
:tags: [hide-input]

from matplotlib.patches import Circle


def draw_pair_potential(c):
    r = np.linspace(0.93, 2.6, 600)                      # distance in units of the molecular diameter
    u = 4 * (r**-12 - r**-6)                             # Lennard-Jones form, in units of the depth of the well
    fig, ax = plt.subplots(figsize=(6, 3.9), constrained_layout=True)
    ax.axhline(0, color=c.gray, lw=0.9)
    ax.fill_between(r, u, 0, where=u > 0, color=c.tint(c.vermilion, 0.78), lw=0)
    ax.fill_between(r, u, 0, where=u < 0, color=c.tint(c.blue, 0.78), lw=0)
    ax.plot(r[u >= 0], u[u >= 0], color=c.vermilion)
    ax.plot(r[u <= 0], u[u <= 0], color=c.blue)
    ax.text(1.06, 1.25, "repulsion at short distance:\nthe molecules cannot overlap,\nwhich gives $b$",
            color=c.vermilion, va="center", ha="left")
    ax.text(1.56, -0.74, "attraction at larger distance:\nthe molecules pull on each other,\nwhich gives $a$",
            color=c.blue, va="center", ha="left")
    ax.set(xlim=(0.8, 2.6), ylim=(-1.25, 2.0), xticks=[1], xticklabels=["molecular\ndiameter"], yticks=[0],
           xlabel="Distance between two molecules", ylabel="Potential energy of the pair")
    return fig


fs.show(draw_pair_potential,
        alt="Potential energy of two molecules against the distance between them. At distances shorter than about one "
            "molecular diameter, the energy rises steeply above zero, labeled repulsion, which gives b. At larger "
            "distances, the energy dips below zero in a shallow well and returns to zero, labeled attraction, which "
            "gives a.")
```

#### Van der Waals Isotherm

```{code-cell} ipython3
:label: fig:vdw-isotherms
:caption: Measured (solid) and van der Waals (dashed) isotherms of CO$_2$ below, at, and above the critical temperature $T_c$. For $T < T_c$, the van der Waals EOS exhibits the characteristic loop associated with phase separation, where the measured isotherm is horizontal; at $T = T_c$ the isotherm has an inflection point; and for $T > T_c$ the liquid–vapor transition disappears.
:tags: [hide-input]

from matplotlib.lines import Line2D

# Van der Waals constants of carbon dioxide and the critical point that they predict
a_m = 3.6551                               # L^2 bar/mol^2
b_m = 0.042816                             # L/mol
R_L_bar = R / 100                          # gas constant in L bar/(mol K)
Tc_vdW = 8 * a_m / (27 * b_m * R_L_bar)    # K
Pc_vdW = a_m / (27 * b_m**2)               # bar
Vc_vdW = 3 * b_m                           # L/mol


def vdw_pressure(T, V):
    """Pressure (bar) of the van der Waals fluid at temperature T (K) and molar volume V (L/mol)."""
    return R_L_bar * T / (V - b_m) - a_m / V**2


# Measured isotherms of carbon dioxide at 0.9, 1.0, and 1.1 times its critical temperature (NIST Chemistry WebBook)
measured = {Tr: pd.read_csv(f"../_static/section-04/isothermal-properties-for-carbon-dioxide-Tr={Tr}.csv")
            for Tr in (0.9, 1.0, 1.1)}


def draw(c):
    colors = {0.9: c.blue, 1.0: c.gray, 1.1: c.vermilion}
    names = {0.9: r"$0.9\,T_c$", 1.0: r"$T_c$", 1.1: r"$1.1\,T_c$"}
    V = np.linspace(0.5 * Vc_vdW, 0.5, 800)
    fig, ax = plt.subplots(figsize=(6, 4), constrained_layout=True)
    for Tr in (0.9, 1.0, 1.1):
        data = measured[Tr]
        V_data, P_data = 1e3 * data["Volume (m3/mol)"].to_numpy(), data["Pressure (bar)"].to_numpy()
        ax.plot(V, vdw_pressure(Tr * Tc_vdW, V), color=colors[Tr], ls=(0, (3.5, 2)), lw=1.5)
        ax.plot(V_data, P_data, color=colors[Tr])
        order = np.argsort(V_data)
        ax.text(0.512, np.interp(0.5, V_data[order], P_data[order]), names[Tr], color=colors[Tr], va="center",
                clip_on=False)
    ax.legend(handles=[Line2D([], [], color=c.ink, label="measured"),
                       Line2D([], [], color=c.ink, ls=(0, (3.5, 2)), lw=1.5, label="van der Waals")],
              loc="upper right", bbox_to_anchor=(0.97, 1.0))
    ax.set(xlabel=r"Molar volume $V_{\mathrm{m}}$ (L/mol)", ylabel="Pressure $P$ (bar)", xlim=(0, 0.5), ylim=(20, 120))
    return fig


fs.show(draw, alt="Pressure against molar volume for carbon dioxide at 0.9, 1.0, and 1.1 times the critical "
                  "temperature. Solid curves are measured and dashed curves are the van der Waals equation. Below "
                  "the critical temperature, the measured isotherm has a horizontal segment where liquid and vapor "
                  "coexist, and the van der Waals isotherm has a loop there. The two sets of curves agree at large "
                  "volume and differ most at small volume.")
```

```{code-cell} ipython3
:label: fig:maxwell-construction
:caption: (a) Van der Waals isotherm of CO$_2$ at $T = 0.9\,T_c$ with the Maxwell construction, which replaces the unphysical loop by a horizontal tie line that defines liquid–vapor coexistence. The two shaded areas are equal. The thin curve joins the ends of the tie lines at all temperatures up to the critical point at $T = T_c$, where the liquid and vapor phases merge. (b) The regions of the diagram. Under the coexistence curve, the color follows the fraction of vapor.
:tags: [hide-input]

# The Maxwell construction follows https://scipython.com/blog/the-maxwell-construction/
from matplotlib.colors import to_rgb
from scipy.optimize import brentq


def vdw_reduced(Tr, Vr):
    """Reduced pressure of the van der Waals fluid at reduced temperature Tr and reduced volume Vr."""
    return 8 * Tr / (3 * Vr - 1) - 3 / Vr**2


def maxwell(Tr):
    """Maxwell construction for an isotherm below the critical temperature.

    Returns the reduced pressure of the tie line and the three reduced volumes at which the isotherm crosses it
    (liquid, middle, vapor). The tie line is placed so that the two areas between it and the loop are equal.
    """
    def crossings(Pr):
        return np.sort(np.roots([3 * Pr, -(Pr + 8 * Tr), 9, -3]).real)

    def area_difference(Pr):
        V_liq, _, V_vap = crossings(Pr)
        area_under_isotherm = 8 * Tr / 3 * np.log((3 * V_vap - 1) / (3 * V_liq - 1)) + 3 / V_vap - 3 / V_liq
        return area_under_isotherm - Pr * (V_vap - V_liq)

    # The tie line lies between the minimum and the maximum of the loop, where dPr/dVr = 0.
    V_min, V_max = np.sort(np.roots([4 * Tr, -9, 6, -1]).real)[1:]
    low, high = max(vdw_reduced(Tr, V_min), 0), vdw_reduced(Tr, V_max)
    span = high - low
    Pr = brentq(area_difference, low + 1e-6 * span, high - 1e-6 * span)
    return Pr, crossings(Pr)


# Coexistence curve: the ends of the tie lines for a series of temperatures up to the critical temperature
Tr_series = 1 - np.linspace(0.45, 0.01, 60) ** 2
ties = [maxwell(Tr) for Tr in Tr_series]
dome_P = np.array([Pr for Pr, _ in ties] + [1.0]) * Pc_vdW
dome_V_liq = np.array([V[0] for _, V in ties] + [1.0]) * Vc_vdW
dome_V_vap = np.array([V[2] for _, V in ties] + [1.0]) * Vc_vdW
dome_V = np.concatenate([dome_V_liq, dome_V_vap[::-1]])
dome = np.concatenate([dome_P, dome_P[::-1]])

Tr_loop = 0.9                              # reduced temperature of the isotherm with the loop (below 1)
Pr_tie, (Vr_liq, Vr_mid, Vr_vap) = maxwell(Tr_loop)
P_tie, V_liq, V_mid, V_vap = Pr_tie * Pc_vdW, Vr_liq * Vc_vdW, Vr_mid * Vc_vdW, Vr_vap * Vc_vdW


def draw(c):
    V = np.linspace(1.02 * b_m, 0.5, 1500)
    P_cold, P_crit = vdw_reduced(Tr_loop, V / Vc_vdW) * Pc_vdW, vdw_reduced(1.0, V / Vc_vdW) * Pc_vdW
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(fs.WIDTH, 3.7), constrained_layout=True, sharex=True, sharey=True)

    # (a) The loop, the tie line, and the two equal areas
    ax1.fill(dome_V, dome, color=c.fill, lw=0)
    ax1.plot(dome_V, dome, color=c.muted, lw=0.9)
    loop = (V >= V_liq) & (V <= V_vap)
    ax1.fill_between(V[loop], P_cold[loop], P_tie, color=c.tint(c.blue, 0.55), lw=0)
    ax1.plot(V, P_crit, color=c.gray, lw=1.6)
    ax1.plot(V, P_cold, color=c.blue)
    ax1.plot([V_liq, V_vap], [P_tie, P_tie], color=c.blue, lw=1.2)
    ax1.plot([V_liq, V_vap], [P_tie, P_tie], "o", mfc=c.ground, mec=c.blue, mew=1.5, ms=6, zorder=5)
    ax1.plot(Vc_vdW, Pc_vdW, "o", mfc=c.ground, mec=c.gray, mew=1.5, ms=6, zorder=5)
    ax1.plot([0, Vc_vdW], [Pc_vdW, Pc_vdW], color=c.gray, lw=0.9, ls=":")
    ax1.plot([Vc_vdW, Vc_vdW], [30, Pc_vdW], color=c.gray, lw=0.9, ls=":")
    ax1.text(0.012, Pc_vdW + 0.8, "$P_c$", color=c.gray, va="bottom")
    ax1.text(Vc_vdW + 0.008, 31, r"$V_{\mathrm{m},c}$", color=c.gray, va="bottom")
    ax1.text(0.30, 68.5, r"$T = T_c$", color=c.gray, va="bottom")
    ax1.text(0.295, 33.6, rf"$T = {Tr_loop:g}\,T_c$", color=c.blue, va="center")
    ax1.annotate("equal areas", xy=(0.5 * (V_liq + V_mid), P_tie - 4.5), xytext=(0.215, 38), color=c.ink,
                 ha="center", va="center", arrowprops=dict(arrowstyle="-", color=c.muted, lw=0.9, shrinkB=2))
    ax1.annotate("", xy=(0.5 * (V_mid + V_vap) - 0.02, P_tie + 1.2), xytext=(0.225, 39.5),
                 arrowprops=dict(arrowstyle="-", color=c.muted, lw=0.9, shrinkA=6, shrinkB=2))
    fs.panel(ax1, "a", "Maxwell construction")

    # (b) The regions of the diagram
    liquid, vapor, fluid = c.tint(c.liquid, 0.7), c.tint(c.gas, 0.7), c.fill
    top = 80.0
    left, right = V <= Vc_vdW, V >= Vc_vdW
    ax2.fill_between(V[left], np.minimum(P_crit[left], top), top, color=fluid, lw=0)
    ax2.fill_between(V[right], Pc_vdW, top, color=fluid, lw=0)
    ax2.fill_between(V[left], np.interp(V[left], dome_V_liq, dome_P, left=0), np.minimum(P_crit[left], top),
                     color=liquid, lw=0)
    ax2.fill_between(V[right], np.interp(V[right], dome_V_vap[::-1], dome_P[::-1], right=0), Pc_vdW, color=vapor, lw=0)
    # Inside the dome, the color follows the fraction of vapor given by the position along the tie line.
    grid_V, grid_P = np.meshgrid(np.linspace(0, 0.5, 500), np.linspace(30, Pc_vdW, 300))
    edge_liq, edge_vap = np.interp(grid_P, dome_P, dome_V_liq), np.interp(grid_P, dome_P, dome_V_vap)
    fraction = np.clip((grid_V - edge_liq) / np.maximum(edge_vap - edge_liq, 1e-9), 0, 1)
    image = np.zeros(grid_V.shape + (4,))
    image[..., :3] = (1 - fraction[..., None]) * np.array(to_rgb(liquid)) + fraction[..., None] * np.array(to_rgb(vapor))
    image[..., 3] = (grid_V >= edge_liq) & (grid_V <= edge_vap)
    ax2.imshow(image, extent=(0, 0.5, 30, Pc_vdW), origin="lower", aspect="auto", interpolation="bilinear", zorder=1)
    ax2.plot(dome_V, dome, color=c.muted, lw=0.9, zorder=2)
    ax2.plot(V[left], P_crit[left], color=c.gray, lw=1.2)
    ax2.plot([Vc_vdW, 0.5], [Pc_vdW, Pc_vdW], color=c.gray, lw=1.2)
    ax2.plot(Vc_vdW, Pc_vdW, "o", mfc=c.ground, mec=c.gray, mew=1.5, ms=6, zorder=5)
    ax2.text(0.063, 60, "Liquid", color=c.ink, rotation=90, ha="center", va="center")
    ax2.text(0.215, 44, "Liquid and vapor", color=c.ink, ha="center", va="center")
    ax2.text(0.42, 62, "Vapor", color=c.ink, ha="center", va="center")
    ax2.text(0.30, 77, "Supercritical fluid", color=c.ink, ha="center", va="top")
    ax2.annotate("critical point", xy=(Vc_vdW, Pc_vdW), xytext=(0.20, 67), color=c.ink, va="center",
                 arrowprops=dict(arrowstyle="-", color=c.muted, lw=0.9, shrinkB=5))
    fs.panel(ax2, "b", "Regions")

    ax1.set(xlim=(0, 0.5), ylim=(30, 80), xlabel=r"Molar volume $V_{\mathrm{m}}$ (L/mol)", ylabel="Pressure $P$ (bar)")
    ax2.set(xlabel=r"Molar volume $V_{\mathrm{m}}$ (L/mol)")
    return fig


fs.show(draw, alt="Two panels of pressure against molar volume for the van der Waals model of carbon dioxide. Panel a: "
                  "the isotherm at 0.9 times the critical temperature has a loop, a horizontal tie line cuts the loop "
                  "into two equal areas, and the critical isotherm passes through the critical point at the top of "
                  "the coexistence curve. Panel b: the same diagram divided into regions. Liquid lies at small "
                  "volume, vapor at large volume, liquid and vapor together under the coexistence curve, and the "
                  "supercritical fluid above the critical temperature and pressure.")
```

Within the loop, part of the van der Waals isotherm has a pressure that increases with volume, which is mechanically unstable and is not observed. The Maxwell construction replaces the loop with a horizontal tie line at the coexistence pressure, placed so that the two areas enclosed between the line and the loop are equal; the ends of the tie line give the molar volumes of the coexisting liquid and vapor.

##### Critical Point

Expanding Equation {eq}`van-der-waals-eos` and rearranging terms gives a cubic equation for $V_\text{m}$:

```{math}
:label: van-der-waals-eos-expanded
V_\text{m}^3 - \left( b_\text{m} + \frac{R T}{P} \right) V_\text{m}^2 + \frac{a_\text{m}}{P} V_\text{m} - \frac{a_\text{m} b_\text{m}}{P} = 0.
```

````{admonition} Complete Derivation of the Cubic Equation in $V_\text{m}$
:class: dropdown

**1. Expand Equation {eq}`van-der-waals-eos` and subtract $RT$ from both sides:**

```{math}
P V_\text{m} - P b_\text{m} + \frac{a_\text{m}}{V_\text{m}} - \frac{a_\text{m} b_\text{m}}{V_\text{m}^2} - R T = 0.
```

**2. Multiply both sides by $V_\text{m}^2 / P$:**

```{math}
V_\text{m}^3 - b_\text{m} V_\text{m}^2 + \frac{a_\text{m}}{P} V_\text{m} - \frac{a_\text{m} b_\text{m}}{P} - \frac{R T}{P} V_\text{m}^2 = 0.
```

**3. Group terms:**

```{math}
V_\text{m}^3 - \left( b_\text{m} + \frac{R T}{P} \right) V_\text{m}^2 + \frac{a_\text{m}}{P} V_\text{m} - \frac{a_\text{m} b_\text{m}}{P} = 0.
```
````

Below $T_c$, Equation {eq}`van-der-waals-eos-expanded` has three real roots for pressures within the van der Waals loop; as $T$ rises to $T_c$, these roots merge into one at the inflection point of the critical isotherm. At the critical point, the cubic therefore has a triple root at $V_{\text{m},c}$, so it can be written as

```{math}
\left( V_\text{m} - V_{\text{m},c} \right)^3
= V_\text{m}^3 - 3 V_{\text{m},c} V_\text{m}^2 + 3 V_{\text{m},c}^2 V_\text{m} - V_{\text{m},c}^3
= 0.
```

Matching coefficients yields

```{math}
V_{\text{m},c} = 3 b_\text{m}, \quad P_{\text{c}} = \frac{a_\text{m}}{27 b_\text{m}^2}, \quad T_{\text{c}} = \frac{8 a_\text{m}}{27 b_\text{m} R}.
```

````{admonition} Complete Derivation of the Critical Point
:class: dropdown

**1. Set coefficients equal and solve for powers of $V_{\text{m},c}$:**

```{math}
V_{\text{m},c} = \frac{1}{3} \left( b_\text{m} + \frac{R T_{\text{c}}}{P_{\text{c}}} \right), \quad
V_{\text{m},c}^2 = \frac{a_\text{m}}{3 P_{\text{c}}}, \quad
V_{\text{m},c}^3 = \frac{a_\text{m} b_\text{m}}{P_{\text{c}}}.
```

**2. Divide $V_{\text{m},c}^3$ by $V_{\text{m},c}^2$ to eliminate $P_{\text{c}}$:**

```{math}
\begin{aligned}
\frac{V_{\text{m},c}^3}{V_{\text{m},c}^2} &= \frac{a_\text{m} b_\text{m}}{P_{\text{c}}} \cdot \frac{3 P_{\text{c}}}{a_\text{m}} \\
V_{\text{m},c} &= 3 b_\text{m}.
\end{aligned}
```

**3. Substitute into $V_{\text{m},c}^3$ and solve for $P_{\text{c}}$:**

```{math}
\begin{aligned}
V_{\text{m},c}^3 &= \frac{a_\text{m} b_\text{m}}{P_{\text{c}}} \\
27 b_\text{m}^3 &= \frac{a_\text{m} b_\text{m}}{P_{\text{c}}} \\
P_{\text{c}} &= \frac{a_\text{m}}{27 b_\text{m}^2}.
\end{aligned}
```

**4. Substitute these results into $V_{\text{m},c}$ and solve for $T_{\text{c}}$:**

```{math}
\begin{aligned}
V_{\text{m},c} &= \frac{1}{3} \left( b_\text{m} + \frac{R T_{\text{c}}}{P_{\text{c}}} \right) \\
3 b_\text{m} &= \frac{1}{3} \left( b_\text{m} + R T_{\text{c}} \cdot \frac{27 b_\text{m}^2}{a_\text{m}} \right) \\
9 b_\text{m} &= b_\text{m} + R T_{\text{c}} \cdot \frac{27 b_\text{m}^2}{a_\text{m}} \\
8 b_\text{m} &= R T_{\text{c}} \cdot \frac{27 b_\text{m}^2}{a_\text{m}} \\
T_{\text{c}} &= \frac{8 a_\text{m}}{27 b_\text{m} R}.
\end{aligned}
```
````

##### Corresponding States

```{code-cell} ipython3
:label: fig:corresponding-states
:caption: Compressibility factor $Z$ of argon, oxygen, water, and carbon dioxide against the reduced pressure at reduced temperatures of (a) 1.0, (b) 1.2, and (c) 2.0. At the same reduced conditions, the four gases show similar trends, consistent with the principle of corresponding states. Water deviates most from the common curve.
:tags: [hide-input]

# Measured isotherms of four substances at three reduced temperatures (NIST Chemistry WebBook)
substances = {"argon": "Ar", "oxygen": r"O$_2$", "water": r"H$_2$O", "carbon-dioxide": r"CO$_2$"}
reduced = {(name, Tr): pd.read_csv(f"../_static/section-04/isothermal-properties-for-{name}-Tr={Tr}.csv")
           for name in substances for Tr in (1.0, 1.2, 2.0)}


def draw(c):
    colors = dict(zip(substances, [c.blue, c.vermilion, c.green, c.purple]))
    fig, axs = plt.subplots(1, 3, figsize=(fs.WIDTH, 3.1), constrained_layout=True, sharex=True, sharey=True)
    for ax, letter, Tr in zip(axs, "abc", (1.0, 1.2, 2.0)):
        ax.axhline(1, color=c.gray, ls="--", lw=1.0)
        for name, label in substances.items():
            data = reduced[name, Tr]
            Z = 1e5 * data["Pressure (bar)"] * data["Volume (m3/mol)"] / (R * data["Temperature (K)"])
            ax.plot(data["Pr"], Z, color=colors[name], lw=1.8, label=label)
        ax.set(xlabel="Reduced pressure $P_r$", xlim=(0, 8), ylim=(0, 1.2), xticks=[0, 2, 4, 6, 8])
        fs.panel(ax, letter, f"$T_r = {Tr}$")
    axs[0].set(ylabel="Compressibility factor $Z$")
    axs[0].legend(loc="lower right", handlelength=1.4)
    return fig


fs.show(draw, alt="Three panels of the compressibility factor against reduced pressure from 0 to 8 for argon, oxygen, "
                  "water, and carbon dioxide, at reduced temperatures of 1.0, 1.2, and 2.0. Within each panel, the "
                  "four curves nearly coincide. At a reduced temperature of 1.0, the curves drop steeply to about 0.2 "
                  "to 0.3 near a reduced pressure of 1 and then rise. At 1.2, the minimum is shallower, and at 2.0, "
                  "the curves stay close to 1.")
```

Introduce reduced variables $V_{\text{m},r}$, $P_r$, and $T_r$ via

```{math}
V_\text{m} = V_{\text{m},r} V_{\text{m},c}, \quad P = P_r P_c, \quad T = T_r T_c,
```

and substitute into Equation {eq}`van-der-waals-eos` to obtain

```{math}
:label: van-der-waals-eos-corresponding-states
\left( P_r + \frac{3}{V_{\text{m},r}^2} \right) \left( V_{\text{m},r} - \frac{1}{3} \right) = \frac{8}{3} T_r.
```

This corresponding-states form shows that all van der Waals fluids obey the same reduced EOS when expressed relative to their critical properties.

```{tip}
Calculate the compressibility factor of a van der Waals fluid at the critical point. Do you notice anything interesting?
```

````{admonition} Complete Derivation of the Corresponding States Equation
:class: dropdown

**1. Substitute the reduced-variable definitions into Equation {eq}`van-der-waals-eos`:**

```{math}
\left( P_r P_c + \frac{a_\text{m}}{V_{\text{m},r}^2 V_{\text{m},c}^2} \right) \left( V_{\text{m},r} V_{\text{m},c} - b_\text{m} \right) = R T_r T_c.
```

**2. Divide both sides by $P_c V_{\text{m},c}$:**

```{math}
\left( P_r + \frac{a_\text{m}}{V_{\text{m},r}^2 V_{\text{m},c}^2 P_c} \right) \left( V_{\text{m},r} - \frac{b_\text{m}}{V_{\text{m},c}} \right) = \frac{R T_c}{P_c V_{\text{m},c}} T_r.
```

**3. Substitute $P_c$, $V_{\text{m},c}$, and $T_c$ in terms of $a_\text{m}$, $b_\text{m}$, and $R$:**

```{math}
\begin{aligned}
\left( P_r + \frac{a_\text{m}}{V_{\text{m},r}^2} \cdot \frac{1}{9 b_\text{m}^2} \cdot \frac{27 b_\text{m}^2}{a_\text{m}} \right) \left( V_{\text{m},r} - \frac{b_\text{m}}{3 b_\text{m}} \right) &= R \cdot \frac{8 a_\text{m}}{27 b_\text{m} R} \cdot \frac{27 b_\text{m}^2}{a_\text{m}} \cdot \frac{1}{3 b_\text{m}} T_r \\
\left( P_r + \frac{3}{V_{\text{m},r}^2} \right) \left( V_{\text{m},r} - \frac{1}{3} \right) &= \frac{8}{3} T_r.
\end{aligned}
```

````

### Computational Studio: Real Gas

In [](#fig:real-gas-studio), you can choose a gas and a temperature and compare the ideal-gas law with the van der Waals EOS and with the Peng–Robinson EOS, a later equation of the same kind.

```{iframe} https://wexlergroup.github.io/chem-4020-5020/computational-studios/real-gas-studio.html
:label: fig:real-gas-studio
:width: 100%
:title: Real gas
:placeholder: ../_static/real-gas-studio/real-gas-studio.png

Compressibility factor and isotherm of a gas according to the ideal-gas law, the van der Waals EOS (Eq. [](#van-der-waals-eos)), and the Peng–Robinson EOS. Panel (a) shows $Z$ against the pressure, and panel (b) shows the pressure against the molar volume on a logarithmic scale. The constants of both EOS follow from the critical temperature $T_c$ and the critical pressure $P_c$ of the gas. Below $T_c$, the thick lines follow the Maxwell construction of [](#fig:maxwell-construction), and the thin lines in panel (b) are the loops that it replaces.
```

You can also open the studio in a new tab: [Real Gas Studio](https://wexlergroup.github.io/chem-4020-5020/computational-studios/real-gas-studio.html).

## Worked Example

### van der Waals critical constants and $Z_c$ for CO$_2$

Use the van der Waals parameters used in the section plots:
$a_\text{m} = 3.6551\ \mathrm{L^2\,bar\,mol^{-2}}$, $b_\text{m}=0.042816\ \mathrm{L\,mol^{-1}}$, and $R=0.08314\ \mathrm{L\,bar\,mol^{-1}\,K^{-1}}$.

The van der Waals critical constants are

```{math}
V_{\text{m},c}=3b_\text{m},\qquad
P_c=\frac{a_\text{m}}{27b_\text{m}^2},\qquad
T_c=\frac{8a_\text{m}}{27b_\text{m}R}.
```

1. **Critical volume**

   ```{math}
   V_{\text{m},c}=3(0.042816)=0.128448\ \mathrm{L/mol}.
   ```

2. **Critical pressure**

   ```{math}
   P_c=\frac{3.6551}{27(0.042816)^2}=73.85\ \mathrm{bar}.
   ```

3. **Critical temperature**

   ```{math}
   T_c=\frac{8(3.6551)}{27(0.042816)(0.08314)}=304.2\ \mathrm{K}.
   ```

4. **Critical compressibility factor**

   ```{math}
   Z_c \equiv \frac{P_cV_{\text{m},c}}{RT_c}
   =\frac{73.85\times 0.12845}{0.08314\times 304.2}
   =0.375=\frac{3}{8}.
   ```

5. **Comparison with experiment**

   The predicted $T_c$ and $P_c$ are within 0.1% of the measured values for CO$_2$, $304.1\ \mathrm{K}$ and $73.8\ \mathrm{bar}$, as expected, because van der Waals parameters are usually determined from the measured critical temperature and pressure. The measured critical molar volume, $0.094\ \mathrm{L/mol}$, and the corresponding $Z_c\approx 0.27$ are both smaller than the van der Waals predictions.

**Result.** The van der Waals model predicts $T_c\approx 304\ \mathrm{K}$, $P_c\approx 74\ \mathrm{bar}$, $V_{\text{m},c}\approx 0.128\ \mathrm{L/mol}$, and $Z_c=3/8$.

## Concept Checks

1. Why can a gas be “coincidentally ideal” (i.e., $Z=1$) even when interactions are not negligible?
2. What physical effects are represented by the $a/V_{\mathrm m}^2$ and $b$ corrections in the van der Waals EOS?
3. At fixed $T$, how does $Z$ typically change as $P$ increases, and why?
4. Why does the reduced-variable form imply similar behavior across different substances at corresponding reduced conditions?

## Key Takeaways

- Nonideality is quantified by $Z=PV/(nRT)$; $Z<1$ suggests net attractions and $Z>1$ suggests net repulsions.
- The van der Waals EOS introduces an attraction parameter $a$ and an excluded-volume parameter $b$.
- Critical constants follow directly from the cubic form of the van der Waals EOS, with $Z_c=3/8$.
- Using reduced variables $(T_r,P_r,V_{\text{m},r})$ yields the corresponding-states form shared by all van der Waals fluids.
