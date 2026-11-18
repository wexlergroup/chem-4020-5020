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
  enumerator: 1.3.%s
---

# 1.3. Ideal Gases

See also: [Course-wide Conventions & Notation](../notation.md)

Run the code: [Notebook for Section 1.3](https://wexlergroup.github.io/chem-4020-5020/lite/notebooks/index.html?path=chapter-01/section-03.ipynb), in which you can change and rerun the Python cells of this section in your browser.

## Overview and Learning Objectives

The ideal-gas model consolidates several empirical gas laws into a single equation of state and is accurate when intermolecular interactions are negligible. In this section, we review the classical gas laws, derive $PV=Nk_{\mathrm{B}}T=nRT$, use the result to estimate microscopic length scales, and see how it motivates an absolute temperature scale.

---

Learning objectives:

- State Boyle’s, Charles’s, Gay–Lussac’s, and Avogadro’s laws and identify the variables held constant in each.
- Derive the ideal-gas equation of state $PV=Nk_{\mathrm{B}}T=nRT$ by combining the gas laws.
- Explain the microscopic assumptions behind ideal-gas behavior and when they are expected to hold.
- Use $V/N=(k_{\mathrm{B}}T/P)$ to estimate number density and mean intermolecular spacing.

## Core Ideas and Derivations

(section:gas-laws)=
### Gas Laws

```{code-cell} ipython3
:label: fig:gas-laws
:caption: The four classical gas laws in chronological order: (a) Boyle's law ([Boyle's original data](https://web.lemoyne.edu/~giunta/classicalcs/boyleform.html)), (b) Charles's law, (c) Gay–Lussac's law, and (d) Avogadro's law. Each curve is calculated for an ideal gas with the two named variables held constant.
:tags: [hide-input]

import sys

import matplotlib.pyplot as plt
import numpy as np
from scipy.constants import R

sys.path.append("..")            # the book folder, which holds figstyle
import figstyle as fs

R_L_bar = R / 100                      # gas constant in L bar/(mol K)
T0, P0, V0 = 290.0, 1.0, 24.53         # the temperature (K), pressure (bar), and volume (L) that are held constant


def draw(c):
    T = np.linspace(273.15, 373.15, 101)
    V = np.linspace(1, 40, 400)
    n = np.linspace(0, 10, 101)        # amount in moles, which is the number of particles in units of N_A
    panels = [
        ("a", "Boyle (1662)", V, R_L_bar * T0 / V, "Volume $V$ (L)", "Pressure $P$ (bar)",
         r"$P = c_{\mathrm{B}}/V$", "constant $T$ and $N$", (0.30, 0.52)),
        ("b", "Charles (1787)", T, R_L_bar * T / P0, "Temperature $T$ (K)", "Volume $V$ (L)",
         r"$V = c_{\mathrm{C}}\,T$", "constant $P$ and $N$", (0.08, 0.80)),
        ("c", "Gay-Lussac (1802)", T, R_L_bar * T / V0, "Temperature $T$ (K)", "Pressure $P$ (bar)",
         r"$P = c_{\mathrm{GL}}\,T$", "constant $V$ and $N$", (0.08, 0.80)),
        ("d", "Avogadro (1811)", n, R_L_bar * T0 / P0 * n, r"Number of particles $N$ ($N_{\mathrm{A}}$)", "Volume $V$ (L)",
         r"$V = c_{\mathrm{A}}\,N$", "constant $P$ and $T$", (0.08, 0.80)),
    ]
    fig, axs = plt.subplots(2, 2, figsize=(fs.WIDTH, 5.4), constrained_layout=True)
    for ax, (letter, title, x, y, xlabel, ylabel, law, held, (tx, ty)) in zip(axs.flat, panels):
        ax.plot(x, y, color=c.blue)
        ax.set(xlabel=xlabel, ylabel=ylabel)
        ax.text(tx, ty, law, transform=ax.transAxes, color=c.blue, fontsize=fs.BASE + 1.5, va="bottom")
        ax.text(tx, ty - 0.03, held, transform=ax.transAxes, color=c.muted, va="top")
        fs.panel(ax, letter, title)
    axs[0, 0].set(xlim=(0, 40), ylim=(0, 25))
    axs[1, 1].set(xlim=(0, 10), ylim=(0, 250))
    return fig


fs.show(draw, alt="Four panels, one for each gas law. Panel a, Boyle: pressure falls as one over the volume at constant "
                  "temperature and number of particles. Panel b, Charles: volume rises in a straight line with "
                  "temperature at constant pressure and number of particles. Panel c, Gay-Lussac: pressure rises in "
                  "a straight line with temperature at constant volume and number of particles. Panel d, Avogadro: "
                  "volume rises in a straight line with the number of particles at constant pressure and temperature.")
```

Written as equations, the four laws are:

````{list-table} Classical gas laws
:header-rows: 1
:name: classical-gas-laws

* - Law
  - Relation
  - Held constant
* - Boyle's law
  - $P = c_{\mathrm{B}} / V$
  - $T$, $N$
* - Charles's law
  - $V = c_{\mathrm{C}}\, T$
  - $P$, $N$
* - Gay–Lussac's law
  - $P = c_{\mathrm{GL}}\, T$
  - $V$, $N$
* - Avogadro's law
  - $V = c_{\mathrm{A}}\, N$
  - $P$, $T$
````

Here, $c_{\mathrm{B}}$, $c_{\mathrm{C}}$, $c_{\mathrm{GL}}$, and $c_{\mathrm{A}}$ are constants whose values depend on the variables held constant. Avogadro's law holds with the same constant $c_{\mathrm{A}}$ for every gas: at a given pressure and temperature, equal volumes of different gases contain equal numbers of particles.

Real gases *approximately* follow these relationships under conditions that favor ideal behavior:

- relatively low pressure,
- relatively low density, and
- sufficiently high temperature.

These conditions collectively minimize intermolecular interactions, allowing the gas to behave *ideally*.

---

### Deriving an Equation of State

In [Section 1.1](section-01.md), we defined an **equation of state** as a mathematical relationship among state variables. Each gas law above relates *two* state variables (among $P, V, T, N$) under conditions where the other two are held constant. To obtain a single equation of state relating all four variables, we can combine these laws using multivariate calculus (as covered in [Math 233](https://math.wustl.edu/calculus-iii)).

#### Total Differential

The total differential $df$ of a function $f$ of $m$ variables $x_1, \ldots, x_m$ is

```{math}
:label: eq-total-differential
df = \sum_{i=1}^m \left( \frac{\partial f}{\partial x_i} \right)_{\{ x_j | j \neq i \}} dx_i,
```

where $\left(\partial f / \partial x_i\right)_{\{ x_j | j \neq i \}}$ denotes the partial derivative of $f$ with respect to $x_i$, holding all other variables constant, and $dx_i$ is an infinitesimal change in $x_i$.

#### Total Differential of Volume

From Boyle's, Charles's, and Avogadro's laws, we can treat $V$ as a function of $P, T,$ and $N$:

```{math}
V = V(P, T, N).
```

Applying Eq. {eq}`eq-total-differential` to $V$ gives

```{math}
dV = \left( \frac{\partial V}{\partial P} \right)_{T,N} dP \;+\; \left( \frac{\partial V}{\partial T} \right)_{P,N} dT \;+\; \left( \frac{\partial V}{\partial N} \right)_{P,T} dN.
```

#### Partial Derivatives via Gas Laws

Each partial derivative follows from the gas law in {numref}`classical-gas-laws` that holds the same two variables constant. Solving Boyle's law for the volume gives $V = c_{\mathrm{B}}/P$, so that

```{math}
\left( \frac{\partial V}{\partial P} \right)_{T,N} = -\frac{c_{\mathrm{B}}}{P^2} = -\frac{V}{P}.
```

Charles's law, $V = c_{\mathrm{C}}\,T$, and Avogadro's law, $V = c_{\mathrm{A}}\,N$, give

```{math}
\left( \frac{\partial V}{\partial T} \right)_{P,N} = c_{\mathrm{C}} = \frac{V}{T}
\qquad \text{and} \qquad
\left( \frac{\partial V}{\partial N} \right)_{P,T} = c_{\mathrm{A}} = \frac{V}{N}.
```

Substituting these three derivatives into the total differential of the volume gives

```{math}
dV \;=\; -\frac{V}{P}\,dP \;+\; \frac{V}{T}\,dT \;+\; \frac{V}{N}\,dN.
```

Dividing by $V$ and rearranging (i.e., using logarithmic differentials) yields

```{math}
d \ln P \;+\; d \ln V \;=\; d \ln T \;+\; d \ln N.
```

#### Integrating the Total Differential

Integrating each term from an initial state $(P_i, V_i, T_i, N_i)$ to a final state $(P_f, V_f, T_f, N_f)$ gives

```{math}
\ln \frac{P_f}{P_i} \;+\; \ln \frac{V_f}{V_i} \;=\; \ln \frac{T_f}{T_i} \;+\; \ln \frac{N_f}{N_i}.
```

Combining the logarithms on each side,

```{math}
\ln \frac{P_f V_f}{P_i V_i} \;=\; \ln \frac{T_f N_f}{T_i N_i},
```

and exponentiating and rearranging gives

```{math}
\frac{P_f V_f}{N_f T_f} \;=\; \frac{P_i V_i}{N_i T_i}.
```

Because the initial and final states are arbitrary, the ratio $PV/(NT)$ must be a constant. Since Avogadro's law holds with the same constant for every gas, the value of this ratio is also the same for every gas. Denoting it by $k_{\mathrm{B}}$, we arrive at the **ideal-gas equation of state**:

```{math}
:label: ideal-gas-equation-of-state
PV \;=\; N k_{\mathrm{B}} T \;=\; nRT,
```

where $n = N / N_{\mathrm{A}}$ is the number of moles and $R = k_{\mathrm{B}} N_{\mathrm{A}}$ is the [molar gas constant](https://physics.nist.gov/cgi-bin/cuu/Value?r) ($R = 8.314\,\mathrm{J\,mol^{-1}\,K^{-1}}$).

### The $P$–$V$–$T$ Surface

For a fixed number of particles, Eq. {eq}`ideal-gas-equation-of-state` assigns one pressure to every pair of volume and temperature, so the states of an ideal gas form a surface. Each gas law of [](#fig:gas-laws) that holds $N$ constant is a cut through this surface: Boyle's law at constant $T$, Charles's law at constant $P$, and Gay–Lussac's law at constant $V$. In [](#fig:pvt-surface), choose the variable to hold constant, and move the state along the corresponding cut.

```{iframe} https://wexlergroup.github.io/chem-4020-5020/computational-studios/pvt-surface.html
:label: fig:pvt-surface
:width: 100%
:title: The P-V-T surface of an ideal gas, with the gas laws as cuts through it
:placeholder: ../_static/pvt-surface/pvt-surface.png

The pressure of one mole of an ideal gas as a surface over the volume and the temperature, from Eq. [](#ideal-gas-equation-of-state). The three curves through the chosen state are the cuts at constant $T$ (Boyle), at constant $P$ (Charles), and at constant $V$ (Gay–Lussac), and the small plots show the same curves in two dimensions.
```

Two things to try:

1. Hold the temperature constant, and halve the volume. Then hold the pressure constant, and halve the volume again. Compare how the state moves on the surface in the two cases.
2. Bring the gas to 2.00 bar and 10.0 L along two different routes, and compare the temperatures that you reach.

---

### Ideal-Gas Assumptions

A gas described by Eq. {eq}`ideal-gas-equation-of-state` is called *ideal* because, under conditions of low pressure, low density, and high temperature, we can adopt the simplifying assumptions of kinetic theory:

1. Each particle’s size is negligible compared to the average distance between particles (point particles).
2. Particles experience no intermolecular forces except during elastic collisions.
3. Collisions conserve total energy and momentum.

Many thermodynamic properties follow cleanly from these assumptions.

---

### Estimating Particle Distances

Rearranging Eq. {eq}`ideal-gas-equation-of-state`, a simple estimate for a typical intermolecular spacing is

```{math}
:label: average-distance-between-particles
\langle d \rangle = \left(\frac{V}{N}\right)^{1/3} \;=\; \left(\frac{k_{\mathrm{B}} T}{P}\right)^{1/3}.
```

This expression shows that increasing temperature or reducing pressure *increases* the typical separation between gas particles.

[](#fig:gas-to-scale) shows this spacing to scale for nitrogen at $300\ \mathrm{K}$ and $1\ \mathrm{bar}$, where Eq. {eq}`average-distance-between-particles` gives $\langle d \rangle = 3.46\ \mathrm{nm}$ (see the worked example below). The molecules are about ten molecular diameters apart, which is why their own volume and their interactions can be neglected. In the liquid, the same molecules touch.

```{code-cell} ipython3
:label: fig:gas-to-scale
:caption: Nitrogen molecules drawn to scale as discs with a diameter of 0.36 nm, in slices 10 nm wide. (a) The gas at 300 K and 1 bar, in a slice that is one mean spacing thick. The bar shows the mean spacing $\langle d \rangle$. (b) The liquid at its boiling point, in a slice that is one molecule thick.
:tags: [hide-input]

from matplotlib.patches import Circle, Rectangle
from scipy.constants import N_A, k as k_B

window = 10.0                                              # width of each slice (nm)
diameter = 0.364                                           # kinetic diameter of N2 (nm)
d_gas = (k_B * 300 / 1e5) ** (1 / 3) * 1e9                 # mean spacing at 300 K and 1 bar (nm), from the equation above
d_liquid = (28.0134e-3 / (806.2 * N_A)) ** (1 / 3) * 1e9   # liquid N2 at its boiling point, 806.2 kg/m^3 (nm)

rng = np.random.default_rng(4020)

# Gas: a slice one mean spacing thick holds (window / d_gas)^2 molecules. They are placed at random.
gas = []
while len(gas) < round((window / d_gas) ** 2):
    point = rng.uniform(0.4, window - 0.4, 2)
    if all(np.hypot(*(point - other)) > 1.5 for other in gas):
        gas.append(point)

# Liquid: one layer of molecules in contact, as a hexagonal arrangement with small random displacements.
spacing = 2 ** (1 / 6) * d_liquid
liquid = [(x + (row % 2) * spacing / 2, row * spacing * np.sqrt(3) / 2)
          for row in range(int(window / (spacing * np.sqrt(3) / 2)) + 1)
          for x in np.arange(0, window, spacing)]
liquid = np.array(liquid) + rng.normal(0, 0.025, (len(liquid), 2))


def draw_to_scale(c):
    fig, axs = plt.subplots(1, 2, figsize=(fs.WIDTH, 3.75), constrained_layout=True)
    for ax, letter, title, points in zip(axs, "ab", ["Gas at 300 K and 1 bar", "Liquid"], [gas, liquid]):
        ax.add_patch(Rectangle((0, 0), window, window, fc="none", ec=c.muted, lw=0.9, zorder=3))
        for x, y in points:
            ax.add_patch(Circle((x, y), diameter / 2, fc=c.blue, ec="none", clip_path=None))
        for patch in ax.patches[1:]:
            patch.set_clip_path(ax.patches[0])
        ax.set(xlim=(-0.1, window + 0.1), ylim=(-1.6, window + 0.1), aspect="equal")
        ax.axis("off")
        fs.panel(ax, letter, title)
    for ax, length, label in [(axs[0], d_gas, rf"$\langle d \rangle$ = {d_gas:.2f} nm"),
                              (axs[1], 1.0, "1 nm")]:
        ax.annotate("", xy=(0, -0.7), xytext=(length, -0.7),
                    arrowprops=dict(arrowstyle="|-|", color=c.ink, lw=1.1, mutation_scale=3, shrinkA=0, shrinkB=0))
        ax.text(length + 0.3, -0.7, label, va="center", ha="left")
    return fig


fs.show(draw_to_scale,
        alt="Two square slices, each 10 nanometers wide, with nitrogen molecules drawn to scale as small discs. In "
            f"the gas at 300 K and 1 bar, the slice holds {len(gas)} widely separated molecules, and a bar marks the "
            f"mean spacing of {d_gas:.2f} nanometers. In the liquid, the slice is filled with molecules in contact.")
```

### Absolute Temperature Scale

```{code-cell} ipython3
:label: fig:absolute-zero
:caption: Molar volume of O$_2$ at 1 Pa against the Celsius temperature (points, from the [NIST Chemistry WebBook](https://webbook.nist.gov/cgi/inchi/InChI%3D1S/O2/c1-2)). The line is an ordinary least-squares fit to the measurements. Extended to zero volume (dashed), it reaches the temperature axis at $-273.15\,{}^{\circ}\mathrm{C}$, which is absolute zero (0 K).
:tags: [hide-input]

import pandas as pd

# Measured molar volume of O2 at a pressure of 1 Pa between 0 and 100 °C (NIST Chemistry WebBook)
data = pd.read_table("../_static/section-03/isobaric-properties-for-oxygen.tsv", sep="\t")
T_C, V_m = data["Temperature (C)"].to_numpy(), data["Volume (m3/mol)"].to_numpy()

# Straight line through the measurements: the temperature at which the line reaches zero volume is the intercept
slope, intercept = np.polyfit(V_m, T_C, 1)


def draw(c):
    fig, ax = plt.subplots(figsize=(6, 3.8), constrained_layout=True)
    ax.plot([intercept, T_C[0]], [0, V_m[0]], color=c.vermilion, ls="--", lw=1.6)
    ax.plot(T_C, V_m, color=c.vermilion, lw=1.6)
    ax.plot(T_C[::5], V_m[::5], "o", color=c.blue, ms=4.5)
    ax.plot(intercept, 0, "o", mfc=c.ground, mec=c.vermilion, mew=1.6, ms=7, clip_on=False, zorder=5)
    label = f"{intercept:.2f} °C".replace("-", "−")
    ax.text(intercept + 24, 40, label, color=c.vermilion, ha="left", va="bottom")
    ax.text(50, 2150, "measured", color=c.blue, ha="center", va="top")
    ax.text(-150, 1250, "straight line,\nextended", color=c.vermilion, ha="right", va="bottom")
    ax.set(xlabel="Temperature (°C)", ylabel=r"Molar volume (m$^3$/mol)", xlim=(-300, 110), ylim=(0, 3300))
    return fig


fs.show(draw, alt="Molar volume of oxygen at 1 pascal against temperature in degrees Celsius. Measured points between "
                  "0 and 100 degrees fall on a straight line. The line, extended to lower temperatures, reaches "
                  f"zero volume at {intercept:.2f} degrees Celsius.".replace("-", "minus "))
```

Charles's law states that, at constant (sufficiently low) pressure, the volume per particle is directly proportional to the gas temperature. If temperature is reported on the Celsius scale ($0\,{}^{\circ}\mathrm{C}$ at the freezing point of water and $100\,{}^{\circ}\mathrm{C}$ at the boiling point), plotting $T(\,{}^{\circ}\mathrm{C})$ versus $v$ gives an approximately linear relationship. Extrapolating this line to the (unphysical) point where the volume would vanish identifies a theoretical lower bound at $-273.15\,{}^{\circ}\mathrm{C}$. Shifting the Celsius scale by $273.15$ degrees places this bound at zero:

```{math}
T(\mathrm{K}) \;=\; T(\,{}^{\circ}\mathrm{C}) \;+\; 273.15.
```

This defines the **Kelvin scale**, an absolute temperature scale that begins at the lowest physically meaningful temperature.

## Worked Example

### Mean intermolecular spacing at 1 bar

Estimate the mean spacing $\langle d\rangle$ between molecules in an ideal gas at $T=300\ \mathrm{K}$ and $P=1.00\ \mathrm{bar}$.

Using

```{math}
\langle d\rangle=\left(\frac{V}{N}\right)^{1/3}=\left(\frac{k_{\mathrm{B}}T}{P}\right)^{1/3},
```

with $k_{\mathrm{B}}=1.38065\times10^{-23}\ \mathrm{J/K}$ and $P=1.00\times10^{5}\ \mathrm{Pa}$:

1. **Compute $k_{\mathrm{B}}T/P$**

   ```{math}
   \frac{k_{\mathrm{B}}T}{P}=\frac{(1.38065\times10^{-23})(300)}{1.00\times10^{5}}
   =4.14\times10^{-26}\ \mathrm{m^3}.
   ```

2. **Take the cube root**

   ```{math}
   \langle d\rangle=(4.14\times10^{-26})^{1/3}
   \approx 3.46\times10^{-9}\ \mathrm{m}
   =3.46\ \mathrm{nm}.
   ```

3. **Compute the number density**

   The number density is the reciprocal of the volume per molecule:

   ```{math}
   \frac{N}{V}=\frac{P}{k_{\mathrm{B}}T}=\frac{1.00\times10^{5}}{(1.38065\times10^{-23})(300)}
   =2.41\times10^{25}\ \mathrm{m^{-3}}.
   ```

**Result.** At $300\ \mathrm{K}$ and $1\ \mathrm{bar}$, an ideal gas contains about $2.4\times10^{25}$ molecules per cubic meter, and the molecules are typically separated by a few nanometers.

## Concept Checks

1. Why does the ratio $PV/(NT)$ have to be constant if the total-differential argument holds for arbitrary initial/final states?
2. Which change (increasing $T$ or increasing $P$) makes a gas *less* ideal, and why?
3. Why must temperature be measured on an absolute (Kelvin) scale in the ideal-gas law?
4. How would $\langle d\rangle$ scale if pressure increased by a factor of 8 at fixed $T$?

## Key Takeaways

- The classical gas laws combine into the ideal-gas equation of state $PV=Nk_{\mathrm{B}}T=nRT$.
- Ideal behavior is favored by **low density/pressure** and **high temperature**, where interactions are small.
- The mean molecular spacing scales as $\langle d\rangle\propto (T/P)^{1/3}$.
- The zero of the Kelvin scale is the temperature at which the extrapolated volume of a gas at low pressure vanishes, $-273.15\,{}^{\circ}\mathrm{C}$.
