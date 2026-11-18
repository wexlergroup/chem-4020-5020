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
  enumerator: 2.5.%s
---


# 2.5. Particle in a Box

See also: [Course-wide Conventions & Notation](../notation.md)

Run the code: [Notebook for Section 2.5](https://wexlergroup.github.io/chem-4020-5020/lite/notebooks/index.html?path=chapter-02/section-05.ipynb), in which you can change and rerun the Python cells of this section in your browser.

## Overview and Learning Objectives

A particle in a box is the simplest quantized model of translation and provides a clean route to the translational partition function. We review the 1D and 3D energy levels, derive $q_{\mathrm{trans}}=V/\Lambda^3$ in the high-temperature/large-volume limit, and use $Q=q^N/N!$ to recover ideal-gas thermodynamic relations.

Learning objectives:

- Write the quantized energy levels for a particle in a 1D box and generalize to a 3D rectangular box and cube.
- Derive the single-particle partition function and its continuum approximation in the limit of small level spacing.
- Define the thermal de Broglie wavelength $\Lambda$ and express $q_{\mathrm{trans}}=V/\Lambda^3$.
- Use $Q=(V/\Lambda^3)^N/N!$ to obtain $U=\tfrac{3}{2}Nk_{\mathrm B}T$, $C_V=\tfrac{3}{2}Nk_{\mathrm B}$, and $P=Nk_{\mathrm B}T/V$.

## Core Ideas and Derivations

### Review of the Particle in a Box

#### Particle in a One-Dimensional Box

[](#fig:particle-in-a-box) shows the first four energy levels and wavefunctions of a particle in a one-dimensional box.

```{code-cell} ipython3
:label: fig:particle-in-a-box
:caption: Energy levels $E_n$ (gray lines) and wavefunctions $\psi_n$ (blue curves) of a particle in a one-dimensional box with walls at $x=\pm L/2$. Each wavefunction is drawn on its energy level, and the energy increases as $n^2$.
:tags: [hide-input]

import sys

import matplotlib.pyplot as plt
import numpy as np

sys.path.append("..")            # the book folder, which holds figstyle
import figstyle as fs


def draw(c):
    fig, ax = plt.subplots(figsize=(4.8, 3.7), constrained_layout=True)
    x = np.linspace(-0.5, 0.5, 200)          # position in units of L
    levels = range(1, 5)

    # Energy levels in units of E_1 = h^2 / (8 m L^2), with each wavefunction drawn on its level
    for n in levels:
        ax.plot([-0.5, 0.5], [n**2, n**2], color=c.gray, lw=1.2)
        ax.plot(x, n**2 + np.sin(n * np.pi * (x + 0.5)), color=c.blue)
        ax.text(0.53, n**2, rf"$\psi_{n}$", color=c.blue, va="center", clip_on=False)

    # Walls and floor of the box
    ax.plot([-0.5, -0.5, 0.5, 0.5], [18, 0, 0, 18], color=c.ink, lw=1.6, solid_capstyle="butt", zorder=3)

    ax.set(xlabel="Position $x$", ylabel="Energy", xlim=(-0.5, 0.5), ylim=(0, 18),
           xticks=[-0.5, 0, 0.5], yticks=[n**2 for n in levels])
    ax.set_xticklabels(["$-L/2$", "0", "$L/2$"])
    ax.set_yticklabels(["$E_1$", "$E_2 = 4E_1$", "$E_3 = 9E_1$", "$E_4 = 16E_1$"])
    ax.tick_params(axis="y", length=0, pad=6)
    for side in ("left", "bottom"):
        ax.spines[side].set_visible(False)
    return fig


fs.show(draw, alt="Box with walls at x equal to minus L over 2 and L over 2. Four horizontal energy levels lie at 1, "
                  "4, 9, and 16 times the lowest energy. On each level, the wavefunction is drawn as a sine curve "
                  "that is zero at both walls, with one half-wave for the first level, two for the second, three for "
                  "the third, and four for the fourth.")
```

For a particle in a one-dimensional box of length $L$, the quantized energy levels are

```{math}
E_n = \frac{h^2}{8 m L^2} \, n^2 \quad \text{for} \; n = 1, 2, 3, \ldots,
```

where $h$ is [Planck’s constant](https://physics.nist.gov/cgi-bin/cuu/Value?h), $m$ is the particle mass, and $n$ is a positive integer.

#### Particle in a Three-Dimensional Box

In three dimensions, the energy levels are

```{math}
E_{n_x, n_y, n_z} = \frac{h^2}{8 m}
\left( \frac{n_x^2}{L_x^2} + \frac{n_y^2}{L_y^2} + \frac{n_z^2}{L_z^2} \right),
```

where $n_x$, $n_y$, and $n_z$ are positive integers, and $L_x$, $L_y$, and $L_z$ are the respective side lengths of a rectangular box.

##### Particle in a Cube

For a cube of side $L$, $L_x = L_y = L_z = L$. The energy simplifies to

```{math}
E_{n_x, n_y, n_z} 
= \frac{h^2}{8mL^2} \left( n_x^2 + n_y^2 + n_z^2 \right)
= \frac{h^2}{8 m V^{2/3}}
\left( n_x^2 + n_y^2 + n_z^2 \right),
```

where $V = L^3$ is the cube volume.

### Partition Function for a Particle in a Cube

The single-particle partition function for a three-dimensional cube is

```{math}
\begin{aligned}
q 
&= \sum_{n_x, n_y, n_z} 
   \exp \left[-\beta E_{n_x, n_y, n_z}\right] \\
&= \sum_{n_x=1}^{\infty} \sum_{n_y=1}^{\infty} \sum_{n_z=1}^{\infty} 
   \exp \left(-\beta \frac{h^2}{8 m L^2} (n_x^2 + n_y^2 + n_z^2)\right).
\end{aligned}
```

Define $\alpha = \beta h^2 / (8 m V^{2/3})$. Then

```{math}
:label: eq:partition-function-cube
q 
= \left( \sum_{n=1}^{\infty} e^{-\alpha n^2} \right)^3.
```

#### Approximation of the Sum

The parameter $\alpha$ is the ratio of $h^2/(8 m V^{2/3})$, the energy that sets the spacing between the levels, to $k_{\mathrm B}T$. For helium at 300 K in a volume of 1.00 L, $\alpha = 2.0\times10^{-19}$ (see the worked example below). When $\alpha \ll 1$ (large $T$ or large $V$), the sum $\sum_{n=1}^{\infty} e^{-\alpha n^2}$ is well approximated by an integral:

```{math}
\sum_{n = 1}^{\infty} e^{-\alpha n^2} 
\approx \int_0^{\infty} e^{-\alpha x^2} \, dx 
= \left(\frac{\pi}{4 \alpha}\right)^{\!\!1/2}.
```

Substituting $\alpha = \beta h^2/(8 m V^{2/3})$ gives

```{math}
\left(\frac{\pi}{4 \alpha}\right)^{\!\!1/2}
= \left(\frac{\pi}{4} \cdot \frac{8 m V^{2/3}}{\beta h^2}\right)^{\!\!1/2}
= \left(\frac{2 \pi m}{\beta h^2}\right)^{\!1/2} V^{1/3}.
```

We often define the thermal de Broglie wavelength $\Lambda$ as

```{math}
\Lambda 
= \left(\frac{h^2}{2 \pi m k_{\mathrm B} T}\right)^{\!1/2}.
```

````{admonition} What is the thermal de Broglie wavelength?
:class: tip

The thermal de Broglie wavelength, $\Lambda$, is a characteristic length scale for the quantum mechanical delocalization of a particle at temperature $T$. When the typical distance between particles is much larger than $\Lambda$, the system behaves classically, and dividing by $N!$ ([Section 2.4](section-04.md)) is the only correction needed for the indistinguishability of the particles. Conversely, when the interparticle spacing is on the order of $\Lambda$ (or smaller), quantum effects become significant, this simple correction fails, and the full quantum treatment of identical particles is required.
````

Hence,

```{math}
\sum_{n = 1}^{\infty} e^{-\alpha n^2} 
\approx \frac{V^{1/3}}{\Lambda}.
```

Cubing this result, as in Eq. {eq}`eq:partition-function-cube`, gives the single-particle partition function. Because the error of the approximation is negligible when $\alpha \ll 1$, the result is written as an equality:

```{math}
q 
= \frac{V}{\Lambda^3}.
```

### Partition Function for $N$ Particles in a Cube

For $N$ identical, non-interacting, indistinguishable particles, [Section 2.4](section-04.md) gives the total partition function $Q$ as

```{math}
:label: eq:partition-function-many-particles-cube
Q 
= \frac{q^N}{N!} 
= \frac{1}{N!}\left(\frac{V}{\Lambda^3}\right)^{\!N}.
```

The factor of $N!$ corrects for particle indistinguishability: permutations of identical particles do not produce new states.

### Ensemble Averages for $N$ Particles in a Cube

#### Natural Logarithm of the Partition Function

Taking the logarithm of $Q$ gives

```{math}
\ln Q 
= N \ln \left(\frac{V}{\Lambda^3}\right) - \ln N!.
```

Often, the Stirling approximation ($\ln N! \approx N \ln N - N$) is used for large $N$. The internal energy, heat capacity, and pressure derived below require only derivatives of $\ln Q$. The entropy is obtained from $\ln Q$ in [Section 4.3](../chapter-04/section-03.md) and evaluated for this gas in [Section 5.1](../chapter-05/section-01.md).

#### Internal Energy

Because $\Lambda = (\beta h^2/2\pi m)^{1/2}$ is proportional to $\beta^{1/2}$, $\ln Q$ depends on $\beta$ as

```{math}
\ln Q = -\frac{3}{2} N \ln \beta + (\text{terms that do not depend on } \beta).
```

The internal energy $U$ then follows from the derivative in [Section 2.3](section-03.md):

```{math}
:label: eq:internal-energy-cube
U 
= -\left(\frac{\partial \ln Q}{\partial \beta}\right)_{N,V}
= \frac{3}{2} \frac{N}{\beta}
= \frac{3}{2} N k_{\mathrm B} T.
```

```{admonition} Interpreting the result
:class: important
Compare $U/N$ with the average translational kinetic energy per particle obtained from kinetic theory in [Section 1.2](../chapter-01/section-02.md), and divide it among the three translational degrees of freedom. The volume does not appear in $U$. Use the expression for $\ln Q$ and the assumptions stated for the particles to explain why.
```

#### Heat Capacity at Constant Volume

From $U = \frac{3}{2}N k_{\mathrm B}T$, it follows that

```{math}
:label: eq:heat-capacity-cube
C_V 
= \left(\frac{\partial U}{\partial T}\right)_{N,V}
= \frac{3}{2} N k_{\mathrm B}.
```

```{important}
This temperature-independent $C_V$ is a hallmark of an ideal (monatomic) gas.
```

#### Pressure

To find the pressure, we use the expression from [Section 2.3](section-03.md),

```{math}
:label: eq:pressure-cube
P 
= k_{\mathrm B} T \left(\frac{\partial \ln Q}{\partial V}\right)_{N,T}
= \frac{N k_{\mathrm B} T}{V},
```

recovering the ideal gas law.

### Translational Partition Function

Because the particle-in-a-box spectrum describes translational motion, the single-particle partition function

```{math}
:label: eq:translational-partition-function
q_{\mathrm{trans}} 
= \frac{V}{\Lambda^3}
```

is called the translational partition function. It underpins the classical description of dilute gases when $\Lambda \ll (V/N)^{1/3}$, the mean spacing between particles, so quantum effects are negligible at ordinary densities and temperatures.

[](#fig:thermal-wavelength) draws the two lengths to scale for helium. In the gas at $300\ \mathrm{K}$ and $1\ \mathrm{bar}$, the conditions of the worked example below, the mean spacing is about 69 times $\Lambda$. In liquid helium at its normal boiling point of $4.22\ \mathrm{K}$, $\Lambda$ is larger than the mean spacing, which is the case that requires the full quantum treatment of identical particles.

```{code-cell} ipython3
:label: fig:thermal-wavelength
:caption: Thermal de Broglie wavelength $\Lambda$ and mean spacing $(V/N)^{1/3}$ of helium atoms, drawn to scale. The atoms are placed on a square grid at the mean spacing, and each disc has the diameter $\Lambda$. The two panels have different scales. (a) Gas at $300\ \mathrm{K}$ and $1\ \mathrm{bar}$, where the mean spacing is 69 times $\Lambda$. Each disc is a barely visible dot, and a gray ring marks its place. (b) Liquid at its normal boiling point of $4.22\ \mathrm{K}$, with a density of $125\ \mathrm{kg\,m^{-3}}$, where $\Lambda$ is 1.13 times the mean spacing and the discs overlap.
:tags: [hide-input]

from matplotlib.patches import Circle
from scipy.constants import h, k, u

m_He = 4.002602 * u              # mass of a helium-4 atom in kg


def thermal_wavelength(T):
    """Thermal de Broglie wavelength of a helium atom in m."""
    return h / np.sqrt(2 * np.pi * m_He * k * T)


T_gas, P_gas = 300, 1.00e5                         # K, Pa
spacing_gas = (k * T_gas / P_gas) ** (1 / 3)       # (V/N)^(1/3) of an ideal gas in m

T_liquid, density_liquid = 4.22, 125               # K, kg/m^3 (liquid helium at its normal boiling point)
spacing_liquid = (m_He / density_liquid) ** (1 / 3)


def draw_atoms(ax, c, spacing, wavelength):
    """Atoms on a grid at the mean spacing, each with a disc of diameter Lambda. Lengths are in units of the spacing."""
    radius = wavelength / (2 * spacing)
    for i in range(4):
        for j in range(3):
            if radius > 0.1:
                ax.add_patch(Circle((i, j), radius, fc=c.tint(c.blue, 0.72), ec=c.blue, lw=1.0, alpha=0.85))
            else:                                  # no outline, which would make the dot larger than Lambda
                ax.add_patch(Circle((i, j), radius, fc=c.blue, ec="none"))
                ax.add_patch(Circle((i, j), 0.11, fc="none", ec=c.muted, lw=0.8))     # a ring to find the dot
    # the mean spacing, between the centers of two neighbors
    for x in (0, 1):
        ax.plot([x, x], [-0.68, -0.92], color=c.muted, lw=0.9)
    ax.annotate("", xy=(0, -0.8), xytext=(1, -0.8),
                arrowprops=dict(arrowstyle="<|-|>", color=c.ink, lw=1.2, mutation_scale=9, shrinkA=0, shrinkB=0))
    ax.text(1.18, -0.8, f"$(V/N)^{{1/3}}$ = {1e9 * spacing:.3g} nm", ha="left", va="center")
    # the thermal wavelength: a dimension line where the disc is large enough, and a leader where it is a dot
    label = f"$\\Lambda$ = {1e9 * wavelength:.3f} nm"
    if radius > 0.1:
        for x in (-radius, radius):
            ax.plot([x, x], [2.68, 2.92], color=c.muted, lw=0.9)
        ax.annotate("", xy=(-radius, 2.8), xytext=(radius, 2.8),
                    arrowprops=dict(arrowstyle="<|-|>", color=c.ink, lw=1.2, mutation_scale=9, shrinkA=0, shrinkB=0))
        ax.text(radius + 0.18, 2.8, label, ha="left", va="center")
    else:
        ax.plot([0.015, 0.3], [2.03, 2.62], color=c.muted, lw=0.9)
        ax.text(0.22, 2.8, label, ha="left", va="center")
    ax.set(xlim=(-0.85, 3.85), ylim=(-1.05, 3.05), aspect="equal")
    ax.axis("off")


def draw(c):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(fs.WIDTH, 3.35), constrained_layout=True)
    draw_atoms(ax1, c, spacing_gas, thermal_wavelength(T_gas))
    draw_atoms(ax2, c, spacing_liquid, thermal_wavelength(T_liquid))
    fs.panel(ax1, "a", "Helium gas at 300 K and 1 bar")
    fs.panel(ax2, "b", "Liquid helium at 4.22 K")
    return fig


fs.show(draw, alt="Two panels, each with twelve helium atoms on a square grid at the mean spacing. In panel a, the gas "
                  "at 300 K and 1 bar, the mean spacing is 3.46 nanometers and the thermal wavelength is 0.050 "
                  "nanometers, so the disc around each atom is a barely visible dot inside a gray ring that marks its place. In panel b, the liquid at 4.22 K, "
                  "the mean spacing is 0.376 nanometers and the thermal wavelength is 0.425 nanometers, so the discs "
                  "of neighboring atoms overlap.")
```

### Computational Studio: Ideal Gas

In [](#fig:ideal-gas-studio), you can change the number of atoms $N$, the volume $V$, and the temperature $T$ of a small sample of a noble gas and follow the quantities that its partition function gives.

```{iframe} https://wexlergroup.github.io/chem-4020-5020/computational-studios/ideal-gas-studio.html
:label: fig:ideal-gas-studio
:width: 100%
:title: Ideal gas
:placeholder: ../_static/ideal-gas-studio/ideal-gas-studio.png

A small sample of an ideal monatomic gas. Panel (a) shows $N$ atoms in a cube of volume $V$, seen from one side. Panel (b) lists the thermal de Broglie wavelength $\Lambda$, the mean distance between atoms, the translational partition function $q_{\mathrm{trans}}$ (Eq. [](#eq:translational-partition-function)), $\ln Q$ (Eq. [](#eq:partition-function-many-particles-cube)), the internal energy (Eq. [](#eq:internal-energy-cube)), the pressure (Eq. [](#eq:pressure-cube)), and the entropy, which [Section 4.3](../chapter-04/section-03.md) obtains from $\ln Q$. Because $N$ is small, $\ln N!$ is evaluated exactly and not with the Stirling approximation.
```

You can also open the studio in a new tab: [Ideal Gas Studio](https://wexlergroup.github.io/chem-4020-5020/computational-studios/ideal-gas-studio.html).

## Worked Example

### Estimating $\Lambda$ and $q_{\mathrm{trans}}$ (helium, 300 K)

For a particle of mass $m$, the thermal de Broglie wavelength is

```{math}
\Lambda=\frac{h}{\sqrt{2\pi m k_{\mathrm B}T}}.
```

Take helium: $m=4.0026\,u=6.65\times10^{-27}\ \mathrm{kg}$, $T=300\ \mathrm{K}$.

1. **Compute $\Lambda$**

   ```{math}
   \Lambda=\frac{6.626\times10^{-34}}{\sqrt{2\pi(6.65\times10^{-27})(1.381\times10^{-23})(300)}}
   \approx 5.04\times10^{-11}\ \mathrm{m}.
   ```

2. **Compute the single-particle translational partition function in $V=1.00\ \mathrm{L}=10^{-3}\ \mathrm{m^3}$**

   ```{math}
   q_{\mathrm{trans}}=\frac{V}{\Lambda^3}
   =\frac{10^{-3}}{(5.04\times10^{-11})^3}
   \approx 7.8\times10^{27}.
   ```

3. **Compute $\alpha$ for the same volume, for which $V^{2/3}=10^{-2}\ \mathrm{m^2}$**

   ```{math}
   \alpha=\frac{h^2}{8mV^{2/3}k_{\mathrm B}T}
   =\frac{(6.626\times10^{-34})^2}{8(6.65\times10^{-27})(10^{-2})(1.381\times10^{-23})(300)}
   \approx 2.0\times10^{-19}.
   ```

4. **Compare the mean spacing between atoms with $\Lambda$ at $P=1.00\ \mathrm{bar}=10^{5}\ \mathrm{Pa}$**

   The ideal gas law gives the number of atoms in the volume,

   ```{math}
   N=\frac{PV}{k_{\mathrm B}T}
   =\frac{(10^{5})(10^{-3})}{(1.381\times10^{-23})(300)}
   \approx 2.41\times10^{22},
   ```

   so the mean spacing between atoms is

   ```{math}
   \left(\frac{V}{N}\right)^{1/3}
   =\left(\frac{10^{-3}}{2.41\times10^{22}}\right)^{1/3}
   \approx 3.46\times10^{-9}\ \mathrm{m}\approx 69\,\Lambda,
   ```

   and the number of thermally accessible translational states per atom is

   ```{math}
   \frac{q_{\mathrm{trans}}}{N}
   =\frac{7.8\times10^{27}}{2.41\times10^{22}}
   \approx 3.2\times10^{5}.
   ```

**Result.** In a macroscopic volume, $\alpha\ll 1$ and $q_{\mathrm{trans}}\gg 1$, consistent with the continuum approximation used to derive $V/\Lambda^3$. At 1 bar, $q_{\mathrm{trans}}$ is also much larger than $N$, which is the condition assumed for $Q=q^N/N!$ in [Section 2.4](section-04.md).

## Concept Checks

1. Why does the sum over quantum numbers become well-approximated by an integral at high $T$ or large $V$?
2. What does the condition "interparticle spacing $\gg \Lambda$" mean physically?
3. Why does translation give $U\propto T$ but not $U\propto V$ for an ideal gas?
4. Where does the $1/N!$ factor enter when connecting microscopic states to macroscopic entropy?

## Key Takeaways

- Quantized translation in a box leads, in the continuum limit, to $q_{\mathrm{trans}}=V/\Lambda^3$.
- $\Lambda$ sets the scale of quantum wavepacket "size"; classical behavior emerges when $\Lambda$ is small compared to typical spacings.
- For an ideal monatomic gas, $Q=(V/\Lambda^3)^N/N!$ reproduces $PV=Nk_{\mathrm B}T$ and $U=\tfrac{3}{2}Nk_{\mathrm B}T$.
- Partition functions connect microscopic spectra to thermodynamic equations of state via derivatives of $\ln Q$.
