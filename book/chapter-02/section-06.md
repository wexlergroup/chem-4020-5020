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
  enumerator: 2.6.%s
---


# 2.6. Harmonic Oscillator

See also: [Course-wide Conventions & Notation](../notation.md)

## Overview and Learning Objectives

The quantum harmonic oscillator is a canonical model for molecular vibrations and lattice modes in solids. In this section we derive the oscillator partition function, compute ensemble averages such as internal energy and heat capacity, and examine how the classical (equipartition) limit emerges when $T\gg \Theta_{\mathrm{vib}}$. We begin with a review of the one-dimensional oscillator, then connect the result to three-dimensional oscillators in the Einstein model.

Learning objectives:

* State the harmonic oscillator energy levels $E_n=\hbar\omega(n+\tfrac12)$ and write the corresponding partition function
* Evaluate the geometric series to obtain $q = e^{-\beta\hbar\omega/2}/(1-e^{-\beta\hbar\omega})$
* Derive the mean energy $U=\tfrac12\hbar\omega+\hbar\omega/(e^{\beta\hbar\omega}-1)$ and interpret the zero-point energy
* Use the vibrational temperature $\Theta_{\mathrm{vib}}=\hbar\omega/k_{\mathrm B}$ to analyze low- and high-temperature limits and connect to the Einstein model

## Core Ideas and Derivations

### Review of the One-Dimensional Harmonic Oscillator

[](#fig:oscillator-levels) shows the potential energy, the first four energy levels, and the corresponding wavefunctions of the one-dimensional harmonic oscillator that is reviewed in this subsection.

```{code-cell} ipython3
:label: fig:oscillator-levels
:caption: Potential energy $V(x)$ (thin curve), energy levels $E_n$ (gray lines), and wavefunctions $\psi_n$ (blue curves) of a one-dimensional harmonic oscillator. Each wavefunction is drawn on its energy level, and the levels are equally spaced by $\hbar\omega$.
:tags: [hide-input]

import sys
from math import factorial, pi

import matplotlib.pyplot as plt
import numpy as np
from scipy.special import eval_hermite

sys.path.append("..")            # the book folder, which holds figstyle
import figstyle as fs


def oscillator_wavefunction(n, x):
    """Wavefunction of level n, with x in units of (hbar / (m omega))^(1/2)."""
    return eval_hermite(n, x) * np.exp(-x**2 / 2) / np.sqrt(2**n * factorial(n) * np.sqrt(pi))


def draw(c):
    fig, ax = plt.subplots(figsize=(4.8, 3.7), constrained_layout=True)
    x = np.linspace(-5, 5, 600)              # displacement in units of (hbar / (m omega))^(1/2)
    levels = range(4)

    # Potential energy in units of hbar omega
    ax.plot(x, x**2 / 2, color=c.ink, lw=1.4)

    # Energy levels, drawn between the classical turning points, with each wavefunction drawn on its level
    for n in levels:
        energy = n + 0.5
        turning_point = np.sqrt(2 * energy)
        ax.plot([-turning_point, turning_point], [energy, energy], color=c.gray, lw=1.2)
        ax.plot(x, energy + 0.57 * oscillator_wavefunction(n, x), color=c.blue)
        ax.text(5.15, energy, rf"$\psi_{n}$", color=c.blue, va="center", clip_on=False)
    ax.text(-3.02, 4.05, "$V(x)$", ha="right", va="center")

    ax.set(xlabel="Displacement $x$", ylabel="Energy", xlim=(-5, 5), ylim=(0, 4.3), xticks=[0],
           yticks=[n + 0.5 for n in levels])
    ax.set_yticklabels([r"$E_0 = \frac{1}{2}\hbar\omega$", r"$E_1 = \frac{3}{2}\hbar\omega$",
                        r"$E_2 = \frac{5}{2}\hbar\omega$", r"$E_3 = \frac{7}{2}\hbar\omega$"])
    ax.tick_params(axis="y", length=0, pad=6)
    ax.spines["left"].set_visible(False)
    return fig


fs.show(draw, alt="Parabolic potential energy against displacement, with four equally spaced horizontal energy "
                  "levels at one half, three halves, five halves, and seven halves of h-bar omega. On each level, "
                  "the wavefunction is drawn. The lowest wavefunction is a single peak at zero displacement, and "
                  "each higher wavefunction has one more node and extends farther from the center.")
```

A one-dimensional harmonic oscillator is a particle of mass $m$ that moves in the potential

```{math}
V(x) = \frac{1}{2} k x^2,
```

where $x$ is the displacement from the equilibrium position and $k$ is the force constant. A larger force constant corresponds to a stiffer spring or, for a molecular vibration, a stiffer bond. In classical mechanics, the energy of the oscillator is the sum of its kinetic and potential energies,

```{math}
E = \frac{p_x^2}{2m} + \frac{1}{2} m \omega^2 x^2,
```

where $p_x$ is the momentum and

```{math}
\omega = \left( \frac{k}{m} \right)^{1/2}
```

is the angular frequency. For the vibration of a diatomic molecule, $x$ is the displacement of the bond length from its equilibrium value, and $m$ is replaced by the reduced mass $\mu$ of the two atoms ([Section 2.7](section-07.md)).

In quantum mechanics, the energy of the oscillator is quantized. For a one-dimensional harmonic oscillator with angular frequency $\omega$, the energy levels are

```{math}
E_n = \hbar \omega \left( n + \frac{1}{2} \right) \quad \text{for} \; n = 0, 1, 2, \ldots,
```

where $\hbar$ is the [reduced Planck constant](https://physics.nist.gov/cgi-bin/cuu/Value?hbar) and $n$ is a nonnegative integer. The lowest level, $n = 0$, has energy $E_0 = \hbar\omega/2$, which is called the **zero-point energy**.

### Partition Function for a Harmonic Oscillator

The (single-oscillator) partition function for a one-dimensional harmonic oscillator is

```{math}
\begin{aligned}
q &= \sum_{n=0}^{\infty} e^{-\beta E_n} \\
&= \sum_{n=0}^{\infty} e^{-\beta \hbar \omega (n + 1/2)} \\
&= e^{-\beta \hbar \omega / 2} \sum_{n=0}^{\infty} e^{-\beta \hbar \omega n}.
\end{aligned}
```

#### Evaluation of the Sum

The term $\sum_{n=0}^{\infty} e^{-\beta \hbar \omega n}$ is a geometric series. For $|x| < 1$,

```{math}
\sum_{n=0}^{\infty} x^n = \frac{1}{1 - x}.
```

Substituting $x = e^{-\beta \hbar \omega}$ gives

```{math}
q = \frac{e^{-\beta \hbar \omega / 2}}{1 - e^{-\beta \hbar \omega}}.
```

````{admonition} Check that $|x| < 1$
:class: dropdown
```{math}
\begin{aligned}
|x| &< 1 \\
\left| e^{-\beta \hbar \omega} \right| &< 1 \\
e^{-\beta \hbar \omega} &< 1 \\
\beta \hbar \omega &> 0 \\
\beta &> 0.
\end{aligned}
```
````

### Ensemble Averages for a Harmonic Oscillator

#### Natural Logarithm of the Partition Function

Taking the natural logarithm of $q$ gives

```{math}
:label: harmonic_oscillator_lnq
\ln q = -\frac{\beta \hbar \omega}{2} - \ln(1 - e^{-\beta \hbar \omega}).
```

#### Internal Energy

For a single oscillator, the mean energy (internal energy) follows from

```{math}
U = -\frac{\partial \ln q}{\partial \beta}
= \frac{\hbar \omega}{2} + \frac{\hbar \omega}{e^{\beta \hbar \omega} - 1}.
```

````{admonition} Derivative of $\ln(1 - e^{-\beta \hbar \omega})$
:class: dropdown
```{math}
\begin{aligned}
\frac{\partial}{\partial \beta} \ln(1 - e^{-\beta \hbar \omega})
&= \frac{1}{1 - e^{-\beta \hbar \omega}}
\frac{\partial}{\partial \beta}\left(1 - e^{-\beta \hbar \omega}\right) \\
&= \frac{\hbar \omega\, e^{-\beta \hbar \omega}}{1 - e^{-\beta \hbar \omega}}
= \frac{\hbar \omega}{e^{\beta \hbar \omega} - 1}.
\end{aligned}
```
````

It is often convenient to define the vibrational temperature $\Theta_{\mathrm{vib}}$ as

```{math}
\Theta_{\mathrm{vib}} = \frac{\hbar \omega}{k_{\mathrm B}}.
```

Using $\beta \hbar \omega = \Theta_{\mathrm{vib}}/T$, the internal energy can be written as

```{math}
U = \frac{\hbar \omega}{2} + \frac{\hbar \omega}{e^{\Theta_{\mathrm{vib}} / T} - 1}
= k_{\mathrm B}\Theta_{\mathrm{vib}}\left[\frac{1}{2}+\frac{1}{e^{\Theta_{\mathrm{vib}}/T}-1}\right].
```

```{code-cell} ipython3
:label: fig:oscillator-energy
:caption: Harmonic oscillator as a function of the reduced temperature $T/\Theta_{\mathrm{vib}}$. (a) Internal energy in units of $\hbar\omega$. The dashed lines are the zero-point energy, which $U$ approaches at low temperature, and $k_{\mathrm B}T$, which it approaches at high temperature. (b) Heat capacity at constant volume in units of $k_{\mathrm B}$, which is derived below. The dashed line is its high-temperature limit.
:tags: [hide-input]

Theta_vib = 805                              # K, for a Cl2 molecule
T = np.linspace(10, 1000, 400)               # K
Tr = T / Theta_vib                           # reduced temperature

U = 0.5 + 1 / (np.exp(1 / Tr) - 1)           # internal energy in units of hbar omega
C_V = Tr**-2 * np.exp(1 / Tr) / (np.exp(1 / Tr) - 1)**2    # heat capacity in units of k_B


def draw(c):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(fs.WIDTH, 3.0), constrained_layout=True)
    dashed = dict(color=c.muted, lw=1.2, ls=(0, (4, 3)))

    ax1.axhline(0.5, **dashed)
    ax1.plot([0, 1.25], [0, 1.25], **dashed)
    ax1.plot(Tr, U, color=c.ink)
    ax1.text(1.25, 0.44, r"zero-point energy, $\hbar\omega/2$", color=c.muted, ha="right", va="top")
    ax1.text(1.02, 0.9, r"$k_{\mathrm{B}}T$", color=c.muted, ha="left", va="top")
    ax1.set(xlabel=r"$T/\Theta_{\mathrm{vib}}$", ylabel=r"$U/\hbar\omega$", xlim=(0, 1.25), ylim=(0, 1.4))
    fs.panel(ax1, "a", "Internal energy")

    ax2.axhline(1, **dashed)
    ax2.plot(Tr, C_V, color=c.ink)
    ax2.text(1.25, 1.03, "high-temperature limit", color=c.muted, ha="right", va="bottom")
    ax2.set(xlabel=r"$T/\Theta_{\mathrm{vib}}$", ylabel=r"$C_V/k_{\mathrm{B}}$", xlim=(0, 1.25), ylim=(0, 1.2))
    fs.panel(ax2, "b", "Heat capacity")
    return fig


fs.show(draw, alt="Two panels against the reduced temperature from 0 to 1.25 for a harmonic oscillator. Panel a: "
                  "the internal energy in units of h-bar omega starts at the zero-point energy of 0.5, stays there up "
                  "to a reduced temperature of about 0.2, and then rises toward the line U equal to k_B T, reaching "
                  "about 1.3 at 1.25. Panel b: the heat capacity in units of k_B is zero up to a reduced "
                  "temperature of about 0.1, rises steeply, and levels off toward a dashed line at 1, reaching about "
                  "0.95 at 1.25.")
```

The dashed lines in panel (a) of [](#fig:oscillator-energy) are the low- and high-temperature limits of $U$. For $T \ll \Theta_{\mathrm{vib}}$, the exponential $e^{\Theta_{\mathrm{vib}}/T}$ is very large, so the second term in the brackets vanishes and

```{math}
U \to \frac{1}{2} k_{\mathrm B}\Theta_{\mathrm{vib}} = \frac{\hbar \omega}{2},
```

which is the zero-point energy. For $T \gg \Theta_{\mathrm{vib}}$, the ratio $x = \Theta_{\mathrm{vib}}/T$ is small, and the expansion $e^x \approx 1 + x + x^2/2$ gives

```{math}
\frac{1}{e^x - 1} \approx \frac{1}{x \left( 1 + x/2 \right)} \approx \frac{1}{x} - \frac{1}{2}.
```

The constant $\tfrac{1}{2}$ in the brackets cancels, and

```{math}
U \to k_{\mathrm B}\Theta_{\mathrm{vib}} \, \frac{T}{\Theta_{\mathrm{vib}}} = k_{\mathrm B} T.
```

This limit can be compared with the equipartition theorem of classical statistical mechanics, which states that each quadratic term in the energy contributes $\tfrac{1}{2}k_{\mathrm B}T$ to the mean energy. For example, the three quadratic terms of the translational kinetic energy, $(p_x^2 + p_y^2 + p_z^2)/2m$, give $\tfrac{3}{2}k_{\mathrm B}T$ ([Section 1.2](../chapter-01/section-02.md)).

The value of $\Theta_{\mathrm{vib}}$ depends on the molecule. {numref}`vibrational-temperatures` lists $\Theta_{\mathrm{vib}}$ for several diatomic molecules, calculated from the harmonic vibrational wavenumber $\tilde{\nu}$ as $\Theta_{\mathrm{vib}} = hc\tilde{\nu}/k_{\mathrm B}$, where $c$ is the speed of light and $hc/k_{\mathrm B} = 1.43878\ \mathrm{cm\,K}$. The plots in this section use the value for $\mathrm{Cl_2}$.

```{list-table} Harmonic vibrational wavenumbers (from the [NIST Chemistry WebBook](https://webbook.nist.gov/chemistry/)) and vibrational temperatures of selected diatomic molecules
:header-rows: 1
:name: vibrational-temperatures

* - Molecule
  - $\tilde{\nu}$ ($\mathrm{cm^{-1}}$)
  - $\Theta_{\mathrm{vib}}$ (K)
* - $\mathrm{H_2}$
  - 4401.21
  - 6332
* - $\mathrm{HCl}$
  - 2990.95
  - 4303
* - $\mathrm{N_2}$
  - 2358.57
  - 3393
* - $\mathrm{CO}$
  - 2169.81
  - 3122
* - $\mathrm{O_2}$
  - 1580.19
  - 2274
* - $\mathrm{Cl_2}$
  - 559.72
  - 805
* - $\mathrm{Br_2}$
  - 325.32
  - 468
* - $\mathrm{I_2}$
  - 214.50
  - 309
```

#### Heat Capacity at Constant Volume

In [Section 2.3](section-03.md), we showed that

```{math}
\sigma_E^2 = \left( \frac{\partial^2 \ln Q}{\partial \beta^2} \right)_{N,V} = k_{\mathrm B} T^2 C_V.
```

For a single oscillator, $Q = q$, and from Equation {eq}`harmonic_oscillator_lnq`, it follows that

```{math}
C_V
= k_{\mathrm B} \left( \frac{\hbar \omega}{k_{\mathrm B} T} \right)^2
\frac{e^{\hbar \omega / k_{\mathrm B} T}}{\left( e^{\hbar \omega / k_{\mathrm B} T} - 1 \right)^2}
= k_{\mathrm B} \left( \frac{\Theta_{\mathrm{vib}}}{T} \right)^2
\frac{e^{\Theta_{\mathrm{vib}} / T}}{\left( e^{\Theta_{\mathrm{vib}} / T} - 1 \right)^2}.
```

````{admonition} Complete Derivation of $C_V$
:class: dropdown
```{math}
\begin{aligned}
\sigma_E^2
&= \frac{\partial^2}{\partial \beta^2}\left( -\frac{\beta \hbar \omega}{2} - \ln(1 - e^{-\beta \hbar \omega}) \right) \\
&= \frac{\partial}{\partial \beta}\left( -\frac{\hbar \omega}{2} - \frac{\hbar \omega}{e^{\beta \hbar \omega} - 1} \right) \\
&= -\hbar \omega \left( \frac{\partial g^{-1}}{\partial g} \frac{\partial g}{\partial \beta} \right)
\quad \text{where} \; g = e^{\beta \hbar \omega} - 1 \\
&= \frac{\hbar \omega}{g^2} \left( \frac{\partial}{\partial \beta} \left( e^{\beta \hbar \omega} - 1 \right) \right) \\
&= \left( \hbar \omega \right)^2 \frac{e^{\beta \hbar \omega}}{\left( e^{\beta \hbar \omega} - 1 \right)^2}.
\end{aligned}
```

Dividing by $k_{\mathrm B} T^2$ and using $\beta = 1/(k_{\mathrm B} T)$ gives the expression for $C_V$ above.
````

Panel (b) of [](#fig:oscillator-energy) shows $C_V$ as a function of the reduced temperature.

```{admonition} Limits of the heat capacity
:class: important
Take the limits $T \ll \Theta_{\mathrm{vib}}$ and $T \gg \Theta_{\mathrm{vib}}$ of the expression for $C_V$, following the steps used for $U$ above, and compare the results with panel (b) of [](#fig:oscillator-energy).
```

### Computational Studio: Harmonic Oscillator

In [](#fig:harmonic-oscillator-studio), you can set the vibrational temperature $\Theta_{\mathrm{vib}}$, for example to a value of {numref}`vibrational-temperatures`, and follow the probabilities of the levels, the internal energy, and the heat capacity at any temperature.

```{iframe} https://wexlergroup.github.io/chem-4020-5020/computational-studios/harmonic-oscillator-studio.html
:label: fig:harmonic-oscillator-studio
:width: 100%
:title: Harmonic oscillator
:placeholder: ../_static/harmonic-oscillator-studio/harmonic-oscillator-studio.png

One harmonic oscillator with the vibrational temperature $\Theta_{\mathrm{vib}}$. Panel (a) shows the probabilities $p_n$ of the levels at the temperature $T$. Panels (b) and (c) show the internal energy in units of $\hbar\omega$ and the heat capacity in units of $k_{\mathrm B}$ as functions of temperature, with the dashed lines of [](#fig:oscillator-energy). The dot marks the temperature $T$ of panel (a). The buttons set $\Theta_{\mathrm{vib}}$ to the values of {numref}`vibrational-temperatures`.
```

You can also open the studio in a new tab: [Harmonic Oscillator Studio](https://wexlergroup.github.io/chem-4020-5020/computational-studios/harmonic-oscillator-studio.html).

### Einstein Model

The Einstein model is a simple model for an atomic crystal, where $N$ identical atoms occupy lattice sites and vibrate about these sites as independent (i.e., non-interacting) three-dimensional harmonic oscillators, all with the same frequency. Because the atoms occupy fixed lattice sites, they can be treated as distinguishable. Each atom vibrates independently along $x$, $y$, and $z$, so its energy is the sum of the energies of three one-dimensional oscillators, and the crystal is a collection of $3N$ independent, distinguishable one-dimensional oscillators. The partition function is therefore a product of $3N$ factors of $q$ ([Section 2.4](section-04.md)):

```{math}
Q = q^{3N} = \left( \frac{e^{-\beta \hbar \omega / 2}}{1 - e^{-\beta \hbar \omega}} \right)^{3N}.
```

Because $\ln Q = 3N \ln q$, the internal energy and the heat capacity are $3N$ times the results for a single oscillator:

```{math}
U = \frac{3}{2} N \hbar \omega + \frac{3N \hbar \omega}{e^{\beta \hbar \omega} - 1},
```

and

```{math}
C_V = 3 N k_{\mathrm B} \left( \frac{\hbar \omega}{k_{\mathrm B} T} \right)^2 \frac{e^{\hbar \omega / k_{\mathrm B} T}}{\left( e^{\hbar \omega / k_{\mathrm B} T} - 1 \right)^2}.
```

#### Dulong–Petit Law

The **Dulong–Petit law** is the empirical observation that the molar heat capacity of many elemental solids at room temperature is close to $3R = 24.9\ \mathrm{J\,mol^{-1}\,K^{-1}}$, where $R$ is the gas constant. [](#fig:dulong-petit) compares the law with molar heat capacities measured at 25 °C. Elements whose heat capacity differs from $3R$ by more than $R$ are labeled.

```{code-cell} ipython3
:label: fig:dulong-petit
:caption: Molar heat capacity at 25 °C of elemental solids (and liquid mercury) as a function of atomic number, in units of the gas constant $R$. The dashed line is the Dulong–Petit law.
:tags: [hide-input]

# https://en.wikipedia.org/wiki/Heat_capacities_of_the_elements_(data_page)

import pandas as pd

df = pd.read_csv("../_static/chapter-02/section-06/dulong_petit.csv")
df = df[df.Source == "use"].copy()
df["Atomic Number"] = df["Atomic Number"].astype(int)
df["Molar (J/mol·K)"] = df["Molar (J/mol·K)"].astype(float)
df["Molar (R)"] = df["Molar (J/mol·K)"] / 8.314
df["Absolute Deviation"] = np.abs(df["Molar (R)"] - 3)

# Remove elemental gases
elemental_gases = ["H", "He", "N", "O", "F", "Ne", "Cl", "Ar", "Kr", "Xe", "Rn"]
df = df[~df.Symbol.isin(elemental_gases)]


def draw(c):
    fig, ax = plt.subplots(figsize=(6.2, 3.6), constrained_layout=True)

    # Dulong-Petit law
    ax.axhline(3, color=c.muted, lw=1.2, ls=(0, (4, 3)))
    ax.text(101.5, 3, "Dulong–Petit\nlaw, $3R$", color=c.muted, va="center", clip_on=False)

    ax.scatter(df["Atomic Number"], df["Molar (R)"], s=20, color=c.blue, edgecolor=c.ground, linewidth=0.5,
               zorder=3)

    # Label the outliers (the two entries for carbon are graphite and diamond)
    outliers = df[df["Absolute Deviation"] > 1]
    for _, row in outliers.iterrows():
        label, offset, ha, va = row["Symbol"], (6, 0), "left", "center"
        if row["Symbol"] == "C":
            allotrope = "graphite" if "graphite" in str(row["Notes"]) else "diamond"
            label = f"C ({allotrope})"
            offset = (6, 3) if allotrope == "graphite" else (6, -4)
        elif row["Symbol"] == "B" or row["Atomic Number"] > 90:
            offset, ha = (-6, 0), "right"
        ax.annotate(label, (row["Atomic Number"], row["Molar (R)"]), xytext=offset, textcoords="offset points",
                    ha=ha, va=va)

    ax.set(xlabel="Atomic number", ylabel="Molar heat capacity, $C/R$", xlim=(0, 100), ylim=(0, 8))
    return fig


fs.show(draw, alt="Molar heat capacity divided by the gas constant against atomic number for 82 elements. Most "
                  "points lie close to the dashed line at 3, which is the Dulong–Petit law. Beryllium, boron, "
                  "graphite, and diamond lie well below the line, between 0.7 and 2. Gadolinium and plutonium lie "
                  "above it, near 4.4, and americium lies near 7.5.")
```

```{admonition} Comparing the data with the Einstein model
:class: important
For one mole of atoms, $N k_{\mathrm B} = R$. Compare the Dulong–Petit law with the heat capacity of the Einstein model. In which temperature range does the model reproduce the law? What does the model suggest about $\Theta_{\mathrm{vib}}$ for the labeled elements that lie below the line (Be, B, and C)?
```

## Worked Example

### Mean vibrational energy at finite temperature

For a harmonic oscillator,

```{math}
U = k_{\mathrm B}\Theta_{\mathrm{vib}}\left[\frac{1}{2}+\frac{1}{e^{\Theta_{\mathrm{vib}}/T}-1}\right].
```

Take $\Theta_{\mathrm{vib}}=805\ \mathrm{K}$ (the value for $\mathrm{Cl_2}$ in {numref}`vibrational-temperatures`) and $T=298\ \mathrm{K}$.

1. **Compute the exponential factor**

   ```{math}
   \frac{\Theta_{\mathrm{vib}}}{T}=\frac{805}{298}=2.70,
   \qquad
   e^{\Theta_{\mathrm{vib}}/T}=e^{2.70}=14.9.
   ```

2. **Evaluate the bracket**

   ```{math}
   \frac{1}{2}+\frac{1}{e^{\Theta_{\mathrm{vib}}/T}-1}
   =\frac{1}{2}+\frac{1}{13.9}=0.5719.
   ```

3. **Compute $U$**

   ```{math}
   U = k_{\mathrm B}(805\ \mathrm{K})(0.5719)=k_{\mathrm B}(460.4\ \mathrm{K}).
   ```

   Per mole, $U_m = R(460.4\ \mathrm{K})=3.83\ \mathrm{kJ\,mol^{-1}}$.

4. **Compute the populations of the two lowest levels**

   The probability of level $n$ is $p_n = e^{-\beta E_n}/q$ ([Section 2.2](section-02.md)). With the partition function derived above, $p_n = e^{-n\Theta_{\mathrm{vib}}/T}\left(1-e^{-\Theta_{\mathrm{vib}}/T}\right)$, so

   ```{math}
   p_0 = 1-\frac{1}{14.9}=0.933,
   \qquad
   p_1 = \frac{1}{14.9}(0.933)=0.063.
   ```

**Result.** At 298 K, $U$ exceeds the zero-point energy, $k_{\mathrm B}(402.5\ \mathrm{K})$, by 14%, and 93% of the oscillators are in the lowest level.

## Concept Checks

1. Why does the harmonic oscillator retain nonzero energy as $T\to 0$?
2. What does the limit $T\gg\Theta_{\mathrm{vib}}$ predict for the heat capacity, and how does that compare to equipartition?
3. How does increasing $\omega$ (stiffer bond) change $\Theta_{\mathrm{vib}}$ and the thermal population of excited levels?
4. Why is vibration often “frozen out” at room temperature for high-frequency modes?

## Key Takeaways

* The harmonic oscillator has evenly spaced energy levels and a closed-form partition function.
* Zero-point energy $\tfrac12\hbar\omega$ persists at $T=0$.
* The vibrational temperature $\Theta_{\mathrm{vib}}$ sets the temperature scale for activating vibrational degrees of freedom.
* In the high-$T$ limit, quantum results approach classical equipartition behavior.
