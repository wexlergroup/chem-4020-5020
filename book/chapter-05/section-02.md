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
  enumerator: 5.2.%s
---


# 5.2. Third Law

See also: [Course-wide Conventions & Notation](../notation.md)

Run the code: [Notebook for Section 5.2](https://wexlergroup.github.io/chem-4020-5020/lite/notebooks/index.html?path=chapter-05/section-02.ipynb), in which you can change and rerun the Python cells of this section in your browser.

## Overview and Learning Objectives

The first and second laws tell us *how* energy and entropy change but leave the entropy scale with an undetermined reference point — the integration constant in $\Delta S = \int \delta q_{\mathrm{rev}}/T$ is never fixed. The third law resolves this by specifying what happens to entropy as $T \to 0$. Using Boltzmann's formula ([Section 4.3](../chapter-04/section-03.md)), we explain why an ideal crystal has $S(0) = 0$, examine how ground-state degeneracy or frozen-in disorder produces a nonzero *residual entropy*, and note two further consequences: heat capacities must vanish as $T \to 0$, and absolute zero itself is unattainable in a finite number of steps.

With the entropy reference point fixed, we can compute **absolute** (third-law) entropies $S^\circ(T)$ for any substance by integrating $C_P/T$ from $0$ to $T$ and adding any phase-transition contributions. These absolute entropies are the basis for tabulated $\Delta_r S^\circ$ and $\Delta_r G^\circ$ values that we will put to use in [Section 5.3](section-03.md).

---

Learning objectives:

- State Planck's formulation of the third law and define what is meant by a pure crystalline substance and an ideal crystal.
- Use $S = k_{\mathrm B}\ln\Omega$ to justify $S(0) = 0$ for an ideal crystal with a unique ground state.
- Compute residual entropy $S_{\mathrm{res}} = k_{\mathrm B}\ln\Omega_0$ (or per mole, $R\ln\Omega_0$) for degenerate or disordered ground states.
- Explain why $C_V$ and $C_P$ must vanish as $T \to 0$ and connect this to the quantum behavior of heat capacities ([Section 2.6](../chapter-02/section-06.md)).
- Describe, qualitatively, why absolute zero is unattainable.

## Core Ideas and Derivations

### Third Law of Thermodynamics

#### Planck's Statement

**Planck's statement of the third law:**

> As $T \to 0\,\mathrm{K}$, the entropy of any **pure crystalline substance** tends to a **constant**.

A key special case is the **ideal (perfect) crystal** — one with no defects (no vacancies, impurities, grain boundaries, or stacking faults). For an ideal crystal, the constant is taken to be zero:

```{math}
:label: eq:third-law
S(0\,\mathrm{K}) = 0 \quad \text{(ideal crystal)}.
```

#### Why the Ideal Crystal Has $S(0) = 0$

Recall from [Section 4.3](../chapter-04/section-03.md) that Boltzmann's formula gives the entropy of an isolated system as

```{math}
S = k_{\mathrm B}\ln\Omega,
```

where $\Omega$ is the number of microstates compatible with the macroscopic constraints. An ideal crystal at $0\,\mathrm{K}$ has a unique ground-state arrangement: every atom sits at its designated lattice site in a single, perfectly ordered configuration. That means $\Omega_0 = 1$, so

```{math}
S(0) = k_{\mathrm B}\ln 1 = 0.
```

The third law, from the statistical-mechanical viewpoint, is the statement that a system with a unique ground state has zero entropy at absolute zero.

#### Real Crystals: Defects and Disorder

Real materials can deviate from this picture. Even at very low temperature, a crystal may contain — or "freeze in" — structural or chemical disorder: vacancies, impurities, grain boundaries, stacking faults, or orientational disorder of molecular units. These imperfections increase the number of distinct microscopic arrangements compatible with the same macroscopic state (i.e., increase $\Omega_0$), leading to a nonzero entropy as $T \to 0$.

---

### Residual Entropy

If a system retains **more than one** accessible microscopic arrangement as $T \to 0$, its entropy approaches a nonzero constant called the **residual entropy**:

```{math}
:label: eq:residual-entropy
S_{\mathrm{res}} = \lim_{T\to 0}\,S(T) = k_{\mathrm B}\ln\Omega_0.
```

#### Example: Orientational Disorder in Solid CO

Solid CO at very low temperature provides a classic example. The CO molecule has a small dipole moment, and the two orientations CO and OC are nearly isoenergetic. If each of $N$ molecules can be frozen in either orientation with roughly equal probability, the ground-state multiplicity is

```{math}
\Omega_0 = 2^N,
```

and the residual entropy is

```{math}
S_{\mathrm{res}} = k_{\mathrm B}\ln(2^N) = Nk_{\mathrm B}\ln 2.
```

Per mole, using $R = N_Ak_{\mathrm B}$:

```{math}
S_{\mathrm{res,m}} = R\ln 2 \approx 5.76\ \mathrm{J\,mol^{-1}\,K^{-1}}.
```

[](#fig:residual-entropy-CO) contrasts the ordered crystal with the disordered one for a small part of a lattice.

% Try this: Change the seed `4020` to another whole number. Each seed draws one of the arrangements that $\Omega_0$ counts. How many different pictures of panel (b) exist for the 36 molecules that are drawn?

```{code-cell} ipython3
:label: fig:residual-entropy-CO
:caption: Schematic of solid CO, drawn in two dimensions. (a) Every molecule has the orientation that the ordered structure assigns to its site, which is a single arrangement. (b) Each molecule is frozen at random in one of its two head-to-tail orientations. The positions of the molecules are the same as in (a), and $N$ molecules have $2^N$ arrangements of this kind.
:tags: [hide-input]

import sys
sys.path.append("..")                              # the book folder, which holds figstyle
import numpy as np
import matplotlib.pyplot as plt
import figstyle as fs

rng = np.random.default_rng(4020)
n = 6                                              # molecules per row
flips = {"ordered": np.ones((n, n)), "disordered": rng.choice([-1, 1], size=(n, n))}
titles = {"ordered": r"Ordered: $\Omega_0 = 1$", "disordered": r"Disordered: $\Omega_0 = 2^N$"}

def draw(c):
    fig, axs = plt.subplots(1, 2, figsize=(fs.WIDTH, 3.6), constrained_layout=True)
    for ax, key, letter in zip(axs, flips, "ab"):
        for i in range(n):
            for j in range(n):
                tilt = np.deg2rad(35 if (i + j) % 2 == 0 else -35)       # herringbone packing
                d = 0.27 * flips[key][i, j] * np.array([np.cos(tilt), np.sin(tilt)])
                ax.plot([j - d[0], j + d[0]], [i - d[1], i + d[1]], "-", color=c.muted, lw=2, zorder=1)
                ax.scatter(j - d[0], i - d[1], s=95, color=c.carbon, zorder=2, lw=0)
                ax.scatter(j + d[0], i + d[1], s=120, color=c.oxygen, zorder=2, lw=0)
        fs.panel(ax, letter, titles[key])
        ax.set_xlim(-0.6, n - 0.4)
        ax.set_ylim(-0.6, n - 0.4)
        ax.set_aspect("equal")
        ax.axis("off")
    axs[0].scatter([], [], s=95, color=c.carbon, label="C")
    axs[0].scatter([], [], s=120, color=c.oxygen, label="O")
    fig.legend(loc="outside lower center", ncols=2)
    return fig

fs.show(draw, alt="Two square lattices of 36 carbon monoxide molecules, each drawn as a small carbon atom bonded to a "
                  "larger oxygen atom. In panel a, every molecule points the way its site prescribes. In panel b, "
                  "about half of the molecules are flipped end for end at random.")
```

[](#fig:solid-CO-interactive) shows one unit cell of the crystal in three dimensions. The molecules sit near the sites of a face-centered cubic lattice, and each one lies along one of the four body diagonals of the cubic cell. Turn the cell, and switch between the ordered structure and a crystal in which each molecule is frozen at random in one of its two orientations.

```{iframe} https://wexlergroup.github.io/chem-4020-5020/computational-studios/solid-carbon-monoxide.html
:label: fig:solid-CO-interactive
:width: 100%
:title: Solid carbon monoxide
:placeholder: ../_static/solid-carbon-monoxide/solid-carbon-monoxide.png

One unit cell of solid CO in its low-temperature phase, with the 14 molecules at the corners and on the faces of the cell. The switch changes between the ordered structure and a crystal in which each molecule is frozen at random in one of its two head-to-tail orientations. A blue outline marks a molecule that is flipped end for end. The atom positions are those of the X-ray structure of @vegard1930.
```

Each press of "Another random arrangement" draws one of the $2^{14} = 16{,}384$ arrangements of the 14 molecules shown. The same count for the $N$ molecules of a crystal gives $\Omega_0 = 2^N$.

```{admonition} Comparison with Experiment
:class: dropdown

The experimentally measured residual entropy of CO is approximately $4.6\ \mathrm{J\,mol^{-1}\,K^{-1}}$, somewhat less than the $R\ln 2 = 5.76\ \mathrm{J\,mol^{-1}\,K^{-1}}$ predicted by the fully random two-orientation model. The discrepancy indicates that CO is not *completely* disordered at low temperature — neighboring molecules have weak orientational correlations that partially reduce the number of accessible configurations. The simple $2^N$ estimate captures the right order of magnitude but overcounts somewhat.
```

```{admonition} Residual Entropy and the Third Law
:class: important

Return to Planck's statement at the start of this section and to the definition of the residual entropy, Eq. {eq}`eq:residual-entropy`. Which values of $\Omega_0$ does the statement allow? Is solid CO a counterexample to Planck's statement, to Eq. {eq}`eq:third-law`, or to neither?
```

The third law's practical utility lies in establishing a universal reference point for entropy: any ideal crystal has $S(0) = 0$.

---

### Consequence: Heat Capacities Vanish as $T \to 0$

The third law has a direct implication for heat capacities. At constant pressure, the entropy at temperature $T$ is related to the entropy at $0\,\mathrm{K}$ by

```{math}
S(T) = S(0) + \int_0^T \frac{C_P(T')}{T'}\,dT'.
```

For this integral to converge (i.e., for $S(T)$ to be finite), $C_P$ must vanish as $T \to 0$. If $C_P$ approached a nonzero constant, the integrand $C_P/T'$ would grow as $1/T'$, and the integral would diverge logarithmically at its lower limit. Therefore:

```{math}
C_P \to 0 \quad \text{as} \quad T \to 0.
```

The same argument applies to $C_V$. This is a nontrivial prediction: it says that the ability of any substance to absorb heat must diminish as absolute zero is approached.

````{admonition} Connection to Module 2
:class: important

[Section 2.6](../chapter-02/section-06.md) derived the heat capacity of the Einstein model of an atomic crystal, which in terms of $\Theta_{\mathrm{vib}} = \hbar\omega/k_{\mathrm B}$ reads

```{math}
C_V = 3Nk_{\mathrm B}\left(\frac{\Theta_{\mathrm{vib}}}{T}\right)^2 \frac{e^{\Theta_{\mathrm{vib}}/T}}{(e^{\Theta_{\mathrm{vib}}/T}-1)^2}.
```

Find its limit as $T \to 0$, and compare it with the classical (equipartition) value $C_V = 3Nk_{\mathrm B}$. Which of the two is compatible with the third law?
````

---

### Consequence: Unattainability of Absolute Zero

A closely related statement — sometimes called the **Nernst form** of the third law — is that no finite sequence of thermodynamic operations can bring a system to exactly $T = 0\,\mathrm{K}$.

The intuitive argument is as follows. Cooling to very low temperature typically alternates two steps, as in adiabatic demagnetization: an isothermal step that removes entropy from the system, and an adiabatic step that lowers the temperature at constant entropy. As $T \to 0$, the entropy of the system approaches the same limiting value whatever the magnetic field (or other control variable), so each successive pair of steps removes less entropy and lowers the temperature by less. The sequence converges but never quite reaches zero. In practice, temperatures below $1\ \mathrm{nK}$ have been achieved in laboratory settings, but $T = 0$ remains a limit, not an endpoint.

[](#fig:adiabatic-demagnetization) shows the construction for a model paramagnet of $N$ independent magnetic moments. Each moment is a two-state system ([Section 2.3](../chapter-02/section-03.md)) whose energy gap is $\varepsilon_0$ when the field is off, because of the interactions between the moments, and $2\varepsilon_0$ when the field is on. The entropy of a two-state system depends on the gap $\varepsilon$ and on the temperature only through the ratio $\varepsilon/(k_{\mathrm B}T)$. An adiabatic step that returns the gap from $2\varepsilon_0$ to $\varepsilon_0$ therefore halves the temperature, and after $n$ cycles that start at $T_0$ the temperature is $T_0/2^n$, which is greater than zero for every finite $n$.

% Try this: Change `ratio` from `2.0` to `3.0`, a stronger field. By what factor does each cycle lower the temperature now?

```{code-cell} ipython3
:label: fig:adiabatic-demagnetization
:caption: Entropy of a model paramagnet with the magnetic field off (blue) and on (orange), and successive cooling cycles between the two curves. Each vertical arrow is an isothermal magnetization, in which entropy leaves the system as heat. Each horizontal arrow is an adiabatic demagnetization at constant entropy. Both curves approach the same entropy as $T \to 0$, so each cycle removes less entropy and lowers the temperature by less than the one before.
:tags: [hide-input]

import sys
sys.path.append("..")                              # the book folder, which holds figstyle
import numpy as np
import matplotlib.pyplot as plt
import figstyle as fs

def S_two_state(x):
    """Entropy of a two-state system in units of k_B, with x = epsilon/(k_B T)."""
    return np.log1p(np.exp(-x)) + x / (np.exp(x) + 1)

# Energy gap of each magnetic moment: eps_0 with the field off, 2 eps_0 with the field on
ratio = 2.0
t = np.linspace(0.01, 1.2, 500)                    # temperature in units of eps_0/k_B

def draw(c):
    fig, ax = plt.subplots(figsize=(6, 4), constrained_layout=True)
    ax.plot(t, S_two_state(1 / t), "-", color=c.blue)
    ax.plot(t, S_two_state(ratio / t), "-", color=c.vermilion)
    ax.text(0.47, 0.42, "field off", color=c.blue, ha="right", fontsize=fs.BASE + 0.5)
    ax.text(0.97, 0.28, "field on", color=c.vermilion, fontsize=fs.BASE + 0.5)

    # Alternate isothermal magnetization (down) and adiabatic demagnetization (left)
    t_n = 0.8
    for n in range(4):
        S_off, S_on = S_two_state(1 / t_n), S_two_state(ratio / t_n)
        step = dict(arrowstyle="-|>", color=c.ink, lw=1.1, shrinkA=0, shrinkB=0)
        ax.annotate("", xy=(t_n, S_on), xytext=(t_n, S_off), arrowprops=step)
        ax.annotate("", xy=(t_n / ratio, S_on), xytext=(t_n, S_on), arrowprops=step)
        t_n /= ratio
    ax.text(0.815, 0.46, "isothermal\nmagnetization", fontsize=fs.BASE - 1, va="center")
    ax.text(0.42, 0.255, "adiabatic\ndemagnetization", fontsize=fs.BASE - 1, ha="left", va="top")

    ax.set_xlabel(r"Temperature, $k_{\mathrm{B}}T/\varepsilon_0$")
    ax.set_ylabel(r"Entropy per magnetic moment, $S/(Nk_{\mathrm{B}})$")
    ax.set_xlim(0, 1.2)
    ax.set_ylim(0, 0.72)
    return fig

fs.show(draw, alt="Entropy per magnetic moment against temperature for a model paramagnet. The curve with the field "
                  "on lies below the curve with the field off, and both go to zero at zero temperature. A staircase "
                  "of vertical and horizontal arrows between the curves starts at a reduced temperature of 0.8 and "
                  "halves the temperature at each step.")
```

---

### Absolute Entropies

With $S(0) = 0$ established for ideal crystals, we can compute the **absolute (third-law) entropy** of any substance at temperature $T$ and pressure $P^\circ$:

```{math}
:label: eq:absolute-entropy
S^\circ(T) = \int_0^T \frac{C_P(T')}{T'}\,dT' + \sum_{\text{transitions}} \frac{\Delta H_{\mathrm{trs}}}{T_{\mathrm{trs}}},
```

where the sum runs over any phase transitions (melting, boiling, etc.) between $0$ and $T$, each contributing $\Delta H_{\mathrm{trs}}/T_{\mathrm{trs}}$ to the entropy.

[](#fig:nitrogen-entropy) applies Eq. {eq}`eq:absolute-entropy` to nitrogen at $1\ \mathrm{atm}$, with the calorimetric measurements of @giauque1933. Nitrogen has two solid phases, $\alpha$ below $35.61\ \mathrm{K}$ and $\beta$ from there to the melting point at $63.14\ \mathrm{K}$, and it boils at $77.32\ \mathrm{K}$. The measurements begin near $15\ \mathrm{K}$ and end at the boiling point. Below $10\ \mathrm{K}$, the authors extrapolated the heat capacity with the Debye model of a crystal, which reduces to the $T^3$ law of Example 1 below at the lowest temperatures. Above the boiling point, the curve is that of the ideal gas with $C_P = \tfrac{7}{2}R$ ([Section 5.3](section-03.md)). The three transitions contribute $6.4$, $11.4$, and $72.1\ \mathrm{J\,mol^{-1}\,K^{-1}}$, the last of which is more than a third of the total. The result at $298.15\ \mathrm{K}$ is $192.3\ \mathrm{J\,mol^{-1}\,K^{-1}}$, which is $0.4\%$ above the value in the [NIST-JANAF Thermochemical Tables](https://janaf.nist.gov/), $191.61\ \mathrm{J\,mol^{-1}\,K^{-1}}$.

% Try this: Take the rotation away from the gas: change `3.5 * R` to `2.5 * R`. By how much does the entropy at 298.15 K fall, and which integral gives this difference?

```{code-cell} ipython3
:label: fig:nitrogen-entropy
:caption: Absolute entropy of nitrogen at $1\ \mathrm{atm}$ from the heat capacities and transition enthalpies measured by @giauque1933. (a) The integrand $C_P/T$ of Eq. [](#eq:absolute-entropy). Symbols are the measurements, and the shaded area under each branch is the entropy gained within that phase. Dashed lines are calculated: the Debye model below $10\ \mathrm{K}$ and the ideal gas above the boiling point. (b) The entropy, which rises continuously within each phase and increases by $\Delta H_{\mathrm{trs}}/T_{\mathrm{trs}}$ (labels, in $\mathrm{J\,mol^{-1}\,K^{-1}}$) at the solid–solid transition, at melting, and at boiling. The step at the boiling point includes the correction of $0.9\ \mathrm{J\,mol^{-1}\,K^{-1}}$ that the authors applied to convert the vapor to an ideal gas.
:tags: [hide-input]

import sys
sys.path.append("..")                              # the book folder, which holds figstyle
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.constants import R
from scipy.integrate import cumulative_trapezoid, quad
import figstyle as fs

# Calorimetric data of Giauque and Clayton (1933) for nitrogen at 1 atm
cal = 4.185                                        # J, the calorie of the paper
df = pd.read_csv("../_static/chapter-05/section-02/nitrogen-heat-capacity.csv", comment="#")
T_trs, T_fus, T_vap, T_end = 35.61, 63.14, 77.32, 298.15      # K
dH = {"transition": 54.71 * cal, "fusion": 172.3 * cal, "vaporization": 1332.9 * cal}   # J/mol
dS_gas = 0.22 * cal                                # J/(mol K), correction of the vapor to the ideal gas

def debye(T, theta=68.0):
    """Debye heat capacity in J/(mol K), used by the paper below 10 K."""
    x = theta / T
    return 9 * R * quad(lambda u: u**4 * np.exp(-u) / np.expm1(-u)**2, 0, x)[0] / x**3

def branch(phase, T_lo, T_hi, C_lo=None):
    """Temperature grid and heat capacity (J/(mol K)) of one phase, by linear interpolation."""
    d = df[df["phase"] == phase].sort_values("T_K")
    t, C = d["T_K"].to_numpy(), cal * d["Cp_cal"].to_numpy()
    lo = np.polyval(np.polyfit(t[:3], C[:3], 1), T_lo) if C_lo is None else C_lo   # linear extrapolation to the ends
    hi = np.polyval(np.polyfit(t[-3:], C[-3:], 1), T_hi)
    keep = (t > T_lo) & (t < T_hi)
    grid = np.linspace(T_lo, T_hi, 400)
    return grid, np.interp(grid, np.r_[T_lo, t[keep], T_hi], np.r_[lo, C[keep], hi])

T_D = np.linspace(0.05, 10, 200)
C_D = np.array([debye(T) for T in T_D])
T_a, C_a = branch("alpha", 10, T_trs, debye(10))
T_b, C_b = branch("beta", T_trs, T_fus)
T_l, C_l = branch("liquid", T_fus, T_vap)
T_g = np.linspace(T_vap, T_end, 400)
C_g = np.full_like(T_g, 3.5 * R)                   # ideal gas: translation and rotation

# Entropy: integrate C_P/T within each phase and add Delta H/T at each transition
segments, jumps, S0 = [], [], 0.0
steps = {"alpha": dH["transition"] / T_trs, "beta": dH["fusion"] / T_fus, "liquid": dH["vaporization"] / T_vap + dS_gas}
for name, T, C in (("debye", T_D, C_D), ("alpha", T_a, C_a), ("beta", T_b, C_b), ("liquid", T_l, C_l), ("gas", T_g, C_g)):
    start = S0 if name != "debye" else C[0] / 3    # T^3 law below the first grid point
    S = start + cumulative_trapezoid(C / T, T, initial=0)
    segments.append((name, T, C, S))
    S0 = S[-1]
    if name in steps:
        jumps.append((T[-1], S0, S0 + steps[name]))
        S0 += steps[name]

def draw(c):
    colors = {"debye": c.solid, "alpha": c.solid, "beta": c.tint(c.solid, 0.38), "liquid": c.liquid, "gas": c.gas}
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(fs.WIDTH, 3.4), constrained_layout=True)

    # (a) The integrand C_P/T; the area under it is the entropy gained within each phase
    for name, T, C, S in segments:
        ax1.fill_between(T, C / T, color=colors[name], alpha=0.22, lw=0)
        ax1.plot(T, C / T, "--" if name in ("debye", "gas") else "-", color=colors[name], lw=0.8 * fs.LINE)
    for phase, marker in (("alpha", "o"), ("beta", "s"), ("liquid", "^")):
        d = df[df["phase"] == phase]
        ax1.plot(d["T_K"], cal * d["Cp_cal"] / d["T_K"], marker, color=colors[phase], ms=4.2, mec=c.ground, mew=0.6, zorder=3)
    for x, label in ((21, r"solid $\alpha$"), (49, r"solid $\beta$"), (70.2, "liquid"), (89, "gas")):
        ax1.text(x, 0.06, label, ha="center", fontsize=fs.BASE - 1, rotation=90 if label == "liquid" else 0)
    ax1.set_xlabel("Temperature (K)")
    ax1.set_ylabel(r"$C_P/T$ (J mol$^{-1}$ K$^{-2}$)")
    ax1.set_xlim(0, 100)
    ax1.set_ylim(0, 1.4)
    fs.panel(ax1, "a", "Integrand")

    # (b) The absolute entropy
    for name, T, C, S in segments:
        ax2.plot(T, S, "-", color=colors[name])
    for (T_x, S_lo, S_hi), label in zip(jumps, ("transition", "fusion", "vaporization")):
        ax2.plot([T_x, T_x], [S_lo, S_hi], ":", color=c.muted, lw=1.2)
        ax2.annotate(f"{label}\n+{dH[label] / T_x:.1f}", (T_x, (S_lo + S_hi) / 2), xytext=(8, 0),
                     textcoords="offset points", va="center", fontsize=fs.BASE - 1)
    ax2.plot(T_end, S0, "o", color=colors["gas"], ms=6.5, mec=c.ground, mew=0.8)
    ax2.annotate(f"{S0:.1f}", (T_end, S0), xytext=(-6, -15), textcoords="offset points", ha="right", fontsize=fs.BASE - 1)
    ax2.set_xlabel("Temperature (K)")
    ax2.set_ylabel(r"$S°$ (J mol$^{-1}$ K$^{-1}$)")
    ax2.set_xlim(0, 310)
    ax2.set_ylim(0, 210)
    fs.panel(ax2, "b", "Absolute entropy")
    return fig

fs.show(draw, alt="Two panels for nitrogen at 1 atm. Panel a: the heat capacity divided by temperature against "
                  "temperature from 0 to 100 K, with measured points for the two solid phases and the liquid and "
                  "shaded areas under each branch. Panel b: the entropy against temperature from 0 to 298 K, which "
                  f"rises within each phase, steps up at the three transitions, and reaches {S0:.1f} joules per mole "
                  "per kelvin.")
```

These absolute entropies are tabulated in standard references ([NIST Chemistry WebBook](https://webbook.nist.gov/chemistry/), NIST-JANAF tables) as $S^\circ(298.15\ \mathrm{K})$. Combined with standard enthalpies of formation, they provide the reaction entropies and Gibbs energies that we will use in [Section 5.3](section-03.md) to analyze the Haber–Bosch process.

---

## Worked Examples

### Example 1: Entropy of a Crystal at Low Temperature

**Problem.** Calorimetric measurements cannot reach $0\,\mathrm{K}$, so the lowest part of the integral in Eq. {eq}`eq:absolute-entropy` is obtained by extrapolation. At sufficiently low temperature, the heat capacity of a nonmetallic crystal is proportional to $T^3$ (the Debye $T^3$ law): $C_P = aT^3$. A crystal has a measured molar heat capacity of $0.24\ \mathrm{J\,mol^{-1}\,K^{-1}}$ at $10\ \mathrm{K}$, and the $T^3$ law holds below this temperature. Find its molar entropy at $10\ \mathrm{K}$, assuming an ideal crystal.

**Solution.**

For an ideal crystal $S(0) = 0$, and no phase transition lies below $10\ \mathrm{K}$, so

```{math}
S(T) = \int_0^T \frac{C_P(T')}{T'}\,dT' = \int_0^T a\,T'^2\,dT' = \frac{aT^3}{3} = \frac{C_P(T)}{3}.
```

At $T = 10\ \mathrm{K}$:

```{math}
S_{\mathrm{m}}(10\ \mathrm{K}) = \frac{0.24\ \mathrm{J\,mol^{-1}\,K^{-1}}}{3} = 0.080\ \mathrm{J\,mol^{-1}\,K^{-1}}.
```

**Result.** Where the $T^3$ law holds, the entropy is one third of the heat capacity at the same temperature, so a single measured value fixes the contribution of the unmeasured range. The integral converges because $C_P$ vanishes as $T \to 0$.

---

### Example 2: Residual Entropy of Ice (Pauling Model)

**Problem.** In Linus Pauling's 1935 model of ice, each oxygen atom is tetrahedrally coordinated to four neighbors, with two hydrogen atoms nearby (covalent O–H bonds) and two farther away (hydrogen bonds). Pauling showed that the number of configurations consistent with the "ice rules" (exactly two H atoms close to each O) gives a ground-state multiplicity of approximately $\Omega_0 \approx (3/2)^N$ for $N$ water molecules. Estimate the molar residual entropy.

**Solution.**

With $\Omega_0 = (3/2)^N$:

```{math}
S_{\mathrm{res}} = k_{\mathrm B}\ln\!\left(\frac{3}{2}\right)^N = Nk_{\mathrm B}\ln\!\left(\frac{3}{2}\right).
```

Per mole:

```{math}
S_{\mathrm{res,m}} = R\ln\!\left(\frac{3}{2}\right) = (8.314)\ln(1.5) = 3.37\ \mathrm{J\,mol^{-1}\,K^{-1}}.
```

**Result.** The predicted residual entropy is $3.37\ \mathrm{J\,mol^{-1}\,K^{-1}}$, in remarkable agreement with the experimental value of approximately $3.41\ \mathrm{J\,mol^{-1}\,K^{-1}}$. Without the ice rules, the $2N$ hydrogen atoms of $N$ molecules, each with two possible positions, would have $2^{2N} = 4^N$ arrangements. The ice rules reduce this number, and the agreement with experiment confirms that the proton disorder in ice is well described by the ice rules even at the lowest accessible temperatures.

```{admonition} Pauling's Count
:class: important

Reconstruct Pauling's estimate. Each hydrogen atom lies between two neighboring oxygen atoms and has two possible positions, close to one oxygen atom or close to the other. (a) Around one oxygen atom, the four neighboring hydrogen atoms have $2^4 = 16$ arrangements. How many of them place exactly two hydrogen atoms close to that oxygen atom, as the ice rules require? (b) Treat the $N$ oxygen atoms as independent, and combine the fraction from (a) with the $4^N$ unconstrained arrangements to obtain $\Omega_0$.
```

## Concept Checks

1. Why does Planck's statement allow $S(0)$ to approach a *constant* without requiring that constant to be zero?
2. What physical features distinguish an ideal crystal from a real crystal, and how do these features affect $\Omega_0$?
3. How can a system have residual entropy without violating the third law?
4. Why does the third law require $C_P \to 0$ as $T \to 0$, and how is this consistent with the quantum heat capacities derived in Section 2.6?
5. In the Pauling ice model, why is $\Omega_0 = (3/2)^N$ rather than $4^N$ or $1$?

## Key Takeaways

- The third law states that $S \to$ a constant as $T \to 0$; for ideal crystals, that constant is zero, fixing the entropy scale.
- Boltzmann's formula explains $S(0) = 0$ when the ground state is unique ($\Omega_0 = 1$).
- Residual entropy ($S_{\mathrm{res}} = k_{\mathrm B}\ln\Omega_0 > 0$) arises from ground-state degeneracy or frozen-in disorder.
- The third law requires heat capacities to vanish at absolute zero, a prediction that can be tested against the heat capacities of [Section 2.6](../chapter-02/section-06.md).
- Absolute (third-law) entropies $S^\circ(T)$ are computed by integrating $C_P/T$ from $0$ to $T$ and provide the standard entropy values used in thermochemistry.
