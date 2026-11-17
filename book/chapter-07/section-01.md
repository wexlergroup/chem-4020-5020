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
  enumerator: 7.1.%s
---


# 7.1. Equilibrium Constant

See also: [Course-wide Conventions & Notation](../notation.md)

The course map shows where this module sits in the course and which relations connect the microscopic and macroscopic descriptions. Select a module to unfold its sections.

<iframe
  src="https://wexlergroup.github.io/chem-4020-5020/computational-studios/course-map.html?module=7&amp;here=7.1"
  title="Course map, opened at Module 7"
  style="width: 100%; height: 830px; border: 0;"
  loading="lazy"
></iframe>

```{image} ../_static/course-map/module-7.png
:class: hidden
:alt: Course map with Module 7 unfolded. Statistical mechanics (Module 2) forms the microscopic strand, and Modules 1 and 3 to 7 form the macroscopic strand. A relation joins the two strands at each module except Module 6, and for Module 7 it is the equilibrium constant in terms of molecular partition functions. Below the map, one card for each of Sections 7.1 and 7.2 lists the microscopic content, the connecting relations, and the macroscopic content of the section.
```

## Overview and Learning Objectives

Module 6 treated one kind of equilibrium — the coexistence of two phases of a single substance — and showed that the condition for equilibrium reduces to equality of chemical potentials, $\mu_\alpha = \mu_\beta$. **Chemical equilibrium** is the multi-species generalization: instead of two phases of one substance, we now have several chemical *species* interconverting through a reaction, and the analogous condition (derived below) is

```{math}
\sum_i \nu_i\,\mu_i = 0.
```

In this section we translate that abstract condition into the quantitative tools students have seen in general chemistry: the extent of reaction $\xi$, the reaction quotient $Q_p$, the equilibrium constant $K$, and the temperature dependence given by the van't Hoff equation. The key ideas are:

- A reaction at constant $T,P$ proceeds in the direction that **decreases** the Gibbs free energy of the system.
- We track progress with the **extent of reaction** $\xi$.
- At constant $T$ and $P$, equilibrium occurs when the Gibbs free energy is minimized, i.e.

  ```{math}
  \left(\frac{\partial G}{\partial \xi}\right)_{T,P}=0.
  ```

This minimization gives the **chemical equilibrium condition** $\Delta_r G=0$, and (for ideal gases) the familiar relationship

```{math}
\Delta_r G = \Delta_r G^{\circ} + RT\ln Q_p,
\qquad
\Delta_r G^{\circ} = -RT\ln K.
```

```{admonition} A note on the symbol $Q$
:class: warning

The symbol $Q$ is overloaded in physical chemistry. In this module we consistently use $Q_p$ (with the pressure subscript) for the ideal-gas **reaction quotient** to distinguish it from the canonical partition function $Q(T,V,N)$ of Modules 2, 4, and 5. See the [course-wide notation page](../notation.md) for the global convention.
```

Learning objectives:

- Use the extent of reaction $\xi$ to relate changes in species amounts via $dn_i=\nu_i\,d\xi$.
- Derive $\Delta_r G=\sum_i \nu_i\mu_i$ and show $(\partial G/\partial\xi)_{T,P}=\Delta_r G$.
- Derive $\Delta_r G=\Delta_r G^\circ + RT\ln Q_p$ for ideal gases and define $Q_p$ as a dimensionless product of pressure ratios.
- Use $\Delta_r G^\circ=-RT\ln K$ and the comparison of $Q_p$ and $K$ to predict reaction direction.
- Apply the van't Hoff equation to predict how $K$ changes with temperature.

## Core Ideas and Derivations

(section:extent-of-reaction)=
### Extent of Reaction $\xi$

Consider a general reaction, written with positive coefficients $a$, $b$, $y$, and $z$:

```{math}
a\,\mathrm{A} + b\,\mathrm{B} \rightleftharpoons y\,\mathrm{Y} + z\,\mathrm{Z}.
```

Each species is given a **stoichiometric coefficient** $\nu_i$ that carries a sign, negative for reactants and positive for products: $\nu_{\mathrm{A}} = -a$, $\nu_{\mathrm{B}} = -b$, $\nu_{\mathrm{Y}} = y$, and $\nu_{\mathrm{Z}} = z$.

A convenient way to describe composition changes is the **extent of reaction** $\xi$ (units of moles). If $n_{i,0}$ are initial amounts, then as the reaction proceeds, the amount of every species is

```{math}
n_i = n_{i,0} + \nu_i\,\xi.
```

Differentiating gives the compact relation

```{math}
:label: eq:dn-from-xi
dn_i = \nu_i\,d\xi.
```

```{admonition} Sign convention for $\nu_i$
:class: note

Equation {eq}`eq:dn-from-xi` requires that the stoichiometric coefficients carry signs: negative for reactants, positive for products. An equivalent formulation keeps all $\nu_i$ positive in the balanced equation and inserts minus signs by hand for reactants. Both conventions yield the same final results; we use the signed form throughout so that expressions like $\Delta_r G \equiv \sum_i \nu_i \mu_i$ and $\Delta \nu \equiv \sum_i \nu_i$ carry their natural meaning.
```

---

(section:reaction-gibbs-energy)=
### Gibbs Free Energy and the Equilibrium Condition

Treat the system Gibbs energy as a function of $T,P$, and composition:

```{math}
G = G(T,P,n_{\mathrm{A}},n_{\mathrm{B}},n_{\mathrm{Y}},n_{\mathrm{Z}},\ldots).
```

Its total differential is

```{math}
dG = \left(\frac{\partial G}{\partial T}\right)_{P,n}\! dT
    + \left(\frac{\partial G}{\partial P}\right)_{T,n}\! dP
    + \sum_i \left(\frac{\partial G}{\partial n_i}\right)_{T,P,n_{j\neq i}} dn_i.
```

The partial derivative $\left(\partial G/\partial n_i\right)_{T,P,n_{j\ne i}}$ was introduced under [Open Systems and the Chemical Potential](#section:chemical-potential) in [Section 6.1](../chapter-06/section-01.md) as the **chemical potential** $\mu_i$, so

```{math}
dG = -S\,dT + V\,dP + \sum_i \mu_i\,dn_i.
```

At constant $T$ and $P$,

```{math}
dG = \sum_i \mu_i\,dn_i.
```

Using Eq. {eq}`eq:dn-from-xi`,

```{math}
dG = \left(\sum_i \nu_i\mu_i\right)d\xi.
```

This motivates the definition of the **Gibbs free energy of reaction**:

```{math}
:label: eq:drG-def
\boxed{\Delta_r G \equiv \sum_i \nu_i\,\mu_i}
```

and therefore

```{math}
\left(\frac{\partial G}{\partial \xi}\right)_{T,P} = \Delta_r G.
```

#### Spontaneity and Equilibrium (Constant $T,P$)

- $\Delta_r G < 0$: reaction proceeds **forward** (toward products).
- $\Delta_r G = 0$: **equilibrium**.
- $\Delta_r G > 0$: reaction proceeds **backward** (toward reactants).

At equilibrium,

```{math}
\boxed{\Delta_r G = 0.}
```

This is the multi-species analog of the phase-equilibrium condition $\mu_\alpha = \mu_\beta$ of [Section 6.1](../chapter-06/section-01.md). Instead of transferring matter between two phases, we now transfer it among several chemical species; equilibrium is again the stationary point of $G$ at fixed $T,P$.

---

(section:reaction-quotient)=
### Reaction Quotient $Q_p$ and the Equilibrium Constant $K$

To connect $\Delta_r G$ to measurable composition variables, we need an expression for the chemical potentials.

#### Ideal-Gas Chemical Potential

[Section 5.3](../chapter-05/section-03.md) showed that for one mole of an ideal gas at constant $T$, integrating $dG = V\,dP$ from $P^\circ$ to $P$ gives $G(T,P) - G^\circ(T) = RT\ln(P/P^\circ)$. For a single species in an ideal-gas *mixture*, the same argument with $P_i$ in place of $P$ yields

```{math}
:label: eq:mu-ideal-gas
\mu_i(T,P_i)=\mu_i^{\circ}(T)+RT\ln\left(\frac{P_i}{P^{\circ}}\right),
```

where $P^{\circ}$ is the standard-state pressure (1 bar, per the [notation conventions](../notation.md)) and $\mu_i^{\circ}(T)$ is the standard chemical potential.

````{admonition} A Numerical Example: Pressure and the Chemical Potential
:class: note

Take a gas with partial pressure $P_i = 2.50\ \mathrm{bar}$ at $T=298\ \mathrm{K}$, with the standard-state pressure $P^\circ = 1.00\ \mathrm{bar}$:

```{math}
\frac{P_i}{P^\circ}=\frac{2.50\ \mathrm{bar}}{1.00\ \mathrm{bar}}=2.50.
```

The pressure contribution to the chemical potential is

```{math}
\mu_i-\mu_i^\circ = RT\ln(2.50) \approx (2.48\ \mathrm{kJ\,mol^{-1}})(0.916) = 2.27\ \mathrm{kJ\,mol^{-1}}.
```
````

```{admonition} Units Inside the Logarithm
:class: important
Repeat the calculation above with both pressures in pascals. Does the result change? Now suppose that the logarithm contained $P_i$ alone, without $P^\circ$: what would happen to $\mu_i-\mu_i^\circ$ when the unit of pressure is changed? Concept Check 2 asks what this means for $Q_p$ and $K$.
```

````{admonition} Later Generalization: Activities
:class: note

For non-ideal systems, the ratio $P_i/P^\circ$ in Eq. {eq}`eq:mu-ideal-gas` is replaced by a more general **activity** $a_i$:

```{math}
\mu_i=\mu_i^\circ + RT\ln a_i.
```

- ideal gas: $a_i = P_i/P^\circ$
- real gas: $a_i = f_i/P^\circ$, where $f_i$ is the **fugacity** (often $f_i = \phi_i P_i$)
- solutions: $a_i$ is built from concentration or mole fraction times an activity coefficient (e.g., $a_i=\gamma_i c_i/c^\circ$ or $a_i=\gamma_i x_i$).

We won't develop activities/fugacities in detail here. Concept Check 5 asks what this replacement does to the derivation that follows.
````

#### Derivation of $\Delta_r G = \Delta_r G^{\circ} + RT\ln Q_p$

Substitute Eq. {eq}`eq:mu-ideal-gas` into Eq. {eq}`eq:drG-def`:

```{math}
\Delta_r G
=\sum_i \nu_i\mu_i^{\circ}(T)
+RT\sum_i\nu_i\ln\left(\frac{P_i}{P^{\circ}}\right).
```

Define the **standard Gibbs energy of reaction**

```{math}
\Delta_r G^{\circ}(T)\equiv \sum_i \nu_i\mu_i^{\circ}(T),
```

and combine the logarithms to define the **reaction quotient** $Q_p$:

```{math}
RT\sum_i\nu_i\ln\left(\frac{P_i}{P^{\circ}}\right)
=RT\ln\left[\prod_i\left(\frac{P_i}{P^{\circ}}\right)^{\nu_i}\right]
\equiv RT\ln Q_p.
```

Therefore,

```{math}
\boxed{\Delta_r G = \Delta_r G^{\circ} + RT\ln Q_p.}
```

#### Equilibrium: $Q_p = K$

At equilibrium, $\Delta_r G=0$, so

```{math}
0 = \Delta_r G^{\circ} + RT\ln K
\qquad\Longrightarrow\qquad
\boxed{\Delta_r G^{\circ} = -RT\ln K,}
```

with

```{math}
\boxed{
K = \exp\!\left(-\frac{\Delta_r G^{\circ}}{RT}\right).
}
```

For an ideal-gas reaction this $K$ is the pressure-based equilibrium constant (the value of $Q_p$ at equilibrium).

#### Using $Q_p$ vs $K$ to Predict Direction

Since $\Delta_r G = RT\ln(Q_p/K)$:

- If $Q_p<K$, then $\Delta_r G<0$ and the reaction proceeds **toward products**.
- If $Q_p>K$, then $\Delta_r G>0$ and the reaction proceeds **toward reactants**.
- If $Q_p=K$, the system is at **equilibrium**.

---

(section:no2-dimerization)=
### Worked Example: Gas-Phase Dimerization of $\mathrm{NO_2}$

Consider the equilibrium

```{math}
2\,\mathrm{NO_2(g)} \rightleftharpoons \mathrm{N_2O_4(g)}.
```

($\mathrm{NO_2}$ is brown; $\mathrm{N_2O_4}$ is colorless.)

#### Thermodynamics and $K$ at 298.15 K

At $298.15\,\mathrm{K}$, the NIST-JANAF tables (see [Computing $K$ from Tabulated Thermochemistry](#section:k-from-tables)) give

```{math}
\begin{aligned}
\Delta_r H^{\circ} &= -57.11\ \mathrm{kJ\,mol^{-1}},\\
\Delta_r S^{\circ} &= -175.69\ \mathrm{J\,mol^{-1}\,K^{-1}},\\
\Delta_r G^{\circ} &= \Delta_r H^{\circ} - T\,\Delta_r S^{\circ} = -4.73\ \mathrm{kJ\,mol^{-1}},
\end{aligned}
```

giving

```{math}
K = \exp\!\left(-\frac{\Delta_r G^{\circ}}{RT}\right)
\approx \exp\!\left(\frac{4.73\times 10^3}{(8.314)(298.15)}\right)
\approx 6.74.
```

#### Equilibrium Composition at a Specified Total Pressure

Suppose we start with 2 mol of $\mathrm{NO_2}$ and 0 mol of $\mathrm{N_2O_4}$. Build an ICE table in terms of $\xi$:

- Initial: $n_{\mathrm{NO_2},0}=2$, $n_{\mathrm{N_2O_4},0}=0$
- Change: $\Delta n_{\mathrm{NO_2}}=-2\xi$, $\Delta n_{\mathrm{N_2O_4}}=+\xi$
- Equilibrium: $n_{\mathrm{NO_2}}=2-2\xi$, $n_{\mathrm{N_2O_4}}=\xi$, $n_{\mathrm{tot}}=2-\xi$

Mole fractions (ideal-gas mixture):

```{math}
y_{\mathrm{NO_2}} = \frac{2-2\xi}{2-\xi},\qquad
y_{\mathrm{N_2O_4}} = \frac{\xi}{2-\xi}.
```

At total pressure $P$, partial pressures are $P_i=y_i P$, and

```{math}
K
=\frac{(P_{\mathrm{N_2O_4}}/P^{\circ})}{(P_{\mathrm{NO_2}}/P^{\circ})^2}
= \frac{\xi(2-\xi)}{(2-2\xi)^2}\,\frac{P^{\circ}}{P}.
```

For $P=P^{\circ}$ and $K=6.74$, numerical root-finding gives

```{math}
\xi_{\mathrm{eq}} \approx 0.81,
```

so

```{math}
n_{\mathrm{NO_2,eq}}\approx 0.38, \qquad n_{\mathrm{N_2O_4,eq}}\approx 0.81.
```

#### $G(\xi)$ Has a Minimum at $\xi_{\mathrm{eq}}$

The equilibrium condition $\Delta_r G = 0$ is the statement that $G(\xi)$ has a stationary point at equilibrium. For this system, with mole fractions $y_i(\xi)$ and $P=P^\circ$, the Gibbs free energy relative to pure reactants is

```{math}
:label: eq:G-xi-NO2
G(\xi) - G(0)
= \xi\,\Delta_r G^{\circ}
+ RT\!\left[(2-2\xi)\ln y_{\mathrm{NO_2}}(\xi) + \xi\ln y_{\mathrm{N_2O_4}}(\xi)\right].
```

The first term is the "bookkeeping" shift in standard chemical potentials, which is linear in $\xi$. The second is the **ideal-mixing** term, $-T\,\Delta_{\mathrm{mix}}S$, where $\Delta_{\mathrm{mix}}S = -R\sum_i n_i\ln y_i$ is the entropy of mixing of the two ideal gases. [](#fig:gibbs-extent) shows their sum.

```{code-cell} ipython3
:label: fig:gibbs-extent
:caption: Gibbs free energy of the reacting mixture relative to pure $\mathrm{NO_2}$, from Eq. [](#eq:G-xi-NO2), for $2\,\mathrm{NO_2(g)} \rightleftharpoons \mathrm{N_2O_4(g)}$ at $298.15\ \mathrm{K}$ and $P = P^\circ$, starting from 2 mol of $\mathrm{NO_2}$. The dashed line marks $\xi_{\mathrm{eq}} = 0.81$, the value computed from $K = 6.74$.
:tags: [hide-input]

import sys

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import brentq

sys.path.append("..")            # the book folder, which holds figstyle
import figstyle as fs

R = 8.314462618          # J mol^-1 K^-1
T = 298.15               # K
dG_std = -4.73e3         # J per mole of reaction
K_eq = np.exp(-dG_std / (R * T))

def G_of_xi(xi):
    # G(xi) - G(0) in J, for 2 mol of NO2 initially and P equal to the standard pressure
    n_NO2, n_N2O4, n_tot = 2 - 2 * xi, xi, 2 - xi
    mix = n_NO2 * np.log(n_NO2 / n_tot) + n_N2O4 * np.log(n_N2O4 / n_tot)
    return xi * dG_std + R * T * mix

xis = np.linspace(1e-4, 1 - 1e-4, 400)

# Numerical root of K = xi (2 - xi) / (2 - 2 xi)^2
xi_eq = brentq(lambda x: x * (2 - x) / (2 - 2 * x)**2 - K_eq, 0.01, 0.99)
G_eq = G_of_xi(xi_eq) / 1e3   # kJ


def draw_extent(c):
    fig, ax = plt.subplots(figsize=(6.0, 4.0), constrained_layout=True)
    ax.axvline(xi_eq, color=c.gray, lw=1.2, ls=(0, (4, 3)))
    ax.plot(xis, G_of_xi(xis) / 1e3, color=c.blue, lw=fs.LINE)
    ax.plot(xi_eq, G_eq, "o", color=c.ink, ms=5.5, zorder=3)
    ax.text(xi_eq - 0.02, -2.5, f"$\\xi_{{\\mathrm{{eq}}}} = {xi_eq:.2f}$", ha="right")
    ax.set(xlim=(0, 1), xlabel="Extent of reaction $\\xi$ (mol)", ylabel="$G(\\xi) - G(0)$ (kJ)")
    return fig


fs.show(draw_extent, alt=f"Gibbs energy of the reacting mixture relative to pure nitrogen dioxide against the extent "
                         f"of reaction from 0 to 1 mole. The curve falls steeply at first, passes through a minimum "
                         f"of {G_eq:.1f} kilojoules at an extent of {xi_eq:.2f} mole, which is marked by a dashed "
                         f"line, and rises again toward pure dinitrogen tetroxide.")
```

The minimum of $G(\xi)$ lies at the $\xi_{\mathrm{eq}}$ computed from $K$, where the slope $(\partial G/\partial\xi)_{T,P} = \Delta_r G$ is zero. The curve is not symmetric, and it is steep at both ends. To see why, evaluate the slope $\Delta_r G = \Delta_r G^\circ + RT\ln Q_p$ in the limits $\xi\to 0$ and $\xi\to 1$, and decide which of the two terms of Eq. {eq}`eq:G-xi-NO2` produces that behavior. Concept Check 4 builds on the result.

#### Interactive Visualization

[](#fig:equilibrium-studio) draws the same curve, $G(\xi) - G(0)$, and lets you change the temperature and the total pressure. Move the red point with the slider for $\xi$. The slope of the red tangent line is $\Delta_r G$, and it is zero at the minimum, where $Q_p = K$. The flask beside the curve shows the mixture at the same extent.

```{iframe} https://wexlergroup.github.io/chem-4020-5020/computational-studios/equilibrium-visualization.html
:label: fig:equilibrium-studio
:width: 100%
:title: Chemical equilibrium
:placeholder: ../_static/equilibrium-visualization/equilibrium-visualization.png

The reacting mixture $2\,\mathrm{NO_2(g)} \rightleftharpoons \mathrm{N_2O_4(g)}$, starting from 2 mol of $\mathrm{NO_2}$, at a temperature and a total pressure that you choose. (a) $G(\xi) - G(0)$ from Eq. [](#eq:G-xi-NO2), with the term $-\xi RT\ln(P/P^\circ)$ added when $P \neq P^\circ$, the tangent at the chosen extent, and $\xi_{\mathrm{eq}}$ (dashed line). (b) The mixture at that extent, with 40 units of $\mathrm{NO_2}$ for the 2 mol (schematic). The brown shade of the flask follows the concentration of $\mathrm{NO_2}$. (c) $\Delta_r G^\circ$, $K$, $Q_p$, the slope $\Delta_r G$, the mole fractions, and $\xi_{\mathrm{eq}}$. $\Delta_r H^\circ$ and $\Delta_r S^\circ$ are the values at $298.15\ \mathrm{K}$, taken as constant, as in [Temperature Dependence: the van't Hoff Equation](#section:vant-hoff).
```

You can also open the studio in a new tab: [Chemical Equilibrium Studio](https://wexlergroup.github.io/chem-4020-5020/computational-studios/equilibrium-visualization.html).

---

(section:vant-hoff)=
### Temperature Dependence: the van't Hoff Equation

The temperature dependence of $K$ follows from the **Gibbs–Helmholtz relation**. Because $G = H - TS$ and $(\partial G/\partial T)_P = -S$ (from $dG = -S\,dT + V\,dP$, [Section 5.1](../chapter-05/section-01.md)),

```{math}
\left(\frac{\partial}{\partial T}\frac{G}{T}\right)_P
= \frac{1}{T}\left(\frac{\partial G}{\partial T}\right)_P - \frac{G}{T^2}
= -\frac{S}{T} - \frac{H - TS}{T^2}
= -\frac{H}{T^2}.
```

Applying this relation to each species in its standard state and summing with the coefficients $\nu_i$ gives

```{math}
\left(\frac{\partial}{\partial T}\frac{\Delta_r G^{\circ}}{T}\right)_P
= -\frac{\Delta_r H^{\circ}}{T^2}.
```

Using $\Delta_r G^{\circ}=-RT\ln K$ gives the **van't Hoff equation**:

```{math}
\boxed{
\left(\frac{\partial \ln K}{\partial T}\right)_P
= \frac{\Delta_r H^{\circ}}{RT^2}.
}
```

An equivalent integrated form between $T_1$ and $T_2$ (assuming $\Delta_r H^{\circ}$ is approximately constant over the range) is

```{math}
\boxed{
\ln\!\left[\frac{K(T_2)}{K(T_1)}\right]
= -\frac{\Delta_r H^{\circ}}{R}\!\left(\frac{1}{T_2}-\frac{1}{T_1}\right).
}
```

#### Example Trend for $2\,\mathrm{NO_2}\rightleftharpoons\mathrm{N_2O_4}$

The integrated form, with $\Delta_r H^{\circ} = -57.11\ \mathrm{kJ\,mol^{-1}}$ and $K(298.15\ \mathrm{K}) = 6.74$, gives the values below. The extent of reaction is for $P=P^{\circ}$.

| $T$ (K) | $K$ | $\xi_{\mathrm{eq}}$ (starting from 2 mol NO$_2$) |
| ---: | ---: | ---: |
| 250 | 570 | 0.98 |
| 298.15 | 6.74 | 0.81 |
| 350 | 0.22 | 0.27 |

The NIST-JANAF tables, which include the temperature dependence of $\Delta_r H^{\circ}$, give $K = 572$ at 250 K and $K = 0.22$ at 350 K, so the assumption of a constant $\Delta_r H^{\circ}$ costs less than 1 % over this range.

Because $\Delta_r H^{\circ}<0$ for dimerization, increasing temperature drives $K$ **down** and shifts equilibrium toward $\mathrm{NO_2}$ (less $\mathrm{N_2O_4}$). This is the quantitative form of Le Châtelier's principle for an exothermic reaction.

The integrated form also says that $\ln K$ plotted against $1/T$ is a straight line with slope $-\Delta_r H^{\circ}/R$, as long as $\Delta_r H^{\circ}$ is constant. [](#fig:vant-hoff) shows such a plot for the exothermic dimerization of $\mathrm{NO_2}$ and for an endothermic reaction, the formation of $\mathrm{NO}$ from $\mathrm{N_2}$ and $\mathrm{O_2}$, for which $\Delta_r H^{\circ} = +180.58\ \mathrm{kJ\,mol^{-1}}$. The line rises toward low temperature for the exothermic reaction and falls for the endothermic one. Each line uses only data at $298.15\ \mathrm{K}$, and it stays within 10 % of the tabulated $K$ over the whole range of its panel, in which $K$ changes by about six orders of magnitude.

```{code-cell} ipython3
:label: fig:vant-hoff
:caption: Plots of $\ln K$ against $1/T$ for (a) the exothermic reaction $2\,\mathrm{NO_2(g)} \rightleftharpoons \mathrm{N_2O_4(g)}$ and (b) the endothermic reaction $\mathrm{N_2(g)} + \mathrm{O_2(g)} \rightleftharpoons 2\,\mathrm{NO(g)}$. The points are from the NIST-JANAF tables. Each line is the integrated van't Hoff equation with the value of $\Delta_r H^{\circ}$ at $298.15\ \mathrm{K}$, drawn through $K(298.15\ \mathrm{K})$. The scale on top gives the temperature.
:tags: [hide-input]

import sys

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

sys.path.append("..")            # the book folder, which holds figstyle
import figstyle as fs

R = 8.314462618e-3               # kJ mol^-1 K^-1
T_ref = 298.15                   # K

# Standard Gibbs free energies of formation (kJ/mol) from the NIST-JANAF tables
T_a = np.array([250, 298.15, 350, 400, 450, 500])                       # K
dfG_NO2 = np.array([48.355, 51.258, 54.445, 57.560, 60.703, 63.867])
dfG_N2O4 = np.array([83.513, 97.787, 113.251, 128.204, 143.163, 158.109])
T_b = np.array([1000, 1200, 1500, 2000, 2500, 3000])                    # K
dfG_NO = np.array([77.775, 75.239, 71.425, 65.060, 58.720, 52.439])

lnK_a = -(dfG_N2O4 - 2 * dfG_NO2) / (R * T_a)            # 2 NO2 = N2O4
lnK_b = -(2 * dfG_NO) / (R * T_b)                        # N2 + O2 = 2 NO

# Integrated van't Hoff equation through K(298.15 K), with the reaction enthalpy at 298.15 K (kJ/mol)
drH_a, lnK_ref_a = 9.079 - 2 * 33.095, -(97.787 - 2 * 51.258) / (R * T_ref)
drH_b, lnK_ref_b = 2 * 90.291, -(2 * 86.600) / (R * T_ref)

def lnK_line(T, drH, lnK_ref):
    return lnK_ref - drH / R * (1 / T - 1 / T_ref)

def reciprocal(x):
    with np.errstate(divide="ignore"):
        return 1 / x


def draw_vant_hoff(c):
    fig, axes = plt.subplots(1, 2, figsize=(fs.WIDTH, 3.5), constrained_layout=True)
    panels = (("a", "Exothermic", T_a, lnK_a, drH_a, lnK_ref_a, [250, 300, 350, 400, 500]),
              ("b", "Endothermic", T_b, lnK_b, drH_b, lnK_ref_b, [1000, 1200, 1500, 2000, 3000]))
    for ax, (letter, title, T, lnK, drH, lnK_ref, ticks) in zip(axes, panels):
        T_fit = np.linspace(T.min(), T.max(), 200)
        ax.plot(1 / T_fit, lnK_line(T_fit, drH, lnK_ref), color=c.blue, lw=fs.LINE)
        ax.plot(1 / T, lnK, "o", color=c.ink, ms=5.5, zorder=3)
        ax.set(xlabel="$1/T$ (K$^{-1}$)", ylabel="$\\ln K$")
        ax.xaxis.set_major_locator(MaxNLocator(4))
        top = ax.secondary_xaxis("top", functions=(reciprocal, reciprocal))
        top.set_xticks(ticks)
        top.set_xlabel("Temperature (K)")
        fs.panel(ax, letter, title)
    axes[0].text(0.0030, -4.6, "$\\Delta_r H° = -57.11$ kJ mol$^{-1}$\nslope $-\\Delta_r H°/R > 0$", ha="left")
    axes[1].text(0.00036, -17.2, "$\\Delta_r H° = +180.58$ kJ mol$^{-1}$\nslope $-\\Delta_r H°/R < 0$", ha="left")
    return fig


fs.show(draw_vant_hoff, alt="Two plots of the natural logarithm of the equilibrium constant against the reciprocal "
                            "of the temperature, each with points from the NIST-JANAF tables on a straight line. "
                            "Panel a, the dimerization of nitrogen dioxide between 250 and 500 kelvin: the line "
                            f"rises from {lnK_a[-1]:.1f} at 500 kelvin to {lnK_a[0]:.1f} at 250 kelvin. Panel b, the "
                            "formation of nitric oxide from nitrogen and oxygen between 1000 and 3000 kelvin: the "
                            f"line falls from {lnK_b[-1]:.1f} at 3000 kelvin to {lnK_b[0]:.1f} at 1000 kelvin.")
```

---

(section:k-from-tables)=
### Computing $K$ from Tabulated Thermochemistry

Under [Reaction Quotient $Q_p$ and the Equilibrium Constant $K$](#section:reaction-quotient) we derived the key link between equilibrium and thermochemistry:

```{math}
K(T) = \exp\!\left(-\frac{\Delta_r G^{\circ}(T)}{RT}\right).
```

So if you can compute $\Delta_r G^{\circ}$ at a temperature of interest, you can immediately compute the equilibrium constant. This subsection is a practical workflow for *getting* $\Delta_r G^{\circ}$ from tabulated thermochemical data.

#### Step 0: Write the Reaction (with Phases) and Choose a Standard State

1. **Balance the reaction** and include phases, e.g. $\mathrm{NH_3(g)}$ vs $\mathrm{NH_3(l)}$.
2. Use a consistent **standard state** across all species (most modern tables use $P^\circ=1\ \mathrm{bar}$ at $T=298.15\ \mathrm{K}$; always check the table header).

#### Step 1: Pull Species Data at the Same Temperature

For each species $i$ in the balanced equation, collect either:

- **Option A (most common at 298.15 K):** $\Delta_f H_i^\circ$ and $S_i^\circ$ (and optionally $\Delta_f G_i^\circ$).
- **Option B (temperature-dependent, e.g. $T\neq 298.15\,\mathrm{K}$):** $H_i^\circ(T)$ and $S_i^\circ(T)$ (or directly $G_i^\circ(T)$ if provided).

Where to get the data:

- **NIST Chemistry WebBook**: fast way to grab $\Delta_f H^\circ(298.15\ \mathrm{K})$ and $S^\circ(298.15\ \mathrm{K})$ for many species (check the phase!).
- **NIST–JANAF Thermochemical Tables**: best when you need $H^\circ(T)$ and $S^\circ(T)$ over a range of temperatures (not just 298.15 K).
- **ATcT (Active Thermochemical Tables)**: best when you care about high-accuracy formation thermochemistry (often includes uncertainties). In practice, many workflows use **ATcT for $\Delta_f H^\circ$** and **JANAF for $S^\circ(T)$**.

#### Step 2: Convert Species Data into Reaction Values

Once you have consistent species data, compute reaction properties using stoichiometric coefficients $\nu_i$ (positive for products, negative for reactants):

```{math}
\Delta_r H^\circ = \sum_i \nu_i\,H_i^\circ
\qquad\text{and}\qquad
\Delta_r S^\circ = \sum_i \nu_i\,S_i^\circ.
```

If your data are **enthalpies of formation**, the same sum applies:

```{math}
\Delta_r H^\circ = \sum_i \nu_i\,\Delta_f H_i^\circ.
```

Because the $\nu_i$ of the reactants are negative, this is the familiar "products minus reactants" form of Hess's law ([Section 3.3](../chapter-03/section-03.md)).

```{admonition} Shortcut for elements
:class: tip
If an element appears in its **standard state**, then by convention $\Delta_f H^\circ = 0$ for that elemental form (e.g., $\mathrm{N_2(g)}$, $\mathrm{H_2(g)}$ at 1 bar). This often simplifies $\Delta_r H^\circ$ calculations dramatically.
```

#### Step 3: Compute $\Delta_r G^\circ$

If you have $\Delta_r H^\circ$ and $\Delta_r S^\circ$ at the same temperature, then

```{math}
\Delta_r G^\circ(T) = \Delta_r H^\circ(T) - T\,\Delta_r S^\circ(T).
```

```{admonition} Unit check
:class: warning
It is very common to have $\Delta H^\circ$ in **kJ/mol** while $S^\circ$ is in **J/(mol·K)**. Convert before combining (e.g., divide entropy by $1000$ to get kJ/(mol·K)).
```

#### Step 4: Compute $K$

Finally,

```{math}
K(T) = \exp\!\left(-\frac{\Delta_r G^\circ(T)}{RT}\right),
\qquad
\log_{10}K(T) = -\frac{\Delta_r G^\circ(T)}{(\ln 10)\,RT}.
```

For an **ideal-gas** reaction, this $K$ corresponds to a pressure-based equilibrium constant with the dimensionless ratio $(P_i/P^\circ)$ inside $Q_p$.

---

#### Worked Example: Ammonia Formation at 298.15 K

For the Haber–Bosch reaction (as written in [Section 5.3](../chapter-05/section-03.md)),

```{math}
\mathrm{N_2(g) + 3H_2(g) \rightarrow 2NH_3(g)}.
```

Using the reaction values quoted in [Section 5.3](../chapter-05/section-03.md) at 298.15 K,

```{math}
\Delta_r H^\circ \approx -92\ \mathrm{kJ\,mol^{-1}},\qquad
\Delta_r S^\circ \approx -198\ \mathrm{J\,mol^{-1}\,K^{-1}},
```

gives

```{math}
\Delta_r G^\circ
\approx -92\ \mathrm{kJ\,mol^{-1}}
- (298.15\ \mathrm{K})\!\left(-0.198\ \mathrm{kJ\,mol^{-1}\,K^{-1}}\right)
\approx -33\ \mathrm{kJ\,mol^{-1}}.
```

Then

```{math}
K(298.15\ \mathrm{K}) = \exp\!\left(-\frac{-33{,}000}{(8.314)(298.15)}\right)
\approx 6\times 10^{5}.
```

**Interpretation.** $K\gg 1$ at room temperature, so the equilibrium strongly favors $\mathrm{NH_3}$ *at the standard state*. ([Section 5.3](../chapter-05/section-03.md) explains why temperature and pressure "knobs" still matter for industrial operation, where the 500 °C kinetic sweet spot drives $K$ down dramatically.)

---

#### Common Pitfalls (and Quick Fixes)

- **Wrong phase:** $\mathrm{H_2O(l)}$ and $\mathrm{H_2O(g)}$ (or $\mathrm{NH_3(l)}$ vs $\mathrm{NH_3(g)}$) have very different $S^\circ$ and $\Delta_f H^\circ$.
- **Mixed standard states:** make sure all species use the same $P^\circ$ convention and the same temperature.
- **Unit mismatches:** always reconcile kJ vs J and "per mole of reaction" vs "per mole of species."
- **Sign mistakes:** products minus reactants (or use $\nu_i$ with the sign convention consistently).

## Concept Checks

1. Why is $\Delta_r G=0$ the correct equilibrium condition at constant $T,P$? What would change if the reaction instead ran at constant $T,V$? (The constant-$T,V$ case is developed in [Section 7.2](section-02.md).)
2. Why must $Q_p$ and $K$ be dimensionless? What role does $P^\circ$ play?
3. How does comparing $Q_p$ to $K$ predict the direction of spontaneous change?
4. Inspecting Eq. {eq}`eq:G-xi-NO2`: which term in $G(\xi)$ would vanish if there were no entropy of mixing? Why does this entropy term prevent the reaction from running to completion even when $\Delta_r G^\circ \ll 0$?
5. What changes in the derivation of $\Delta_r G = \Delta_r G^\circ + RT\ln Q_p$ if the mixture is non-ideal and activities replace partial-pressure ratios?

## Key Takeaways

- Extent of reaction $\xi$ provides a compact way to track composition changes via stoichiometry: $dn_i = \nu_i\,d\xi$.
- At constant $T,P$, equilibrium corresponds to minimizing $G$, giving $\Delta_r G = \sum_i\nu_i\mu_i = 0$ — the multi-species analog of the phase-equilibrium condition of [Section 6.1](../chapter-06/section-01.md).
- For ideal gases, $\Delta_r G = \Delta_r G^\circ + RT\ln Q_p$ and $\Delta_r G^\circ = -RT\ln K$.
- The sign of $\Delta_r G$ (equivalently the comparison $Q_p$ vs $K$) determines spontaneous direction.
- The temperature dependence of $K$ is governed by van't Hoff: $\partial\ln K/\partial T = \Delta_r H^\circ/(RT^2)$.
- Thermochemical tables (NIST WebBook, JANAF, ATcT) provide the $\Delta_f H^\circ$ and $S^\circ$ values needed to compute $\Delta_r G^\circ$ and hence $K$.
