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


# 7.1. Equilibrium Constant

[Course-wide Conventions & Notation](../notation.md)

## Overview and Learning Objectives

Chapter 6 treated one kind of equilibrium — the coexistence of two phases of a single substance — and showed that the condition for equilibrium reduces to equality of chemical potentials, $\mu_\alpha = \mu_\beta$. **Chemical equilibrium** is the multi-species generalization: instead of two phases of one substance, we now have several chemical *species* interconverting through a reaction, and the analogous condition (derived below) is

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

The symbol $Q$ is overloaded in physical chemistry. In this chapter we consistently use $Q_p$ (with the pressure subscript) for the ideal-gas **reaction quotient** to distinguish it from the canonical partition function $Q(T,V,N)$ of Chapters 2, 4, and 5. See the [course-wide notation page](../notation.md) for the global convention.
```

Learning objectives:

- Use the extent of reaction $\xi$ to relate changes in species amounts via $dn_i=\nu_i\,d\xi$.
- Derive $\Delta_r G=\sum_i \nu_i\mu_i$ and show $(\partial G/\partial\xi)_{T,P}=\Delta_r G$.
- Derive $\Delta_r G=\Delta_r G^\circ + RT\ln Q_p$ for ideal gases and define $Q_p$ as a dimensionless product of pressure ratios.
- Use $\Delta_r G^\circ=-RT\ln K$ and the comparison of $Q_p$ and $K$ to predict reaction direction.
- Apply the van't Hoff equation to predict how $K$ changes with temperature.

## Core Ideas and Derivations

### 7.1.1 Extent of reaction $\xi$

Consider a general reaction written with stoichiometric coefficients:

```{math}
\nu_A A + \nu_B B \rightleftharpoons \nu_Y Y + \nu_Z Z.
```

A convenient way to describe composition changes is the **extent of reaction** $\xi$ (units of moles). If $n_{i,0}$ are initial amounts, then as the reaction proceeds,

```{math}
n_A = n_{A,0} - \nu_A\,\xi,\qquad
n_B = n_{B,0} - \nu_B\,\xi,\qquad
n_Y = n_{Y,0} + \nu_Y\,\xi,\qquad
n_Z = n_{Z,0} + \nu_Z\,\xi.
```

Differentiating gives the compact relation

```{math}
:label: eq:dn-from-xi
dn_i = \nu_i\,d\xi
```

if we adopt the **signed stoichiometry convention**: $\nu_i<0$ for reactants and $\nu_i>0$ for products.

```{admonition} Sign convention for $\nu_i$
:class: note

Equation {eq}`eq:dn-from-xi` requires that the stoichiometric coefficients carry signs: negative for reactants, positive for products. An equivalent formulation keeps all $\nu_i$ positive in the balanced equation and inserts minus signs by hand for reactants. Both conventions yield the same final results; we use the signed form throughout so that expressions like $\Delta_r G \equiv \sum_i \nu_i \mu_i$ and $\Delta \nu \equiv \sum_i \nu_i$ carry their natural meaning.
```

---

### 7.1.2 Gibbs free energy and the equilibrium condition

Treat the system Gibbs energy as a function of $T,P$, and composition:

```{math}
G = G(T,P,n_A,n_B,n_Y,n_Z,\ldots).
```

Its total differential is

```{math}
dG = \left(\frac{\partial G}{\partial T}\right)_{P,n}\! dT
    + \left(\frac{\partial G}{\partial P}\right)_{T,n}\! dP
    + \sum_i \left(\frac{\partial G}{\partial n_i}\right)_{T,P,n_{j\neq i}} dn_i.
```

The partial derivative $\left(\partial G/\partial n_i\right)_{T,P,n_{j\ne i}}$ was introduced in § 6.1.2 as the **chemical potential** $\mu_i$, so

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

#### Spontaneity and equilibrium (constant $T,P$)

- $\Delta_r G < 0$: reaction proceeds **forward** (toward products).
- $\Delta_r G = 0$: **equilibrium**.
- $\Delta_r G > 0$: reaction proceeds **backward** (toward reactants).

At equilibrium,

```{math}
\boxed{\Delta_r G = 0.}
```

This is the multi-species analog of the phase-equilibrium condition $\mu_\alpha = \mu_\beta$ from § 6.1.2. Instead of transferring matter between two phases, we now transfer it among several chemical species; equilibrium is again the stationary point of $G$ at fixed $T,P$.

---

### 7.1.3 Reaction quotient $Q_p$ and the equilibrium constant $K$

To connect $\Delta_r G$ to measurable composition variables, we need an expression for the chemical potentials.

#### Ideal-gas chemical potential

Section 5.3 showed that for one mole of an ideal gas at constant $T$, integrating $dG = V\,dP$ from $P^\circ$ to $P$ gives $G(T,P) - G^\circ(T) = RT\ln(P/P^\circ)$. For a single species in an ideal-gas *mixture*, the same argument with $P_i$ in place of $P$ yields

```{math}
:label: eq:mu-ideal-gas
\mu_i(T,P_i)=\mu_i^{\circ}(T)+RT\ln\left(\frac{P_i}{P^{\circ}}\right),
```

where $P^{\circ}$ is the standard-state pressure (1 bar, per the [notation conventions](../notation.md)) and $\mu_i^{\circ}(T)$ is the standard chemical potential.

````{admonition} Why is $K$ unitless if it's built from pressures?
:class: note

Any logarithm in thermodynamics must take a **dimensionless** argument. That's why we never write $\ln P$ (pressure has units), and instead write $\ln(P/P^\circ)$, where $P^\circ$ is the **standard-state pressure** (1 bar here).

**Numerical example.** Take a gas with partial pressure $P_i = 2.50\ \text{bar}$ at $T=298\ \text{K}$, with $P^\circ = 1.00\ \text{bar}$:

```{math}
\frac{P_i}{P^\circ}=\frac{2.50\ \text{bar}}{1.00\ \text{bar}}=2.50 \quad \text{(dimensionless)}.
```

The pressure contribution to the chemical potential is

```{math}
\mu_i-\mu_i^\circ = RT\ln(2.50) \approx (2.48\ \text{kJ mol}^{-1})(0.916) = 2.27\ \text{kJ mol}^{-1}.
```

**Why this matters:** the ratio makes the result **unit-independent**. Computing in pascals instead of bar gives the same ratio $2.50$ and therefore the same $\ln(2.50)$. Without the ratio, changing units would (incorrectly) change the number inside the logarithm.

**Later generalization.** For non-ideal systems, we keep the same *structure* but replace the ratio $P_i/P^\circ$ with a more general **activity** $a_i$:

```{math}
\mu_i=\mu_i^\circ + RT\ln a_i,\qquad Q=\prod_i a_i^{\nu_i}.
```

- ideal gas: $a_i = P_i/P^\circ$
- real gas: $a_i = f_i/P^\circ$, where $f_i$ is the **fugacity** (often $f_i = \phi_i P_i$)
- solutions: $a_i$ is built from concentration or mole fraction times an activity coefficient (e.g., $a_i=\gamma_i c_i/c^\circ$ or $a_i=\gamma_i x_i$).

We won't develop activities/fugacities in detail here, but this is where they will "plug in" later.
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

For an ideal-gas reaction this $K$ is the pressure-based equilibrium constant (and equals the limit of $Q_p$ at equilibrium).

#### Using $Q_p$ vs $K$ to predict direction

Since $\Delta_r G = RT\ln(Q_p/K)$:

- If $Q_p<K$, then $\Delta_r G<0$ and the reaction proceeds **toward products**.
- If $Q_p>K$, then $\Delta_r G>0$ and the reaction proceeds **toward reactants**.
- If $Q_p=K$, the system is at **equilibrium**.

---

### 7.1.4 Worked example: gas-phase dimerization of $\mathrm{NO_2}$

Consider the equilibrium

```{math}
2\,\mathrm{NO_2(g)} \rightleftharpoons \mathrm{N_2O_4(g)}.
```

($\mathrm{NO_2}$ is brown; $\mathrm{N_2O_4}$ is colorless.)

#### Thermodynamics and $K$ at 298.15 K

At $298.15\,\mathrm{K}$ (values from § 5.3-style tabulated data):

```{math}
\Delta_r H^{\circ} = -57.1\ \mathrm{kJ\,mol^{-1}},\qquad
\Delta_r S^{\circ} = -175.7\ \mathrm{J\,K^{-1}\,mol^{-1}},\qquad
\Delta_r G^{\circ} = -4.74\ \mathrm{kJ\,mol^{-1}},
```

giving

```{math}
K = \exp\!\left(-\frac{\Delta_r G^{\circ}}{RT}\right)
\approx \exp\!\left(\frac{4.74\times 10^3}{(8.314)(298.15)}\right)
\approx 6.74.
```

#### Equilibrium composition at a specified total pressure

Suppose we start with 2 mol of $\mathrm{NO_2}$ and 0 mol of $\mathrm{N_2O_4}$. Build an ICE table in terms of $\xi$:

- Initial: $n_{\mathrm{NO_2},0}=2$, $n_{\mathrm{N_2O_4},0}=0$
- Change: $\Delta n_{\mathrm{NO_2}}=-2\xi$, $\Delta n_{\mathrm{N_2O_4}}=+\xi$
- Equilibrium: $n_{\mathrm{NO_2}}=2-2\xi$, $n_{\mathrm{N_2O_4}}=\xi$, $n_{\text{tot}}=2-\xi$

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
\xi_{eq} \approx 0.81,
```

so

```{math}
n_{\mathrm{NO_2},eq}\approx 0.38, \qquad n_{\mathrm{N_2O_4},eq}\approx 0.81.
```

#### $G(\xi)$ has a minimum at $\xi_{\mathrm{eq}}$

The equilibrium condition $\Delta_r G = 0$ is the statement that $G(\xi)$ has a stationary point at equilibrium. For this system, with mole fractions $y_i(\xi)$ and $P=P^\circ$, the Gibbs free energy relative to pure reactants is

```{math}
:label: eq:G-xi-NO2
G(\xi) - G(0)
= \xi\,\Delta_r G^{\circ}
+ RT\!\left[(2-2\xi)\ln y_{\mathrm{NO_2}}(\xi) + \xi\ln y_{\mathrm{N_2O_4}}(\xi)\right].
```

The first term is the "bookkeeping" shift in standard chemical potentials; the second is the **ideal-mixing entropy** that prevents the reaction from running to completion. Plotting:

```{code-cell} ipython3
import numpy as np
import matplotlib.pyplot as plt

R = 8.314462618          # J mol^-1 K^-1
T = 298.15
dG_std = -4.74e3         # J per mole of reaction
K_eq = np.exp(-dG_std / (R*T))

def G_of_xi(xi):
    n_NO2 = 2 - 2*xi
    n_N2O4 = xi
    n_tot = 2 - xi
    mix = 0.0
    if n_NO2 > 0:
        mix += n_NO2 * np.log(n_NO2 / n_tot)
    if n_N2O4 > 0:
        mix += n_N2O4 * np.log(n_N2O4 / n_tot)
    return xi*dG_std + R*T*mix   # J (per initial 2 mol NO2)

xis = np.linspace(1e-4, 1 - 1e-4, 400)
G_vals = np.array([G_of_xi(x) for x in xis]) / 1e3   # kJ

# Analytical xi_eq from K = xi(2-xi)/(2-2xi)^2
from scipy.optimize import brentq
xi_eq = brentq(lambda x: x*(2-x)/(2-2*x)**2 - K_eq, 0.01, 0.99)

fig, ax = plt.subplots(figsize=(6.0, 4.0))
ax.plot(xis, G_vals, lw=2)
ax.axvline(xi_eq, color="k", ls="--", lw=1)
ax.plot([xi_eq], [G_of_xi(xi_eq)/1e3], "o", color="tab:red", zorder=5)
ax.annotate(fr"$\xi_{{\mathrm{{eq}}}} \approx {xi_eq:.2f}$",
            xy=(xi_eq, G_of_xi(xi_eq)/1e3),
            xytext=(xi_eq - 0.35, G_of_xi(xi_eq)/1e3 + 1.0),
            arrowprops=dict(arrowstyle="->", lw=1))
ax.set_xlabel(r"Extent of reaction $\xi$ (mol)")
ax.set_ylabel(r"$G(\xi) - G(0)$ (kJ)")
ax.set_title(r"$2\,\mathrm{NO_2}(g) \rightleftharpoons \mathrm{N_2O_4}(g)$, " 
             r"$T=298.15$ K, $P=P^\circ$")
fig.subplots_adjust(left=0.14, right=0.96, top=0.90, bottom=0.14)
plt.show()
```

The minimum of $G(\xi)$ lies precisely at the $\xi_{eq}$ computed from $K$. Note that the curve is asymmetric: the left branch is driven by $\Delta_r G^\circ < 0$ (pushing $\xi$ upward from pure reactants), while the right branch is pulled back up by the $RT\sum n_i \ln y_i$ mixing term as the mixture becomes dominated by $\mathrm{N_2O_4}$. The "tug-of-war" between these two terms is *the* reason chemical equilibrium is generically an interior minimum rather than either pure reactants or pure products.

#### Interactive visualization

<div style="width: 100%; border: 1px solid #cbd5e1; border-radius: 12px; overflow: hidden; box-shadow: 0 4px 12px rgba(0,0,0,0.08);">
  <iframe
    src="https://chem-4020-5020-auqq.vercel.app/"
    title="Equilibrium Visualization"
    style="width: 100%; height: 1200px; border: 0;"
    loading="lazy"
  ></iframe>
</div>

---

### 7.1.5 Temperature dependence: the van't Hoff equation

The temperature dependence of $K$ follows from the Gibbs–Helmholtz relation:

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

#### Example trend for $2\,\mathrm{NO_2}\rightleftharpoons\mathrm{N_2O_4}$

At $P=P^{\circ}$:

| $T$ (K) | $K$ | $\xi_{eq}$ (starting from 2 mol NO$_2$) |
| ---: | ---: | ---: |
| 250 | 205.18 | 0.97 |
| 298.15 | 6.74 | 0.81 |
| 350 | 0.17 | 0.23 |

Because $\Delta_r H^{\circ}<0$ for dimerization, increasing temperature drives $K$ **down** and shifts equilibrium toward $\mathrm{NO_2}$ (less $\mathrm{N_2O_4}$). This is the quantitative form of Le Châtelier's principle for an exothermic reaction.

---

### 7.1.6 Computing $K$ from tabulated thermochemistry

In § 7.1.3 we derived the key link between equilibrium and thermochemistry:

```{math}
K(T) = \exp\!\left(-\frac{\Delta_r G^{\circ}(T)}{RT}\right).
```

So if you can compute $\Delta_r G^{\circ}$ at a temperature of interest, you can immediately compute the equilibrium constant. This subsection is a practical workflow for *getting* $\Delta_r G^{\circ}$ from tabulated thermochemical data.

#### Step 0: Write the reaction (with phases) and choose a standard state

1. **Balance the reaction** and include phases, e.g. $\mathrm{NH_3(g)}$ vs $\mathrm{NH_3(\ell)}$.
2. Use a consistent **standard state** across all species (most modern tables use $P^\circ=1\ \mathrm{bar}$ at $T=298.15\ \mathrm{K}$; always check the table header).

#### Step 1: Pull species data at the same temperature

For each species $i$ in the balanced equation, collect either:

- **Option A (most common at 298.15 K):** $\Delta_f H_i^\circ$ and $S_i^\circ$ (and optionally $\Delta_f G_i^\circ$).
- **Option B (temperature-dependent, e.g. $T\neq 298.15\,\mathrm{K}$):** $H_i^\circ(T)$ and $S_i^\circ(T)$ (or directly $G_i^\circ(T)$ if provided).

Where to get the data:

- **NIST Chemistry WebBook**: fast way to grab $\Delta_f H^\circ(298.15\ \mathrm{K})$ and $S^\circ(298.15\ \mathrm{K})$ for many species (check the phase!).
- **NIST–JANAF Thermochemical Tables**: best when you need $H^\circ(T)$ and $S^\circ(T)$ over a range of temperatures (not just 298.15 K).
- **ATcT (Active Thermochemical Tables)**: best when you care about high-accuracy formation thermochemistry (often includes uncertainties). In practice, many workflows use **ATcT for $\Delta_f H^\circ$** and **JANAF for $S^\circ(T)$**.

#### Step 2: Convert species data into reaction values

Once you have consistent species data, compute reaction properties using stoichiometric coefficients $\nu_i$ (positive for products, negative for reactants):

```{math}
\Delta_r H^\circ = \sum_i \nu_i\,H_i^\circ
\qquad\text{and}\qquad
\Delta_r S^\circ = \sum_i \nu_i\,S_i^\circ.
```

If your data are **enthalpies of formation**, this becomes the familiar "products minus reactants" form:

```{math}
\Delta_r H^\circ
= \sum_{\text{products}} \nu_p\,\Delta_f H_p^\circ
-\sum_{\text{reactants}} \nu_r\,\Delta_f H_r^\circ.
```

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

#### Worked example: ammonia formation at 298.15 K

For the Haber–Bosch reaction (as written in § 5.3),

```{math}
\mathrm{N_2(g) + 3H_2(g) \rightarrow 2NH_3(g)}.
```

Using the reaction values quoted in § 5.3 at 298.15 K,

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

**Interpretation.** $K\gg 1$ at room temperature, so the equilibrium strongly favors $\mathrm{NH_3}$ *at the standard state*. (Section 5.3 explains why temperature and pressure "knobs" still matter for industrial operation, where the 500 °C kinetic sweet spot drives $K$ down dramatically.)

---

#### Common pitfalls (and quick fixes)

- **Wrong phase:** $\mathrm{H_2O(\ell)}$ and $\mathrm{H_2O(g)}$ (or $\mathrm{NH_3(\ell)}$ vs $\mathrm{NH_3(g)}$) have very different $S^\circ$ and $\Delta_f H^\circ$.
- **Mixed standard states:** make sure all species use the same $P^\circ$ convention and the same temperature.
- **Unit mismatches:** always reconcile kJ vs J and "per mole of reaction" vs "per mole of species."
- **Sign mistakes:** products minus reactants (or use $\nu_i$ with the sign convention consistently).

## Concept Checks

1. Why is $\Delta_r G=0$ the correct equilibrium condition at constant $T,P$? What would change if the reaction instead ran at constant $T,V$? (The constant-$T,V$ case is developed in § 7.2.)
2. Why must $Q_p$ and $K$ be dimensionless? What role does $P^\circ$ play?
3. How does comparing $Q_p$ to $K$ predict the direction of spontaneous change?
4. Inspecting Eq. {eq}`eq:G-xi-NO2`: which term in $G(\xi)$ would vanish if the mixture were not ideal (i.e., if there were no entropy of mixing)? Why does this entropy term prevent the reaction from running to completion even when $\Delta_r G^\circ \ll 0$?
5. What changes in the derivation of $\Delta_r G = \Delta_r G^\circ + RT\ln Q_p$ if the mixture is non-ideal and activities replace partial-pressure ratios?

## Key Takeaways

- Extent of reaction $\xi$ provides a compact way to track composition changes via stoichiometry: $dn_i = \nu_i\,d\xi$.
- At constant $T,P$, equilibrium corresponds to minimizing $G$, giving $\Delta_r G = \sum_i\nu_i\mu_i = 0$ — the multi-species analog of the phase-equilibrium condition of § 6.1.2.
- For ideal gases, $\Delta_r G = \Delta_r G^\circ + RT\ln Q_p$ and $\Delta_r G^\circ = -RT\ln K$.
- The sign of $\Delta_r G$ (equivalently the comparison $Q_p$ vs $K$) determines spontaneous direction.
- The temperature dependence of $K$ is governed by van't Hoff: $\partial\ln K/\partial T = \Delta_r H^\circ/(RT^2)$.
- Thermochemical tables (NIST WebBook, JANAF, ATcT) provide the $\Delta_f H^\circ$ and $S^\circ$ values needed to compute $\Delta_r G^\circ$ and hence $K$.
