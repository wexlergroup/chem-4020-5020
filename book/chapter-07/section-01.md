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

Chemical equilibrium is reached when a reacting system can no longer lower its Gibbs free energy at the specified constraints. This section introduces the extent of reaction, derives the equilibrium condition $\Delta_r G=0$ at constant $T,P$, and connects $\Delta_r G$ to the reaction quotient $Q$ and equilibrium constant $K$.

In this section we connect **Gibbs free energy** to (i) whether a reaction is spontaneous and (ii) how far it proceeds before reaching **chemical equilibrium**.

Key ideas:

- A reaction proceeds in the direction that **decreases** the Gibbs free energy of the system.
- We track progress using the **extent of reaction** $\xi$.
- At constant $T$ and $P$, equilibrium occurs when the Gibbs free energy is **minimized**, i.e. when

```{math}
\left(\frac{\partial G}{\partial \xi}\right)_{T,P}=0.
```

This leads to the **chemical equilibrium condition** $\Delta_r G=0$, and (for ideal gases) the familiar relationship

```{admonition} Notational warning
:class: warning

The symbol $Q$ is overloaded in this course.

- **Here (equilibrium thermodynamics):** $Q$ (or $Q_p$) is the **reaction quotient**,
  a *dimensionless* product of activities.
  For ideal gases: $Q_p=\prod_i(P_i/P^{\circ})^{\nu_i}$.

- **Elsewhere (statistical mechanics):** $Q$ (often $Z$) can denote a **partition function**.
  In Section 7.2 we use $\mathcal{Q}$ for partition functions.
```

```{math}
\Delta_r G = \Delta_r G^{\circ} + RT\ln Q,
\qquad
\Delta_r G^{\circ} = -RT\ln K.
```

---

Learning objectives:

- Use the extent of reaction $\xi$ to relate changes in species amounts via $dn_i=\nu_i\,d\xi$.
- Derive $\Delta_r G=\sum_i \nu_i\mu_i$ and show $(\partial G/\partial\xi)_{T,P}=\Delta_r G$.
- Derive $\Delta_r G=\Delta_r G^\circ + RT\ln Q$ for ideal gases and define $Q$ as a dimensionless product of activities.
- Use $\Delta_r G^\circ=-RT\ln K$ and the comparison of $Q$ and $K$ to predict reaction direction.

## Core Ideas and Derivations

### 7.1.1 Extent of reaction $\xi$

Consider a general reaction written with stoichiometric coefficients:

```{math}
\nu_A A + \nu_B B \rightleftharpoons \nu_Y Y + \nu_Z Z.
```

A convenient way to describe composition changes is with the **extent of reaction** $\xi$ (units of moles). If $n_{i,0}$ are initial amounts, then as the reaction proceeds,

```{math}
n_A = n_{A,0} - \nu_A\,\xi,\qquad
n_B = n_{B,0} - \nu_B\,\xi,\qquad
n_Y = n_{Y,0} + \nu_Y\,\xi,\qquad
n_Z = n_{Z,0} + \nu_Z\,\xi.
```

Differentiating gives the compact relation

```{math}
dn_i = \nu_i\,d\xi
```

**if** we adopt the sign convention that $\nu_i<0$ for reactants and $\nu_i>0$ for products. (Equivalently: keep $\nu_i$ positive in the balanced equation and insert minus signs for reactants; both conventions lead to the same final results.)

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

The partial derivative $\left(\frac{\partial G}{\partial n_i}\right)_{T,P}$ is the **chemical potential** $\mu_i$, so

```{math}
dG = -S\,dT + V\,dP + \sum_i \mu_i\,dn_i.
```

At constant $T$ and $P$,

```{math}
dG = \sum_i \mu_i\,dn_i.
```

Using $dn_i = \nu_i d\xi$ gives

```{math}
dG = \left(\sum_i \nu_i\mu_i\right)d\xi.
```

This motivates the definition of the **Gibbs free energy of reaction**:

```{math}
\boxed{
\Delta_r G \equiv \sum_i \nu_i\,\mu_i
}
```

and therefore

```{math}
\left(\frac{\partial G}{\partial \xi}\right)_{T,P} = \Delta_r G.
```

#### Spontaneity and equilibrium (constant $T,P$)

- $\Delta_r G < 0$: reaction proceeds **forward** (toward products)
- $\Delta_r G = 0$: **equilibrium**
- $\Delta_r G > 0$: reaction proceeds **backward** (toward reactants)

At equilibrium,

```{math}
\boxed{\Delta_r G = 0.}
```

---

### 7.1.3 Reaction quotient $Q$ and the equilibrium constant $K$

```{admonition} Notational warning
:class: warning

In this section, $Q$ (or $Q_p$) means the **reaction quotient**, not a partition function.
```

To connect $\Delta_r G$ to measurable composition variables, we need an expression for the chemical potentials.

#### Ideal-gas chemical potential

For an ideal gas,

```{math}
\mu_i(T,P_i)=\mu_i^{\circ}(T)+RT\ln\left(\frac{P_i}{P^{\circ}}\right),
```

where $P^{\circ}$ is the standard-state pressure (often 1 bar), and $\mu_i^{\circ}(T)$ is the standard chemical potential.

````{admonition} Why is $K_P$ unitless if it's built from pressures?
:class: note

Any logarithm in thermodynamics must take a **dimensionless** argument. That's why we never write $\ln P$ (pressure has units), and instead write $\ln(P/P^\circ)$, where $P^\circ$ is the **standard-state pressure** (often $1\ \text{bar}$).

**Numerical example (units cancel):**
Take a gas with partial pressure $P_i = 2.50\ \text{bar}$ at $T=298\ \text{K}$, with $P^\circ = 1.00\ \text{bar}$.

```{math}
\frac{P_i}{P^\circ}=\frac{2.50\ \text{bar}}{1.00\ \text{bar}}=2.50 \quad \text{(dimensionless)}
```

So the pressure contribution to chemical potential is

```{math}
\mu_i-\mu_i^\circ = RT\ln\!\left(\frac{P_i}{P^\circ}\right)
= RT\ln(2.50).
```

At $298\ \text{K}$, $RT \approx 2.48\ \text{kJ mol}^{-1}$, and $\ln(2.50)=0.916$, so

```{math}
\mu_i-\mu_i^\circ \approx (2.48)(0.916)=2.27\ \text{kJ mol}^{-1}.
```

**Why this matters:** the ratio makes the result **unit-independent**.
If you compute the same ratio in pascals,
```{math}
\frac{2.50\times 10^{5}\ \text{Pa}}{1.00\times 10^{5}\ \text{Pa}}=2.50,
```
you get the same $\ln(2.50)$. Without the ratio, changing units would (incorrectly) change the number inside the logarithm.

**Connection to $Q$ and $K$:**
This same idea is why the reaction quotient is written as
```{math}
Q=\prod_i\left(\frac{P_i}{P^\circ}\right)^{\nu_i},
```
so that $Q$ (and $K$) are always **dimensionless**, which is required by $\Delta_r G=\Delta_r G^\circ + RT\ln Q$.

**Later generalization (flag only):**
For non-ideal systems, we keep the same *structure* but replace the simple ratio $P_i/P^\circ$ with a more general **activity** $a_i$:
```{math}
\mu_i=\mu_i^\circ + RT\ln a_i,\qquad Q=\prod_i a_i^{\nu_i}.
```

* ideal gas: $a_i = P_i/P^\circ$
* real gas: $a_i = f_i/P^\circ$ where $f_i$ is the **fugacity** (often $f_i=\phi_i P_i$)
* solutions: $a_i$ is built from concentration/mole fraction times an activity coefficient (e.g., $a_i=\gamma_i,c_i/c^\circ$ or $a_i=\gamma_i x_i$).

(We won't develop activities/fugacities in detail yet, but this is where they will "plug in" later.)
````

#### Derivation of $\Delta_r G = \Delta_r G^{\circ} + RT\ln Q$

Substitute $\mu_i(T,P_i)$ into $\Delta_r G=\sum_i\nu_i\mu_i$:

```{math}
\Delta_r G
=\sum_i \nu_i\mu_i^{\circ}(T)
+RT\sum_i\nu_i\ln\left(\frac{P_i}{P^{\circ}}\right).
```

Define the **standard Gibbs energy of reaction**

```{math}
\Delta_r G^{\circ}(T)\equiv \sum_i \nu_i\mu_i^{\circ}(T),
```

and combine the logarithms to define the **reaction quotient** $Q$ (here, a pressure-based form $Q_p$):

```{math}
RT\sum_i\nu_i\ln\left(\frac{P_i}{P^{\circ}}\right)
=RT\ln\left[\prod_i\left(\frac{P_i}{P^{\circ}}\right)^{\nu_i}\right]
\equiv RT\ln Q_p.
```

Therefore,

```{math}
\boxed{\Delta_r G = \Delta_r G^{\circ} + RT\ln Q_p.}
```

#### Equilibrium: $Q = K$

At equilibrium $\Delta_r G=0$, so

```{math}
0 = \Delta_r G^{\circ} + RT\ln K
\qquad\Longrightarrow\qquad
\boxed{\Delta_r G^{\circ} = -RT\ln K.}
```

For the ideal-gas, pressure-based constant this is

```{math}
\boxed{
K_p = \exp\left(-\frac{\Delta_r G^{\circ}}{RT}\right).
}
```

#### Using $Q$ vs $K$ to predict direction

Since $\Delta_r G = RT\ln(Q/K)$:

- If $Q<K$, then $\Delta_r G<0$ and the reaction proceeds **toward products**.
- If $Q>K$, then $\Delta_r G>0$ and the reaction proceeds **toward reactants**.
- If $Q=K$, the system is at **equilibrium**.

---

### 7.1.4 Worked example: gas-phase dimerization of $\mathrm{NO_2}$

Consider the equilibrium

```{math}
2\,\mathrm{NO_2(g)} \rightleftharpoons \mathrm{N_2O_4(g)}.
```

($\mathrm{NO_2}$ is brown; $\mathrm{N_2O_4}$ is colorless.)

#### Thermodynamics and $K_p$ at 298.15 K

At $298.15\,\mathrm{K}$ (values from lecture notes):

```{math}
\Delta_r H^{\circ} = -57.1\ \mathrm{kJ\,mol^{-1}},\qquad
\Delta_r S^{\circ} = -175.7\ \mathrm{J\,K^{-1}\,mol^{-1}},\qquad
\Delta_r G^{\circ} = -4.74\ \mathrm{kJ\,mol^{-1}}.
```

Then

```{math}
K_p = \exp\left(-\frac{\Delta_r G^{\circ}}{RT}\right)
\approx \exp\left(-\frac{-4.74\times10^3}{(8.314)(298.15)}\right)
\approx 6.74.
```

#### Equilibrium composition at a specified total pressure

Suppose we start with **2 mol of** $\mathrm{NO_2}$ and **0 mol of** $\mathrm{N_2O_4}$.

Use an ICE table in terms of $\xi$:

- Initial: $n_{\mathrm{NO_2},0}=2$, $n_{\mathrm{N_2O_4},0}=0$
- Change: $\Delta n_{\mathrm{NO_2}}=-2\xi$, $\Delta n_{\mathrm{N_2O_4}}=+\xi$
- Equilibrium:

  ```{math}
  n_{\mathrm{NO_2}}=2-2\xi,\qquad
  n_{\mathrm{N_2O_4}}=\xi.
  ```

Total moles:

```{math}
n_{\text{tot}} = (2-2\xi)+\xi = 2-\xi.
```

Mole fractions (ideal gas mixture):

```{math}
y_{\mathrm{NO_2}} = \frac{2-2\xi}{2-\xi},\qquad
y_{\mathrm{N_2O_4}} = \frac{\xi}{2-\xi}.
```

Partial pressures in terms of the total pressure $P$:

```{math}
P_{\mathrm{NO_2}} = y_{\mathrm{NO_2}}P
=\frac{2-2\xi}{2-\xi}P,
\qquad
P_{\mathrm{N_2O_4}} = y_{\mathrm{N_2O_4}}P
=\frac{\xi}{2-\xi}P.
```

For this reaction,

```{math}
K_p
=\frac{\left(P_{\mathrm{N_2O_4}}/P^{\circ}\right)}{\left(P_{\mathrm{NO_2}}/P^{\circ}\right)^2}.
```

Substituting the expressions above and simplifying gives a single equation for $\xi$:

```{math}
\boxed{
K_p
=\frac{\xi(2-\xi)}{(2-2\xi)^2}\left(\frac{P^{\circ}}{P}\right).
}
```

If $P=P^{\circ}$ and $K_p=6.74$, solving yields

```{math}
\xi_{eq} \approx 0.81.
```

So at equilibrium (when $P=P^{\circ}$):

```{math}
n_{\mathrm{NO_2},eq} = 2-2\xi_{eq} \approx 0.38,\qquad
n_{\mathrm{N_2O_4},eq} = \xi_{eq} \approx 0.81.
```

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
\ln\left[\frac{K(T_2)}{K(T_1)}\right]
=
-\frac{\Delta_r H^{\circ}}{R}\left(\frac{1}{T_2}-\frac{1}{T_1}\right).
}
```

#### Example trend for $2\,\mathrm{NO_2}\rightleftharpoons\mathrm{N_2O_4}$

From the lecture notes, for $P=P^{\circ}$:

| $T$ (K) | $K_p$ | $\xi_{eq}$ (starting from 2 mol NO$_2$) |
| ---: | ---: | ---: |
| 250 | 205.18 | 0.97 |
| 298.15 | 6.74 | 0.81 |
| 350 | 0.17 | 0.23 |

Because $\Delta_r H^{\circ}<0$ for dimerization, increasing temperature drives $K_p$ **down** and shifts equilibrium toward $\mathrm{NO_2}$ (less $\mathrm{N_2O_4}$).

---

### 7.1.6 Visualizing equilibrium as a minimum of $G(\xi)$

If we write the Gibbs free energy of the reacting mixture as a function of extent, $G(\xi)$, then:

- the **slope** is $\left(\partial G/\partial \xi\right)_{T,P}=\Delta_r G$
- the **equilibrium extent** $\xi_{eq}$ is located where the slope is zero
- equilibrium corresponds to a **minimum** in $G(\xi)$ at fixed $T,P$

For the $\mathrm{NO_2/N_2O_4}$ system, the notes sketch $G$ vs. $\xi$ with a clear minimum at $\xi_{eq}$: the system "rolls downhill" in $G$ until it reaches that minimum.

---

### 7.1.7 Data Skills: Using NIST/JANAF/ATcT to compute $\Delta_r H^\circ$, $\Delta_r S^\circ$, $\Delta_r G^\circ$, and $K$

In Section 7.1.3 we derived the key link between equilibrium and thermochemistry:

```{math}
\Delta_r G^\circ(T) = -RT\ln K(T).
```

So if you can compute $\Delta_r G^\circ$ at a temperature of interest, you can immediately compute the equilibrium constant:

```{math}
K(T) = \exp\!\left(-\frac{\Delta_r G^\circ(T)}{RT}\right).
```

This mini-section is a practical workflow for *getting* $\Delta_r G^\circ$ from tabulated thermochemical data.

#### Step 0: Write the reaction (with phases) and choose a standard state

1. **Balance the reaction** and include phases, e.g. $\mathrm{NH_3(g)}$ vs $\mathrm{NH_3(\ell)}$.
2. Use a consistent **standard state** across all species (most modern tables use $P^\circ=1\ \mathrm{bar}$ at $T=298.15\ \mathrm{K}$; always check the table header).

#### Step 1: Pull species data at the same temperature

For each species $i$ in the balanced equation, collect *either*:

- **Option A (most common at 298.15 K):**  
  $\Delta_f H_i^\circ$ and $S_i^\circ$ (and optionally $\Delta_f G_i^\circ$)

- **Option B (temperature-dependent, e.g. $T\neq 298.15$):**  
  $H_i^\circ(T)$ and $S_i^\circ(T)$ (or directly $G_i^\circ(T)$ if provided)

Where to get the data:

- **NIST Chemistry WebBook**: fast way to grab $\Delta_f H^\circ(298.15\ \mathrm{K})$ and $S^\circ(298.15\ \mathrm{K})$ for many species (check the phase!).  
- **NIST/JANAF Thermochemical Tables**: best when you need $H^\circ(T)$ and $S^\circ(T)$ over a range of temperatures (not just 298.15 K).  
- **ATcT (Active Thermochemical Tables)**: best when you care about *high-accuracy* formation thermochemistry (often includes uncertainties). In practice, many workflows use **ATcT for $\Delta_f H^\circ$** and **JANAF for $S^\circ(T)$**.

(Links to these sources are collected in Appendix Y.)

#### Step 2: Convert species data into reaction values

Once you have consistent species data, compute reaction properties using stoichiometric coefficients $\nu_i$
(products positive, reactants negative):

```{math}
\Delta_r H^\circ = \sum_i \nu_i\,H_i^\circ
\qquad\text{and}\qquad
\Delta_r S^\circ = \sum_i \nu_i\,S_i^\circ.
```

If your data are **enthalpies of formation**, this becomes the familiar “products minus reactants” form:

```{math}
\Delta_r H^\circ
= \sum_{\text{products}} \nu_p\,\Delta_f H_p^\circ
-\sum_{\text{reactants}} \nu_r\,\Delta_f H_r^\circ.
```

> **Shortcut for elements:** if an element appears in its **standard state**, then by convention  
> $\Delta_f H^\circ = 0$ for that elemental form (e.g., $\mathrm{N_2(g)}$, $\mathrm{H_2(g)}$ at 1 bar).  
> This often simplifies $\Delta_r H^\circ$ calculations dramatically.

#### Step 3: Compute $\Delta_r G^\circ$

If you have $\Delta_r H^\circ$ and $\Delta_r S^\circ$ at the same temperature, then

```{math}
\Delta_r G^\circ(T) = \Delta_r H^\circ(T) - T\,\Delta_r S^\circ(T).
```

**Unit check:** it is very common to have $\Delta H^\circ$ in **kJ/mol** while $S^\circ$ is in **J/(mol·K)**.  
Convert before combining (e.g., divide entropy by $1000$ to get kJ/(mol·K)).

#### Step 4: Compute $K$

Finally,

```{math}
K(T) = \exp\!\left(-\frac{\Delta_r G^\circ(T)}{RT}\right),
\qquad
\log_{10}K(T) = -\frac{\Delta_r G^\circ(T)}{2.303\,RT}.
```

For an **ideal-gas** reaction, this $K$ corresponds to a pressure-based equilibrium constant $K_p$ with the
dimensionless ratio $(P_i/P^\circ)$ inside $Q$.

---

#### Worked example: ammonia formation at 298.15 K

For the Haber–Bosch reaction (as written in Section 5.3),

```{math}
\mathrm{N_2(g) + 3H_2(g) \rightarrow 2NH_3(g)}.
```

A typical 298.15 K workflow is:

1. Pull $\Delta_f H^\circ$ and $S^\circ$ for $\mathrm{N_2(g)}$, $\mathrm{H_2(g)}$, and $\mathrm{NH_3(g)}$ from NIST (or JANAF/ATcT).
2. Compute $\Delta_r H^\circ$ and $\Delta_r S^\circ$ using stoichiometry.
3. Compute $\Delta_r G^\circ = \Delta_r H^\circ - T\Delta_r S^\circ$.
4. Compute $K_p = \exp\!\left(-\Delta_r G^\circ/RT\right)$.

Using the **reaction values quoted in Section 5.3** at 298.15 K,

```{math}
\Delta_r H^\circ \approx -92\ \mathrm{kJ\,mol^{-1}},\qquad
\Delta_r S^\circ \approx -198\ \mathrm{J\,mol^{-1}\,K^{-1}},
```

gives

```{math}
\Delta_r G^\circ
\approx -92\ \mathrm{kJ\,mol^{-1}}
- (298.15\ \mathrm{K})\left(-0.198\ \mathrm{kJ\,mol^{-1}\,K^{-1}}\right)
\approx -33\ \mathrm{kJ\,mol^{-1}}.
```

Then

```{math}
K_p(298.15\ \mathrm{K})
= \exp\!\left(-\frac{\Delta_r G^\circ}{RT}\right)
\approx \exp\!\left(\frac{33{,}000}{(8.314)(298.15)}\right)
\approx 10^{5}\text{–}10^{6}.
```

**Interpretation:** $K_p \gg 1$ at room temperature, so the equilibrium strongly favors $\mathrm{NH_3}$ *at the standard state*.
(Section 5.3 explains why temperature and pressure “knobs” still matter for industrial operation.)

---

#### Common pitfalls (and quick fixes)

- **Wrong phase:** $\mathrm{H_2O(\ell)}$ and $\mathrm{H_2O(g)}$ (or $\mathrm{NH_3(\ell)}$ vs $\mathrm{NH_3(g)}$) have very different $S^\circ$ and $\Delta_f H^\circ$.
- **Mixed standard states:** make sure all species use the same $P^\circ$ convention and the same temperature.
- **Unit mismatches:** always reconcile kJ vs J and “per mole of reaction” vs “per mole of species.”
- **Sign mistakes:** products minus reactants (or use $\nu_i$ with the sign convention consistently).

## Worked Example

### NO$_2$ dimerization: solving for $\xi_{eq}$ at $P=P^\circ$

Reaction:

```{math}
2\,\mathrm{NO_2(g)} \rightleftharpoons \mathrm{N_2O_4(g)}.
```

Start with $n_{\mathrm{NO_2},0}=2$, $n_{\mathrm{N_2O_4},0}=0$. At extent $\xi$:

```{math}
n_{\mathrm{NO_2}}=2-2\xi,\qquad n_{\mathrm{N_2O_4}}=\xi,\qquad n_{\text{tot}}=2-\xi.
```

Mole fractions:

```{math}
y_{\mathrm{NO_2}}=\frac{2-2\xi}{2-\xi},\qquad y_{\mathrm{N_2O_4}}=\frac{\xi}{2-\xi}.
```

At total pressure $P=P^\circ$, partial pressures are $P_i=y_iP^\circ$, and the ideal-gas equilibrium constant is

```{math}
K_p=\frac{(P_{\mathrm{N_2O_4}}/P^\circ)}{(P_{\mathrm{NO_2}}/P^\circ)^2}
=\frac{y_{\mathrm{N_2O_4}}}{y_{\mathrm{NO_2}}^2}.
```

Therefore,

```{math}
K_p=\frac{\xi/(2-\xi)}{\left[(2-2\xi)/(2-\xi)\right]^2}
=\frac{\xi(2-\xi)}{(2-2\xi)^2}.
```

Given $K_p=6.74$, solve $\xi(2-\xi)=6.74(2-2\xi)^2$, yielding $\xi_{eq}\approx 0.81$.

**Result.** $n_{\mathrm{NO_2},eq}\approx 2-2(0.81)=0.38$ and $n_{\mathrm{N_2O_4},eq}\approx 0.81$ at $P=P^\circ$.

## Concept Checks

1. Why is $\Delta_r G=0$ the correct equilibrium condition at constant $T,P$ but not necessarily at constant $T,V$?
2. Why must $Q$ and $K$ be dimensionless? What role does $P^\circ$ play?
3. How does comparing $Q$ to $K$ predict the direction of spontaneous change?
4. What changes in the derivation if the mixture is non-ideal and activities replace partial-pressure ratios?

## Key Takeaways

- Extent of reaction $\xi$ provides a compact way to track composition changes via stoichiometry.
- At constant $T,P$, equilibrium corresponds to minimizing $G$, giving $\Delta_r G=0$.
- For ideal gases, $\Delta_r G=\Delta_r G^\circ+RT\ln Q$ and $\Delta_r G^\circ=-RT\ln K$.
- The sign of $\Delta_r G$ (equivalently $Q$ vs $K$) determines spontaneous direction.
