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

## Overview

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

## 7.1.1 Extent of reaction $\xi$

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

## 7.1.2 Gibbs free energy and the equilibrium condition

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

### Spontaneity and equilibrium (constant $T,P$)

- $\Delta_r G < 0$: reaction proceeds **forward** (toward products)
- $\Delta_r G = 0$: **equilibrium**
- $\Delta_r G > 0$: reaction proceeds **backward** (toward reactants)

At equilibrium,

```{math}
\boxed{\Delta_r G = 0.}
```

---

## 7.1.3 Reaction quotient $Q$ and the equilibrium constant $K$

```{admonition} Notational warning
:class: warning

In this section, $Q$ (or $Q_p$) means the **reaction quotient**, not a partition function.
```

To connect $\Delta_r G$ to measurable composition variables, we need an expression for the chemical potentials.

### Ideal-gas chemical potential

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

### Derivation of $\Delta_r G = \Delta_r G^{\circ} + RT\ln Q$

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

### Equilibrium: $Q = K$

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

### Using $Q$ vs $K$ to predict direction

Since $\Delta_r G = RT\ln(Q/K)$:

- If $Q<K$, then $\Delta_r G<0$ and the reaction proceeds **toward products**.
- If $Q>K$, then $\Delta_r G>0$ and the reaction proceeds **toward reactants**.
- If $Q=K$, the system is at **equilibrium**.

---

## 7.1.4 Worked example: gas-phase dimerization of $\mathrm{NO_2}$

Consider the equilibrium

```{math}
2\,\mathrm{NO_2(g)} \rightleftharpoons \mathrm{N_2O_4(g)}.
```

($\mathrm{NO_2}$ is brown; $\mathrm{N_2O_4}$ is colorless.)

### Thermodynamics and $K_p$ at 298.15 K

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

### Equilibrium composition at a specified total pressure

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

---

## 7.1.5 Temperature dependence: the van't Hoff equation

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

### Example trend for $2\,\mathrm{NO_2}\rightleftharpoons\mathrm{N_2O_4}$

From the lecture notes, for $P=P^{\circ}$:

| $T$ (K) | $K_p$ | $\xi_{eq}$ (starting from 2 mol NO$_2$) |
| ---: | ---: | ---: |
| 250 | 205.18 | 0.97 |
| 298.15 | 6.74 | 0.81 |
| 350 | 0.17 | 0.23 |

Because $\Delta_r H^{\circ}<0$ for dimerization, increasing temperature drives $K_p$ **down** and shifts equilibrium toward $\mathrm{NO_2}$ (less $\mathrm{N_2O_4}$).

---

## 7.1.6 Visualizing equilibrium as a minimum of $G(\xi)$

If we write the Gibbs free energy of the reacting mixture as a function of extent, $G(\xi)$, then:

- the **slope** is $\left(\partial G/\partial \xi\right)_{T,P}=\Delta_r G$
- the **equilibrium extent** $\xi_{eq}$ is located where the slope is zero
- equilibrium corresponds to a **minimum** in $G(\xi)$ at fixed $T,P$

For the $\mathrm{NO_2/N_2O_4}$ system, the notes sketch $G$ vs. $\xi$ with a clear minimum at $\xi_{eq}$: the system "rolls downhill" in $G$ until it reaches that minimum.

---

## Key takeaways

1. Use $\xi$ to relate **composition** to reaction progress.
2. At constant $T,P$, equilibrium occurs when $\Delta_r G = 0$.
3. For ideal gases, $\Delta_r G = \Delta_r G^{\circ} + RT\ln Q_p$.
4. At equilibrium, $Q_p=K_p$ and $\Delta_r G^{\circ}=-RT\ln K_p$.
5. Comparing $Q$ to $K$ predicts the direction of spontaneous change.
6. The van't Hoff equation links $K$ to temperature through $\Delta_r H^{\circ}$.
