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


# 6.2. Clausius-Clapeyron

[Course-wide Conventions & Notation](../notation.md)

## Overview and Learning Objectives

The Clausius–Clapeyron equation is a practical approximation of the Clapeyron equation for liquid–vapor equilibrium when the vapor is ideal and dominates the molar volume change. It predicts how saturation vapor pressure changes with temperature and provides a convenient way to estimate $\Delta H_{vap}$ from vapor-pressure data.

In the previous section we introduced **phase equilibrium** using Gibbs free energy / chemical potential equality.
In this lecture we push that idea further:

- **Chemical potential** $\mu_i$ is the *intensive* variable that controls **matter exchange** in an **open system** (a system where $n_i$ can change).
- The phase-equilibrium condition $\mu_\alpha = \mu_\beta$ (or equivalently $g_\alpha=g_\beta$ for a pure substance) lets us derive the **Clapeyron / Clausius–Clapeyron equations**.
- Those equations tell us the **slope of phase boundaries** and give a practical way to estimate **enthalpies of phase transition** (especially $\Delta H_{vap}$) from vapor-pressure data.

(See the “open system” sketch and the Clausius–Clapeyron derivation on pages 1–2 of the lecture notes, and the water example on page 3.)

---

Learning objectives:

- State the phase-equilibrium condition $\mu_\alpha=\mu_\beta$ and derive the Clapeyron slope $dP/dT=\Delta H_{tr}/(T\Delta v)$.
- Apply the approximations $v_g\gg v_l$ and $v_g=RT/P$ to obtain $d\ln P/dT=\Delta H_{vap}/(RT^2)$.
- Integrate the Clausius–Clapeyron equation under constant $\Delta H_{vap}$ to relate two $(T,P)$ points.
- Use vapor-pressure data to estimate $\Delta H_{vap}$ and interpret trends (e.g., hydrogen bonding).

## Core Ideas and Derivations

### 6.1.1 Chemical potential and open systems

#### What makes a system “open”?

An **open system** can exchange **energy** and **matter** with its surroundings across a boundary.

- Energy exchange: heat (often visualized via temperature $T$)
- Matter exchange: particles / moles (controlled by chemical potential $\mu$)

The notes emphasize: **$\mu$ controls matter exchange across a boundary.**

#### Fundamental differentials including composition

For a multicomponent system (sum over *components* $i$):

```{math}
dU = T\,dS - P\,dV + \sum_i \mu_i\,dN_i
```

and

```{math}
dG = -S\,dT + V\,dP + \sum_i \mu_i\,dN_i .
```

Here the chemical potential is defined by the partial derivative

```{math}
\boxed{\mu_i \equiv \left(\frac{\partial G}{\partial N_i}\right)_{T,P,N_{j\neq i}} }.
```

> “Components” means the *minimum* number of chemical species needed to specify the composition of all phases in the system.

#### Microscopic interpretation (intuition)

Because $\mu_i$ is “how much $G$ changes when you add a little of $i$,” you can think of it as:

- the **driving force** for particles to move between regions/phases, and
- an “energy cost” (in the Gibbs sense) of adding particles.

Matter tends to flow from **higher $\mu$** to **lower $\mu$** until equilibrium is reached.

---

### 6.1.2 Phase equilibrium as “no driving force for transfer”

Consider two phases, $\alpha$ and $\beta$, separated by a **matter-permeable boundary**. Let an amount $dn$ be transferred from $\alpha$ to $\beta$.

Write the Gibbs change due to transfer:

```{math}
dG = \mu_\alpha\,dn_\alpha + \mu_\beta\,dn_\beta.
```

If the *total* amount of matter in the combined system is conserved, then

```{math}
dn_\alpha + dn_\beta = 0 \quad \Rightarrow \quad dn_\alpha = -dn_\beta.
```

Substitute:

```{math}
dG = -\mu_\alpha\,dn_\beta + \mu_\beta\,dn_\beta
= (\mu_\beta-\mu_\alpha)\,dn_\beta.
```

At **phase equilibrium**, transfer should not lower $G$ in either direction, so the only way for this to hold for arbitrary small transfers is:

```{math}
\boxed{\mu_\alpha = \mu_\beta.}
```

This is the key equilibrium condition we will now differentiate to get the slope of coexistence curves.

---

### 6.1.3 Clapeyron equation (general slope of a coexistence curve)

Along a phase boundary, $\mu_\alpha(T,P)=\mu_\beta(T,P)$.  
Moving an infinitesimal amount *along the boundary* keeps them equal, so:

```{math}
d\mu_\alpha = d\mu_\beta.
```

For a one-component phase, the differential of chemical potential can be written in molar form:

```{math}
d\mu = -\bar{s}\,dT + \bar{v}\,dP,
```

where $\bar{s}=S/n$ and $\bar{v}=V/n$ are **molar entropy** and **molar volume**.

Apply this to both phases:

```{math}
d\mu_\alpha = -\bar{s}_\alpha\,dT + \bar{v}_\alpha\,dP,
```

```{math}
d\mu_\beta = -\bar{s}_\beta\,dT + \bar{v}_\beta\,dP.
```

Set them equal and rearrange:

```{math}
-\bar{s}_\alpha\,dT + \bar{v}_\alpha\,dP
=
-\bar{s}_\beta\,dT + \bar{v}_\beta\,dP
```

```{math}
\Rightarrow \boxed{\frac{dP}{dT} = \frac{\bar{s}_\beta-\bar{s}_\alpha}{\bar{v}_\beta-\bar{v}_\alpha}}
= \boxed{\frac{\Delta \bar{s}}{\Delta \bar{v}}}.
```

This is the **Clapeyron equation** (the notes label this step as the “Clausius equation”).

Using $\Delta \bar{h} = T\,\Delta\bar{s}$ at a phase transition (latent heat relation), we can also write:

```{math}
\boxed{\frac{dP}{dT} = \frac{\Delta \bar{h}_{tr}}{T\,\Delta \bar{v}} }.
```

---

### 6.1.4 Clausius–Clapeyron (liquid–vapor, ideal-gas approximation)

For **vaporization** ($l \rightleftharpoons g$), usually

- $\bar{v}_g \gg \bar{v}_l$ so $\Delta \bar{v} \approx \bar{v}_g$,
- and if the vapor behaves ideally, $\bar{v}_g = RT/P$.

Start from Clapeyron:

```{math}
\frac{dP}{dT} = \frac{\Delta\bar{h}_{vap}}{T\,\Delta\bar{v}}
\approx
\frac{\Delta\bar{h}_{vap}}{T\,(RT/P)}.
```

So

```{math}
\frac{dP}{dT} \approx \frac{P\,\Delta\bar{h}_{vap}}{R\,T^2}
\quad\Rightarrow\quad
\boxed{\frac{d\ln P}{dT} = \frac{\Delta\bar{h}_{vap}}{R\,T^2}}.
```

This is the **Clausius–Clapeyron equation** in differential form.

#### Integrated form (constant $\Delta H_{vap}$ approximation)

If $\Delta\bar{h}_{vap}$ is approximately constant over the temperature interval, integrate:

```{math}
\boxed{\ln\left(\frac{P_2}{P_1}\right)
=
-\frac{\Delta\bar{h}_{vap}}{R}\left(\frac{1}{T_2}-\frac{1}{T_1}\right)}
=
\frac{\Delta\bar{h}_{vap}}{R}\left(\frac{1}{T_1}-\frac{1}{T_2}\right).
```

A common algebraic rearrangement (used on the board) is

```{math}
\ln\left(\frac{P_2}{P_1}\right)
=
\frac{\Delta\bar{h}_{vap}}{R}
\left(\frac{T_2-T_1}{T_1T_2}\right).
```

---

### 6.1.5 Worked example: estimating $\Delta H_{vap}$ of water from vapor pressure

The board example uses two points on the vapor-pressure curve of water:

- At $T_1 = 363.2\,\text{K}$, $P_1 \approx 5.2\times 10^2\,\text{torr}$ (measured vapor pressure).
- At the **normal boiling point**, $T_2 = 373.2\,\text{K}$ and $P_2 = 760\,\text{torr}$ (1 atm).

Use the integrated Clausius–Clapeyron equation:

```{math}
\ln\left(\frac{P_2}{P_1}\right)
=
\frac{\Delta H_{vap}}{R}\left(\frac{1}{T_1}-\frac{1}{T_2}\right),
```

so

```{math}
\boxed{
\Delta H_{vap}
=
R\,\frac{\ln(P_2/P_1)}{(1/T_1 - 1/T_2)}
}.
```

Plugging in the board values gives an estimate on the order of

```{math}
\Delta H_{vap}(\text{H}_2\text{O}) \approx 4.1\times 10^1\ \text{kJ mol}^{-1}
```

(the board computation reports **$\approx 40.8\ \text{kJ mol}^{-1}$**).

#### Physical interpretation (hydrogen bonding)

The notes add a molecular-level interpretation:

- Vaporizing water requires disrupting hydrogen-bonding interactions present in the liquid.
- If one H-bond is roughly $\sim 20\ \text{kJ mol}^{-1}$, then

  ```{math}
  \Delta H_{vap} \sim 40\ \text{kJ mol}^{-1}
  ```

  corresponds (very roughly) to “breaking” about **two** H-bond interactions per molecule on average to go from $\text{H}_2\text{O}(l)$ to $\text{H}_2\text{O}(g)$.

---

## Worked Example

### Estimating $\Delta H_{vap}$ from two vapor-pressure points

Using the integrated Clausius–Clapeyron form

```{math}
\ln\left(\frac{P_2}{P_1}\right)=\frac{\Delta H_{vap}}{R}\left(\frac{1}{T_1}-\frac{1}{T_2}\right),
```

estimate $\Delta H_{vap}$ for water given:

- $T_1=363.2\ \mathrm{K}$, $P_1=5.2\times10^2\ \mathrm{torr}$
- $T_2=373.2\ \mathrm{K}$, $P_2=760\ \mathrm{torr}$

Solve for $\Delta H_{vap}$:

```{math}
\Delta H_{vap}=R\,\frac{\ln(P_2/P_1)}{(1/T_1-1/T_2)}.
```

1. **Compute the logarithm**

   ```{math}
   \ln\left(\frac{P_2}{P_1}\right)=\ln\left(\frac{760}{520}\right)=0.379.
   ```

2. **Compute the temperature factor**

   ```{math}
   \left(\frac{1}{T_1}-\frac{1}{T_2}\right)=\frac{1}{363.2}-\frac{1}{373.2}
   =7.38\times10^{-5}\ \mathrm{K^{-1}}.
   ```

3. **Compute $\Delta H_{vap}$**

   ```{math}
   \Delta H_{vap}=(8.314)\frac{0.379}{7.38\times10^{-5}}
   =4.27\times10^{4}\ \mathrm{J\,mol^{-1}}=42.7\ \mathrm{kJ\,mol^{-1}}.
   ```

**Result.** The two-point estimate gives $\Delta H_{vap}\sim 4\times10^{1}\ \mathrm{kJ\,mol^{-1}}$, consistent with tabulated values (and sensitive to the chosen data points/approximations).

## Concept Checks

1. Why does $d\ln P/dT$ depend on $1/T^2$ rather than simply on $1/T$?
2. What physical effect would cause a Clausius–Clapeyron plot to deviate from a straight line?
3. Why is it usually safe to approximate $\Delta v\approx v_g$ for vaporization but not for melting of many solids?
4. What is the sign of $dP/dT$ for vaporization, and what does it imply about the phase boundary slope?

## Key Takeaways

- Clapeyron relates coexistence-curve slope to latent heat and volume change: $dP/dT=\Delta H_{tr}/(T\Delta v)$.
- For vaporization with ideal vapor, Clausius–Clapeyron gives $d\ln P/dT=\Delta H_{vap}/(RT^2)$.
- The integrated form extracts $\Delta H_{vap}$ from two vapor-pressure points (or from a linear fit of $\ln P$ vs $1/T$).
- Deviations from linearity signal changing $\Delta H_{vap}(T)$ and/or non-ideal behavior.
