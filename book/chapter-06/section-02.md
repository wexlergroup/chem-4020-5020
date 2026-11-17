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
  enumerator: 6.2.%s
---


# 6.2. Clapeyron and Clausius–Clapeyron

See also: [Course-wide Conventions & Notation](../notation.md)

## Overview and Learning Objectives

In [Section 6.1](section-01.md) we showed that two phases of a pure substance are in equilibrium at a given $(T,P)$ if and only if their chemical potentials are equal:

```{math}
:label: eq:phase-equilibrium-recap
\mu_\alpha(T,P) = \mu_\beta(T,P).
```

Because this is a single scalar constraint on two variables, it picks out a *curve* in the $(T,P)$ plane — the coexistence curve we saw in the phase diagram. What it does *not* immediately give us is the **slope** $dP/dT$ of that curve, which is what we would need to predict, for instance, how a boiling point shifts with pressure or how a vapor pressure changes with temperature.

In this section we differentiate the equilibrium condition along the coexistence curve to obtain the **Clapeyron equation**, a general expression for $dP/dT$ that applies to any first-order phase transition. We then specialize it to the liquid–vapor case under two further approximations — the gas molar volume dominates, and the vapor is ideal — to obtain the **Clausius–Clapeyron equation**, which is the working tool for analyzing vapor-pressure data.

Learning objectives:

- Differentiate $\mu_\alpha(T,P) = \mu_\beta(T,P)$ along a coexistence curve to derive the Clapeyron equation $dP/dT = \Delta\bar s/\Delta\bar v = \Delta\bar h_{\mathrm{trs}}/(T\,\Delta\bar v)$.
- Apply the Clapeyron equation to ice → water and explain physically why the solid–liquid coexistence curve for water has a *negative* slope, in contrast to most substances.
- Apply the approximations $\bar v_{\mathrm{g}} \gg \bar v_{\mathrm{l}}$ and $\bar v_{\mathrm{g}} = RT/P$ to obtain the differential Clausius–Clapeyron equation $d\ln P / dT = \Delta\bar h_{\mathrm{vap}}/(RT^2)$.
- Integrate the Clausius–Clapeyron equation under the constant-$\Delta\bar h_{\mathrm{vap}}$ approximation and use a two-point dataset to estimate $\Delta\bar h_{\mathrm{vap}}$ from vapor-pressure measurements.
- Diagnose the failure of the constant-$\Delta\bar h_{\mathrm{vap}}$ assumption from residuals of a $\ln P$-vs-$1/T$ fit, and use a local-slope estimate to follow $\Delta\bar h_{\mathrm{vap}}(T)$ where the vapor is nearly ideal.

## Core Ideas and Derivations

(section:clapeyron)=
### The Clapeyron Equation

Pick a point $(T,P)$ on a coexistence curve between phases $\alpha$ and $\beta$, and consider an infinitesimal displacement $(dT,dP)$ that stays *on* the curve. By construction, the equilibrium condition continues to hold along the curve, so $\mu_\alpha$ and $\mu_\beta$ change by equal amounts:

```{math}
d\mu_\alpha = d\mu_\beta.
```

For a one-component, single-phase system, the Gibbs differential per mole reads

```{math}
:label: eq:dmu
d\mu = -\bar s\,dT + \bar v\,dP,
```

which is just the molar form of $dG = -S\,dT + V\,dP$ specialized to a pure substance (so that $\mu = \bar g$, as established under [Open Systems and the Chemical Potential](#section:chemical-potential) in [Section 6.1](section-01.md)). Apply Eq. {eq}`eq:dmu` to each phase along the coexistence curve:

```{math}
-\bar s_\alpha\,dT + \bar v_\alpha\,dP = -\bar s_\beta\,dT + \bar v_\beta\,dP.
```

Collecting $dT$ and $dP$ terms,

```{math}
(\bar v_\beta - \bar v_\alpha)\,dP = (\bar s_\beta - \bar s_\alpha)\,dT,
```

and using the entropy–latent-heat relation $\Delta\bar h_{\mathrm{trs}} = T\,\Delta\bar s_{\mathrm{trs}}$ derived under [Latent Heat and Entropy Jump at a First-Order Transition](#section:latent-heat) in [Section 6.1](section-01.md), we obtain the **Clapeyron equation**:

```{math}
:label: eq:clapeyron
\boxed{\frac{dP}{dT} = \frac{\Delta\bar s_{\mathrm{trs}}}{\Delta\bar v_{\mathrm{trs}}} = \frac{\Delta\bar h_{\mathrm{trs}}}{T\,\Delta\bar v_{\mathrm{trs}}}}.
```

Three things to notice about Eq. {eq}`eq:clapeyron`:

1. **It is exact.** No approximation has been made beyond the assumption that the transition is first-order (so $\Delta\bar v$ and $\Delta\bar h$ are well-defined finite jumps) and that both phases are pure.
2. **The sign of $dP/dT$ is set by the sign of $\Delta\bar v$,** since $\Delta\bar h$ and $T$ are positive for any "normal" transition (the higher-temperature phase is the higher-enthalpy one). Most substances have $\bar v_{\mathrm{l}} > \bar v_{\mathrm{s}}$ at the melting line and $\bar v_{\mathrm{g}} \gg \bar v_{\mathrm{l}}$ at the boiling line, so both slopes are positive — an increase in pressure raises both the melting and the boiling temperature.
3. **The magnitude of $dP/dT$ is set by $\Delta\bar v$.** Solid–liquid transitions have very small $\Delta\bar v$ and consequently very steep coexistence curves; liquid–vapor transitions have very large $\Delta\bar v$ and consequently shallow coexistence curves. This is the reason why melting points are nearly insensitive to ordinary pressure changes while boiling points are very sensitive.

### Worked Example: Melting of Ice and the Negative-Slope Coexistence Curve

Water is the textbook exception to "$\bar v_{\mathrm{l}} > \bar v_{\mathrm{s}}$": ice floats on liquid water because $\bar v_{\mathrm{s}} > \bar v_{\mathrm{l}}$ near the normal melting point. By Eq. {eq}`eq:clapeyron`, this forces the *slope* of the ice–water coexistence curve to be **negative**.

**Inputs (at $T_{\mathrm{fus}} = 273.15$ K, $P = 1$ bar):**

- $\Delta\bar h_{\mathrm{fus}} = 6.010\ \mathrm{kJ\,mol^{-1}}$
- $\rho_{\mathrm{s}} = 0.917\ \mathrm{g\,cm^{-3}}$, $\rho_{\mathrm{l}} = 0.9998\ \mathrm{g\,cm^{-3}}$
- $M(\mathrm{H_2O}) = 18.015\ \mathrm{g\,mol^{-1}}$

**Molar volumes** of the two phases:

```{math}
\bar v_{\mathrm{s}} = \frac{M}{\rho_{\mathrm{s}}} = \frac{18.015}{0.917}\ \mathrm{cm^3\,mol^{-1}} = 19.65\ \mathrm{cm^3\,mol^{-1}},
```

```{math}
\bar v_{\mathrm{l}} = \frac{M}{\rho_{\mathrm{l}}} = \frac{18.015}{0.9998}\ \mathrm{cm^3\,mol^{-1}} = 18.02\ \mathrm{cm^3\,mol^{-1}}.
```

The molar volume **decreases** on melting:

```{math}
\Delta\bar v_{\mathrm{fus}} = \bar v_{\mathrm{l}} - \bar v_{\mathrm{s}} = -1.63\ \mathrm{cm^3\,mol^{-1}} = -1.63\times 10^{-6}\ \mathrm{m^3\,mol^{-1}}.
```

**Clapeyron slope.** Insert into Eq. {eq}`eq:clapeyron`:

```{math}
\begin{aligned}
\frac{dP}{dT} &= \frac{\Delta\bar h_{\mathrm{fus}}}{T_{\mathrm{fus}}\,\Delta\bar v_{\mathrm{fus}}}
= \frac{6010\ \mathrm{J\,mol^{-1}}}{(273.15\ \mathrm{K})(-1.63\times 10^{-6}\ \mathrm{m^3\,mol^{-1}})} \\
&\approx -1.35\times 10^{7}\ \mathrm{Pa\,K^{-1}}
\approx -135\ \mathrm{bar\,K^{-1}}.
\end{aligned}
```

**Physical interpretation.** Inverting the slope gives $dT/dP \approx -7.4\ \mathrm{mK\,bar^{-1}}$: increasing the pressure on ice by one bar lowers its melting point by about seven thousandths of a kelvin. Even at 100 bar — well above any pressure encountered in everyday situations — the melting point depression is only about 0.74 K. This is enough to let the ice–water phase boundary bend perceptibly to the *left* on the phase diagram, but it is far too small to account for the popular "pressure melting" explanation of ice skating, where contact pressures of a few hundred bar would depress $T_{\mathrm{fus}}$ by at most a degree or two (skating works mainly because of friction-generated heating and a thin pre-existing surface liquid layer, not pressure melting).

(section:clausius-clapeyron)=
### The Clausius–Clapeyron Equation (Liquid–Vapor with Ideal Vapor)

For the liquid–vapor coexistence curve, the molar-volume change is overwhelmingly dominated by the gas:

```{math}
\bar v_{\mathrm{g}} \gg \bar v_{\mathrm{l}} \quad\Longrightarrow\quad \Delta\bar v_{\mathrm{vap}} \approx \bar v_{\mathrm{g}}.
```

(For water at the normal boiling point, $\bar v_{\mathrm{g}} \approx 3.0\times 10^{-2}\ \mathrm{m^3\,mol^{-1}}$ versus $\bar v_{\mathrm{l}} \approx 1.9\times 10^{-5}\ \mathrm{m^3\,mol^{-1}}$ — a ratio of more than $10^3$, so dropping $\bar v_{\mathrm{l}}$ introduces an error of less than 0.1%.)

If we further approximate the vapor as an ideal gas at the saturation pressure,

```{math}
\bar v_{\mathrm{g}} = \frac{RT}{P_{\mathrm{sat}}},
```

then the Clapeyron equation {eq}`eq:clapeyron` becomes

```{math}
\frac{dP_{\mathrm{sat}}}{dT} \approx \frac{\Delta\bar h_{\mathrm{vap}}}{T\,\bar v_{\mathrm{g}}} = \frac{\Delta\bar h_{\mathrm{vap}}}{T\,(RT/P_{\mathrm{sat}})} = \frac{P_{\mathrm{sat}}\,\Delta\bar h_{\mathrm{vap}}}{RT^2}.
```

Dividing both sides by $P_{\mathrm{sat}}$ and recognizing $d\ln P_{\mathrm{sat}} = dP_{\mathrm{sat}}/P_{\mathrm{sat}}$ gives the **Clausius–Clapeyron equation** in differential form:

```{math}
:label: eq:CC-differential
\boxed{\frac{d\ln P_{\mathrm{sat}}}{dT} = \frac{\Delta\bar h_{\mathrm{vap}}}{RT^2}}.
```

The quadratic $T^{-2}$ dependence, rather than the $T^{-1}$ that one might naively expect, has its origin in the derivation: two of its steps each bring in a factor of $T$. Concept Check 1 asks for them.

Without the two approximations, the Clapeyron equation {eq}`eq:clapeyron` for the liquid–vapor curve can be written in the same form, and it remains exact:

```{math}
:label: eq:CC-exact
\frac{d\ln P_{\mathrm{sat}}}{dT} = \frac{\Delta\bar h_{\mathrm{vap}}}{RT^2\,\Delta Z},
\qquad
\Delta Z \equiv \frac{P_{\mathrm{sat}}\,\Delta\bar v_{\mathrm{vap}}}{RT}.
```

Here $\Delta Z$ is the difference between the compressibility factors $Z = P\bar v/RT$ ([Section 1.4](../chapter-01/section-04.md)) of the vapor and the liquid. The Clausius–Clapeyron equation sets $\Delta Z = 1$. For water, $\Delta Z = 0.999$ at 20 °C, $0.984$ at 100 °C, and $0.924$ at 180 °C, so the approximation is excellent at low vapor pressure and weakens as the vapor becomes denser.

#### Integrated Form (Constant $\Delta \bar h_{\mathrm{vap}}$ Approximation)

If $\Delta\bar h_{\mathrm{vap}}$ varies only weakly with $T$ over the range of interest, it can be pulled outside the integral:

```{math}
\int_{T_1}^{T_2} d\ln P_{\mathrm{sat}}
= \frac{\Delta\bar h_{\mathrm{vap}}}{R}\int_{T_1}^{T_2}\frac{dT}{T^2},
```

so that

```{math}
:label: eq:CC-integrated
\boxed{\ln\!\left(\frac{P_2}{P_1}\right)
= -\frac{\Delta\bar h_{\mathrm{vap}}}{R}\left(\frac{1}{T_2}-\frac{1}{T_1}\right)
= \frac{\Delta\bar h_{\mathrm{vap}}}{R}\left(\frac{1}{T_1}-\frac{1}{T_2}\right)}.
```

Two algebraically equivalent rearrangements of Eq. {eq}`eq:CC-integrated` are useful for different purposes:

- For **two known $(T,P)$ points**, solve for the unknown $\Delta\bar h_{\mathrm{vap}}$:

  ```{math}
  :label: eq:CC-two-point
  \Delta\bar h_{\mathrm{vap}} = R\,\frac{\ln(P_2/P_1)}{(1/T_1 - 1/T_2)}.
  ```

- For **fitting a dataset**, recognize that treating one of the two points of Eq. {eq}`eq:CC-integrated` as a fixed reference gives

  ```{math}
  \ln P_{\mathrm{sat}}(T) = -\frac{\Delta\bar h_{\mathrm{vap}}}{R}\,\frac{1}{T} + \text{const.}
  ```

  A plot of $\ln P_{\mathrm{sat}}$ versus $1/T$ should be a straight line with slope $-\Delta\bar h_{\mathrm{vap}}/R$ — provided the constant-$\Delta\bar h_{\mathrm{vap}}$ approximation actually holds.

### Worked Example: Estimating $\Delta\bar h_{\mathrm{vap}}$ of Water from Two Vapor-Pressure Points

Use the integrated Clausius–Clapeyron form {eq}`eq:CC-two-point` to estimate $\Delta\bar h_{\mathrm{vap}}$ for water given the following two points on its saturation curve:

- $T_1 = 363.15\ \mathrm{K}$ (90 °C), $P_1 = 70.18\ \mathrm{kPa}$
- $T_2 = 373.15\ \mathrm{K}$ (100 °C), $P_2 = 101.42\ \mathrm{kPa}$

Both points are taken from the saturation tables for water of the NIST Chemistry WebBook.

**Step 1.** Compute the logarithm of the pressure ratio:

```{math}
\ln\!\left(\frac{P_2}{P_1}\right) = \ln\!\left(\frac{101.42}{70.18}\right) = 0.3682.
```

**Step 2.** Compute the inverse-temperature factor:

```{math}
\frac{1}{T_1} - \frac{1}{T_2}
= \frac{1}{363.15\ \mathrm{K}} - \frac{1}{373.15\ \mathrm{K}}
= 7.380\times 10^{-5}\ \mathrm{K^{-1}}.
```

**Step 3.** Combine using Eq. {eq}`eq:CC-two-point`:

```{math}
\begin{aligned}
\Delta\bar h_{\mathrm{vap}}
&= (8.314\ \mathrm{J\,mol^{-1}\,K^{-1}})\,\frac{0.3682}{7.380\times 10^{-5}\ \mathrm{K^{-1}}} \\
&\approx 4.15\times 10^{4}\ \mathrm{J\,mol^{-1}}
= 41.5\ \mathrm{kJ\,mol^{-1}}.
\end{aligned}
```

**Discussion.** The two-point estimate of $41.5\ \mathrm{kJ\,mol^{-1}}$ sits 2 % above the calorimetric value $\Delta\bar h_{\mathrm{vap}}(\mathrm{H_2O}, 100\ ^\circ\mathrm{C}) = 40.65\ \mathrm{kJ\,mol^{-1}}$. The discrepancy is *not* just measurement noise, and it has two sources. First, $\Delta\bar h_{\mathrm{vap}}$ is itself slightly temperature-dependent (it decreases with rising $T$ and vanishes at the critical point), so an estimate built from data points spanning 90–100 °C reflects an *average* over that interval rather than the value at the upper endpoint. That average is $40.9\ \mathrm{kJ\,mol^{-1}}$ ($\Delta\bar h_{\mathrm{vap}}$ is $41.12\ \mathrm{kJ\,mol^{-1}}$ at 90 °C). Second, the remaining 1.5 % comes from the approximations of the Clausius–Clapeyron equation itself. At 95 °C, $\Delta Z = 0.986$ in Eq. {eq}`eq:CC-exact`, so the slope of $\ln P_{\mathrm{sat}}$ is larger than $\Delta\bar h_{\mathrm{vap}}/RT^2$ by that factor: $40.9/0.986 = 41.5\ \mathrm{kJ\,mol^{-1}}$. We will return to both effects in the mini-lab below.

#### A Molecular Reading of the Result

Why is $\Delta\bar h_{\mathrm{vap}}$ for water of order 40 kJ mol$^{-1}$? The molecular interpretation gives a satisfying back-of-the-envelope answer that ties the macroscopic measurement to the intermolecular forces of [Section 1.4](../chapter-01/section-04.md). In liquid water, each molecule takes part in nearly four hydrogen bonds. Each hydrogen bond is shared between two molecules, so there are about two per molecule, and a single hydrogen bond is worth $\sim 20\ \mathrm{kJ\,mol^{-1}}$. Vaporizing one mole of water to an effectively non-interacting gas therefore costs roughly

```{math}
\Delta\bar h_{\mathrm{vap}}(\mathrm{H_2O}) \sim 2\times 20\ \mathrm{kJ\,mol^{-1}} = 40\ \mathrm{kJ\,mol^{-1}},
```

in striking agreement with experiment. The estimate is rough — it omits dispersion interactions and treats each H-bond as an isolated bond — but it captures the right order of magnitude and identifies the dominant physics. The same picture bears on the large value of $\Delta\bar s_{\mathrm{vap}}$ for water and its deviation from Trouton's rule ([Trouton's Rule and Water as an Exception](#section:trouton-rule) in [Section 6.1](section-01.md)): on vaporization, water takes up more entropy as well as more enthalpy than an "average" liquid, and Concept Check 4 of Section 6.1 asks what that says about the liquid.

## Mini-Lab: Vapor Pressure → $\Delta\bar h_{\mathrm{vap}}$, Linearity Limits, and $\Delta\bar h_{\mathrm{vap}}(T)$

### Goals

1. Use a **Clausius–Clapeyron plot** to estimate an *average* $\Delta\bar h_{\mathrm{vap}}$ from a vapor-pressure dataset.
2. Diagnose **when the straight-line model breaks down** by examining residuals.
3. Estimate a **temperature-dependent** $\Delta\bar h_{\mathrm{vap}}(T)$ from local slopes of the same data, and see where the estimate departs from the reference values.

:::{important}
This whole mini-lab assumes the **liquid–vapor approximation** developed under [The Clausius–Clapeyron Equation (Liquid–Vapor with Ideal Vapor)](#section:clausius-clapeyron): the gas is ideal and $\bar v_{\mathrm{g}} \gg \bar v_{\mathrm{l}}$. Deviations from a perfect straight line in $\ln P$ versus $1/T$ can come from any of (i) real-gas behavior at high $P$, (ii) non-negligible liquid molar volume, and (iii) the temperature dependence of $\Delta\bar h_{\mathrm{vap}}$ itself.
:::

### Dataset: Saturation Vapor Pressure of Water

The values are the reference saturation pressures of water (the IAPWS-95 formulation, as tabulated in the NIST Chemistry WebBook). Use them as if they were experimental measurements; they are smooth and consistent, so they are good for analysis practice rather than for assessing experimental uncertainty.

```{list-table}
:header-rows: 1
:name: tab:water-vap-pressure-mini-lab

* - $T$ (°C)
  - $P_\mathrm{sat}$ (kPa)
* - 20
  - 2.3393
* - 30
  - 4.2470
* - 40
  - 7.3849
* - 50
  - 12.352
* - 60
  - 19.946
* - 70
  - 31.201
* - 80
  - 47.414
* - 90
  - 70.182
* - 100
  - 101.42
* - 120
  - 198.67
* - 140
  - 361.54
* - 160
  - 618.23
* - 180
  - 1002.8
```

### Part A. Global Clausius–Clapeyron Fit (the "Straight-Line" Model)

**Task.** Plot $y = \ln(P_{\mathrm{sat}}/1\ \mathrm{bar})$ versus $x = 1/T$ (with $T$ in K), fit a line, and extract $\Delta\bar h_{\mathrm{vap}}$ from the slope.

```{code-cell} ipython3
import io
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats

R = 8.314462618  # J mol^-1 K^-1

data = """T_C,P_kPa
20,2.3393
30,4.2470
40,7.3849
50,12.352
60,19.946
70,31.201
80,47.414
90,70.182
100,101.42
120,198.67
140,361.54
160,618.23
180,1002.8
"""

df = pd.read_csv(io.StringIO(data))
df["T_K"] = df["T_C"] + 273.15
df["P_bar"] = df["P_kPa"] / 100.0  # 100 kPa = 1 bar
df["invT"] = 1.0 / df["T_K"]
df["lnP"] = np.log(df["P_bar"])    # ln(P/1 bar)

# Linear regression: lnP = m*(1/T) + b
m, b, r, p, se = stats.linregress(df["invT"], df["lnP"])
DeltaHvap = -m * R  # J/mol

print(f"slope m            = {m:10.2f} K")
print(f"intercept b        = {b:10.3f}")
print(f"R^2                = {r**2:10.6f}")
print(f"DeltaHvap (global) = {DeltaHvap/1000:10.2f} kJ/mol")

# Plot
x = df["invT"].to_numpy()
y = df["lnP"].to_numpy()
x_fit = np.linspace(x.min(), x.max(), 200)
y_fit = m*x_fit + b

plt.figure(figsize=(6,4))
plt.plot(x, y, "o", label="data")
plt.plot(x_fit, y_fit, "-", label="linear fit")
plt.xlabel(r"$1/T\ (\mathrm{K^{-1}})$")
plt.ylabel(r"$\ln(P_\mathrm{sat}/1\,\mathrm{bar})$")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()
```

**Interpretation prompts (one or two sentences each):**

- What units should $\Delta\bar h_{\mathrm{vap}}$ have, and how do they emerge from the slope formula?
- Why is the slope of $\ln P_{\mathrm{sat}}$ versus $1/T$ negative?
- Does the line *look* straight over the full temperature range?

### Part B. Linearity Limits: Residuals Beat "It Looks Straight"

A Clausius–Clapeyron plot can look linear by eye even when the constant-$\Delta\bar h_{\mathrm{vap}}$ assumption is quietly failing. A quick diagnostic is a **residual plot**.

```{code-cell} ipython3
df["lnP_fit"] = m*df["invT"] + b
df["residual"] = df["lnP"] - df["lnP_fit"]

plt.figure(figsize=(6,4))
plt.axhline(0, linewidth=1)
plt.plot(df["invT"], df["residual"], "o-")
plt.xlabel(r"$1/T\ (\mathrm{K^{-1}})$")
plt.ylabel("residual = data − fit (in ln units)")
plt.grid(True)
plt.tight_layout()
plt.show()
```

**What to look for:**

- **Random scatter around zero** → the linear model is adequate over this range.
- **A systematic curve** (e.g. negative at both ends, positive in the middle, or vice versa) → the straight-line model is missing physics, most likely because $\Delta\bar h_{\mathrm{vap}}(T)$ is changing with $T$ across the dataset.

### Part C. Two-Range Fits: Do You Get the Same $\Delta\bar h_{\mathrm{vap}}$?

If $\Delta\bar h_{\mathrm{vap}}$ were truly constant, then fitting two separate temperature subranges should give the same slope within uncertainty. If they disagree systematically, that is hard evidence that $\Delta\bar h_{\mathrm{vap}}$ depends on $T$.

```{code-cell} ipython3
def fit_subset(df_sub):
    m, b, r, p, se = stats.linregress(df_sub["invT"], df_sub["lnP"])
    return -m * R / 1000, r**2

for name, lo, hi in (("Low-T ", 20, 90), ("High-T", 120, 180)):
    dH, r2 = fit_subset(df[df["T_C"].between(lo, hi)])
    print(f"{name} fit ({lo}-{hi} °C): DeltaHvap = {dH:.2f} kJ/mol, R^2 = {r2:.6f}")
```

### Part D. The Local Slope and $\Delta\bar h_{\mathrm{vap}}(T)$

From the differential form {eq}`eq:CC-differential`,

```{math}
\Delta\bar h_{\mathrm{vap}}(T) \approx R\,T^2\left(\frac{d\ln P_{\mathrm{sat}}}{dT}\right),
```

so we can estimate $\Delta\bar h_{\mathrm{vap}}(T)$ from finite differences of the $(T,\ln P_{\mathrm{sat}})$ data. The code below applies the two-point formula {eq}`eq:CC-two-point` to each pair of neighboring points and assigns the result to the midpoint of the interval. It then compares the result with the reference enthalpy of vaporization at the same temperature (NIST Chemistry WebBook).

```{code-cell} ipython3
T = df["T_K"].to_numpy()
lnP = df["lnP"].to_numpy()

# two-point estimate for each pair of neighboring points
T_mid = 0.5 * (T[1:] + T[:-1])
local_slope = -R * np.diff(lnP) / np.diff(1.0 / T) / 1000   # kJ/mol

# reference enthalpy of vaporization at the midpoints (kJ/mol)
dH_ref = [43.99, 43.56, 43.13, 42.69, 42.25, 41.81,
          41.35, 40.89, 40.17, 39.16, 38.08, 36.91]

local = pd.DataFrame({
    "T_mid_C": T_mid - 273.15,
    "local_slope_kJmol": local_slope.round(2),
    "DeltaHvap_ref_kJmol": dH_ref,
})

plt.figure(figsize=(6,4))
plt.plot(local["T_mid_C"], local["local_slope_kJmol"], "o-", label="local slope")
plt.plot(local["T_mid_C"], local["DeltaHvap_ref_kJmol"], "s--", color="0.4",
         label=r"reference $\Delta\bar h_\mathrm{vap}$")
plt.xlabel(r"$T\ (^\circ\mathrm{C})$")
plt.ylabel("kJ/mol")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()

local
```

**Interpretation.** The local slope decreases with increasing temperature, as $\Delta\bar h_{\mathrm{vap}}$ does, but the two separate as the temperature rises: at 170 °C the local slope gives $39.5\ \mathrm{kJ\,mol^{-1}}$ and the enthalpy of vaporization is $36.9\ \mathrm{kJ\,mol^{-1}}$. The reason is that the local slope is not $\Delta\bar h_{\mathrm{vap}}$ itself. By Eq. {eq}`eq:CC-exact`, $RT^2\,d\ln P_{\mathrm{sat}}/dT = \Delta\bar h_{\mathrm{vap}}/\Delta Z$, and $\Delta Z$ falls below 1 as the vapor becomes denser and less ideal. Physically, as $T$ rises toward the critical point, the liquid and vapor phases become more and more similar in density and structure, so both the molar enthalpy difference and the molar volume difference between them shrink. At the critical point itself, $\Delta\bar h_{\mathrm{vap}}\to 0$ and $\Delta\bar v_{\mathrm{vap}}\to 0$ together: there is no longer a meaningful "vaporization" because the two phases have merged.

```{dropdown} Check Your Work (Expected Ballpark Results)
Using the dataset above, typical results are:

- Global fit (all 13 points): $\Delta\bar h_{\mathrm{vap}}\approx 42\ \mathrm{kJ\,mol^{-1}}$ (an average over the whole 20–180 °C range).
- Low-$T$ fit (20–90 °C): $\Delta\bar h_{\mathrm{vap}}\approx 43\ \mathrm{kJ\,mol^{-1}}$.
- High-$T$ fit (120–180 °C): $\Delta\bar h_{\mathrm{vap}}\approx 40\ \mathrm{kJ\,mol^{-1}}$.
- Local slope: trends downward across the dataset, from about $44\ \mathrm{kJ\,mol^{-1}}$ at 25 °C to about $39.5\ \mathrm{kJ\,mol^{-1}}$ at 170 °C, while the reference $\Delta\bar h_{\mathrm{vap}}$ falls from $44.0$ to $36.9\ \mathrm{kJ\,mol^{-1}}$.

If your numbers differ noticeably, check that you used (1) **kelvin**, (2) **natural log** (not $\log_{10}$), and (3) a **dimensionless** pressure inside the log (e.g. $P/1\ \mathrm{bar}$).
```

### Reflection Questions

1. Over what temperature window does the Clausius–Clapeyron straight-line model look "good enough" *and* have residuals that are roughly random?
2. List two physical reasons why $\ln P$ versus $1/T$ might curve over a wider range.
3. If you needed a better model than a single $\Delta\bar h_{\mathrm{vap}}$, what would you do? (Examples: fit two ranges, use a published vapor-pressure correlation such as Antoine's equation, or include the temperature dependence of $\Delta\bar h_{\mathrm{vap}}$ via heat-capacity corrections.)

## Computational Studio: Clausius–Clapeyron

Interact with real vapor-pressure data for ethanol. The studio transforms the dataset to visualize the Clausius–Clapeyron linearization ($\ln P$ versus $1/T$), lets you dynamically adjust the regression bounds to explore linearity limits, and reports the corresponding $\Delta\bar h_{\mathrm{vap}}$ extracted from the slope.

You can open the studio in a new tab: [Clausius–Clapeyron Studio](https://wexlergroup.github.io/chem-4020-5020/computational-studios/clausius-clapeyron-studio.html).

## Concept Checks

1. Why does $d\ln P/dT$ in the Clausius–Clapeyron equation depend on $1/T^2$ rather than simply on $1/T$? Trace the two factors of $T$ back to the steps of the derivation.
2. What is the sign of $dP/dT$ for the **ice → water** transition, and what specific feature of the molar volumes is responsible for that sign? Contrast with the corresponding transition for a "normal" substance such as benzene.
3. Why is it usually safe to approximate $\Delta\bar v \approx \bar v_{\mathrm{g}}$ for vaporization but **not** for melting?
4. The local-slope analysis in Part D of the mini-lab shows $\Delta\bar h_{\mathrm{vap}}(T)$ decreasing as $T$ rises. Explain qualitatively what would happen to $d\ln P/dT$ — and to the apparent slope of a $\ln P$-versus-$1/T$ plot — if you were to extend the dataset all the way up to the critical point.
5. A student fits a Clausius–Clapeyron line to vapor-pressure data and extracts $\Delta\bar h_{\mathrm{vap}} = 38\ \mathrm{kJ\,mol^{-1}}$ for a polar liquid that should give around 45 kJ mol$^{-1}$ from calorimetry. List two distinct error sources that could be responsible.

## Key Takeaways

- Differentiating the equilibrium condition $\mu_\alpha = \mu_\beta$ along a coexistence curve produces the **Clapeyron equation** $dP/dT = \Delta\bar s_{\mathrm{trs}}/\Delta\bar v_{\mathrm{trs}} = \Delta\bar h_{\mathrm{trs}}/(T\,\Delta\bar v_{\mathrm{trs}})$, which is exact for any first-order transition between pure phases.
- The **sign** of $dP/dT$ is set by the sign of $\Delta\bar v$. For water, $\bar v_{\mathrm{s}} > \bar v_{\mathrm{l}}$ at the melting line, so the ice–water coexistence curve has a *negative* slope ($\sim -135$ bar K$^{-1}$).
- For liquid–vapor equilibrium with $\bar v_{\mathrm{g}} \gg \bar v_{\mathrm{l}}$ and an ideal vapor, the Clapeyron equation specializes to the **Clausius–Clapeyron equation** $d\ln P_{\mathrm{sat}}/dT = \Delta\bar h_{\mathrm{vap}}/(RT^2)$.
- Under the constant-$\Delta\bar h_{\mathrm{vap}}$ approximation, the integrated Clausius–Clapeyron equation gives $\ln(P_2/P_1) = -(\Delta\bar h_{\mathrm{vap}}/R)(1/T_2 - 1/T_1)$, which extracts $\Delta\bar h_{\mathrm{vap}}$ from two vapor-pressure points or from a linear fit of $\ln P_{\mathrm{sat}}$ versus $1/T$.
- Deviations from a straight line in $\ln P_{\mathrm{sat}}$ versus $1/T$ — best diagnosed from residuals, not by eye — signal the breakdown of the constant-$\Delta\bar h_{\mathrm{vap}}$ approximation. The local slope $RT^2\,d\ln P_{\mathrm{sat}}/dT$ equals $\Delta\bar h_{\mathrm{vap}}/\Delta Z$ with $\Delta Z = P_{\mathrm{sat}}\,\Delta\bar v_{\mathrm{vap}}/RT$: it follows $\Delta\bar h_{\mathrm{vap}}(T)$ closely at low vapor pressure and departs from it as the vapor becomes denser.
- For water, $\Delta\bar h_{\mathrm{vap}} \approx 40\ \mathrm{kJ\,mol^{-1}}$ is consistent with the molecular picture of breaking roughly two hydrogen bonds per molecule on going from liquid to vapor — the macroscopic shadow of the intermolecular forces introduced in [Section 1.4](../chapter-01/section-04.md).
