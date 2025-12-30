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

# 6.1. Phase Diagrams

## Overview

This section applies **thermodynamic equilibrium conditions** to **phases of matter** and shows how those conditions generate the familiar **phase diagram** (e.g., a $P$–$T$ diagram for a pure substance).

From the lecture hand-written notes (Apr 16):

1. **Phase equilibrium condition:** $G_\alpha = G_\beta$ (shorthand for “equal molar Gibbs free energy / equal chemical potential”) at a phase transition.
2. **First-order transition signatures:** $H_\alpha \neq H_\beta$ (**latent heat**) and $S_\alpha \neq S_\beta$ (**entropy jump**) at phase equilibrium.
3. **Open systems:** the **chemical potential** $\mu$ is the natural quantity to use.
4. **Coexistence curves:** the **(Clausius–)Clapeyron equation** provides the slope(s) of phase boundaries and lets us draw phase diagrams.

## 6.1.1 A quick tour of a $P$–$T$ phase diagram (one-component)

A schematic $P$–$T$ phase diagram for a single component contains regions where one phase is stable (solid, liquid, gas) separated by **coexistence curves** where two phases are in equilibrium.

Key landmarks:

- **Triple point (tp):** all three phases coexist (solid–liquid–gas) at a single $(T,P)$.
- **Critical point (cp):** the end of the liquid–gas coexistence curve. Beyond the critical point the system is a **supercritical fluid (scf)** with no distinct liquid–gas boundary.

The coexistence curves are **phase-transition lines**: moving across a line changes which phase is thermodynamically stable at that $(T,P)$.

## 6.1.2 Thermodynamic criterion for phase equilibrium

### Gibbs free energy and chemical potential

For a closed, fixed-composition system,

```{math}
dG = V\,dP - S\,dT.
```

For an **open** system (composition can change),

```{math}
dG = V\,dP - S\,dT + \sum_i \mu_i\,dn_i,
```

where the **chemical potential** is

```{math}
\mu_i \equiv \left(\frac{\partial G}{\partial n_i}\right)_{T,P,n_{j\neq i}}.
```

For a **one-component** system, the chemical potential equals the **molar Gibbs free energy**:

```{math}
\mu = g \equiv \frac{G}{n}.
```

### Phase equilibrium condition ($\mu_\alpha = \mu_\beta$)

At equilibrium between two phases $\alpha$ and $\beta$ of the _same substance_,

```{math}
\mu_\alpha(T,P) = \mu_\beta(T,P),
```

i.e. the molar Gibbs free energies are equal at the phase boundary.

A useful way to remember the stability criterion:

- At fixed $(T,P)$ the **stable phase** is the one with the **lowest** $g$ (or $\mu$).
- A **phase transition** occurs where the $g$–surfaces for two phases intersect.

## 6.1.3 What “$G_\alpha = G_\beta$” looks like on graphs

The lecture notes emphasize that phase transitions correspond to intersections of Gibbs free-energy curves, and that the **slopes** of those curves encode entropy and volume.

### (A) $G$ vs. $T$ at fixed pressure

At constant pressure,

```{math}
\left(\frac{\partial G}{\partial T}\right)_P = -S.
```

So on a $G$–$T$ plot:

- The **slope** is negative.
- The phase with larger entropy has a **more negative** slope.

Typically,

```{math}
S_g > S_l > S_s,
```

so the gas line is steepest (most negative slope), then liquid, then solid.

At the transition temperature $T_{tr}$ (for that pressure), the two curves cross:

```{math}
G_\alpha(T_{tr}) = G_\beta(T_{tr}),
```

and the system switches which phase gives the lowest $G$.

### (B) $G$ vs. $P$ at fixed temperature

At constant temperature,

```{math}
\left(\frac{\partial G}{\partial P}\right)_T = V.
```

So on a $G$–$P$ plot:

- The **slope** is the (molar) volume.
- The phase with larger molar volume rises more steeply with pressure.

Typically,

```{math}
V_g > V_l > V_s,
```

so the gas line has the largest slope, then liquid, then solid. Intersections of these lines correspond to the pressures where phase coexistence occurs at that temperature.

## 6.1.4 Latent heat and entropy jump at phase equilibrium

At many phase transitions (fusion, vaporization, sublimation) the first derivatives of $G$ are discontinuous. The lecture notes highlight two consequences:

- **Latent heat:** the enthalpy changes discontinuously across the boundary.
- **Entropy jump:** the entropy changes discontinuously across the boundary.

### Definitions at the transition temperature

For a transition $\alpha \rightarrow \beta$ occurring at $(T_{tr},P_{tr})$:

- Latent heat (molar enthalpy change):

  ```{math}
  \Delta H_{tr} \equiv H_\beta(T_{tr}) - H_\alpha(T_{tr}).
  ```

- Entropy jump:

  ```{math}
  \Delta S_{tr} \equiv S_\beta(T_{tr}) - S_\alpha(T_{tr}).
  ```

Because $\Delta G=0$ at coexistence and $\Delta G = \Delta H - T\Delta S$, we obtain the key relation

```{math}
\boxed{\Delta H_{tr} = T_{tr}\,\Delta S_{tr}}.
```

### Common special cases

- **Fusion (melting) at $T_{fus}$:**

  ```{math}
  \Delta H_{fus} = H_l(T_{fus}) - H_s(T_{fus}), \qquad
  \Delta S_{fus} = S_l(T_{fus}) - S_s(T_{fus}).
  ```

- **Vaporization at $T_{vap}$ (boiling point at the chosen pressure):**

  ```{math}
  \Delta H_{vap} = H_g(T_{vap}) - H_l(T_{vap}), \qquad
  \Delta S_{vap} = S_g(T_{vap}) - S_l(T_{vap}).
  ```

#### Trouton's rule (rule-of-thumb)

The notes mention **Trouton's rule (1884)**: for many liquids at their _normal_ boiling point,

```{math}
\Delta S_{vap} \approx 85.9\;\mathrm{J\,mol^{-1}\,K^{-1}} \approx 10.3\,R.
```

This is a useful “order-of-magnitude” estimate; real substances can deviate (especially when there is strong association / hydrogen bonding, or when the liquid structure is unusually ordered).

#### Worked example: using tabulated $C_P(T)$ + an entropy jump (water at 1 bar)

**Goal.** Estimate $\Delta S$ and $\Delta H$ when **1 mol** of water is heated _reversibly_ at constant pressure $p=1\ \mathrm{bar}$ from **298.15 K** (liquid water) to **800 K** (steam), crossing the liquid–vapor boundary.

Within a _single phase_ at constant pressure:

```{math}
\Delta S = \int_{T_1}^{T_2} \frac{C_P(T)}{T}\,dT,
\qquad
\Delta H = \int_{T_1}^{T_2} C_P(T)\,dT.
```

At a first-order phase transition (here, vaporization) the entropy has a _finite jump_:

```{math}
\Delta S_{vap}(T_{tr}) = \frac{\Delta H_{vap}(T_{tr})}{T_{tr}}.
```

We'll use a **coarse** set of $C_P^\circ(T)$ values from the **NIST–JANAF thermochemical tables** for $\mathrm{H_2O(l)}$ and $\mathrm{H_2O(g)}$ at the standard-state pressure $p^\circ=0.1\ \mathrm{MPa}=1\ \mathrm{bar}$. In that dataset, the liquid–vapor transition occurs at

```{math}
T_{tr}=372.780\ \mathrm{K}\quad\text{(very close to 373.15 K at 1 atm)}.
```

##### Data used (excerpt)

| Phase | $T$ (K) | $C_P^\circ$ (J mol$^{-1}$ K$^{-1}$) |
| --- | ---: | ---: |
| $\mathrm{H_2O(l)}$ | 298.15 | 75.351 |
| | 300 | 75.349 |
| | 320 | 75.344 |
| | 340 | 75.388 |
| | 360 | 75.679 |
| | 372.780 | 75.962 |
| $\mathrm{H_2O(g)}$ | 300 | 33.596 |
| | 400 | 34.262 |
| | 500 | 35.226 |
| | 600 | 36.325 |
| | 700 | 37.495 |
| | 800 | 38.721 |

For the gas-phase heating integral, we linearly interpolate $C_P^\circ$ to $T_{tr}$ using the 300–400 K values.

To estimate the **latent heat** at the transition we use the tabulated standard enthalpies of formation (same JANAF source) and interpolate them to $T_{tr}$:

```{math}
\Delta H_{vap}(T_{tr})
\approx
\Delta_f H^\circ_{g}(T_{tr})-\Delta_f H^\circ_{l}(T_{tr}).
```

```{code-cell} ipython3
import numpy as np

# ---- Cp(T) data (J/mol/K) ----
T_tr = 372.780  # K (liquid <-> vapor at 1 bar in the JANAF tables)

T_liq = np.array([298.15, 300, 320, 340, 360, T_tr])
Cp_liq = np.array([75.351, 75.349, 75.344, 75.388, 75.679, 75.962])

T_gas_tab = np.array([300, 400, 500, 600, 700, 800])
Cp_gas_tab = np.array([33.596, 34.262, 35.226, 36.325, 37.495, 38.721])

# Interpolate Cp_gas at T_tr (between 300 and 400 K)
Cp_gas_tr = Cp_gas_tab[0] + (T_tr - T_gas_tab[0])/(T_gas_tab[1] - T_gas_tab[0])*(Cp_gas_tab[1] - Cp_gas_tab[0])

T_gas = np.concatenate([[T_tr], T_gas_tab[1:]])
Cp_gas = np.concatenate([[Cp_gas_tr], Cp_gas_tab[1:]])

# ---- Heating contributions (trapezoid rule) ----
dS_liq = np.trapezoid(Cp_liq/T_liq, T_liq)          # J/mol/K
dH_liq = np.trapezoid(Cp_liq, T_liq)/1000           # kJ/mol

dS_gas = np.trapezoid(Cp_gas/T_gas, T_gas)          # J/mol/K
dH_gas = np.trapezoid(Cp_gas, T_gas)/1000           # kJ/mol

# ---- Vaporization jump from JANAF delta_f H° (linear interpolation to T_tr) ----
# liquid delta_f H° at 360 and 400 K (kJ/mol)
Hf_liq_360, Hf_liq_400 = -283.874, -282.591
Hf_liq_tr = Hf_liq_360 + (T_tr-360)/(400-360)*(Hf_liq_400 - Hf_liq_360)

# gas delta_f H° at 300 and 400 K (kJ/mol)
Hf_gas_300, Hf_gas_400 = -241.844, -242.846
Hf_gas_tr = Hf_gas_300 + (T_tr-300)/(400-300)*(Hf_gas_400 - Hf_gas_300)

dH_vap = Hf_gas_tr - Hf_liq_tr                 # kJ/mol
dS_vap = dH_vap*1000/T_tr                      # J/mol/K

print(f"Heating (liquid, 298.15 -> {T_tr:.3f} K):  ΔH = {dH_liq:6.3f} kJ/mol,  ΔS = {dS_liq:6.2f} J/mol/K")
print(f"Vaporization at {T_tr:.3f} K (1 bar):        ΔH = {dH_vap:6.2f} kJ/mol,  ΔS = {dS_vap:6.2f} J/mol/K")
print(f"Heating (gas,    {T_tr:.3f} -> 800 K):       ΔH = {dH_gas:6.2f} kJ/mol,  ΔS = {dS_gas:6.2f} J/mol/K")

dH_total = dH_liq + dH_vap + dH_gas
dS_total = dS_liq + dS_vap + dS_gas
print(f"\nTOTAL (298.15 K liquid -> 800 K steam, 1 bar): ΔH ≈ {dH_total:6.2f} kJ/mol,  ΔS ≈ {dS_total:6.1f} J/mol/K")
```

**What to notice.**

- The entropy change is dominated by the **phase transition**: $\Delta S_{vap}\approx 1.1\times 10^2\ \mathrm{J\,mol^{-1}\,K^{-1}}$.
- Water is a classic _exception_ to Trouton's rule (hydrogen bonding makes the liquid unusually ordered, increasing $\Delta S_{vap}$).

```{code-cell} ipython3
import matplotlib.pyplot as plt

# Cumulative entropy relative to 298.15 K along the path (trapezoid segments)
S_liq_cum = np.concatenate([[0.0], np.cumsum(0.5*(Cp_liq[1:]/T_liq[1:] + Cp_liq[:-1]/T_liq[:-1]) * np.diff(T_liq))])
S_gas_cum = np.concatenate([[0.0], np.cumsum(0.5*(Cp_gas[1:]/T_gas[1:] + Cp_gas[:-1]/T_gas[:-1]) * np.diff(T_gas))])

T_path = np.concatenate([T_liq, T_gas[1:]])
S_path = np.concatenate([S_liq_cum, S_liq_cum[-1] + dS_vap + S_gas_cum[1:]])

plt.figure()
plt.plot(T_path, S_path)
plt.axvline(T_tr, linestyle="--")
plt.xlabel("Temperature (K)")
plt.ylabel(r"$\Delta S$ from 298.15 K (J mol$^{-1}$ K$^{-1}$)")
plt.title("Heating water at 1 bar: entropy rise + vaporization jump")
plt.show()
```

## 6.1.5 Clapeyron and Clausius–Clapeyron equations (phase-boundary slopes)

The coexistence curves in a phase diagram come from enforcing $\mu_\alpha(T,P)=\mu_\beta(T,P)$ _along the boundary_.

For a one-component phase, the differential of the chemical potential is

```{math}
d\mu = -s\,dT + v\,dP,
```

where $s$ and $v$ are **molar** entropy and **molar** volume.

Along the coexistence curve,

```{math}
d\mu_\alpha = d\mu_\beta
```

so

```{math}
-s_\alpha\,dT + v_\alpha\,dP = -s_\beta\,dT + v_\beta\,dP.
```

Rearranging gives the **Clapeyron equation**

```{math}
\boxed{\frac{dP}{dT} = \frac{\Delta s}{\Delta v} = \frac{\Delta H_{tr}}{T\,\Delta v}},
```

where $\Delta s = s_\beta-s_\alpha$, $\Delta v=v_\beta-v_\alpha$, and we used $\Delta H_{tr}=T\Delta s$.

### Clausius–Clapeyron (liquid–vapor approximation)

For liquid–vapor equilibrium:

- $v_g \gg v_l$ so $\Delta v \approx v_g$
- if the vapor is ideal: $v_g = \dfrac{RT}{P}$

Insert these into Clapeyron:

```{math}
\frac{dP}{dT} \approx \frac{\Delta H_{vap}}{T\,(RT/P)} = \frac{P\,\Delta H_{vap}}{R T^2}.
```

Equivalently,

```{math}
\boxed{\frac{d\ln P_{sat}}{dT} = \frac{\Delta H_{vap}}{R T^2}}.
```

If $\Delta H_{vap}$ is approximately constant over the temperature range of interest, integration gives

```{math}
\ln P_{sat}(T) \approx -\frac{\Delta H_{vap}}{R}\,\frac{1}{T} + C,
```

which is the working form often used to fit vapor-pressure data and sketch the liquid–gas coexistence curve.

## Key takeaways

- Phase boundaries are where **chemical potentials match**: $\mu_\alpha=\mu_\beta$.
- Intersections of $G(T)$ or $G(P)$ curves visualize phase transitions:
  - $(\partial G/\partial T)_P=-S$  (entropy controls slope in $G$–$T$)
  - $(\partial G/\partial P)_T=V$  (volume controls slope in $G$–$P$)
- First-order transitions show **latent heat** and **entropy jumps**, related by $\Delta H = T\Delta S$.
- The **Clapeyron** (and **Clausius–Clapeyron**) equations determine **phase-boundary slopes** in a $P$–$T$ diagram.
