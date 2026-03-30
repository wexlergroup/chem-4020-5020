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

# Review

## 1. Checklist of Key Concepts

### Section 3.1: Conservation of Energy

1. **First Law of Thermodynamics**
   - The internal energy change of a system equals the heat absorbed plus the work done on it:

   ```{math}
   \Delta U = q + w.
   ```

   - In differential form: $dU = \delta q + \delta w$, or equivalently $\delta q = dU - \delta w$.

2. **State Functions vs. Path Functions**
   - $U$ is a **state function**: $\Delta U$ depends only on the initial and final states, and $dU$ is an exact differential ($\oint dU = 0$).
   - $q$ and $w$ are **path functions**: their values depend on the process path, and $\delta q$, $\delta w$ are inexact differentials.

3. **Sign Convention** (Chemistry Convention)
   - $q > 0$: heat absorbed by the system.
   - $w > 0$: work done **on** the system.
   - For $PV$ work: $\delta w = -P\,dV$, so compression ($dV < 0$) gives $w > 0$.

4. **Generalized Work**
   - Many work modes share a generalized force $\times$ generalized displacement structure: $PV$ work ($-P\,dV$), surface work ($\gamma\,dA$), elastic work ($k\,dl$), electrical work ($\mathcal{E}\,dq_{\mathrm{el}}$), and chemical work ($\mu\,dN$).

5. **Free Expansion**
   - Expansion into vacuum: $P_{\mathrm{ext}} = 0$, so $w = 0$. If the container is also insulated, $q = 0$ and $\Delta U = 0$. For an ideal gas, $\Delta T = 0$ as well.
   - Free expansion is the extreme case of an irreversible process: the gas changes state while doing no work and exchanging no heat.

---

### Section 3.2: Applications of the First Law

1. **Thermodynamic Processes**
   - **Quasi-static**: carried out infinitesimally slowly; the system remains near equilibrium at all times.
   - **Reversible**: quasi-static *and* free of dissipative effects (friction, turbulence, viscous drag, etc.); can be reversed with no net change to system or surroundings.
   - **Irreversible**: any real process that violates one or more conditions for reversibility.

2. **Reversible Work Bounds**
   - For a given expansion ($\Delta V > 0$), a reversible process does the **maximum** work on the surroundings ($|w_{\mathrm{rev}}| \geq |w_{\mathrm{irrev}}|$).
   - For a given compression, a reversible process requires the **minimum** work input.

3. **Five-Step First Law Workflow**
   1. Choose two independent variables (e.g., $V, T$ or $P, T$).
   2. Rewrite the First Law in terms of those variables.
   3. Apply process constraints (isothermal, isochoric, isobaric, adiabatic).
   4. Specify the equation of state (e.g., ideal gas).
   5. Integrate to find $q$, $w$, $\Delta U$.

4. **First Law with $V$ and $T$ as Independent Variables**
   - General form:

   ```{math}
   \delta q \;=\; C_V \, dT \;+\; \left[\left(\frac{\partial U}{\partial V}\right)_T + P\right] dV.
   ```

   - For an ideal gas, $(\partial U/\partial V)_T = 0$, simplifying to $\delta q = C_V\,dT + P\,dV$.

5. **Common Ideal-Gas Processes**
   - **Isochoric** ($dV = 0$): $w = 0$, $q = C_V\,\Delta T$, $\Delta U = C_V\,\Delta T$.
   - **Isothermal** ($dT = 0$): $\Delta U = 0$ (ideal gas), $q = -w = Nk_{\mathrm{B}}T\ln(V_2/V_1)$.
   - **Adiabatic** ($\delta q = 0$): $\Delta U = w = C_V\,\Delta T$, and $TV^{\gamma - 1} = \text{const}$ (equivalently $PV^\gamma = \text{const}$), where $\gamma = C_P/C_V$.

6. **Microscopic Interpretation of the First Law**
   - From $U = \sum_i p_i E_i$:

   ```{math}
   dU = \underbrace{\sum_i E_i\, dp_i}_{\delta q} + \underbrace{\sum_i p_i\, dE_i}_{\delta w}.
   ```

   - **Heat** changes *which states are occupied* (probabilities shift; energy levels fixed).
   - **Work** changes *the energies of the states* (energy levels shift; probabilities fixed).

---

### Section 3.3: Enthalpy

1. **Definition and Motivation**
   - **Enthalpy**: $H = U + PV$.
   - At constant pressure with $PV$-only work: $\delta q_P = dH$, so $q_P = \Delta H$.
   - Defines a state function that plays the same role at constant $P$ that $U$ plays at constant $V$.

2. **Heat Capacity at Constant Pressure**

   ```{math}
   C_P = \left(\frac{\partial H}{\partial T}\right)_P
   = \left(\frac{\partial U}{\partial T}\right)_P + P\left(\frac{\partial V}{\partial T}\right)_P.
   ```

   - For an ideal gas: $C_P - C_V = Nk_{\mathrm{B}}$ (or $nR$ per mole), so $C_P > C_V$.

3. **When $q_P = \Delta H$ Fails**
   - The relation holds only when $PV$ work is the sole form of work. Non-$PV$ work at constant pressure (e.g., electrical work in an electrochemical cell) breaks the equivalence.

4. **Standard States and Formation Enthalpies**
   - Standard pressure: $P^\circ = 1\,\text{bar}$; reference temperature typically 298.15 K.
   - Standard enthalpy of formation, $\Delta H_f^\circ$: enthalpy change when 1 mol of compound is formed from elements in their standard states.
   - By convention, $\Delta H_f^\circ = 0$ for elements in their standard states.

5. **Hess's Law**
   - Because $H$ is a state function, enthalpy changes are path-independent and additive:

   ```{math}
   \Delta H_{\mathrm{rxn}}^\circ
   = \sum_{p} \nu_p \,\Delta H_{f,p}^\circ
   - \sum_{r} \nu_r \,\Delta H_{f,r}^\circ.
   ```

6. **First Law Toolkit Summary**

   | Constraint | Relevant state function | Key relation |
   | :---: | :---: | :---: |
   | Constant $V$ | $U$ | $q_V = \Delta U$ |
   | Constant $P$ | $H = U + PV$ | $q_P = \Delta H$ |

---

## 2. Checklist of Most Important Equations

Below is a unified list of the major equations from Sections 3.1–3.3.

A. **First Law of Thermodynamics**

```{math}
\Delta U = q + w
\qquad\text{(finite)},
\qquad
dU = \delta q + \delta w
\qquad\text{(differential)}.
```

- **Applicability**: any closed system. $q > 0$ = heat absorbed by the system; $w > 0$ = work done on the system.

---

B. **$PV$ Work**

```{math}
\delta w = -P_{\mathrm{ext}}\,dV
\qquad\Longrightarrow\qquad
w = -\int_{V_1}^{V_2} P_{\mathrm{ext}}\,dV.
```

- **Applicability**: any expansion or compression. For reversible processes, $P_{\mathrm{ext}} = P$ (system pressure). For irreversible processes, $P_{\mathrm{ext}}$ is determined by the surroundings.

---

C. **First Law in $(V, T)$ Variables** (for $PV$-only work)

```{math}
\delta q = C_V\,dT + \left[\left(\frac{\partial U}{\partial V}\right)_T + P\right] dV.
```

- For an ideal gas, $(\partial U/\partial V)_T = 0$, giving $\delta q = C_V\,dT + P\,dV$.

---

D. **Isothermal Reversible Work** (ideal gas)

```{math}
w = -nRT\ln\!\left(\frac{V_2}{V_1}\right),
\qquad
q = -w = nRT\ln\!\left(\frac{V_2}{V_1}\right).
```

- $\Delta U = 0$ for an ideal-gas isothermal process.

---

E. **Adiabatic Relations** (ideal gas, reversible)

```{math}
TV^{\gamma - 1} = \text{const},
\qquad
PV^{\gamma} = \text{const},
\qquad
\gamma = \frac{C_P}{C_V}.
```

- For a monatomic ideal gas, $\gamma = 5/3$ and the $TV$ relation becomes $TV^{2/3} = \text{const}$.

---

F. **Isochoric Process** ($dV = 0$)

```{math}
q_V = \Delta U = \int_{T_1}^{T_2} C_V\,dT.
```

---

G. **Enthalpy**

```{math}
H = U + PV.
```

- At constant pressure with $PV$-only work: $q_P = \Delta H$.

---

H. **Heat Capacity at Constant Pressure**

```{math}
C_P = \left(\frac{\partial H}{\partial T}\right)_P,
\qquad
\Delta H = \int_{T_1}^{T_2} C_P\,dT.
```

- For an ideal gas: $C_P - C_V = Nk_{\mathrm{B}}$ (per particle) or $C_P - C_V = nR$ (per mole).

---

I. **Hess's Law**

```{math}
\Delta H_{\mathrm{rxn}}^\circ
= \sum_{p} \nu_p\,\Delta H_{f,p}^\circ
- \sum_{r} \nu_r\,\Delta H_{f,r}^\circ.
```

- **Applicability**: any reaction at standard conditions, using tabulated formation enthalpies.

---

J. **Microscopic First Law**

```{math}
dU = \underbrace{\sum_i E_i\, dp_i}_{\delta q}
\;+\; \underbrace{\sum_i p_i\, dE_i}_{\delta w}.
```

- **Applicability**: closed system with $PV$-only work. Heat redistributes probabilities; work shifts energy levels.
