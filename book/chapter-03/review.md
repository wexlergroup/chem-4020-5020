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

The course map shows where this module sits in the course and which relations connect the microscopic and macroscopic descriptions. Select a module to unfold its sections.

<iframe
  src="https://wexlergroup.github.io/chem-4020-5020/computational-studios/course-map.html?module=3&amp;here=review"
  title="Course map, opened at Module 3"
  style="width: 100%; height: 830px; border: 0;"
  loading="lazy"
></iframe>

```{image} ../_static/course-map/module-3.png
:class: hidden
:alt: Course map with Module 3 unfolded. Statistical mechanics (Module 2) forms the microscopic strand, and Modules 1 and 3 to 7 form the macroscopic strand. A relation joins the two strands at each module except Module 6, and for Module 3 it is the internal energy as the ensemble average of the energy. Below the map, one card for each of Sections 3.1 to 3.3 lists the microscopic content, the connecting relations, and the macroscopic content of the section.
```

## 1. Checklist of Key Concepts

### Section 3.1: Conservation of Energy

See [Section 3.1](section-01.md).

1. **Mechanical Equivalent of Heat**
   - Joule found that a definite quantity of mechanical work always produces the same quantity of heat: about $4.19\ \mathrm{J}$ of work raises the temperature of $1\ \mathrm{g}$ of water by $1\,\text{°C}$.

2. **First Law of Thermodynamics**
   - The internal energy change of a closed system equals the heat absorbed plus the work done on it:

   ```{math}
   \Delta U = q + w.
   ```

   - In differential form: $dU = \delta q + \delta w$, or equivalently $\delta q = dU - \delta w$.

3. **State Functions vs. Path Functions**
   - $U$ is a **state function**: $\Delta U$ depends only on the initial and final states, and $dU$ is an exact differential ($\oint dU = 0$).
   - $q$ and $w$ are **path functions**: their values depend on the process path, and $\delta q$, $\delta w$ are inexact differentials.

4. **Sign Convention** (Chemistry Convention)
   - $q > 0$: heat absorbed by the system.
   - $w > 0$: work done **on** the system.
   - Some textbooks define $w$ as the work done **by** the system, so the convention of a source should be checked before one of its formulas is used.
   - For $PV$ work: $\delta w = -P_{\mathrm{ext}}\,dV$, so compression ($dV < 0$) gives $w > 0$. For a reversible process, $P_{\mathrm{ext}} = P$.

5. **Generalized Work**
   - Many work modes share a generalized force $\times$ generalized displacement structure: $PV$ work ($-P_{\mathrm{ext}}\,dV$), surface work ($\gamma\,dA$), elastic work ($f\,dl$), electrical work ($\mathcal{E}\,dq_{\mathrm{el}}$), and chemical work ($\mu\,dN$).

6. **Free Expansion**
   - Expansion into vacuum: $P_{\mathrm{ext}} = 0$, so $w = 0$. If the container is also insulated, $q = 0$ and $\Delta U = 0$. For an ideal gas, $\Delta T = 0$ as well.
   - Free expansion is the extreme case of an irreversible process: the gas changes state while doing no work and exchanging no heat.

---

### Section 3.2: Applications of the First Law

See [Section 3.2](section-02.md).

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
   4. Specify the equation of state and an expression for $U$ (e.g., ideal gas).
   5. Integrate to find $q$, $w$, $\Delta U$.

4. **First Law with $V$ and $T$ as Independent Variables**
   - General form (reversible process, $PV$-only work):

   ```{math}
   \delta q \;=\; C_V \, dT \;+\; \left[\left(\frac{\partial U}{\partial V}\right)_T + P\right] dV.
   ```

   - $C_V = \left(\frac{\partial U}{\partial T}\right)_V$ is the heat capacity at constant volume.
   - For an ideal gas, $U$ depends only on $T$ ([Section 2.5](../chapter-02/section-05.md)), so $(\partial U/\partial V)_T = 0$, simplifying to $\delta q = C_V\,dT + P\,dV$.

5. **Common Ideal-Gas Processes**
   - **Isochoric** ($dV = 0$): $w = 0$, $q = C_V\,\Delta T$, $\Delta U = C_V\,\Delta T$.
   - **Isothermal** ($dT = 0$): $\Delta U = 0$ (ideal gas) and, for a reversible process, $q = -w = Nk_{\mathrm{B}}T\ln(V_2/V_1)$.
   - **Adiabatic** ($\delta q = 0$): $\Delta U = w = C_V\,\Delta T$ and, for a reversible process, $TV^{\gamma - 1} = \text{const}$ (equivalently $PV^\gamma = \text{const}$), where $\gamma = C_P/C_V$.

6. **Microscopic Interpretation of the First Law**
   - From $U = \sum_i p_i E_i$, for a reversible process in a closed system with $PV$-only work:

   ```{math}
   dU = \underbrace{\sum_i E_i\, dp_i}_{\delta q} + \underbrace{\sum_i p_i\, dE_i}_{\delta w}.
   ```

   - **Heat** changes *which states are occupied* (probabilities shift; energy levels fixed).
   - **Work** changes *the energies of the states* (energy levels shift; probabilities fixed).

---

### Section 3.3: Enthalpy

See [Section 3.3](section-03.md).

1. **Definition and Motivation**
   - **Enthalpy**: $H = U + PV$.
   - At constant pressure with $PV$-only work: $\delta q_P = dH$, so $q_P = \Delta H$.
   - Defines a state function that plays the same role at constant $P$ that $U$ plays at constant $V$.

2. **Heat Capacity at Constant Pressure**

   ```{math}
   C_P = \left(\frac{\partial H}{\partial T}\right)_P
   = \left(\frac{\partial U}{\partial T}\right)_P + P\left(\frac{\partial V}{\partial T}\right)_P.
   ```

   - Compare with the heat capacity at constant volume, $C_V = \left(\frac{\partial U}{\partial T}\right)_V$ ([Section 3.2](section-02.md)).

3. **Condition for $q_P = \Delta H$**
   - The relation holds only when $PV$ work is the sole form of work. The table of work terms in [Section 3.1](section-01.md) lists other kinds of work.

4. **Measuring Enthalpy Changes**
   - For a change of temperature at constant pressure with no reaction or phase change, $\Delta H = \int_{T_i}^{T_f} C_P\,dT \approx C_P\,\Delta T$.
   - **Calorimetry**: a process exchanges heat with a known mass $m$ of water (heat capacity $C_P$ per gram), whose temperature changes by $\Delta T$. The heat absorbed by the process at constant pressure is $q_P = -m\,C_P\,\Delta T = \Delta H$.

5. **Standard States and Formation Enthalpies**
   - Standard pressure: $P^\circ = 1\,\text{bar}$; reference temperature typically 298.15 K.
   - The standard state of an element is its most stable form at $P^\circ$. Phosphorus is the usual exception: its reference state is white phosphorus.
   - Standard enthalpy of reaction: $\Delta_r H^\circ = \sum_p \nu_p H_p^\circ - \sum_r \nu_r H_r^\circ$, where $H_p^\circ$ and $H_r^\circ$ are standard molar enthalpies. Only differences of enthalpy can be measured.
   - Standard enthalpy of formation, $\Delta_f H^\circ$: enthalpy change when 1 mol of a compound in its standard state is formed from its elements in their standard states.
   - By convention, $\Delta_f H^\circ = 0$ for elements in their standard states.

6. **Hess's Law**
   - Because $H$ is a state function, enthalpy changes are path-independent and additive:

   ```{math}
   \Delta_r H^\circ
   = \sum_{p} \nu_p \,\Delta_f H_p^\circ
   - \sum_{r} \nu_r \,\Delta_f H_r^\circ.
   ```

7. **First Law Toolkit Summary**

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

C. **First Law in $(V, T)$ Variables** (for a reversible process with $PV$-only work)

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

- $C_V = \left(\frac{\partial U}{\partial T}\right)_V$ is the heat capacity at constant volume.

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

- **Applicability**: a change of temperature at constant pressure with no reaction or phase change.

---

I. **Hess's Law**

```{math}
\Delta_r H^\circ
= \sum_{p} \nu_p\,\Delta_f H_p^\circ
- \sum_{r} \nu_r\,\Delta_f H_r^\circ.
```

- **Applicability**: any reaction at standard conditions, using tabulated formation enthalpies.

---

J. **Microscopic First Law**

```{math}
dU = \underbrace{\sum_i E_i\, dp_i}_{\delta q}
\;+\; \underbrace{\sum_i p_i\, dE_i}_{\delta w}.
```

- **Applicability**: reversible process in a closed system with $PV$-only work. Heat redistributes probabilities; work shifts energy levels.

## 3. Practice Problems

Each problem is a multi-part question of exam length. Work through it on paper before you open the solution.

### Problem 1. Two paths between the same states

Draws on [Section 3.1](section-01.md) and [Section 3.2](section-02.md).

One mole of nitrogen, treated as an ideal gas with $C_{V,m}=\tfrac{5}{2}R$, starts in state A at $T_A=300\ \mathrm{K}$ and $P_A=1.00\ \mathrm{bar}$. It ends in state B at $P_B=3.00\ \mathrm{bar}$ with half the initial volume, $V_B=V_A/2$. Two reversible paths connect the states:

- **Path 1.** The gas is compressed at constant pressure $P_A$ from $V_A$ to $V_B$ and then heated at constant volume $V_B$ until the pressure is $P_B$.
- **Path 2.** The gas is heated at constant volume $V_A$ until the pressure is $P_B$ and then compressed at constant pressure $P_B$ from $V_A$ to $V_B$.

Use $R=8.3145\ \mathrm{J\,mol^{-1}\,K^{-1}}$.

**(a)** Find the temperature of state B and of the intermediate state of each path.

**(b)** Find $w$, $\Delta U$, and $q$ for each step of Path 1.

**(c)** Find $w$, $\Delta U$, and $q$ for each step of Path 2.

**(d)** Add the steps to obtain $w$, $q$, and $\Delta U$ for each path, in joules. Then calculate $\Delta U$ from the temperatures of A and B alone.

````{admonition} Solution
:class: dropdown

Every answer is written first as a multiple of $nRT_A=(1.00\ \mathrm{mol})(8.3145\ \mathrm{J\,mol^{-1}\,K^{-1}})(300\ \mathrm{K})=2494\ \mathrm{J}$.

1. **Temperatures**

   For a fixed amount of an ideal gas, $T$ is proportional to $PV$. State B has three times the pressure and half the volume of state A, so $T_B=\tfrac{3}{2}T_A=450\ \mathrm{K}$. The intermediate state of Path 1 has $P_A$ and $V_A/2$, so $T_1=\tfrac{1}{2}T_A=150\ \mathrm{K}$. The intermediate state of Path 2 has $3P_A$ and $V_A$, so $T_2=3T_A=900\ \mathrm{K}$.

2. **Relations for the two kinds of step**

   The work done on the gas in a reversible step is $w=-\int P\,dV$, and for an ideal gas $\Delta U=C_V\,\Delta T=\tfrac{5}{2}nR\,\Delta T$ for any step. The heat then follows from the First Law, $q=\Delta U-w$.

   - At constant volume, $w=0$ and $q=\Delta U$.
   - At constant pressure, $w=-P\,\Delta V=-nR\,\Delta T$.

3. **Path 1, compression at constant pressure ($300\ \mathrm{K}\to150\ \mathrm{K}$)**

   ```{math}
   \begin{aligned}
   w&=-nR\,(T_1-T_A)=+\tfrac{1}{2}\,nRT_A,\\
   \Delta U&=\tfrac{5}{2}nR\,(T_1-T_A)=-\tfrac{5}{4}\,nRT_A,\\
   q&=\Delta U-w=-\tfrac{7}{4}\,nRT_A.
   \end{aligned}
   ```

4. **Path 1, heating at constant volume ($150\ \mathrm{K}\to450\ \mathrm{K}$)**

   ```{math}
   w=0,
   \qquad
   q=\Delta U=\tfrac{5}{2}nR\,(T_B-T_1)=+\tfrac{5}{2}\,nRT_A.
   ```

5. **Path 2, heating at constant volume ($300\ \mathrm{K}\to900\ \mathrm{K}$)**

   ```{math}
   w=0,
   \qquad
   q=\Delta U=\tfrac{5}{2}nR\,(T_2-T_A)=+5\,nRT_A.
   ```

6. **Path 2, compression at constant pressure ($900\ \mathrm{K}\to450\ \mathrm{K}$)**

   ```{math}
   \begin{aligned}
   w&=-nR\,(T_B-T_2)=+\tfrac{3}{2}\,nRT_A,\\
   \Delta U&=\tfrac{5}{2}nR\,(T_B-T_2)=-\tfrac{15}{4}\,nRT_A,\\
   q&=\Delta U-w=-\tfrac{21}{4}\,nRT_A.
   \end{aligned}
   ```

7. **Totals for each path**

   | | $w/(nRT_A)$ | $q/(nRT_A)$ | $\Delta U/(nRT_A)$ | $w$ (J) | $q$ (J) | $\Delta U$ (J) |
   |---|---:|---:|---:|---:|---:|---:|
   | Path 1 | $+\tfrac{1}{2}$ | $+\tfrac{3}{4}$ | $+\tfrac{5}{4}$ | +1247 | +1871 | +3118 |
   | Path 2 | $+\tfrac{3}{2}$ | $-\tfrac{1}{4}$ | $+\tfrac{5}{4}$ | +3742 | −624 | +3118 |

8. **$\Delta U$ from the end states**

   ```{math}
   \Delta U=\tfrac{5}{2}nR\,(T_B-T_A)=\tfrac{5}{2}(8.3145\ \mathrm{J\,K^{-1}})(150\ \mathrm{K})=3118\ \mathrm{J},
   ```

   the value obtained by adding $q$ and $w$ along either path.

**Result.** $T_B=450\ \mathrm{K}$, with intermediate temperatures of $150\ \mathrm{K}$ (Path 1) and $900\ \mathrm{K}$ (Path 2). Path 1 has $w=+1247\ \mathrm{J}$ and $q=+1871\ \mathrm{J}$; Path 2 has $w=+3742\ \mathrm{J}$ and $q=-624\ \mathrm{J}$. Both give $\Delta U=+3118\ \mathrm{J}$.
````

### Problem 2. Adiabatic compression in a diesel cylinder

Draws on [Section 3.2](section-02.md) and [Section 3.3](section-03.md).

In the compression stroke of a diesel engine, the air in the cylinder is compressed to $1/18$ of its initial volume before the fuel is injected. Treat the air as a diatomic ideal gas with $C_{V,m}=\tfrac{5}{2}R$, and treat the compression as reversible and adiabatic. The air starts at $T_1=300\ \mathrm{K}$ and $P_1=1.00\ \mathrm{bar}$. Use $R=8.3145\ \mathrm{J\,mol^{-1}\,K^{-1}}$.

**(a)** Follow the workflow of Section 3.2 to show that $T\,V^{2/5}$ is constant along the path.

**(b)** Find the final temperature $T_2$ and the final pressure $P_2$.

**(c)** Find $w$, $\Delta U$, and $\Delta H$ per mole of air.

**(d)** Diesel fuel ignites in hot air without a spark at about $500\ \mathrm{K}$. Does the compression alone bring the air to this temperature?

````{admonition} Solution
:class: dropdown

1. **First Law for an adiabatic step**

   For an ideal gas in a reversible process with only $PV$ work, $\delta q=C_V\,dT+P\,dV$. Setting $\delta q=0$ and using $P=nRT/V$ and $C_V=\tfrac{5}{2}nR$,

   ```{math}
   \frac{5}{2}\,nR\,dT=-\frac{nRT}{V}\,dV
   \quad\Longrightarrow\quad
   \frac{5}{2}\,\frac{dT}{T}=-\frac{dV}{V}.
   ```

2. **Integration**

   ```{math}
   \frac{5}{2}\ln\frac{T_2}{T_1}=-\ln\frac{V_2}{V_1}
   \quad\Longrightarrow\quad
   \frac{T_2}{T_1}=\left(\frac{V_1}{V_2}\right)^{2/5},
   ```

   so $T_1V_1^{2/5}=T_2V_2^{2/5}$.

3. **Final temperature**

   ```{math}
   T_2=T_1\left(\frac{V_1}{V_2}\right)^{2/5}=(300\ \mathrm{K})(18)^{0.4}=(300\ \mathrm{K})(3.178)=953\ \mathrm{K}.
   ```

4. **Final pressure**

   From $PV=nRT$ at the two states,

   ```{math}
   P_2=P_1\,\frac{V_1}{V_2}\,\frac{T_2}{T_1}=(1.00\ \mathrm{bar})(18)(3.178)=57.2\ \mathrm{bar}.
   ```

5. **Work and internal energy**

   Because $q=0$, $w=\Delta U$. With $\Delta T=953.3-300=653.3\ \mathrm{K}$,

   ```{math}
   w=\Delta U=C_{V,m}\,\Delta T=\tfrac{5}{2}(8.3145\ \mathrm{J\,mol^{-1}\,K^{-1}})(653.3\ \mathrm{K})=13.58\ \mathrm{kJ\,mol^{-1}}.
   ```

6. **Enthalpy**

   For an ideal gas, $\Delta H=\Delta U+\Delta(PV)=\Delta U+R\,\Delta T$ per mole, so

   ```{math}
   \begin{aligned}
   \Delta H&=13.58\ \mathrm{kJ\,mol^{-1}}+(8.3145\ \mathrm{J\,mol^{-1}\,K^{-1}})(653.3\ \mathrm{K})\\
   &=13.58+5.43=19.01\ \mathrm{kJ\,mol^{-1}}.
   \end{aligned}
   ```

7. **Ignition**

   The air reaches $953\ \mathrm{K}$, well above $500\ \mathrm{K}$, so fuel injected at the end of the stroke ignites without a spark.

**Result.** $T_2=953\ \mathrm{K}$ and $P_2=57.2\ \mathrm{bar}$. Per mole of air, $w=\Delta U=+13.58\ \mathrm{kJ}$ and $\Delta H=+19.01\ \mathrm{kJ}$. The compression heats the air past the ignition temperature of the fuel.
````

### Problem 3. A reaction enthalpy at the process temperature

Draws on [Section 3.3](section-03.md).

In the contact process for sulfuric acid, sulfur dioxide is oxidized over a catalyst near $700\ \mathrm{K}$:

```{math}
\mathrm{SO_2}(g)+\tfrac{1}{2}\,\mathrm{O_2}(g)\;\longrightarrow\;\mathrm{SO_3}(g).
```

The table lists data from the [NIST-JANAF tables](https://janaf.nist.gov/). Treat each heat capacity as constant between $298.15\ \mathrm{K}$ and $700\ \mathrm{K}$ at its value at $500\ \mathrm{K}$.

| Species | $\Delta_f H^\circ(298.15\ \mathrm{K})$ ($\mathrm{kJ\,mol^{-1}}$) | $C_{P,m}(500\ \mathrm{K})$ ($\mathrm{J\,mol^{-1}\,K^{-1}}$) | $\Delta_f H^\circ(700\ \mathrm{K})$ ($\mathrm{kJ\,mol^{-1}}$) |
|---|---:|---:|---:|
| $\mathrm{SO_2}(g)$ | −296.842 | 46.576 | −306.291 |
| $\mathrm{O_2}(g)$ | 0 | 31.091 | 0 |
| $\mathrm{SO_3}(g)$ | −395.765 | 63.100 | −405.014 |

**(a)** Find $\Delta_r H^\circ$ at $298.15\ \mathrm{K}$.

**(b)** Find $\Delta_r H^\circ$ at $700\ \mathrm{K}$ from a cycle of three steps at constant pressure: cool the reactants from $700\ \mathrm{K}$ to $298.15\ \mathrm{K}$, carry out the reaction at $298.15\ \mathrm{K}$, and heat the product to $700\ \mathrm{K}$.

**(c)** Compare your result with the value that follows from the formation enthalpies at $700\ \mathrm{K}$.

**(d)** How much heat is released when $1.00\ \mathrm{kg}$ of $\mathrm{SO_3}$ ($M=80.063\ \mathrm{g\,mol^{-1}}$) forms at $700\ \mathrm{K}$ and constant pressure?

````{admonition} Solution
:class: dropdown

1. **Reaction enthalpy at 298.15 K**

   From the formation enthalpies,

   ```{math}
   \begin{aligned}
   \Delta_r H^\circ(298.15\ \mathrm{K})&=(-395.765)-\left[(-296.842)+\tfrac{1}{2}(0)\right]\\
   &=-98.923\ \mathrm{kJ\,mol^{-1}}.
   \end{aligned}
   ```

2. **Cooling the reactants**

   The reactants are one mole of $\mathrm{SO_2}$ and half a mole of $\mathrm{O_2}$, with a combined heat capacity of $46.576+\tfrac{1}{2}(31.091)=62.122\ \mathrm{J\,K^{-1}}$. With $\Delta T=298.15-700=-401.85\ \mathrm{K}$,

   ```{math}
   \Delta H_1=(62.122\ \mathrm{J\,K^{-1}})(-401.85\ \mathrm{K})=-24.96\ \mathrm{kJ}.
   ```

3. **Reaction at 298.15 K**

   ```{math}
   \Delta H_2=-98.92\ \mathrm{kJ}.
   ```

4. **Heating the product**

   ```{math}
   \Delta H_3=(63.100\ \mathrm{J\,K^{-1}})(+401.85\ \mathrm{K})=+25.36\ \mathrm{kJ}.
   ```

5. **Reaction enthalpy at 700 K**

   The three steps take the reactants at $700\ \mathrm{K}$ to the product at $700\ \mathrm{K}$, so

   ```{math}
   \begin{aligned}
   \Delta_r H^\circ(700\ \mathrm{K})&=\Delta H_1+\Delta H_2+\Delta H_3\\
   &=-24.96-98.92+25.36=-98.53\ \mathrm{kJ\,mol^{-1}}.
   \end{aligned}
   ```

6. **Comparison with the tabulated value**

   The formation enthalpies at $700\ \mathrm{K}$ give $\Delta_r H^\circ(700\ \mathrm{K})=(-405.014)-(-306.291)=-98.723\ \mathrm{kJ\,mol^{-1}}$. The estimate with constant heat capacities differs from it by $0.19\ \mathrm{kJ\,mol^{-1}}$, or $0.2\%$.

7. **Heat released per kilogram**

   One kilogram of $\mathrm{SO_3}$ is $1000/80.063=12.49\ \mathrm{mol}$. At constant pressure the heat absorbed is $q_P=\Delta H$, so

   ```{math}
   q_P=(12.49\ \mathrm{mol})(-98.53\ \mathrm{kJ\,mol^{-1}})=-1.23\times10^{3}\ \mathrm{kJ}.
   ```

   The reaction releases $1.23\ \mathrm{MJ}$ of heat per kilogram of $\mathrm{SO_3}$.

**Result.** $\Delta_r H^\circ=-98.92\ \mathrm{kJ\,mol^{-1}}$ at $298.15\ \mathrm{K}$ and $-98.53\ \mathrm{kJ\,mol^{-1}}$ at $700\ \mathrm{K}$ from the cycle, against the tabulated $-98.72\ \mathrm{kJ\,mol^{-1}}$. Forming $1.00\ \mathrm{kg}$ of $\mathrm{SO_3}$ at $700\ \mathrm{K}$ releases $1.23\ \mathrm{MJ}$.
````

### Problem 4. $\Delta U$ and $\Delta H$ from a bomb calorimeter

Draws on [Section 3.1](section-01.md) and [Section 3.3](section-03.md).

A bomb calorimeter holds its contents at constant volume. It is calibrated by burning $0.9000\ \mathrm{g}$ of benzoic acid, whose energy of combustion at constant volume is $-26.434\ \mathrm{kJ\,g^{-1}}$; the temperature of the calorimeter rises by $2.379\ \mathrm{K}$. Then $1.500\ \mathrm{g}$ of urea ($M=60.056\ \mathrm{g\,mol^{-1}}$) is burned in the same calorimeter near $298\ \mathrm{K}$, and the temperature rises by $1.582\ \mathrm{K}$. The combustion reaction is

```{math}
\mathrm{CO(NH_2)_2}(s)+\tfrac{3}{2}\,\mathrm{O_2}(g)\;\longrightarrow\;\mathrm{CO_2}(g)+\mathrm{N_2}(g)+2\,\mathrm{H_2O}(l).
```

Use $R=8.3145\ \mathrm{J\,mol^{-1}\,K^{-1}}$, $\Delta_f H^\circ(\mathrm{CO_2},g)=-393.51\ \mathrm{kJ\,mol^{-1}}$, and $\Delta_f H^\circ(\mathrm{H_2O},l)=-285.83\ \mathrm{kJ\,mol^{-1}}$.

**(a)** Find the heat capacity of the calorimeter, $C_{\mathrm{cal}}$.

**(b)** Find $\Delta U$ for the combustion of one mole of urea.

**(c)** Find $\Delta H$ for the combustion of one mole of urea at $298.15\ \mathrm{K}$. Treat the gases as ideal and neglect the volumes of the solid and the liquid.

**(d)** Find $\Delta_f H^\circ$ of solid urea, and compare it with the value in the [NIST Chemistry WebBook](https://webbook.nist.gov/chemistry/), $-333.11\ \mathrm{kJ\,mol^{-1}}$.

````{admonition} Solution
:class: dropdown

1. **Heat capacity of the calorimeter**

   The benzoic acid releases $(0.9000\ \mathrm{g})(26.434\ \mathrm{kJ\,g^{-1}})=23.79\ \mathrm{kJ}$, which the calorimeter absorbs. Therefore

   ```{math}
   C_{\mathrm{cal}}=\frac{23.79\ \mathrm{kJ}}{2.379\ \mathrm{K}}=10.00\ \mathrm{kJ\,K^{-1}}.
   ```

2. **Heat of the urea combustion at constant volume**

   The calorimeter absorbs $(10.00\ \mathrm{kJ\,K^{-1}})(1.582\ \mathrm{K})=15.82\ \mathrm{kJ}$, so the reacting system gives up this much heat: $q_V=-15.82\ \mathrm{kJ}$.

3. **$\Delta U$ per mole**

   The amount of urea is $n=1.500/60.056=0.02498\ \mathrm{mol}$. At constant volume, $q_V=\Delta U$, so

   ```{math}
   \Delta U=\frac{-15.82\ \mathrm{kJ}}{0.02498\ \mathrm{mol}}=-633.3\ \mathrm{kJ\,mol^{-1}}.
   ```

4. **$\Delta H$ per mole**

   From $H=U+PV$, $\Delta H=\Delta U+\Delta(PV)$. With ideal gases and negligible volumes of the solid and the liquid, $\Delta(PV)=\Delta n_{\mathrm{gas}}\,RT$ at constant temperature. The reaction consumes $\tfrac{3}{2}$ mol of gas and forms 2 mol, so $\Delta n_{\mathrm{gas}}=+\tfrac{1}{2}$ and

   ```{math}
   \begin{aligned}
   \Delta H&=\Delta U+\Delta n_{\mathrm{gas}}\,RT\\
   &=-633.3\ \mathrm{kJ\,mol^{-1}}+\tfrac{1}{2}(8.3145\times10^{-3}\ \mathrm{kJ\,mol^{-1}\,K^{-1}})(298.15\ \mathrm{K})\\
   &=-633.3+1.24=-632.1\ \mathrm{kJ\,mol^{-1}}.
   \end{aligned}
   ```

5. **Formation enthalpy of urea**

   For the combustion reaction, $\Delta_r H^\circ=\Delta_f H^\circ(\mathrm{CO_2})+2\,\Delta_f H^\circ(\mathrm{H_2O},l)-\Delta_f H^\circ(\text{urea})$, with zero for $\mathrm{N_2}$ and $\mathrm{O_2}$. Solving for the formation enthalpy of urea,

   ```{math}
   \Delta_f H^\circ(\text{urea})=(-393.51)+2(-285.83)-(-632.1)=-333.1\ \mathrm{kJ\,mol^{-1}},
   ```

   which agrees with the NIST value, $-333.11\ \mathrm{kJ\,mol^{-1}}$.

**Result.** $C_{\mathrm{cal}}=10.00\ \mathrm{kJ\,K^{-1}}$. For the combustion of urea, $\Delta U=-633.3\ \mathrm{kJ\,mol^{-1}}$ and $\Delta H=-632.1\ \mathrm{kJ\,mol^{-1}}$, which gives $\Delta_f H^\circ(\text{urea})=-333.1\ \mathrm{kJ\,mol^{-1}}$.
````
