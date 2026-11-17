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

## Checklist of Key Concepts

### Section 1.1. Course Introduction

- **Thermodynamic Systems**
  - **System vs. Surroundings**
  - **Types**: isolated (no exchange of energy or matter), closed (exchanges energy but not matter), open (exchanges energy and matter)
- **State Variables and Equilibrium**
  - **State variables** (e.g., $P$, $V$, $T$) and other **state functions** (e.g., internal energy $U$)
  - **Path functions** (e.g., heat $q$, work $w$) depend on the path taken between states
  - **Equilibrium** implies no net change in macroscopic properties over time
- **Energy, Heat, and Work**
  - **Kinetic vs. Potential energy**
  - **Heat**: energy transfer due to temperature difference
  - **Work**: energy transfer when a force acts through a distance

### Section 1.2. Kinetic Theory

- **Kinetic Theory Assumptions**
  1. Large number of particles (identical, random motion)
  2. Negligible particle size compared to the average distance between particles
  3. Collisions are elastic (no energy loss)
  4. No long-range interparticle forces (except during collisions)
  5. Newton’s laws (classical mechanics) apply
- **Pressure Origin**
  - Microscopic collisions of particles with container walls
  - Force per unit area results from changes in particle momentum
- **Equipartition Theorem**
  - Average translational kinetic energy per particle $\langle E_{\mathrm{kin}} \rangle = \frac{3}{2} k_{\mathrm{B}} T$
  - Relates $T$ to the mean-square speed

### Section 1.3. Ideal Gases

- **Classical Gas Laws**
  - **Boyle** ($P \propto 1/V$ at constant $T,\,N$)
  - **Charles** ($V \propto T$ at constant $P,\,N$)
  - **Gay-Lussac** ($P \propto T$ at constant $V,\,N$)
  - **Avogadro** ($V \propto N$ at constant $P,\,T$)
- **Ideal Gas Equation**
  - $PV = N k_{\mathrm{B}} T = n R T$
  - Assumes point-like, non-interacting particles
- **Absolute Temperature Scale**
  - Kelvin scale derived by extrapolating volumes to zero at $-273.15\,\text{°C}$

### Section 1.4. Real Gases

- **Deviations from Ideal Behavior**
  - Become significant at high $P$ (high density) or low $T$
  - **Compressibility factor** $Z = \frac{PV}{nRT}$ indicates deviation ($Z=1$ ideal; $Z<1$ net attraction; $Z>1$ net repulsion)
- **van der Waals Equation**
  - Introduces parameters $a$ (accounts for attractions) and $b$ (excluded volume)
  - Shows the characteristic van der Waals “loop” below the critical temperature (associated with liquid–vapor coexistence)
- **Critical Point and Corresponding States**
  - Critical conditions $(T_c, P_c, V_{m,c})$ mark the end of the liquid–vapor coexistence line
  - Reduced variables $T_r = T/T_c,\; P_r = P/P_c,\; V_{m,r} = V_m / V_{m,c}$ collapse different substances onto similar curves (principle of corresponding states)

---

## Checklist of Most Important Equations

1. **Ideal Gas Law**

   $$
   PV = N k_{\mathrm{B}} T = n R T.
   $$
   - **Applicability**: low pressures, low densities, and relatively high temperatures (particles effectively non-interacting).

2. **Pressure from Kinetic Theory**

   $$
   P = \frac{N m \langle v^2 \rangle}{3V}.
   $$
   - **Applicability**: idealized gas of point particles with elastic collisions; derived under kinetic-theory assumptions.

3. **Average Kinetic Energy / Equipartition**

   $$
   \langle E_{\mathrm{kin}} \rangle
   \;=\; \frac{3}{2} k_{\mathrm{B}} T
   \quad\Longleftrightarrow\quad
   \frac{1}{2} m \langle v^2 \rangle
   = \frac{3}{2} k_{\mathrm{B}} T.
   $$
   - **Applicability**: classical (high-temperature) regime, ignoring quantum effects and molecular internal modes.

4. **Root-Mean-Square (rms) Speed**

   $$
   v_{\mathrm{rms}}
   = \sqrt{\frac{3 k_{\mathrm{B}} T}{m}}.
   $$
   - **Applicability**: same assumptions as equipartition (kinetic theory of gases).

5. **Compressibility Factor**

   $$
   Z = \frac{PV}{nRT}.
   $$
   - **Interpretation**:
     - $Z = 1$: ideal gas
     - $Z < 1$: net attractive forces
     - $Z > 1$: net repulsive forces
   - **Applicability**: any real gas to quantify deviation from ideality.

6. **van der Waals Equation of State** (in molar form)

   $$
   \left(P + \frac{a_m}{V_m^2}\right)\;\left(V_m - b_m\right)\;=\;R\,T.
   $$
   - **Applicability**: real gases at moderate deviations from ideality; fails under extreme conditions (very high $P$, near liquefaction, etc.).

7. **Critical-Point Relationships (van der Waals)**

   $$
   V_{m,c} = 3\,b_m,\quad
   P_c = \frac{a_m}{27\,b_m^2},\quad
   T_c = \frac{8\,a_m}{27\,b_m\,R}.
   $$
   - **Interpretation**: defines the critical temperature $T_c$, critical pressure $P_c$, and critical molar volume $V_{m,c}$ for a van der Waals fluid.

8. **Corresponding States** (reduced variables)

   $$
   \left(P_r + \frac{3}{V_{m,r}^2}\right)\;\left(V_{m,r} - \frac{1}{3}\right)\;=\;\frac{8}{3}\,T_r
   \quad
   \text{where}
   \quad
   P_r=\frac{P}{P_c},\quad
   T_r=\frac{T}{T_c},\quad
   V_{m,r}=\frac{V_m}{V_{m,c}}.
   $$
   - **Applicability**: near or above critical conditions for fluids that approximate van der Waals behavior; demonstrates “universal” behavior across substances when scaled by critical properties.

## Practice Problems

Each problem is a multi-part question of exam length. Work through it on paper before you open the solution.

### Problem 1. Identifying a noble gas from its molecular speed

Draws on [Section 1.2](section-02.md) and [Section 1.3](section-03.md).

A sample of a noble gas contains $n=0.250\ \mathrm{mol}$ in a container of volume $V=5.00\ \mathrm{L}$ at $P=1.50\ \mathrm{bar}$. The rms speed of its atoms is $v_{\mathrm{rms}}=668\ \mathrm{m/s}$. Treat the gas as ideal.

**(a)** Use the kinetic-theory expression for the pressure to find the molar mass $M$ of the gas, and identify the gas.

**(b)** Find the temperature of the sample in kelvin and in degrees Celsius.

**(c)** Find the total translational kinetic energy $E_{\mathrm{kin}}$ of the sample and show that it equals $\tfrac{3}{2}PV$. Would a sample of argon with the same $P$ and $V$ have a different $E_{\mathrm{kin}}$?

````{admonition} Solution
:class: dropdown

Use $R=8.3145\ \mathrm{J\,mol^{-1}\,K^{-1}}$, $P=1.50\times10^{5}\ \mathrm{Pa}$, and $V=5.00\times10^{-3}\ \mathrm{m^3}$, so that $PV=750\ \mathrm{J}$.

1. **Molar mass from the pressure**

   With $N=nN_{\mathrm{A}}$, $M=N_{\mathrm{A}}m$, and $\langle v^2\rangle=v_{\mathrm{rms}}^2$, the kinetic-theory pressure $P=Nm\langle v^2\rangle/(3V)$ becomes $P=nMv_{\mathrm{rms}}^2/(3V)$. Solving for $M$,

   ```{math}
   M=\frac{3PV}{n\,v_{\mathrm{rms}}^2}
   =\frac{3(750\ \mathrm{J})}{(0.250\ \mathrm{mol})(668\ \mathrm{m/s})^2}
   =0.02017\ \mathrm{kg/mol}
   =20.2\ \mathrm{g/mol}.
   ```

   The noble gas with this molar mass is neon ($20.18\ \mathrm{g/mol}$).

2. **Temperature**

   ```{math}
   T=\frac{PV}{nR}
   =\frac{750\ \mathrm{J}}{(0.250\ \mathrm{mol})(8.3145\ \mathrm{J\,mol^{-1}\,K^{-1}})}
   =360.8\ \mathrm{K}.
   ```

   On the Celsius scale, $T({}^{\circ}\mathrm{C})=360.8-273.15=87.7$.

3. **Total translational kinetic energy**

   ```{math}
   \begin{aligned}
   E_{\mathrm{kin}}&=\frac{1}{2}Nm\langle v^2\rangle
   =\frac{1}{2}nMv_{\mathrm{rms}}^2\\
   &=\frac{1}{2}(0.250\ \mathrm{mol})(0.02017\ \mathrm{kg/mol})(668\ \mathrm{m/s})^2
   =1125\ \mathrm{J}.
   \end{aligned}
   ```

4. **Comparison with $\tfrac{3}{2}PV$**

   Substituting $nMv_{\mathrm{rms}}^2=3PV$ from step 1 gives

   ```{math}
   E_{\mathrm{kin}}=\frac{1}{2}\,(3PV)=\frac{3}{2}PV=\frac{3}{2}(750\ \mathrm{J})=1125\ \mathrm{J},
   ```

   which is also $\tfrac{3}{2}nRT$. The expression $\tfrac{3}{2}PV$ contains neither the mass of the atoms nor their number, so a sample of argon with the same $P$ and $V$ has the same $E_{\mathrm{kin}}$.

**Result.** The gas is neon at $360.8\ \mathrm{K}$ ($87.7\,{}^{\circ}\mathrm{C}$), and $E_{\mathrm{kin}}=\tfrac{3}{2}PV=1125\ \mathrm{J}$ for any ideal gas with this pressure and volume.
````

### Problem 2. A small change of state from the total differential

Draws on [Section 1.3](section-03.md) and [Section 1.4](section-04.md).

Carbon dioxide is held at $T=320\ \mathrm{K}$, above its critical temperature, with molar volume $V_\text{m}=0.300\ \mathrm{L\,mol^{-1}}$. Model it with the van der Waals equation of state solved for the pressure,

```{math}
P=\frac{RT}{V_\text{m}-b_\text{m}}-\frac{a_\text{m}}{V_\text{m}^2},
```

with $a_\text{m}=3.6551\ \mathrm{L^2\,bar\,mol^{-2}}$ and $b_\text{m}=0.042816\ \mathrm{L\,mol^{-1}}$ (the values used in Section 1.4) and $R=0.083145\ \mathrm{L\,bar\,mol^{-1}\,K^{-1}}$.

**(a)** Write the total differential $dP$ in terms of $dT$ and $dV_\text{m}$.

**(b)** The gas is heated from $320\ \mathrm{K}$ to $330\ \mathrm{K}$ at constant $V_\text{m}$. Integrate $dP$ to find the change in pressure, and compare it with the change for an ideal gas with the same $V_\text{m}$.

**(c)** Evaluate $(\partial P/\partial V_\text{m})_T$ at the initial state. What does its sign say about the isotherm at this state?

````{admonition} Solution
:class: dropdown

1. **Partial derivatives**

   ```{math}
   \left(\frac{\partial P}{\partial T}\right)_{V_\text{m}}=\frac{R}{V_\text{m}-b_\text{m}},
   \qquad
   \left(\frac{\partial P}{\partial V_\text{m}}\right)_{T}=-\frac{RT}{(V_\text{m}-b_\text{m})^2}+\frac{2a_\text{m}}{V_\text{m}^3}.
   ```

2. **Total differential**

   ```{math}
   dP=\frac{R}{V_\text{m}-b_\text{m}}\,dT
   +\left[-\frac{RT}{(V_\text{m}-b_\text{m})^2}+\frac{2a_\text{m}}{V_\text{m}^3}\right]dV_\text{m}.
   ```

3. **Heating at constant $V_\text{m}$**

   With $dV_\text{m}=0$, only the first term remains, and its coefficient does not depend on $T$. Integrating from $T_i=320\ \mathrm{K}$ to $T_f=330\ \mathrm{K}$,

   ```{math}
   \Delta P=\frac{R\,(T_f-T_i)}{V_\text{m}-b_\text{m}}
   =\frac{(0.083145)(10)}{0.300-0.042816}
   =\frac{0.83145}{0.257184}
   =3.23\ \mathrm{bar}.
   ```

   As a check, the equation of state gives $P_i=103.45-40.61=62.84\ \mathrm{bar}$ and $P_f=106.69-40.61=66.07\ \mathrm{bar}$, which differ by $3.23\ \mathrm{bar}$.

4. **Comparison with an ideal gas**

   For an ideal gas, $P=RT/V_\text{m}$, so at constant $V_\text{m}$

   ```{math}
   \Delta P_{\text{ideal}}=\frac{R\,(T_f-T_i)}{V_\text{m}}=\frac{0.83145}{0.300}=2.77\ \mathrm{bar}.
   ```

   The van der Waals value is larger by the factor $V_\text{m}/(V_\text{m}-b_\text{m})=1.17$. The term $a_\text{m}/V_\text{m}^2$ does not depend on $T$, so it drops out of the change at constant $V_\text{m}$.

5. **Slope of the isotherm at the initial state**

   ```{math}
   \begin{aligned}
   \left(\frac{\partial P}{\partial V_\text{m}}\right)_{T}
   &=-\frac{(0.083145)(320)}{(0.257184)^2}+\frac{2(3.6551)}{(0.300)^3}\\
   &=-402.3+270.7
   =-131.5\ \mathrm{bar\,mol\,L^{-1}}.
   \end{aligned}
   ```

   The slope is negative, so along the $320\ \mathrm{K}$ isotherm the pressure falls as the molar volume increases at this state, and a small compression raises the pressure.

**Result.** Heating by $10\ \mathrm{K}$ at constant $V_\text{m}$ raises the van der Waals pressure by $3.23\ \mathrm{bar}$ (from $62.84$ to $66.07\ \mathrm{bar}$), compared with $2.77\ \mathrm{bar}$ for an ideal gas, and $(\partial P/\partial V_\text{m})_T=-131.5\ \mathrm{bar\,mol\,L^{-1}}$ at the initial state.
````

### Problem 3. Ethane against the van der Waals model

Draws on [Section 1.4](section-04.md).

At $T=350\ \mathrm{K}$ and $P=55.0\ \mathrm{bar}$, the measured molar volume of ethane is $V_\text{m}=0.3761\ \mathrm{L\,mol^{-1}}$ ([NIST Chemistry WebBook](https://webbook.nist.gov/chemistry/fluid/)). The van der Waals parameters of ethane, determined from its measured critical temperature and pressure, are $a_\text{m}=5.580\ \mathrm{L^2\,bar\,mol^{-2}}$ and $b_\text{m}=0.06513\ \mathrm{L\,mol^{-1}}$. Use $R=0.083145\ \mathrm{L\,bar\,mol^{-1}\,K^{-1}}$.

**(a)** Compute the compressibility factor $Z$ from the measured values. Do attractive or repulsive interactions dominate at this state?

**(b)** Compute the pressure that the van der Waals equation predicts at the measured $T$ and $V_\text{m}$, its percent error against the measured pressure, and the corresponding $Z$.

**(c)** Compute the van der Waals critical constants of ethane, then $T_r$ and $V_{\text{m},r}$ for this state. Obtain $P_r$ from the reduced equation of state and compare it with $P/P_c$ from part (b).

````{admonition} Solution
:class: dropdown

1. **Measured compressibility factor**

   ```{math}
   Z=\frac{PV_\text{m}}{RT}
   =\frac{(55.0)(0.3761)}{(0.083145)(350)}
   =\frac{20.69}{29.10}
   =0.711.
   ```

   Because $Z<1$, attractive interactions dominate at this state.

2. **van der Waals pressure**

   ```{math}
   \begin{aligned}
   P&=\frac{RT}{V_\text{m}-b_\text{m}}-\frac{a_\text{m}}{V_\text{m}^2}
   =\frac{29.10}{0.3761-0.06513}-\frac{5.580}{(0.3761)^2}\\
   &=93.58-39.45
   =54.13\ \mathrm{bar}.
   \end{aligned}
   ```

   The percent error is $(54.13-55.0)/55.0=-1.6\%$, and the corresponding compressibility factor is

   ```{math}
   Z_{\text{vdW}}=\frac{(54.13)(0.3761)}{29.10}=0.700.
   ```

3. **van der Waals critical constants**

   ```{math}
   \begin{aligned}
   V_{\text{m},c}&=3b_\text{m}=3(0.06513)=0.1954\ \mathrm{L\,mol^{-1}},\\
   P_c&=\frac{a_\text{m}}{27b_\text{m}^2}=\frac{5.580}{27(0.06513)^2}=48.72\ \mathrm{bar},\\
   T_c&=\frac{8a_\text{m}}{27b_\text{m}R}=\frac{8(5.580)}{27(0.06513)(0.083145)}=305.3\ \mathrm{K}.
   \end{aligned}
   ```

4. **Reduced variables**

   ```{math}
   T_r=\frac{350}{305.3}=1.146,
   \qquad
   V_{\text{m},r}=\frac{0.3761}{0.1954}=1.925.
   ```

5. **Reduced pressure**

   Solving the reduced equation of state for $P_r$,

   ```{math}
   \begin{aligned}
   P_r&=\frac{8T_r}{3V_{\text{m},r}-1}-\frac{3}{V_{\text{m},r}^2}
   =\frac{8(1.146)}{3(1.925)-1}-\frac{3}{(1.925)^2}\\
   &=1.921-0.810
   =1.111.
   \end{aligned}
   ```

   From part (b), $P/P_c=54.13/48.72=1.111$, in agreement.

**Result.** The measured state has $Z=0.711$. The van der Waals equation gives $P=54.13\ \mathrm{bar}$ ($1.6\%$ below the measured $55.0\ \mathrm{bar}$) and $Z=0.700$, and the state corresponds to $T_r=1.146$, $V_{\text{m},r}=1.925$, and $P_r=1.111$.
````

### Problem 4. Absolute zero from an isobar

Draws on [Section 1.1](section-01.md) and [Section 1.3](section-03.md).

The table lists the molar volume of neon at $P=0.0100\ \mathrm{bar}$ at three temperatures ([NIST Chemistry WebBook](https://webbook.nist.gov/chemistry/fluid/)).

| $T\ ({}^{\circ}\mathrm{C})$ | $V_\text{m}\ (\mathrm{m^3\,mol^{-1}})$ |
|---:|---:|
| $-100.00$ | 1.4397 |
| 0.00 | 2.2711 |
| 100.00 | 3.1026 |

**(a)** Find the slope and the intercept of the straight line $V_\text{m}=c_0+c_1\,T({}^{\circ}\mathrm{C})$ through the first and last rows, confirm that the middle row lies on this line, and find the Celsius temperature at which the extrapolated volume vanishes.

**(b)** Use the slope and the pressure to obtain the molar gas constant $R$.

**(c)** Find the number density $N/V$ of the gas at $0.00\,{}^{\circ}\mathrm{C}$.

````{admonition} Solution
:class: dropdown

1. **Slope**

   ```{math}
   c_1=\frac{3.1026-1.4397}{100.00-(-100.00)}
   =\frac{1.6629}{200.00}
   =8.3145\times10^{-3}\ \mathrm{m^3\,mol^{-1}\,{}^{\circ}C^{-1}}.
   ```

2. **Intercept and check**

   At $T({}^{\circ}\mathrm{C})=0$, halfway between the two rows, the line gives

   ```{math}
   c_0=\frac{1.4397+3.1026}{2}=2.2712\ \mathrm{m^3\,mol^{-1}},
   ```

   which matches the middle row ($2.2711$) to within one unit in the last digit, so the three points lie on one line.

3. **Extrapolation to zero volume**

   Setting $V_\text{m}=0$,

   ```{math}
   T({}^{\circ}\mathrm{C})=-\frac{c_0}{c_1}=-\frac{2.2712}{8.3145\times10^{-3}}=-273.16.
   ```

   The extrapolated volume vanishes at $-273.16\,{}^{\circ}\mathrm{C}$, within $0.01\,{}^{\circ}\mathrm{C}$ of the accepted $-273.15\,{}^{\circ}\mathrm{C}$.

4. **Gas constant**

   A temperature step of one degree Celsius equals a step of one kelvin, and $V_\text{m}=RT/P$ at constant $P$ has slope $R/P$. With $P=0.0100\ \mathrm{bar}=1.00\times10^{3}\ \mathrm{Pa}$,

   ```{math}
   R=c_1P=(8.3145\times10^{-3})(1.00\times10^{3})
   =8.31\ \mathrm{J\,mol^{-1}\,K^{-1}}.
   ```

5. **Number density at $0.00\,{}^{\circ}\mathrm{C}$**

   One mole contains $N_{\mathrm{A}}$ atoms in the volume $V_\text{m}$:

   ```{math}
   \frac{N}{V}=\frac{N_{\mathrm{A}}}{V_\text{m}}
   =\frac{6.022\times10^{23}\ \mathrm{mol^{-1}}}{2.2711\ \mathrm{m^3\,mol^{-1}}}
   =2.65\times10^{23}\ \mathrm{m^{-3}}.
   ```

**Result.** The line through the data has slope $8.3145\times10^{-3}\ \mathrm{m^3\,mol^{-1}\,K^{-1}}$ and reaches zero volume at $-273.16\,{}^{\circ}\mathrm{C}$. The slope gives $R=8.31\ \mathrm{J\,mol^{-1}\,K^{-1}}$, and the gas contains $2.65\times10^{23}$ atoms per cubic meter at $0.00\,{}^{\circ}\mathrm{C}$.
````
