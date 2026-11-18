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
  src="https://wexlergroup.github.io/chem-4020-5020/computational-studios/course-map.html?module=2&amp;here=review"
  title="Course map, opened at Module 2"
  style="width: 100%; height: 830px; border: 0;"
  loading="lazy"
></iframe>

```{image} ../_static/course-map/module-2.png
:class: hidden
:alt: Course map with Module 2 unfolded. Statistical mechanics (Module 2) forms the microscopic strand, and Modules 1 and 3 to 7 form the macroscopic strand. A relation joins the two strands at each module except Module 6. Module 2 is the microscopic strand itself. Below the map, one card for each of Sections 2.1 to 2.8 lists the microscopic content, the connecting relations, and the macroscopic content of the section.
```

## 1. Checklist of Key Concepts

### Section 2.1: Introduction to Statistical Mechanics

See [Section 2.1](section-01.md).

1. **Microscopic–Macroscopic Connection**
   - Macroscopic (thermodynamic) properties can be understood as **statistical averages** of microscopic properties.
   - “Expected value” (ensemble average) is the central idea:

   ```{math}
   \langle X \rangle \;=\; \sum_i X_i\,p_i.
   ```

2. **Arithmetic Average vs. Expected Value**
   - Arithmetic average of a dataset of $n$ values:

   ```{math}
   \bar{X} \;=\; \frac{1}{n}\sum_{j=1}^n X_j.
   ```

   - Expected value:

   ```{math}
   \langle X \rangle \;=\; \sum_i X_i\,p_i,
   \quad p_i = \text{probability of microstate } i.
   ```

3. **Microstates and Ensembles**
   - **Microcanonical ensemble** $(N, V, E)$: system is isolated, fixed total energy $E$.
   - **Canonical ensemble** $(N, V, T)$: system in thermal contact with a reservoir at temperature $T$.
   - **Grand canonical ensemble** $(\mu, V, T)$: system can exchange both energy and particles with a reservoir.

4. **Fundamental Postulate (Microcanonical)**
   - For an isolated system, **each accessible microstate is equally probable**: $p_i = 1/M$, where $M$ is the number of accessible microstates.

---

### Section 2.2: Canonical Ensemble

See [Section 2.2](section-02.md).

1. **Closed System**
   - Exchanges energy (heat) with surroundings; no exchange of matter.

2. **Boltzmann Factor and Partition Function**
   - Probability of microstate $i$:

   ```{math}
   p_i \;=\; \frac{e^{-\beta E_i}}{Q},
   \quad \beta \;=\;\frac{1}{k_{\mathrm B}\,T}.
   ```

   - **Partition function** $Q$:

   ```{math}
   Q \;=\; \sum_{i} e^{-\beta E_i},
   ```

   which normalizes probabilities.

3. **Two-State System**
   - A simple example with energies $E_1$ and $E_2$.
   - Partition function:

   ```{math}
   Q \;=\; e^{-\beta E_1} + e^{-\beta E_2}.
   ```

   - Probabilities (with $\Delta E = E_2 - E_1 > 0$):

   ```{math}
   p_1 = \frac{1}{1 + e^{-\beta\,\Delta E}},
   \quad
   p_2 = 1 - p_1.
   ```

4. **Interpretation of $Q$**
   - $Q$ is like an “effective count” of accessible microstates.
   - At low $T$, only the lowest energy states matter; at high $T$, many states are accessible.

---

### Section 2.3: Ensemble Averages

See [Section 2.3](section-03.md).

1. **Internal Energy**
   - $U \equiv \langle E \rangle$.
   - In the canonical ensemble:

   ```{math}
   U
   \;=\;
   \frac{1}{Q} \sum_i E_i\, e^{-\beta E_i}
   \;=\;
   -\,\bigl(\tfrac{\partial \ln Q}{\partial \beta}\bigr)_{N,V}.
   ```

2. **Heat Capacity $(C_V)$**
   - Measures how internal energy changes with temperature:

   ```{math}
   C_V = \bigl(\tfrac{\partial U}{\partial T}\bigr)_{N,V}.
   ```

   - Also related to **energy fluctuations**:

   ```{math}
   \sigma_E^2
   \;=\;
   \langle (E - \langle E\rangle)^2\rangle
   \;=\;
   k_{\mathrm B}\,T^2\,C_V.
   ```

3. **Pressure** (brief introduction)
   - In the canonical ensemble:

   ```{math}
   P
   \;=\;
   k_{\mathrm B}\,T\,
   \bigl(\tfrac{\partial\ln Q}{\partial V}\bigr)_{N,T}.
   ```

---

### Section 2.4: Molecular Partition Functions

See [Section 2.4](section-04.md).

1. **Many Identical and Independent Particles**
   - For $N$ independent **distinguishable** particles with one-particle partition function $q$, the total partition function is

   ```{math}
   Q \;=\; q^N.
   ```

   - For independent **indistinguishable** particles, when the thermally accessible one-particle states far outnumber the particles,

   ```{math}
   Q \;=\; \frac{q^N}{N!}.
   ```

2. **Molecular Partition Function**
   - Within the Born–Oppenheimer, rigid-rotor, and harmonic-oscillator approximations, a molecule’s energy decomposes into translational, rotational, vibrational, and electronic contributions:

   ```{math}
   \varepsilon = \varepsilon_{\mathrm{trans}} + \varepsilon_{\mathrm{rot}} + \varepsilon_{\mathrm{vib}} + \varepsilon_{\mathrm{elec}}.
   ```

   - Thus,

   ```{math}
   q
   \;=\;
   q_{\mathrm{trans}}\,
   q_{\mathrm{rot}}\,
   q_{\mathrm{vib}}\,
   q_{\mathrm{elec}}.
   ```

---

### Section 2.5: Particle in a Box

See [Section 2.5](section-05.md).

1. **Quantum Levels**
   - For a 1D box $(0 \le x \le L)$,

   ```{math}
   E_n = \frac{h^2}{8mL^2}\,n^2,\;\; n=1,2,3,\dots
   ```

   - In 3D (a rectangular box of sides $L_x, L_y, L_z$):

   ```{math}
   E_{n_x,n_y,n_z}
   \;=\;
   \frac{h^2}{8m}
   \left(\frac{n_x^2}{L_x^2} + \frac{n_y^2}{L_y^2} + \frac{n_z^2}{L_z^2}\right).
   ```

2. **Partition Function (3D Box)**
   - For a cubic box of volume $V = L^3$:

   ```{math}
   q
   \;=\;
   \sum_{n_x=1}^\infty \sum_{n_y=1}^\infty \sum_{n_z=1}^\infty
   \exp\!\Bigl[
    -\beta \frac{h^2}{8mL^2}\,(n_x^2 + n_y^2 + n_z^2)
    \Bigr].
    ```

   - When $\alpha = h^2/(8mV^{2/3}k_{\mathrm B}T) \ll 1$ (large $T$ or large $V$), we approximate the sum by an integral and obtain the **classical** partition function:

   ```{math}
   q_{\mathrm{trans}}
   \;=\;
   \frac{V}{\Lambda^3},
   \quad
   \Lambda
   \;=\;
   \sqrt{\frac{h^2}{2\pi m k_{\mathrm B} T}}\,\,
   \text{(thermal de Broglie wavelength).}
   ```

3. **Many-Particle System (Ideal Gas)**
   - For $N$ indistinguishable, non-interacting particles, when $\Lambda \ll (V/N)^{1/3}$, the mean spacing between particles:

   ```{math}
   Q
   \;=\;
   \frac{q^N}{N!}
   \;=\;
   \frac{1}{N!}
   \biggl(\frac{V}{\Lambda^3}\biggr)^{\!N}.
   ```

   - Leads directly to the ideal gas law and to:

   ```{math}
   U \;=\; \frac{3}{2}\,N\,k_{\mathrm B}T,
   \quad
   C_V \;=\; \frac{3}{2}\,N\,k_{\mathrm B},
   \quad
   P = \frac{N\,k_{\mathrm B}\,T}{V}.
   ```

---

### Section 2.6: Harmonic Oscillator

See [Section 2.6](section-06.md).

1. **Quantum Harmonic Oscillator**
   - Energy levels for a 1D harmonic oscillator of frequency $\omega$:

   ```{math}
   E_n
   \;=\;
   \hbar\omega\,\Bigl(n + \tfrac{1}{2}\Bigr),
   \quad
   n=0,1,2,\dots
   ```

   - The lowest level, $n = 0$, has the **zero-point energy** $E_0 = \hbar\omega/2$.

2. **Partition Function**
   - Summing over these energy levels:

   ```{math}
   q
   \;=\;
   \sum_{n=0}^{\infty} e^{-\beta E_n}
   \;=\;
   \frac{e^{-\frac{1}{2}\,\beta\,\hbar\omega}}{1 - e^{-\beta\,\hbar\omega}}.
   ```

   - Define the **vibrational temperature** $\Theta_{\mathrm{vib}} = \hbar\omega/k_{\mathrm B}$, so that $\beta\hbar\omega = \Theta_{\mathrm{vib}}/T$.

3. **Ensemble Averages**
   - **Internal Energy**:

   ```{math}
   U
   \;=\;
   \frac{\hbar\omega}{2}
   \;+\;
   \frac{\hbar\omega}{\,e^{\beta\,\hbar\omega}-1\,}.
   ```

   - **Heat Capacity**:

   ```{math}
   C_V
   \;=\;
   k_{\mathrm B}
   \Bigl(\frac{\hbar\omega}{k_{\mathrm B} T}\Bigr)^{2}\,
   \frac{\,e^{\hbar\omega/(k_{\mathrm B}T)}\,}
        {\bigl(e^{\hbar\omega/(k_{\mathrm B}T)} - 1\bigr)^{2}}.
   ```

   - At low temperature ($T \ll \Theta_{\mathrm{vib}}$), $U \to \hbar\omega/2$. At high temperature ($T \gg \Theta_{\mathrm{vib}}$), $U \to k_{\mathrm B}T$. For the limits of $C_V$, see the prompt "Limits of the heat capacity" in [Section 2.6](section-06.md).

4. **Equipartition Theorem**
   - In classical statistical mechanics, each quadratic term in the energy contributes $\tfrac{1}{2}k_{\mathrm B}T$ to the mean energy.

5. **Einstein Model and Dulong–Petit Law**
   - The Einstein model treats a crystal of $N$ atoms as $3N$ independent, distinguishable one-dimensional harmonic oscillators with the same $\omega$, so that

   ```{math}
   Q \;=\; q^{3N},
   ```

     and $U$ and $C_V$ are $3N$ times the results for a single oscillator.
   - The Dulong–Petit law is the empirical observation that the molar heat capacity of many elemental solids at room temperature is close to $3R$.

---

### Section 2.7: Linear Rigid Rotor

See [Section 2.7](section-07.md).

1. **Energy Levels**
   - For a rigid, linear rotor of moment of inertia $I$:

   ```{math}
   E_J
   \;=\;
   \frac{\hbar^2}{2I}\,J\,(J+1),
   \quad
   J=0,1,2,\dots
   ```

   - Each level $E_J$ has degeneracy $g_J = 2J + 1$.

2. **Rotational Partition Function**
   - Exact form:

   ```{math}
   q_{\mathrm{rot}}
   \;=\;
   \sum_{J=0}^\infty
   (2J + 1)\, e^{-\beta \,\frac{\hbar^2}{2I}\,J\,(J+1)}.
   ```

   - Define the **rotational temperature** $\Theta_{\mathrm{rot}} = \frac{\hbar^2}{2k_{\mathrm B}I}$, so that $E_J = k_{\mathrm B}\Theta_{\mathrm{rot}}\,J(J+1)$.
   - **High-temperature limit** ($T \gg \Theta_{\mathrm{rot}}$) gives

   ```{math}
   q_{\mathrm{rot}}
   \;\approx\;
   \frac{T}{\Theta_{\mathrm{rot}}}
   \quad (\text{for a heteronuclear diatomic}).
   ```

   - For a **homonuclear** diatomic, include a symmetry number $\sigma=2$, so

   ```{math}
   q_{\mathrm{rot}}
   \;\approx\;
   \frac{T}{\sigma\,\Theta_{\mathrm{rot}}}.
   ```

3. **Ensemble Averages** (High-$T$ Approximation)
   - **Internal Energy** ($U_{\mathrm{rot}}$):

   ```{math}
   U_{\mathrm{rot}} \;\approx\; k_{\mathrm B}\,T.
   ```

     (To compare this result with the equipartition theorem, see the prompt "Comparing with equipartition" in [Section 2.7](section-07.md).)
   - **Heat Capacity**:

   ```{math}
   C_V^{(\mathrm{rot})} \;\approx\; k_{\mathrm B}.
   ```

---

### Section 2.8: Molecular Statistical Mechanics

See [Section 2.8](section-08.md).

1. **Combining All Degrees of Freedom**
   - For a general molecule, the total one-molecule partition function is

   ```{math}
   q
   \;=\;
   q_{\mathrm{trans}}\,
   q_{\mathrm{rot}}\,
   q_{\mathrm{vib}}\,
   q_{\mathrm{elec}}.
   ```

   - For a polyatomic molecule, $q_{\mathrm{rot}}$ contains up to three rotational temperatures, $\Theta_{\mathrm{rot},A}$, $\Theta_{\mathrm{rot},B}$, and $\Theta_{\mathrm{rot},C}$, and $q_{\mathrm{vib}}$ contains one vibrational temperature $\Theta_{\mathrm{vib},j}$ for each normal mode.

2. **Rotation of Nonlinear Molecules**
   - A nonlinear molecule has three principal moments of inertia, $I_A$, $I_B$, and $I_C$, each with a rotational temperature $\Theta_{\mathrm{rot},i} = \hbar^2/(2k_{\mathrm B}I_i)$.
   - In a **spherical top**, all three are equal; in a **symmetric top**, two are equal; in an **asymmetric top**, all three are different.
   - High-temperature limit ($T \gg \Theta_{\mathrm{rot},i}$):

   ```{math}
   q_{\mathrm{rot}}
   \;\approx\;
   \frac{\sqrt{\pi}}{\sigma}\,
   \sqrt{\frac{T^3}{\Theta_{\mathrm{rot},A}\,\Theta_{\mathrm{rot},B}\,\Theta_{\mathrm{rot},C}}}.
   ```

3. **Symmetry Number**
   - The symmetry number $\sigma$ is the number of indistinguishable orientations of the molecule that can be reached by rotating it as a rigid body, including the original orientation (for example, 2 for $\mathrm{H_2O}$, 3 for $\mathrm{NH_3}$, and 12 for $\mathrm{CH_4}$).

4. **Vibration of Polyatomic Molecules**
   - A molecule with $N$ atoms has $3N-6$ normal modes if it is nonlinear and $3N-5$ if it is linear.
   - In the harmonic approximation,

   ```{math}
   q_{\mathrm{vib}}
   \;=\;
   \prod_{j} \frac{e^{-\Theta_{\mathrm{vib},j}/(2T)}}{1 - e^{-\Theta_{\mathrm{vib},j}/T}},
   ```

     where the product runs over the normal modes.

5. **Electronic Partition Function**
   - With energies measured from the ground level ($\varepsilon_1 = 0$),

   ```{math}
   q_{\mathrm{elec}}
   \;=\;
   \sum_{i} g_i\,e^{-\beta\varepsilon_i}
   \;=\;
   g_1 + g_2\,e^{-\beta\varepsilon_2} + \cdots,
   ```

     where $g_i$ is the degeneracy of the $i$-th electronic level. If every term after the first is negligible, $q_{\mathrm{elec}} \approx g_1$.

6. **Summary Table**
   - The summary table in [Section 2.8](section-08.md) collects $q_{\mathrm{trans}}$, $q_{\mathrm{rot}}$, $q_{\mathrm{vib}}$, and $q_{\mathrm{elec}}$ for atoms, diatomic molecules, and linear and nonlinear polyatomic molecules.

---

## 2. Checklist of Most Important Equations

Below is a **unified list** of the major equations from Sections 2.1–2.8.

A. **Expected Value of a General Variable**

```{math}
\langle X \rangle
\;=\;
\sum_i X_i\,p_i.
```

- **Applicability**: general definition in statistical mechanics/probability theory.

---

B. **Microcanonical Ensemble (Fundamental Postulate)**

```{math}
p_i
\;=\;
\frac{1}{M}
\quad(\text{for all accessible microstates}).
```

- **Applicability**: isolated system, $(N, V, E)$, with $M$ accessible microstates.

---

C. **Canonical Ensemble Probability**

```{math}
p_i
\;=\;
\frac{\,e^{-\beta\,E_i}\,}{Q},
\quad
\beta = \frac{1}{k_{\mathrm B}\,T}.
```

- **Applicability**: closed system in thermal contact at $T$.

---

D. **Canonical Partition Function**

```{math}
Q
\;=\;
\sum_i e^{-\beta\,E_i}.
```

- **Applicability**: $(N, V, T)$ ensemble; sum/integral over all microstates.

---

E. **Internal Energy (Canonical)**

```{math}
U
\;=\;
\langle E\rangle
\;=\;
\frac{1}{Q}\,\sum_i E_i\,e^{-\beta E_i}
\;=\;
-\bigl(\tfrac{\partial \ln Q}{\partial \beta}\bigr)_{N,V}.
```

---

F. **Heat Capacity at Constant Volume**

```{math}
C_V
\;=\;
\bigl(\tfrac{\partial U}{\partial T}\bigr)_{N,V}.
```

Also,

```{math}
\sigma_E^2
\;=\;
k_{\mathrm B}\,T^2\,C_V.
```

---

G. **Pressure (Canonical)**

```{math}
P
\;=\;
k_{\mathrm B}\,T \,\bigl(\tfrac{\partial \ln Q}{\partial V}\bigr)_{N,T}.
```

---

H. **Molecular Systems: Indistinguishability Factor**

```{math}
Q
\;=\;
\frac{q^N}{N!},
```

- **Applicability**: identical, independent, indistinguishable particles whose thermally accessible one-particle states far outnumber the particles (e.g., gases at ordinary temperatures and pressures). For independent distinguishable particles, $Q = q^N$.

---

I. **Translational Partition Function** (Particle in a 3D box at high $T$)

```{math}
q_{\mathrm{trans}}
\;=\;
\frac{V}{\Lambda^3},
\quad
\Lambda
=
\sqrt{\frac{h^2}{2\pi m k_{\mathrm B}T}}.
```

---

J. **Harmonic Oscillator Partition Function (1D)**

```{math}
q_{\mathrm{HO}}
\;=\;
\frac{\,e^{-\tfrac{1}{2}\beta\hbar\omega}\,}{1 - e^{-\beta\,\hbar\omega}}.
```

- **Internal Energy**:

```{math}
U_{\mathrm{HO}}
\;=\;
\frac{\hbar\omega}{2}
+
\frac{\hbar\omega}{\,e^{\beta\,\hbar\omega}-1\,}.
```

- **Heat Capacity** (with $\Theta_{\mathrm{vib}} = \hbar\omega/k_{\mathrm B}$):

```{math}
C_V^{(\mathrm{HO})}
\;=\;
k_{\mathrm B}
\Bigl(\frac{\Theta_{\mathrm{vib}}}{T}\Bigr)^{2}\,
\frac{e^{\Theta_{\mathrm{vib}}/T}}{\bigl(e^{\Theta_{\mathrm{vib}}/T} - 1\bigr)^{2}}.
```

---

K. **Rotational Partition Function (Linear Rotor, High $T$)**

```{math}
q_{\mathrm{rot}}
\;\approx\;
\frac{T}{\sigma\,\Theta_{\mathrm{rot}}},
\quad
\Theta_{\mathrm{rot}} = \frac{\hbar^2}{2k_{\mathrm B}\,I}.
```

- **Internal Energy**:

```{math}
U_{\mathrm{rot}} \;\approx\; k_{\mathrm B}\,T.
```

- **Heat Capacity**:

```{math}
C_V^{(\mathrm{rot})} \;\approx\; k_{\mathrm B}.
```

---

L. **Polyatomic Molecules**

- **General form** (when the energy of one molecule is a sum of translational, rotational, vibrational, and electronic parts):

```{math}
q
\;=\;
q_{\mathrm{trans}}\;
q_{\mathrm{rot}}\;
q_{\mathrm{vib}}\;
q_{\mathrm{elec}}.
```

- **Rotation of a nonlinear molecule** (high $T$):

```{math}
q_{\mathrm{rot}}
\;\approx\;
\frac{\sqrt{\pi}}{\sigma}\,
\sqrt{\frac{T^3}{\Theta_{\mathrm{rot},A}\,\Theta_{\mathrm{rot},B}\,\Theta_{\mathrm{rot},C}}}.
```

- **Vibration** (product over the $3N-6$ normal modes of a nonlinear molecule or the $3N-5$ normal modes of a linear molecule with $N$ atoms):

```{math}
q_{\mathrm{vib}}
\;=\;
\prod_{j} \frac{e^{-\Theta_{\mathrm{vib},j}/(2T)}}{1 - e^{-\Theta_{\mathrm{vib},j}/T}}.
```

- **Electronic** (when every term after the first is negligible):

```{math}
q_{\mathrm{elec}}
\;\approx\;
g_1.
```

## 3. Practice Problems

Each problem is a multi-part question of exam length. Work through it on paper before you open the solution.

### Problem 1. Fine-structure levels of the silicon atom

Draws on [Section 2.2](section-02.md), [Section 2.3](section-03.md), and [Section 2.8](section-08.md).

The lowest electronic levels of the silicon atom are three closely spaced levels. The table lists their degeneracies and their energies above the lowest level as wavenumbers ([NIST Atomic Spectra Database](https://www.nist.gov/pml/atomic-spectra-database)).

| Level $i$ | $g_i$ | $\varepsilon_i/hc\ (\mathrm{cm^{-1}})$ |
|---:|---:|---:|
| 1 | 1 | 0.000 |
| 2 | 3 | 77.115 |
| 3 | 5 | 223.157 |

The fourth level lies at $6298.850\ \mathrm{cm^{-1}}$ and is left out of this problem. Use $hc/k_{\mathrm B}=1.43878\ \mathrm{cm\,K}$ and $R=8.3145\ \mathrm{J\,mol^{-1}\,K^{-1}}$.

**(a)** Convert the energies of levels 2 and 3 to temperatures, $\varepsilon_i/k_{\mathrm B}$. Then find $q_{\mathrm{elec}}$ and the probabilities $p_1$, $p_2$, and $p_3$ of the three levels at $T=298.15\ \mathrm{K}$.

**(b)** Repeat the calculation at $T=1000\ \mathrm{K}$, and find the molar electronic energy $U_m$, measured from the lowest level, at both temperatures.

**(c)** Find, to the nearest kelvin, the temperature at which half of the atoms are in the lowest level.

````{admonition} Solution
:class: dropdown

1. **Level energies as temperatures**

   ```{math}
   \begin{aligned}
   \frac{\varepsilon_2}{k_{\mathrm B}}&=(1.43878\ \mathrm{cm\,K})(77.115\ \mathrm{cm^{-1}})=110.95\ \mathrm{K},\\
   \frac{\varepsilon_3}{k_{\mathrm B}}&=(1.43878\ \mathrm{cm\,K})(223.157\ \mathrm{cm^{-1}})=321.07\ \mathrm{K}.
   \end{aligned}
   ```

2. **Partition function at 298.15 K**

   With $\varepsilon_1=0$, the electronic partition function of [Section 2.8](section-08.md) is $q_{\mathrm{elec}}=g_1+g_2\,e^{-\varepsilon_2/k_{\mathrm B}T}+g_3\,e^{-\varepsilon_3/k_{\mathrm B}T}$. The two exponentials are $e^{-110.95/298.15}=0.6893$ and $e^{-321.07/298.15}=0.3407$, so

   ```{math}
   q_{\mathrm{elec}}=1+3(0.6893)+5(0.3407)=1+2.068+1.703=4.771.
   ```

3. **Probabilities at 298.15 K**

   The probability of level $i$ is its term in the sum divided by $q_{\mathrm{elec}}$:

   ```{math}
   p_1=\frac{1}{4.771}=0.210,
   \qquad
   p_2=\frac{2.068}{4.771}=0.433,
   \qquad
   p_3=\frac{1.703}{4.771}=0.357.
   ```

4. **Partition function and probabilities at 1000 K**

   The exponentials are $e^{-110.95/1000}=0.8950$ and $e^{-321.07/1000}=0.7254$, so

   ```{math}
   q_{\mathrm{elec}}=1+3(0.8950)+5(0.7254)=1+2.685+3.627=7.312,
   ```

   ```{math}
   p_1=\frac{1}{7.312}=0.137,
   \qquad
   p_2=\frac{2.685}{7.312}=0.367,
   \qquad
   p_3=\frac{3.627}{7.312}=0.496.
   ```

5. **Molar electronic energy**

   The mean energy of one atom is $\sum_i p_i\,\varepsilon_i$ ([Section 2.3](section-03.md)), so $U_m=N_{\mathrm A}\sum_i p_i\,\varepsilon_i=R\sum_i p_i\,(\varepsilon_i/k_{\mathrm B})$. Carrying four figures in the probabilities,

   ```{math}
   \begin{aligned}
   U_m(298.15\ \mathrm{K})&=R\,[(0.4334)(110.95\ \mathrm{K})+(0.3570)(321.07\ \mathrm{K})]\\
   &=R\,(162.7\ \mathrm{K})=1.35\ \mathrm{kJ\,mol^{-1}},\\
   U_m(1000\ \mathrm{K})&=R\,[(0.3672)(110.95\ \mathrm{K})+(0.4960)(321.07\ \mathrm{K})]\\
   &=R\,(200.0\ \mathrm{K})=1.66\ \mathrm{kJ\,mol^{-1}}.
   \end{aligned}
   ```

6. **Temperature at which $p_1=\tfrac{1}{2}$**

   Because $p_1=1/q_{\mathrm{elec}}$, the condition is $q_{\mathrm{elec}}=2$, that is,

   ```{math}
   3\,e^{-110.95/T}+5\,e^{-321.07/T}=1 .
   ```

   This equation has no closed-form solution, so try values of $T$. The left side is $0.750+0.090=0.840$ at $80\ \mathrm{K}$, $0.874+0.141=1.016$ at $90\ \mathrm{K}$, and $0.862+0.136=0.998$ at $89\ \mathrm{K}$. The left side increases with $T$, so the solution lies just above $89\ \mathrm{K}$.

**Result.** At $298.15\ \mathrm{K}$, $q_{\mathrm{elec}}=4.771$ and $(p_1,p_2,p_3)=(0.210,\,0.433,\,0.357)$; at $1000\ \mathrm{K}$, $q_{\mathrm{elec}}=7.312$ and $(p_1,p_2,p_3)=(0.137,\,0.367,\,0.496)$. The molar electronic energy is $1.35\ \mathrm{kJ\,mol^{-1}}$ and $1.66\ \mathrm{kJ\,mol^{-1}}$ at the two temperatures, and half of the atoms are in the lowest level at $T\approx 89\ \mathrm{K}$.
````

### Problem 2. Vibration of bromine

Draws on [Section 2.6](section-06.md).

The harmonic vibrational wavenumber of $\mathrm{Br_2}$ is $\tilde{\nu}=325.32\ \mathrm{cm^{-1}}$ ([NIST Chemistry WebBook](https://webbook.nist.gov/chemistry/)). Treat the vibration as a one-dimensional harmonic oscillator, with the energy measured from the bottom of the potential well as in Section 2.6. Use $hc/k_{\mathrm B}=1.43878\ \mathrm{cm\,K}$ and $R=8.3145\ \mathrm{J\,mol^{-1}\,K^{-1}}$.

**(a)** Find the vibrational temperature $\Theta_{\mathrm{vib}}$.

**(b)** At $T=298.15\ \mathrm{K}$, find the partition function $q$, the probabilities $p_0$, $p_1$, and $p_2$ of the three lowest levels, and the fraction of the molecules that are in excited levels.

**(c)** Find the molar vibrational energy in excess of the zero-point energy, $U_m-\tfrac{1}{2}R\,\Theta_{\mathrm{vib}}$, at $298.15\ \mathrm{K}$ and at $1000\ \mathrm{K}$.

**(d)** Find the molar vibrational heat capacity $C_{V,m}$ at the same two temperatures, as a multiple of $R$ and in $\mathrm{J\,mol^{-1}\,K^{-1}}$.

````{admonition} Solution
:class: dropdown

1. **Vibrational temperature**

   ```{math}
   \Theta_{\mathrm{vib}}=\frac{hc\,\tilde{\nu}}{k_{\mathrm B}}
   =(1.43878\ \mathrm{cm\,K})(325.32\ \mathrm{cm^{-1}})=468\ \mathrm{K}.
   ```

2. **Partition function at 298.15 K**

   Here $\Theta_{\mathrm{vib}}/T=468/298.15=1.570$, so $e^{-\Theta_{\mathrm{vib}}/T}=0.2081$ and $e^{-\Theta_{\mathrm{vib}}/(2T)}=0.4562$. Then

   ```{math}
   q=\frac{e^{-\Theta_{\mathrm{vib}}/(2T)}}{1-e^{-\Theta_{\mathrm{vib}}/T}}
   =\frac{0.4562}{1-0.2081}=0.576.
   ```

3. **Probabilities of the levels at 298.15 K**

   The probability of level $n$ is $p_n=e^{-\beta E_n}/q=e^{-n\Theta_{\mathrm{vib}}/T}\left(1-e^{-\Theta_{\mathrm{vib}}/T}\right)$, so

   ```{math}
   \begin{aligned}
   p_0&=1-0.2081=0.792,\\
   p_1&=(0.2081)(0.7919)=0.165,\\
   p_2&=(0.2081)^2(0.7919)=0.034.
   \end{aligned}
   ```

   The fraction of the molecules in excited levels is $1-p_0=0.208$.

4. **Vibrational energy in excess of the zero-point energy**

   From the internal energy of Section 2.6, $U_m-\tfrac{1}{2}R\,\Theta_{\mathrm{vib}}=R\,\Theta_{\mathrm{vib}}/(e^{\Theta_{\mathrm{vib}}/T}-1)$, with $R\,\Theta_{\mathrm{vib}}=3891\ \mathrm{J\,mol^{-1}}$. At $298.15\ \mathrm{K}$, $e^{\Theta_{\mathrm{vib}}/T}=4.805$; at $1000\ \mathrm{K}$, $\Theta_{\mathrm{vib}}/T=0.468$ and $e^{\Theta_{\mathrm{vib}}/T}=1.5968$. Therefore

   ```{math}
   \begin{aligned}
   U_m-\tfrac{1}{2}R\,\Theta_{\mathrm{vib}}&=\frac{3891\ \mathrm{J\,mol^{-1}}}{4.805-1}=1.02\ \mathrm{kJ\,mol^{-1}}
   &&(298.15\ \mathrm{K}),\\
   U_m-\tfrac{1}{2}R\,\Theta_{\mathrm{vib}}&=\frac{3891\ \mathrm{J\,mol^{-1}}}{1.5968-1}=6.52\ \mathrm{kJ\,mol^{-1}}
   &&(1000\ \mathrm{K}).
   \end{aligned}
   ```

5. **Heat capacity**

   Per mole, the heat capacity of Section 2.6 is $C_{V,m}=R\,(\Theta_{\mathrm{vib}}/T)^2\,e^{\Theta_{\mathrm{vib}}/T}/(e^{\Theta_{\mathrm{vib}}/T}-1)^2$. Therefore

   ```{math}
   \begin{aligned}
   C_{V,m}&=R\,\frac{(1.570)^2(4.805)}{(3.805)^2}=0.818\,R=6.80\ \mathrm{J\,mol^{-1}\,K^{-1}}
   &&(298.15\ \mathrm{K}),\\
   C_{V,m}&=R\,\frac{(0.468)^2(1.5968)}{(0.5968)^2}=0.982\,R=8.16\ \mathrm{J\,mol^{-1}\,K^{-1}}
   &&(1000\ \mathrm{K}).
   \end{aligned}
   ```

**Result.** $\Theta_{\mathrm{vib}}=468\ \mathrm{K}$. At $298.15\ \mathrm{K}$, $q=0.576$, $(p_0,p_1,p_2)=(0.792,\,0.165,\,0.034)$, and $20.8\%$ of the molecules are in excited levels. The vibrational energy in excess of the zero-point energy is $1.02\ \mathrm{kJ\,mol^{-1}}$ at $298.15\ \mathrm{K}$ and $6.52\ \mathrm{kJ\,mol^{-1}}$ at $1000\ \mathrm{K}$, and $C_{V,m}$ is $0.818\,R$ and $0.982\,R$.
````

### Problem 3. The partition function of carbon dioxide

Draws on [Section 2.5](section-05.md), [Section 2.6](section-06.md), [Section 2.7](section-07.md), and [Section 2.8](section-08.md).

Carbon dioxide is a linear molecule with molar mass $M=44.0095\ \mathrm{g\,mol^{-1}}$, rotational constant $\tilde{B}=0.39021\ \mathrm{cm^{-1}}$, and symmetry number $\sigma=2$. Its four normal modes have the wavenumbers $1333\ \mathrm{cm^{-1}}$ (symmetric stretch), $667\ \mathrm{cm^{-1}}$ (bend, two modes), and $2349\ \mathrm{cm^{-1}}$ (antisymmetric stretch). The molecular constants are from the [NIST Computational Chemistry Comparison and Benchmark Database](https://cccbdb.nist.gov/). The ground electronic level is not degenerate, and the excited electronic levels can be ignored, so $q_{\mathrm{elec}}=1$. Use $hc/k_{\mathrm B}=1.43878\ \mathrm{cm\,K}$ and $R=8.3145\ \mathrm{J\,mol^{-1}\,K^{-1}}$.

**(a)** Find the thermal de Broglie wavelength $\Lambda$ at $T=298.15\ \mathrm{K}$ and the translational partition function $q_{\mathrm{trans}}$ of one molecule in the molar volume of an ideal gas at $298.15\ \mathrm{K}$ and $1\ \mathrm{bar}$.

**(b)** Find $\Theta_{\mathrm{rot}}$ and $q_{\mathrm{rot}}$ at $298.15\ \mathrm{K}$.

**(c)** Find the vibrational temperatures of the four normal modes, $q_{\mathrm{vib}}$ at $298.15\ \mathrm{K}$, and the molecular partition function $q$.

**(d)** Find $C_{V,m}/R$ at $298.15\ \mathrm{K}$ as the sum of the translational, rotational, and vibrational contributions. Compare it with $C_{V,m}=28.815\ \mathrm{J\,mol^{-1}\,K^{-1}}$, the value that corresponds to the ideal-gas $C_{P,m}=37.129\ \mathrm{J\,mol^{-1}\,K^{-1}}$ of the [NIST-JANAF tables](https://janaf.nist.gov/) (for an ideal gas the two differ by $R$).

````{admonition} Solution
:class: dropdown

1. **Thermal de Broglie wavelength**

   The mass of one molecule is $m=M/N_{\mathrm A}=7.308\times10^{-26}\ \mathrm{kg}$, so

   ```{math}
   \begin{aligned}
   \Lambda&=\frac{h}{\sqrt{2\pi m k_{\mathrm B}T}}
   =\frac{6.626\times10^{-34}}{\sqrt{2\pi(7.308\times10^{-26})(1.381\times10^{-23})(298.15)}}\\
   &=1.524\times10^{-11}\ \mathrm{m}.
   \end{aligned}
   ```

2. **Translational partition function**

   The molar volume is $V=RT/P=(8.3145)(298.15)/(1.00\times10^{5})=2.479\times10^{-2}\ \mathrm{m^3}$, so

   ```{math}
   q_{\mathrm{trans}}=\frac{V}{\Lambda^3}
   =\frac{2.479\times10^{-2}}{(1.524\times10^{-11})^3}
   =7.00\times10^{30}.
   ```

3. **Rotational partition function**

   ```{math}
   \begin{aligned}
   \Theta_{\mathrm{rot}}&=\frac{hc\,\tilde{B}}{k_{\mathrm B}}=(1.43878)(0.39021)=0.5614\ \mathrm{K},\\
   q_{\mathrm{rot}}&=\frac{T}{\sigma\,\Theta_{\mathrm{rot}}}=\frac{298.15}{2(0.5614)}=265.5.
   \end{aligned}
   ```

4. **Vibrational temperatures**

   With $\Theta_{\mathrm{vib},j}=hc\,\tilde{\nu}_j/k_{\mathrm B}$, the four modes have $\Theta_{\mathrm{vib},j}=1918\ \mathrm{K}$, $959.7\ \mathrm{K}$, $959.7\ \mathrm{K}$, and $3380\ \mathrm{K}$.

5. **Vibrational partition function**

   Each mode contributes the factor $e^{-\Theta_{\mathrm{vib},j}/(2T)}/(1-e^{-\Theta_{\mathrm{vib},j}/T})$:

   | Mode | $\Theta_{\mathrm{vib},j}/T$ | $e^{-\Theta_{\mathrm{vib},j}/(2T)}$ | $1-e^{-\Theta_{\mathrm{vib},j}/T}$ | Factor |
   |---|---:|---:|---:|---:|
   | Symmetric stretch | 6.433 | 0.04010 | 0.9984 | 0.04016 |
   | Bend (each of two) | 3.219 | 0.2000 | 0.9600 | 0.2083 |
   | Antisymmetric stretch | 11.337 | 0.003454 | 1.0000 | 0.003454 |

   ```{math}
   q_{\mathrm{vib}}=(0.04016)(0.2083)^2(0.003454)=6.02\times10^{-6}.
   ```

6. **Molecular partition function**

   ```{math}
   \begin{aligned}
   q&=q_{\mathrm{trans}}\,q_{\mathrm{rot}}\,q_{\mathrm{vib}}\,q_{\mathrm{elec}}\\
   &=(7.00\times10^{30})(265.5)(6.02\times10^{-6})(1)
   =1.12\times10^{28}.
   \end{aligned}
   ```

7. **Heat capacity**

   Translation contributes $\tfrac{3}{2}R$ ([Section 2.5](section-05.md)), and the rotation of a linear molecule contributes $R$ ([Section 2.7](section-07.md)). Each normal mode contributes $R\,x_j^2\,e^{x_j}/(e^{x_j}-1)^2$ with $x_j=\Theta_{\mathrm{vib},j}/T$ ([Section 2.6](section-06.md)), which is $0.0667\,R$ for the symmetric stretch, $0.4497\,R$ for each bending mode, and $0.0015\,R$ for the antisymmetric stretch. The sum is

   ```{math}
   \frac{C_{V,m}}{R}=\frac{3}{2}+1+\left[0.0667+2(0.4497)+0.0015\right]=2.5+0.968=3.468,
   ```

   so $C_{V,m}=28.83\ \mathrm{J\,mol^{-1}\,K^{-1}}$. The NIST-JANAF value, $28.815\ \mathrm{J\,mol^{-1}\,K^{-1}}$, corresponds to $C_{V,m}/R=3.466$. The two agree to within $0.1\%$.

**Result.** At $298.15\ \mathrm{K}$ and $1\ \mathrm{bar}$, $\Lambda=1.524\times10^{-11}\ \mathrm{m}$, $q_{\mathrm{trans}}=7.00\times10^{30}$, $q_{\mathrm{rot}}=265.5$, $q_{\mathrm{vib}}=6.02\times10^{-6}$, and $q=1.12\times10^{28}$. The model gives $C_{V,m}=3.468\,R$, against $3.466\,R$ from the NIST-JANAF tables.
````

### Problem 4. A gas confined to a line

Draws on [Section 2.3](section-03.md), [Section 2.4](section-04.md), and [Section 2.5](section-05.md).

$N$ atoms of mass $m$ move along a narrow channel of length $L$ and do not interact with one another. Each atom is a particle in a one-dimensional box, with the energy levels $E_n=h^2n^2/(8mL^2)$, $n=1,2,3,\ldots$ ([Section 2.5](section-05.md)).

**(a)** Show that the partition function of one atom is $q=L/\Lambda$. Evaluate the sum over $n$ with the same integral that Section 2.5 uses for the cube.

**(b)** Write the partition function $Q$ of the $N$ indistinguishable atoms, and find $U$ and $C_V$.

**(c)** The gas pushes on each end of the channel with a force $f$, and the work done on the gas when the length changes is $-f\,dL$. The force therefore takes the place of the pressure, and $L$ takes the place of $V$, in the pressure formula of [Section 2.3](section-03.md): $f=k_{\mathrm B}T\,(\partial\ln Q/\partial L)_{N,T}$. Find the equation of state that relates $f$, $L$, $N$, and $T$.

**(d)** Evaluate $\Lambda$, $q$, $U$, and $f$ for $N=20$ argon atoms ($m=39.948\,u=6.634\times10^{-26}\ \mathrm{kg}$) in a channel of length $L=10.0\ \mathrm{nm}$ at $T=300\ \mathrm{K}$.

````{admonition} Solution
:class: dropdown

1. **Partition function of one atom**

   With $\alpha=\beta h^2/(8mL^2)$, the sum over the levels is $q=\sum_{n=1}^{\infty}e^{-\alpha n^2}$. Replacing the sum with the integral of Section 2.5,

   ```{math}
   q=\int_0^{\infty}e^{-\alpha x^2}\,dx
   =\left(\frac{\pi}{4\alpha}\right)^{1/2}
   =\left(\frac{2\pi m}{\beta h^2}\right)^{1/2}L
   =\frac{L}{\Lambda},
   ```

   where $\Lambda=\left(h^2/2\pi m k_{\mathrm B}T\right)^{1/2}$ is the thermal de Broglie wavelength.

2. **Partition function of $N$ atoms**

   For independent, indistinguishable atoms ([Section 2.4](section-04.md)),

   ```{math}
   Q=\frac{q^N}{N!}=\frac{1}{N!}\left(\frac{L}{\Lambda}\right)^{N},
   \qquad
   \ln Q=N\ln L-N\ln\Lambda-\ln N!.
   ```

3. **Internal energy and heat capacity**

   Because $\Lambda$ is proportional to $\beta^{1/2}$, $\ln Q=-\tfrac{1}{2}N\ln\beta$ plus terms that do not depend on $\beta$. Then

   ```{math}
   U=-\left(\frac{\partial\ln Q}{\partial\beta}\right)_{N,L}=\frac{N}{2\beta}=\frac{1}{2}Nk_{\mathrm B}T,
   \qquad
   C_V=\left(\frac{\partial U}{\partial T}\right)_{N,L}=\frac{1}{2}Nk_{\mathrm B}.
   ```

4. **Equation of state**

   Only the term $N\ln L$ in $\ln Q$ depends on $L$, so

   ```{math}
   f=k_{\mathrm B}T\left(\frac{\partial\ln Q}{\partial L}\right)_{N,T}=\frac{Nk_{\mathrm B}T}{L},
   \qquad\text{that is,}\qquad
   fL=Nk_{\mathrm B}T.
   ```

5. **Numerical values for argon**

   ```{math}
   \begin{aligned}
   \Lambda&=\frac{6.626\times10^{-34}}{\sqrt{2\pi(6.634\times10^{-26})(1.381\times10^{-23})(300)}}
   =1.595\times10^{-11}\ \mathrm{m},\\
   q&=\frac{L}{\Lambda}=\frac{1.00\times10^{-8}}{1.595\times10^{-11}}=627.
   \end{aligned}
   ```

   With $k_{\mathrm B}T=4.142\times10^{-21}\ \mathrm{J}$,

   ```{math}
   \begin{aligned}
   U&=\frac{1}{2}(20)(4.142\times10^{-21}\ \mathrm{J})=4.14\times10^{-20}\ \mathrm{J},\\
   f&=\frac{(20)(4.142\times10^{-21}\ \mathrm{J})}{1.00\times10^{-8}\ \mathrm{m}}=8.28\times10^{-12}\ \mathrm{N}.
   \end{aligned}
   ```

**Result.** $q=L/\Lambda$, $Q=(L/\Lambda)^N/N!$, $U=\tfrac{1}{2}Nk_{\mathrm B}T$, $C_V=\tfrac{1}{2}Nk_{\mathrm B}$, and $fL=Nk_{\mathrm B}T$. For 20 argon atoms in a $10.0\ \mathrm{nm}$ channel at $300\ \mathrm{K}$, $\Lambda=1.595\times10^{-11}\ \mathrm{m}$, $q=627$, $U=4.14\times10^{-20}\ \mathrm{J}$, and $f=8.28\ \mathrm{pN}$.
````

## 4. Further Reading

**Textbooks**

- @mcquarrie1997, Chapter 17 (the Boltzmann factor and partition functions) and Chapter 18 (partition functions and ideal gases). Chapters 3 and 5 review the particle in a box, the harmonic oscillator, and the rigid rotor.
- @dill2011, Chapter 1 (probability), Chapter 5 (entropy and the Boltzmann law), Chapter 10 (the Boltzmann distribution law), and Chapter 11 (the statistical mechanics of simple gases and solids).
- @hill1986, Chapters 1 and 2 (ensembles and thermodynamics), Chapter 5 (monatomic crystals), and Chapter 8 (the ideal diatomic gas).

**Original papers**

- @petit1819 report that the heat capacity per mole of atoms is nearly the same for the solid elements they measured, which is the Dulong–Petit law of [Section 2.6](section-06.md). An [English translation](https://web.lemoyne.edu/giunta/PETIT.html) appeared in *Annals of Philosophy* **14**, 189–198 (1819).
- @einstein1907 applies Planck's quantization of energy to the vibrations of the atoms of a crystal and obtains the heat capacity of the Einstein model of [Section 2.6](section-06.md).
