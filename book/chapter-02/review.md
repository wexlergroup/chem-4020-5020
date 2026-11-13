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
