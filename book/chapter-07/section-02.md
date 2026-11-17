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
  enumerator: 7.2.%s
---


# 7.2. Equilibrium Constants from Microscopic Properties

See also: [Course-wide Conventions & Notation](../notation.md)

## Overview and Learning Objectives

[Section 7.1](section-01.md) derived the equilibrium constant from macroscopic thermochemistry: given $\Delta_r H^\circ$ and $\Delta_r S^\circ$ (e.g., from NIST/JANAF/ATcT tables), one computes $\Delta_r G^\circ$ and then $K = \exp(-\Delta_r G^\circ/RT)$.

In this section we derive the same $K$ from **microscopic** information — molecular masses, rotational and vibrational constants, electronic degeneracies, and bond dissociation energies. This completes the statistical-mechanical bridge built in Modules 2, 4, and 5: Module 4 connected the canonical partition function to the Helmholtz free energy, $A = -k_{\mathrm B}T\ln Q$; Module 5 showed that $A$ is the master potential at constant $T,V$, generating $S$, $U$, $P$, and — as we now develop — the chemical potentials $\mu_i$. Once we have $\mu_i$ in terms of the single-molecule partition functions $q_i$, the equilibrium condition $\sum_i\nu_i\mu_i = 0$ gives $K$ directly in terms of molecular properties.

The payoff: for any gas-phase reaction, we can predict $K(T)$ from spectroscopic data alone, and we can *see* which molecular features (masses, moments of inertia, bond strengths) drive the equilibrium.

```{admonition} Notation note
:class: note

In this section, partition functions are $Q$ (many-particle, canonical) and $q$ (single-particle), consistent with Modules 2, 4, and 5 and the [course-wide conventions](../notation.md). The reaction quotient $Q_p$ of [Section 7.1](section-01.md) does not appear in this section, so no notational collision arises. All partition functions here are of the single- or many-particle kind.
```

Learning objectives:

- Write the equilibrium condition at fixed $T,V$ as $\sum_i \nu_i\mu_i=0$ from minimizing the Helmholtz free energy.
- Use $Q = \prod_i q_i^{N_i}/N_i!$ for an ideal-gas mixture to derive $\mu_i$ from $A = -k_{\mathrm B}T\ln Q$.
- Express equilibrium constants in terms of single-molecule partition functions, and relate $K_p$ and $K_c$ through $\Delta\nu$.
- Identify how translational, rotational, vibrational, and electronic factors contribute to $K$ for a specific reaction, and compute $K(T)$ numerically from molecular constants.

## Core Ideas and Derivations

(section:helmholtz-equilibrium)=
### The Right Potential for Equilibrium at Constant $T,V$

Consider a general gas-phase reaction

```{math}
a\,\mathrm{A(g)} + b\,\mathrm{B(g)} \rightleftharpoons y\,\mathrm{Y(g)} + z\,\mathrm{Z(g)},
```

under the **ideal-gas** assumption used throughout this module, with the signed stoichiometric coefficients $\nu_i$ of [Section 7.1](section-01.md).

[Section 5.1](../chapter-05/section-01.md) established that at constant $T$ and $V$, the relevant thermodynamic potential is the **Helmholtz free energy** $A$. When the amounts of the species can change, its differential acquires a term $\mu_i\,dn_i$ for each species, as the differential of $G$ does in [Section 6.1](../chapter-06/section-01.md):

```{math}
dA = -S\,dT - P\,dV + \sum_i \mu_i\,dn_i.
```

At constant $T$ and $V$,

```{math}
dA = \sum_i \mu_i\,dn_i = \left(\sum_i \nu_i\,\mu_i\right)d\xi,
```

where the second equality uses $dn_i = \nu_i\,d\xi$ from [Section 7.1](section-01.md) (under [Extent of Reaction $\xi$](#section:extent-of-reaction)). Equilibrium at fixed $T,V$ occurs at a minimum of $A$, so

```{math}
\left(\frac{\partial A}{\partial \xi}\right)_{T,V} = 0
\qquad\Longrightarrow\qquad
\boxed{\sum_i \nu_i\,\mu_i = 0.}
```

This is the **chemical equilibrium condition** written in terms of chemical potentials. It has the same form as the condition of [Section 7.1](section-01.md). What differs is the potential that is stationary: $G$ at fixed $T,P$, and $A$ at fixed $T,V$. The chemical potential can be taken from either potential (the second expression follows from $G = A + PV$):

```{math}
\mu_i = \left(\frac{\partial G}{\partial n_i}\right)_{T,P,n_{j\neq i}}
= \left(\frac{\partial A}{\partial n_i}\right)_{T,V,n_{j\neq i}}.
```

Concept Check 1 asks why the two conditions agree, and Concept Check 1 of [Section 7.1](section-01.md) asks what changes when the volume is held fixed instead of the pressure.

---

(section:mixture-partition-function)=
### Partition Function of an Ideal-Gas Mixture

From [Section 2.4](../chapter-02/section-04.md) (factorization over independent subsystems, and the factor $1/N!$ for indistinguishable particles), the canonical partition function of a mixture of ideal gases factorizes:

```{math}
:label: eq:Q-mixture
\boxed{
Q(T,V,\{N_i\}) = \prod_i Q_i(T,V,N_i),
\qquad
Q_i(T,V,N_i) = \frac{q_i(T,V)^{N_i}}{N_i!},
}
```

where $q_i$ is the **single-molecule** partition function for species $i$. Each species is treated independently because ideal-gas molecules don't interact.

---

(section:mu-from-partition-function)=
### Chemical Potentials from the Partition Function

[Section 4.3](../chapter-04/section-03.md) established

```{math}
A = -k_{\mathrm B}T\ln Q,
```

and [Section 5.1](../chapter-05/section-01.md) developed $A$ as the "master potential" at fixed $T,V$, with $S$, $U$, $P$ all following from its derivatives. The chemical potential $\mu_i$ is the natural extension of that derivative machinery to composition:

```{math}
\mu_i^{\mathrm{(per\ molecule)}}
\equiv \left(\frac{\partial A}{\partial N_i}\right)_{T,V,\{N_{j\neq i}\}}
= -k_{\mathrm B}T\left(\frac{\partial \ln Q}{\partial N_i}\right)_{T,V,\{N_{j\neq i}\}}.
```

Using Eq. {eq}`eq:Q-mixture`,

```{math}
\ln Q = \sum_i \left(N_i\ln q_i - \ln N_i!\right),
```

and applying Stirling's approximation in the thermodynamic limit ($\ln N_i! \approx N_i\ln N_i - N_i$, as introduced in [Section 2.5](../chapter-02/section-05.md)),

```{math}
\left(\frac{\partial \ln Q}{\partial N_i}\right)_{T,V,\{N_{j\neq i}\}}
= \ln q_i - \ln N_i = \ln\!\left(\frac{q_i}{N_i}\right).
```

So the per-molecule chemical potential is

```{math}
\mu_i^{\mathrm{(per\ molecule)}} = -k_{\mathrm B}T\ln\!\left(\frac{q_i}{N_i}\right).
```

Multiplying by Avogadro's number $N_{\mathrm A}$ to express $\mu_i$ on a **per-mole** basis (as used throughout [Section 7.1](section-01.md)),

```{math}
:label: eq:mu-from-q
\boxed{
\mu_i = -RT\ln\!\left(\frac{q_i}{N_i}\right) = RT\ln\!\left(\frac{N_i}{q_i}\right),
}
```

where $R = N_{\mathrm A}k_{\mathrm B}$. Equation {eq}`eq:mu-from-q` is the microscopic counterpart of the thermodynamic expression $\mu_i = \mu_i^\circ(T) + RT\ln(P_i/P^\circ)$ from [Section 7.1](section-01.md): both express the chemical potential of a gas-phase species, one in terms of partial pressure and a reference $\mu_i^\circ(T)$, the other in terms of microscopic counting through $q_i/N_i$.

---

(section:k-from-partition-functions)=
### Equilibrium Constants from Partition Functions

Inserting Eq. {eq}`eq:mu-from-q` into the equilibrium condition $\sum_i \nu_i\mu_i=0$:

```{math}
\sum_i \nu_i\ln\!\left(\frac{q_i}{N_i}\right)=0
\qquad\Longrightarrow\qquad
\prod_i\left(\frac{q_i}{N_i}\right)^{\nu_i}=1.
```

Equivalently,

```{math}
:label: eq:K-partition-numbers
\boxed{
\prod_i N_i^{\nu_i} = \prod_i q_i^{\nu_i}.
}
```

#### Relating to $K_c$ and $K_p$

Define a number concentration $c_i = N_i/V$. Dividing Eq. {eq}`eq:K-partition-numbers` by $V^{\Delta\nu}$ (with $\Delta\nu \equiv \sum_i \nu_i$) gives

```{math}
K_c \equiv \prod_i c_i^{\nu_i} = \prod_i\!\left(\frac{q_i}{V}\right)^{\nu_i}.
```

For an ideal gas, $P_i = c_i\, k_{\mathrm B}T$, so the pressure-based equilibrium constant is

```{math}
:label: eq:Kp-Kc
\boxed{
K_p \equiv \prod_i\left(\frac{P_i}{P^\circ}\right)^{\nu_i} = K_c\,\left(\frac{k_{\mathrm B}T}{P^\circ}\right)^{\Delta\nu}.
}
```

The prefactor $(k_{\mathrm B}T/P^\circ)^{\Delta\nu}$ converts between "molecules per unit volume" and "pressure ratios relative to $P^\circ$." When $\Delta\nu = 0$, this prefactor is unity and $K_p = K_c$ is dimensionless directly.

So, **once you can compute** the single-molecule partition functions $q_i$, you can compute $K_c$ (and then $K_p$).

#### The Simplest Case: $\mathrm{A} \rightleftharpoons \mathrm{B}$

For an isomerization, $\mathrm{A} \rightleftharpoons \mathrm{B}$, Eq. {eq}`eq:K-partition-numbers` reads

```{math}
:label: eq:K-isomerization
K = \frac{N_{\mathrm{B}}}{N_{\mathrm{A}}} = \frac{q_{\mathrm{B}}}{q_{\mathrm{A}}}.
```

Both partition functions are sums over the levels of one molecule, and both measure the energy from the same zero. A molecule that can convert between the two forms therefore has a single set of levels, of which some belong to A and the others to B. At equilibrium every level is populated according to one Boltzmann distribution, and $K$ is the total population of the levels of B divided by the total population of the levels of A.

[](#fig:K-populations) shows a model in which the lowest level of B lies above the lowest level of A, while the levels of B lie closer together. At the lower temperature almost all molecules occupy the lowest levels, which belong to A, and $K < 1$. At the higher temperature the many levels of B together hold more molecules than the few levels of A, and $K > 1$. The energy of the lowest level favors A, and the number of accessible levels favors B.

```{code-cell} ipython3
:label: fig:K-populations
:caption: Populations of the levels of a model molecule with two forms, A and B, at (a) $k_{\mathrm{B}}T = \varepsilon$ and (b) $k_{\mathrm{B}}T = 4\varepsilon$. The levels of A lie at $0, 2\varepsilon, 4\varepsilon, \dots$ and those of B at $3\varepsilon, 3.5\varepsilon, 4\varepsilon, \dots$. The length of each bar is the population of its level relative to the lowest level, and the dashed curve is the Boltzmann factor $e^{-E/k_{\mathrm{B}}T}$, which is the same for both forms. $K$ is the total population of the levels of B divided by that of the levels of A, Eq. [](#eq:K-isomerization). The sums include the levels above $12\varepsilon$, which are not drawn.
:tags: [hide-input]

import sys

import numpy as np
import matplotlib.pyplot as plt

sys.path.append("..")            # the book folder, which holds figstyle
import figstyle as fs

# Model: energies in units of epsilon, temperatures as t = k_B T / epsilon
spacing_A, spacing_B, offset_B = 2.0, 0.5, 3.0
E_top = 12.0
levels_A = np.arange(0.0, E_top + 1e-9, spacing_A)
levels_B = np.arange(offset_B, E_top + 1e-9, spacing_B)

def K_model(t):
    # q_B / q_A with the sums over all levels (geometric series)
    q_A = 1 / (1 - np.exp(-spacing_A / t))
    q_B = np.exp(-offset_B / t) / (1 - np.exp(-spacing_B / t))
    return q_B / q_A


def draw_populations(c):
    fig, axes = plt.subplots(1, 2, figsize=(fs.WIDTH, 3.9), sharey=True, constrained_layout=True)
    columns = (("A", 0.2, levels_A, c.blue), ("B", 1.95, levels_B, c.vermilion))
    for ax, letter, t, title in zip(axes, "ab", (1.0, 4.0),
                                    ("$k_{\\mathrm{B}}T = \\varepsilon$", "$k_{\\mathrm{B}}T = 4\\varepsilon$")):
        for name, x0, levels, color in columns:
            ax.hlines(levels, x0, x0 + 0.22, color=c.muted, lw=1.0)
            ax.barh(levels, np.exp(-levels / t), left=x0 + 0.30, height=0.3, color=color, lw=0)
            E = np.linspace(levels[0], E_top + 0.3, 200)
            ax.plot(x0 + 0.30 + np.exp(-E / t), E, color=c.gray, lw=1.2, ls=(0, (4, 3)))
            ax.text(x0 + 0.11, -0.75, name, color=color, ha="center", va="top", fontweight="bold")
        ax.text(1.85, 13.2, f"$K = N_{{\\mathrm{{B}}}}/N_{{\\mathrm{{A}}}} = {K_model(t):.2f}$", ha="center")
        ax.set(xlim=(0, 3.4), ylim=(-0.6, 13.9), xticks=[])
        ax.spines["bottom"].set_visible(False)
        fs.panel(ax, letter, title)
    axes[0].set_ylabel("Energy, $E/\\varepsilon$")
    return fig


fs.show(draw_populations, alt="Two panels, each with two ladders of energy levels on one energy axis. The levels of "
                              "form A start at zero and are 2 units apart. The levels of form B start at 3 units "
                              "and are half a unit apart. A bar on each level shows its population, which falls "
                              "with energy along the same exponential curve for both forms. In panel a, at the "
                              "lower temperature, nearly all the population is on the lowest level of A, and K is "
                              f"{K_model(1.0):.2f}. In panel b, at four times that temperature, many levels of B "
                              f"are populated, and K is {K_model(4.0):.2f}.")
```

In [](#fig:equilibrium-levels-interactive), you can change the temperature and the energy of the lowest level of B in the same model, and follow $q_{\mathrm{A}}$, $q_{\mathrm{B}}$, and $K$.

```{iframe} https://wexlergroup.github.io/chem-4020-5020/computational-studios/equilibrium-levels.html
:label: fig:equilibrium-levels-interactive
:width: 100%
:title: Equilibrium from energy levels
:placeholder: ../_static/equilibrium-levels/equilibrium-levels.png

Populations of the levels of a molecule with two forms, A and B, and the equilibrium constant $K = q_{\mathrm{B}}/q_{\mathrm{A}}$. Panel (a) shows the levels of both forms on one energy scale, with a bar for the population of each level. Panel (b) shows $K$ against temperature for the chosen energy of the lowest level of B, on a logarithmic scale, with the present temperature marked. The dashed line marks $K = 1$.
```

Two things to try:

1. Start from the levels of [](#fig:K-populations), and raise the temperature until $K = 1$. Compare the two partition functions at this temperature, and find the levels that hold most of the molecules of each form.
2. Set the lowest level of B to the same energy as the lowest level of A, and raise the temperature from its lowest to its highest value. Neither form is now lower in energy. Follow how $K$ changes, and relate the change to the number of levels of each form that lie within $k_{\mathrm{B}}T$ of the lowest level.

---

(section:hi-formation)=
### Example: $\mathrm{HI(g)}$ Formation

Consider

```{math}
\mathrm{H_2(g) + I_2(g) \rightleftharpoons 2HI(g)}.
```

Here,

```{math}
\Delta\nu = 2 - 1 - 1 = 0,
```

so $K_p = K_c$ and the equilibrium constant is

```{math}
:label: eq:K-HI-ratio
\boxed{
K = \frac{q_{\mathrm{HI}}^2}{q_{\mathrm{H_2}}\,q_{\mathrm{I_2}}}.
}
```

The $\Delta\nu = 0$ choice is pedagogically convenient: it makes volume factors cancel exactly (see the Worked Example below), so $K$ depends only on intrinsic molecular properties (masses, rotational and vibrational constants, bond strengths). With $\Delta\nu \neq 0$ reactions, one would need the $(k_{\mathrm B}T/P^\circ)^{\Delta\nu}$ factor from Eq. {eq}`eq:Kp-Kc` to convert between $K_c$ and $K_p$.

---

(section:diatomic-partition-function)=
### Single-Molecule Partition Function for a Diatomic Ideal Gas

Within the rigid-rotor/harmonic-oscillator approximations developed in Module 2 (see [Section 2.8](../chapter-02/section-08.md) for the summary table), the single-molecule partition function for a diatomic gas factors as

```{math}
q = q_{\mathrm{trans}}\,q_{\mathrm{rot}}\,q_{\mathrm{vib}}\,q_{\mathrm{elec}}.
```

````{admonition} Energy-zero convention
:class: note

A partition function depends on *where* one places the zero of energy, and the partition functions of different molecules can be combined only if they share one zero. In this section the zero of molecular energy lies at the **dissociated atoms at rest**. The ground vibrational level, $v=0$, of a bound diatomic molecule then sits at the energy $-D_0$, where $D_0$ is the bond dissociation energy measured from $v=0$. (The vibrational quantum number $v$ is the $n$ of [Section 2.6](../chapter-02/section-06.md).)

With this choice the bond energy enters the electronic factor as $e^{+D_0/k_{\mathrm B}T}$ (per molecule) or $e^{+D_0/RT}$ (per mole), and the vibrational energies are measured from the level $v=0$, so that level $v$ lies $v\,k_{\mathrm B}\Theta_{\mathrm{vib}}$ above it. The vibrational factor is therefore

```{math}
q_{\mathrm{vib}} = \sum_{v=0}^{\infty} e^{-v\,\Theta_{\mathrm{vib}}/T} = \frac{1}{1-e^{-\Theta_{\mathrm{vib}}/T}}.
```

This is not the $q_{\mathrm{vib}}$ of [Section 2.6](../chapter-02/section-06.md), which measures the energy from the bottom of the potential well. Either zero may be used, provided that $q_{\mathrm{vib}}$ and the bond energy refer to the same one. Concept Check 3 asks you to show this for $K$.
````

Putting the four factors together (cf. [Section 2.8](../chapter-02/section-08.md)):

```{math}
:label: eq:q-diatomic
\boxed{
q_{\mathrm{diatomic}}
=
\underbrace{\frac{V}{\Lambda^3}}_{q_{\mathrm{trans}}}
\;\underbrace{\frac{T}{\sigma\,\Theta_{\mathrm{rot}}}}_{q_{\mathrm{rot}}}
\;\underbrace{\frac{1}{1-e^{-\Theta_{\mathrm{vib}}/T}}}_{q_{\mathrm{vib}}}
\;\underbrace{g_1\,e^{D_0/(RT)}}_{q_{\mathrm{elec}}}.
}
```

where the thermal de Broglie wavelength is

```{math}
\Lambda \equiv \sqrt{\frac{h^2}{2\pi m k_{\mathrm B}T}}.
```

#### Notes on Parameters

- $m$: molecular mass.
- $\sigma$: rotational symmetry number — $\sigma=2$ for homonuclear diatomics (H$_2$, I$_2$), $\sigma=1$ for heteronuclear (HI). See [Section 2.7](../chapter-02/section-07.md).
- $\Theta_{\mathrm{rot}} \equiv \hbar^2/(2 I k_{\mathrm B})$: characteristic rotational temperature (from [Section 2.7](../chapter-02/section-07.md)).
- $\Theta_{\mathrm{vib}} \equiv h\nu/k_{\mathrm B}$: characteristic vibrational temperature (from [Section 2.6](../chapter-02/section-06.md)).
- $g_1$: ground-state electronic degeneracy. For the species in the HI reaction, all three ground states are nondegenerate singlets ($^1\Sigma_g^+$ for $\mathrm{H_2}$ and $\mathrm{I_2}$, $^1\Sigma^+$ for HI), so $g_1 = 1$ for each.
- $D_0$: bond dissociation energy referenced from the $v=0$ level (per-mole units when combined with $RT$, per-molecule when combined with $k_{\mathrm B}T$).

---

(section:k-hi-molecular)=
### $K$ for $\mathrm{H_2 + I_2 \rightleftharpoons 2HI}$ in Terms of Molecular Parameters

Substituting Eq. {eq}`eq:q-diatomic` into $K = q_{\mathrm{HI}}^2/(q_{\mathrm{H_2}}\,q_{\mathrm{I_2}})$, and using the fact that volume cancels when $\Delta\nu=0$ (see the Worked Example), yields

```{math}
:label: eq:K-HI-full
\boxed{
\begin{aligned}
K ={}&
\left(\frac{m_{\mathrm{HI}}^2}{m_{\mathrm{H_2}}\,m_{\mathrm{I_2}}}\right)^{3/2}
\;
\frac{\sigma_{\mathrm{H_2}}\sigma_{\mathrm{I_2}}}{\sigma_{\mathrm{HI}}^2}\,
\frac{\Theta_{\mathrm{rot}}^{\mathrm{H_2}}\,\Theta_{\mathrm{rot}}^{\mathrm{I_2}}}{\left(\Theta_{\mathrm{rot}}^{\mathrm{HI}}\right)^2}
\\[1ex]
&\times
\frac{\bigl(1-e^{-\Theta_{\mathrm{vib}}^{\mathrm{H_2}}/T}\bigr)\bigl(1-e^{-\Theta_{\mathrm{vib}}^{\mathrm{I_2}}/T}\bigr)}{\bigl(1-e^{-\Theta_{\mathrm{vib}}^{\mathrm{HI}}/T}\bigr)^2}
\;
\exp\!\left(\frac{2D_0^{\mathrm{HI}}-D_0^{\mathrm{H_2}}-D_0^{\mathrm{I_2}}}{RT}\right).
\end{aligned}
}
```

With $\sigma_{\mathrm{H_2}}=\sigma_{\mathrm{I_2}}=2$ and $\sigma_{\mathrm{HI}}=1$, the symmetry-number prefactor is $\sigma_{\mathrm{H_2}}\sigma_{\mathrm{I_2}}/\sigma_{\mathrm{HI}}^2 = 4$. With the energy zero at the separated atoms, the zero-point energies are contained in the $D_0$ values, and the vibrational block holds only the $(1-e^{-\Theta_{\mathrm{vib}}/T})$ "hot-mode" factors.

Equation {eq}`eq:K-HI-full` makes the microscopic physics very explicit:

- **Translation** contributes the $m^{3/2}$ mass dependence.
- **Rotation** contributes the symmetry-number and $\Theta_{\mathrm{rot}}$ factors.
- **Vibration** contributes the $(1-e^{-\Theta_{\mathrm{vib}}/T})$ factors that encode how many vibrational quanta are thermally accessible.
- **Bond strengths** enter through the dissociation energies $D_0$ in the exponential; the numerator $2D_0^{\mathrm{HI}}-D_0^{\mathrm{H_2}}-D_0^{\mathrm{I_2}}$ is essentially the negative of the reaction internal energy at 0 K.

---

(section:k-hi-numerical)=
### Numerical Evaluation: $K(700\,\mathrm{K})$ from Spectroscopic Data

To see the microscopic-to-macroscopic bridge at work, we evaluate Eq. {eq}`eq:K-HI-full` at $T = 700\ \mathrm{K}$ using spectroscopic constants. The rotational and vibrational temperatures of $\mathrm{H_2}$ and $\mathrm{I_2}$ are those of the tables in [Section 2.6](../chapter-02/section-06.md) and [Section 2.7](../chapter-02/section-07.md). The values for HI are computed in the same way, from $\tilde{\nu} = 2309.01\ \mathrm{cm^{-1}}$ and $\tilde{B} = 6.4264\ \mathrm{cm^{-1}}$ ([NIST Chemistry WebBook](https://webbook.nist.gov/chemistry/)).

| Species | $M$ ($\mathrm{g\,mol^{-1}}$) | $\Theta_{\mathrm{rot}}$ (K) | $\Theta_{\mathrm{vib}}$ (K) | $D_0$ ($\mathrm{kJ\,mol^{-1}}$) | $\sigma$ | $g_1$ |
| :-- | --: | --: | --: | --: | --: | --: |
| $\mathrm{H_2}$ | 2.016 | 87.55 | 6332 | 432.07 | 2 | 1 |
| $\mathrm{I_2}$ | 253.808 | 0.05377 | 309 | 148.81 | 2 | 1 |
| $\mathrm{HI}$ | 127.912 | 9.25 | 3322 | 294.67 | 1 | 1 |

```{code-cell} ipython3
import numpy as np

R = 8.314462618      # J mol^-1 K^-1

# Spectroscopic / thermochemical constants
mol = {
    "H2": dict(M=2.016e-3,   Theta_rot=87.55,   Theta_vib=6332.0, D0=432.07e3, sigma=2),
    "I2": dict(M=253.808e-3, Theta_rot=0.05377, Theta_vib=309.0,  D0=148.81e3, sigma=2),
    "HI": dict(M=127.912e-3, Theta_rot=9.25,    Theta_vib=3322.0, D0=294.67e3, sigma=1),
}
T = 700.0  # K

# Mass factor (translational)
mass_factor = (mol["HI"]["M"]**2 / (mol["H2"]["M"] * mol["I2"]["M"]))**1.5

# Rotational factor
rot_factor = (mol["H2"]["sigma"] * mol["I2"]["sigma"] / mol["HI"]["sigma"]**2) \
             * (mol["H2"]["Theta_rot"] * mol["I2"]["Theta_rot"]
                / mol["HI"]["Theta_rot"]**2)

# Vibrational factor: "hot-mode" (1 - exp(-Th_vib/T)) factors
def vfac(Tvib): return 1 - np.exp(-Tvib/T)
vib_factor = (vfac(mol["H2"]["Theta_vib"]) * vfac(mol["I2"]["Theta_vib"])
              / vfac(mol["HI"]["Theta_vib"])**2)

# Bond-energy factor: exp((2 D0_HI - D0_H2 - D0_I2)/RT)
dD0 = 2*mol["HI"]["D0"] - mol["H2"]["D0"] - mol["I2"]["D0"]  # J/mol
bond_factor = np.exp(dD0 / (R*T))

K = mass_factor * rot_factor * vib_factor * bond_factor

print(f"T = {T:.0f} K")
print(f"  Mass (translational) factor:       {mass_factor:10.3f}")
print(f"  Rotational factor:                 {rot_factor:10.4f}")
print(f"  Vibrational factor:                {vib_factor:10.4f}")
print(f"  Bond-energy factor exp(dD0/RT):    {bond_factor:10.4f}")
print(f"  2 D0(HI) - D0(H2) - D0(I2):        {dD0/1e3:+.2f} kJ/mol")
print()
print(f"  K(700 K) from molecular constants: {K:.1f}")
print(f"  K(700 K) from the NIST-JANAF tables:  56.8")
print(f"  Measured K (Bodenstein, ~700 K):      ~54")
```

**Discussion.** The predicted $K \approx 62$ compares with $K = 56.8$ from the NIST-JANAF tables at 700 K and with a measured value near 54. Agreement at the 10 to 15 % level is about what one should expect from the rigid-rotor/harmonic-oscillator approximation. The largest single correction comes from the rotation of $\mathrm{H_2}$. With $\Theta_{\mathrm{rot}}/T = 0.125$, the high-temperature form $q_{\mathrm{rot}} = T/(\sigma\Theta_{\mathrm{rot}})$ of [Section 2.7](../chapter-02/section-07.md) lies about 4 % below the sum over the rotational levels, and correcting it lowers $K$ by the same 4 %. Corrections for anharmonicity, vibration–rotation coupling, and centrifugal distortion (none of which are in our model) each shift $K$ by a few percent or less.

More importantly, the decomposition makes the physics transparent:

- The **translational mass factor** ($\approx 181$) is the largest single contribution. The ratio $m_{\mathrm{HI}}^2/(m_{\mathrm{H_2}}m_{\mathrm{I_2}}) = (m_{\mathrm{H}}+m_{\mathrm{I}})^2/(4 m_{\mathrm{H}} m_{\mathrm{I}})$ is driven far from unity when the two atomic masses differ by a lot — for H and I, $m_{\mathrm{I}}/m_{\mathrm{H}} \approx 126$, so $m_{\mathrm{HI}}^2 \gg 4\,m_{\mathrm{H}} m_{\mathrm{I}}$. If instead the atoms had similar masses (e.g., $\mathrm{D_2+X_2}\rightleftharpoons 2\mathrm{DX}$ with $m_{\mathrm{D}}\sim m_{\mathrm{X}}$), this factor would be close to 1.
- The **rotational factor** ($\approx 0.22$) is well below unity because I$_2$ is unusually heavy and has an extraordinarily small $\Theta_{\mathrm{rot}} \approx 0.05$ K; the corresponding $q_{\mathrm{rot}}(\mathrm{I_2})$ is very large, which favors keeping mass in I$_2$ rather than redistributing it into HI. The factor-of-4 symmetry-number prefactor (from $\sigma_{\mathrm{H_2}}\sigma_{\mathrm{I_2}}/\sigma_{\mathrm{HI}}^2$) nudges equilibrium toward products but doesn't overcome this.
- The **vibrational factor** ($\approx 0.36$) is also below unity, and it comes almost entirely from I$_2$. The stiff stretches of H$_2$ and HI are still nearly frozen out at 700 K, so their factors $(1-e^{-\Theta_{\mathrm{vib}}/T})$ are $0.9999$ and $0.991$. The low-frequency mode of I$_2$ is partially excited, and its factor, which stands in the numerator, is $0.357$.
- The **bond-energy factor** ($\approx 4.3$) contains the "chemistry": $2D_0^{\mathrm{HI}} - D_0^{\mathrm{H_2}} - D_0^{\mathrm{I_2}} \approx +8.5\ \mathrm{kJ\,mol^{-1}}$ says the two HI bonds together are slightly stronger than one H–H plus one I–I bond, a modest exothermicity that at 700 K gives just over a factor of 4 in $K$.

The calculation can be repeated at other temperatures. [](#fig:K-HI) compares Eq. {eq}`eq:K-HI-full` with the NIST-JANAF tables between 300 and 1500 K and shows how each of the four factors changes with temperature.

```{code-cell} ipython3
:label: fig:K-HI
:caption: (a) Equilibrium constant of $\mathrm{H_2 + I_2 \rightleftharpoons 2HI}$ from Eq. [](#eq:K-HI-full) (line) and from the NIST-JANAF tables (points). The formula lies 8 to 13 % above the tabulated values over the whole range. (b) The four factors of Eq. [](#eq:K-HI-full). In this model the mass and rotational factors do not depend on temperature. The bond-energy factor falls from 30 at 300 K to 2.0 at 1500 K, and the vibrational factor from 0.64 to 0.23.
:tags: [hide-input]

import sys

import matplotlib.pyplot as plt

sys.path.append("..")            # the book folder, which holds figstyle
import figstyle as fs

def hot(Theta_vib, T):
    return 1 - np.exp(-Theta_vib / T)

def factors(T):
    # The four factors of the molecular formula for K at temperature T (K)
    vib = (hot(mol["H2"]["Theta_vib"], T) * hot(mol["I2"]["Theta_vib"], T)
           / hot(mol["HI"]["Theta_vib"], T)**2)
    bond = np.exp(dD0 / (R * T))
    return mass_factor * np.ones_like(T), rot_factor * np.ones_like(T), vib, bond

# K from the NIST-JANAF tables: exp(-[2 dfG(HI) - dfG(I2, g)] / RT)
T_ref = np.array([300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500])
K_ref = np.array([673.9, 245.9, 128.2, 80.57, 56.81, 43.25, 34.74, 29.05, 25.04, 22.09,
                  19.86, 18.11, 16.71])

T_plot = np.linspace(300, 1500, 300)
f_mass, f_rot, f_vib, f_bond = factors(T_plot)



def draw_constant(c):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(fs.WIDTH, 3.4), constrained_layout=True)

    # (a) The molecular formula against the tables
    ax1.semilogy(T_plot, f_mass * f_rot * f_vib * f_bond, color=c.blue, lw=fs.LINE)
    ax1.semilogy(T_ref, K_ref, "o", color=c.gray, ms=5.5)
    ax1.text(520, 230, "molecular formula", color=c.blue)
    ax1.text(330, 30, "NIST-JANAF", color=c.muted)
    ax1.set(xlabel="Temperature (K)", ylabel="$K$")
    ax1.set_yticks([20, 50, 100, 200, 500], ["20", "50", "100", "200", "500"])
    ax1.minorticks_off()
    fs.panel(ax1, "a", "Equilibrium constant")

    # (b) The four factors of the formula
    ax2.axhline(1, color=c.gray, lw=1.2, ls=(0, (4, 3)))
    for factor, label, color, y in ((f_mass, "mass", c.blue, 220), (f_bond, "bond energy", c.purple, 5.2),
                                    (f_vib, "vibrational", c.vermilion, 0.47), (f_rot, "rotational", c.green, 0.135)):
        ax2.semilogy(T_plot, factor, color=color, lw=fs.LINE)
        ax2.text(900, y, label, color=color)
    ax2.set(ylim=(0.1, 400), xlabel="Temperature (K)", ylabel="Factor")
    fs.panel(ax2, "b", "The four factors")
    return fig


fs.show(draw_constant, alt="Two panels for the equilibrium of hydrogen, iodine, and hydrogen iodide between 300 and "
                           "1500 K. Panel a plots the equilibrium constant on a logarithmic axis. It falls from about "
                           "700 at 300 K to about 18 at 1500 K, and the curve from the molecular formula lies slightly "
                           "above the points from the NIST-JANAF tables. Panel b plots the four factors of the "
                           "formula. The mass factor, near 180, and the rotational factor, near 0.14, do not change "
                           "with temperature. The bond-energy factor falls from about 30 toward 2, and the "
                           "vibrational factor falls from near 1 to about 0.3.")
```

The lesson: every digit of $K$ can be traced back to a specific molecular property, and the relative sizes of the four factors tell us which part of molecular structure matters most for this particular equilibrium. A student who computed $\Delta_r G^\circ$ from thermochemical tables ([Computing $K$ from Tabulated Thermochemistry](#section:k-from-tables) in [Section 7.1](section-01.md)) and a student who computed $K$ from spectroscopic constants (this subsection) would arrive at nearly the same number by very different routes — and each route teaches something different about *why* the equilibrium sits where it does.

---

## Worked Example

### Why Volume Cancels When $\Delta\nu=0$: $\mathrm{H_2+I_2\rightleftharpoons 2HI}$

Writing each single-molecule partition function as $q_i = q_{\mathrm{trans},i}\,q_{\mathrm{int},i}$ with $q_{\mathrm{trans},i}=V/\Lambda_i^3$ (and $q_{\mathrm{int},i}$ collecting the rotational, vibrational, and electronic factors, all of which are $V$-independent):

```{math}
K = \frac{q_{\mathrm{HI}}^2}{q_{\mathrm{H_2}}\,q_{\mathrm{I_2}}}
= \frac{(V/\Lambda_{\mathrm{HI}}^3)^2\,q_{\mathrm{int,HI}}^2}{(V/\Lambda_{\mathrm{H_2}}^3)(V/\Lambda_{\mathrm{I_2}}^3)\,q_{\mathrm{int,H_2}}\,q_{\mathrm{int,I_2}}}.
```

Because $\Delta\nu = 2-1-1 = 0$, the volume factors cancel:

```{math}
\frac{V^2}{V\cdot V}=1,
```

leaving

```{math}
K=\left(\frac{\Lambda_{\mathrm{H_2}}^3\Lambda_{\mathrm{I_2}}^3}{\Lambda_{\mathrm{HI}}^6}\right)\,
\frac{q_{\mathrm{int,HI}}^2}{q_{\mathrm{int,H_2}}\,q_{\mathrm{int,I_2}}}.
```

Since $\Lambda\propto m^{-1/2}$, the translational contribution collapses to the mass factor $\left(m_{\mathrm{HI}}^2/(m_{\mathrm{H_2}}m_{\mathrm{I_2}})\right)^{3/2}$ seen in Eq. {eq}`eq:K-HI-full`.

**Result.** For reactions with $\Delta\nu=0$, $K$ is independent of volume (and pressure) in the ideal-gas model; microscopic physics enters through masses and internal partition functions only. For reactions with $\Delta\nu\ne 0$, the uncancelled factors of $V$ reappear as the $(k_{\mathrm B}T/P^\circ)^{\Delta\nu}$ conversion in Eq. {eq}`eq:Kp-Kc`.

## Concept Checks

1. Why is the Helmholtz free energy the relevant potential for equilibrium at constant $T,V$, and why does the resulting condition $\sum_i\nu_i\mu_i=0$ match the condition derived in [Section 7.1](section-01.md) at constant $T,P$?
2. Where does Stirling's approximation enter the derivation of $\mu_i = -RT\ln(q_i/N_i)$, and what physical limit makes it reasonable?
3. In Eq. {eq}`eq:q-diatomic`, how would the vibrational factor change if we chose the energy zero at the bottom of the potential well rather than at $v=0$? Show that $K$ is unchanged by this rechoice (as it must be).
4. Vibrational contributions to $K$ often matter most at high temperatures. Why? Reason from the temperature dependence of $(1-e^{-\Theta_{\mathrm{vib}}/T})$.
5. Suppose all three species in a $\Delta\nu=0$ reaction had identical masses, moments of inertia, vibrational frequencies, and $g_1$ values, but different $D_0$ values. What would $K$ reduce to, and what does this limit teach you about when "bond counting" arguments are quantitatively trustworthy?

## Key Takeaways

- Microscopic partition functions provide a route to equilibrium constants by linking $A = -k_{\mathrm B}T\ln Q$ ([Section 4.3](../chapter-04/section-03.md)) to chemical potentials via $\mu_i = (\partial A/\partial N_i)_{T,V}$.
- For ideal-gas mixtures, $\mu_i = -RT\ln(q_i/N_i)$ (per-mole), and the equilibrium condition $\sum_i\nu_i\mu_i=0$ delivers $K$ as a ratio of single-molecule partition functions.
- Translational, rotational, vibrational, and electronic factors contribute in separable ways under the common rigid-rotor/harmonic-oscillator approximations of Module 2.
- For $\Delta\nu=0$ gas reactions, volume (and hence pressure) dependence cancels in the ideal-gas equilibrium constant; for $\Delta\nu\ne 0$, the $(k_{\mathrm B}T/P^\circ)^{\Delta\nu}$ factor converts between $K_c$ and $K_p$.
- For the $\mathrm{H_2+I_2\rightleftharpoons 2HI}$ benchmark, the four factors in Eq. {eq}`eq:K-HI-full` predict $K(700\,\mathrm{K})\approx 62$, within 10 to 15 % of the tabulated value (56.8) and of the measured value (about 54) — a concrete illustration of how molecular structure alone determines the position of equilibrium.
