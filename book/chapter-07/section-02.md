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


# 7.2. Reaction Spontaneity from Microscopic Properties

[Course-wide Conventions & Notation](../notation.md)

## Overview

In Section 7.1 we described equilibrium and spontaneity using macroscopic thermodynamics
(e.g., $\Delta_r G = \Delta_r G^\circ + RT\ln Q$ and $\Delta_r G=0$ at equilibrium).

In this section we connect those macroscopic ideas to **microscopic** information by using
the **canonical partition function** of an ideal-gas mixture. The key goal is:

- use statistical mechanics to compute the **equilibrium constant** ($K$) from molecular properties  
  (translation, rotation, vibration, electronic structure),
- and therefore predict whether a reaction is product-favored or reactant-favored at a given $T$.

```{admonition} Notational warning
:class: warning

- In **thermodynamics**, $Q$ denotes the **reaction quotient** (Section 7.1).
- In **statistical mechanics**, $Q$ (often $Z$) is also used for the **canonical partition function**.
  In this section we write $\mathcal{Q}$ for partition functions to avoid confusion.
```

---

## 7.2.1 The right potential for equilibrium at constant $T,V$

Consider a general gas-phase reaction

```{math}
\nu_A A(g) + \nu_B B(g) \rightleftharpoons \nu_Y Y(g) + \nu_Z Z(g),
```

and assume **ideal-gas behavior** (as in the lecture notes).

For problems at **constant temperature and volume**, the natural thermodynamic potential is the
**Helmholtz free energy** $A$, whose differential is

```{math}
dA = -S\,dT - P\,dV + \sum_i \mu_i\,dn_i.
```

At constant $T$ and $V$,

```{math}
dA = \sum_i \mu_i\,dn_i.
```

If the reaction progress is described by an extent of reaction $\xi$, then (with the usual sign convention)
$dn_i = \nu_i\,d\xi$, giving

```{math}
dA = \left(\sum_i \nu_i\,\mu_i\right)d\xi.
```

Equilibrium at fixed $T,V$ occurs at a minimum of $A$, so

```{math}
\left(\frac{\partial A}{\partial \xi}\right)_{T,V} = 0
\quad\Longrightarrow\quad
\boxed{\sum_i \nu_i\,\mu_i = 0.}
```

This is the **chemical equilibrium condition** written in terms of chemical potentials.

---

## 7.2.2 Partition function of a mixture of ideal gases

For a mixture of ideal gases, the canonical partition function factorizes:

```{math}
\boxed{
\mathcal{Q}(T,V,\{N_i\}) = \prod_i \mathcal{Q}_i(T,V,N_i)
}
\qquad
\text{where } \{N_i\}=\{N_A,N_B,N_Y,N_Z,\dots\}.
```

For each species $i$, an ideal gas of **indistinguishable** particles has

```{math}
\boxed{
\mathcal{Q}_i(T,V,N_i)=\frac{q_i(T,V)^{N_i}}{N_i!},
}
```

where $q_i$ is the **single-particle partition function** for species $i$.

---

## 7.2.3 Connecting $\mathcal{Q}$ to chemical potentials

Statistical mechanics connects Helmholtz free energy to the partition function:

```{math}
\boxed{
A = -RT\ln \mathcal{Q}.
}
```

The chemical potential follows from

```{math}
\mu_i = \left(\frac{\partial A}{\partial N_i}\right)_{T,V,\{N_{j\neq i}\}}
= -RT\left(\frac{\partial \ln\mathcal{Q}}{\partial N_i}\right)_{T,V,\{N_{j\neq i}\}}.
```

Using the factorized form,

```{math}
\ln\mathcal{Q} = \sum_i \left(N_i\ln q_i - \ln N_i!\right).
```

Applying Stirling’s approximation (valid for large $N_i$),

```{math}
\ln N_i! \approx N_i\ln N_i - N_i,
```

one finds

```{math}
\boxed{
\mu_i = -RT\ln\left(\frac{q_i}{N_i}\right)
= RT\ln\left(\frac{N_i}{q_i}\right).
}
```

---

## 7.2.4 Equilibrium constants from partition functions

Insert $\mu_i = -RT\ln(q_i/N_i)$ into the equilibrium condition $\sum_i \nu_i\mu_i=0$:

```{math}
\sum_i \nu_i\ln\left(\frac{q_i}{N_i}\right)=0
\quad\Longrightarrow\quad
\prod_i\left(\frac{q_i}{N_i}\right)^{\nu_i}=1.
```

For the reaction $\nu_A A + \nu_B B \rightleftharpoons \nu_Y Y + \nu_Z Z$, this becomes

```{math}
\ln\left[
\frac{
\left(q_Y/N_Y\right)^{\nu_Y}\,\left(q_Z/N_Z\right)^{\nu_Z}
}{
\left(q_A/N_A\right)^{\nu_A}\,\left(q_B/N_B\right)^{\nu_B}
}
\right]=0,
```

or equivalently,

```{math}
\boxed{
\frac{N_Y^{\nu_Y}N_Z^{\nu_Z}}{N_A^{\nu_A}N_B^{\nu_B}}
=
\frac{q_Y^{\nu_Y}q_Z^{\nu_Z}}{q_A^{\nu_A}q_B^{\nu_B}}.
}
```

### Relating to $K_c$ and $K_p$

Define number concentrations $c_i = N_i/V$. Then

```{math}
K_c \equiv \frac{c_Y^{\nu_Y}c_Z^{\nu_Z}}{c_A^{\nu_A}c_B^{\nu_B}}.
```

For ideal gases, the pressure-based and concentration-based equilibrium constants satisfy

```{math}
\boxed{
K_p = K_c\,(RT)^{\Delta\nu},
}
\qquad
\Delta\nu \equiv (\nu_Y+\nu_Z)-(\nu_A+\nu_B).
```

So, **once you can compute** the single-particle partition functions $q_i$,
you can compute $K_c$ (and then $K_p$).

---

## 7.2.5 Example: spontaneity of $\mathrm{HI(g)}$ formation

Consider

```{math}
\mathrm{H_2(g) + I_2(g) \rightleftharpoons 2HI(g)}.
```

Here,

```{math}
\Delta\nu = 2 - 1 - 1 = 0,
```

so

```{math}
\boxed{K_p = K_c.}
```

From the partition-function result,

```{math}
\boxed{
K = \frac{c_{HI}^2}{c_{H_2}c_{I_2}}
= \frac{q_{HI}^2}{q_{H_2}\,q_{I_2}}.
}
```

This explicitly shows how equilibrium (and therefore spontaneity trends) can be predicted from
microscopic molecular properties encoded in the $q_i$.

---

## 7.2.6 Single-particle partition function for a diatomic ideal gas

For a diatomic molecule, the single-particle partition function is typically written as a product of
translational and internal contributions:

```{math}
q = q_{\mathrm{trans}}\,q_{\mathrm{rot}}\,q_{\mathrm{vib}}\,q_{\mathrm{el}}.
```

In the (common) approximations used in the lecture, this can be written compactly as

```{math}
\boxed{
q_{\text{diatomic}}
=
\frac{V}{\Lambda^3}
\;
\frac{T}{\sigma\,\theta_{\mathrm{rot}}}
\;
\frac{e^{-\theta_{\mathrm{vib}}/(2T)}}{1-e^{-\theta_{\mathrm{vib}}/T}}
\;
g_{\mathrm{el}}\,e^{D_0/(RT)}.
}
```

The brackets in the lecture notes correspond to:

- **translational:** $q_{\mathrm{trans}}=V/\Lambda^3$
- **rotational:** $q_{\mathrm{rot}}\approx T/(\sigma\theta_{\mathrm{rot}})$
- **vibrational:** $q_{\mathrm{vib}}=e^{-\theta_{\mathrm{vib}}/(2T)}/(1-e^{-\theta_{\mathrm{vib}}/T})$
- **electronic:** $q_{\mathrm{el}}\approx g_{\mathrm{el}}e^{D_0/(RT)}$

### Thermal de Broglie wavelength

```{math}
\boxed{
\Lambda \equiv \left(\frac{h^2}{2\pi m k_B T}\right)^{1/2}
}
```

### Notes on parameters

- $m$: molecular mass  
- $\sigma$: rotational symmetry number  
  ($\sigma=2$ for homonuclear diatomics like $\mathrm{H_2}$ and $\mathrm{I_2}$; $\sigma=1$ for heteronuclear $\mathrm{HI}$)
- $\theta_{\mathrm{rot}}$: characteristic rotational temperature
- $\theta_{\mathrm{vib}}$: characteristic vibrational temperature
- $g_{\mathrm{el}}$: electronic degeneracy factor
- $D_0$: dissociation energy referenced from $v=0$

> Depending on the choice of energy zero, the vibrational factor may also be written as
> $q_{\mathrm{vib}} = 1/(1-e^{-\theta_{\mathrm{vib}}/T})$ if the zero-point energy is absorbed into $D_0$.
> The final expression for $K$ below is consistent with the lecture form.

---

## 7.2.7 $K$ for $\mathrm{H_2 + I_2 \rightleftharpoons 2HI}$ in terms of molecular parameters

Using $K = q_{HI}^2/(q_{H_2}q_{I_2})$ and the diatomic-gas forms above, the lecture notes give

```{math}
\boxed{
K
=
\left(\frac{m_{HI}^2}{m_{H_2}m_{I_2}}\right)^{3/2}
\,
\frac{4\,\theta_{\mathrm{rot}}^{H_2}\,\theta_{\mathrm{rot}}^{I_2}}{\left(\theta_{\mathrm{rot}}^{HI}\right)^2}
\,
\frac{\left(1-e^{-\theta_{\mathrm{vib}}^{H_2}/T}\right)\left(1-e^{-\theta_{\mathrm{vib}}^{I_2}/T}\right)}{\left(1-e^{-\theta_{\mathrm{vib}}^{HI}/T}\right)^2}
\,
\exp\left(\frac{2D_0^{HI}-D_0^{H_2}-D_0^{I_2}}{RT}\right).
}
```

This expression makes the microscopic physics very explicit:

- **translation** contributes the $m^{3/2}$ mass dependence,
- **rotation** contributes symmetry-number and $\theta_{\mathrm{rot}}$ factors (the factor of 4 comes from $\sigma_{H_2}=\sigma_{I_2}=2$ vs. $\sigma_{HI}=1$),
- **vibration** contributes the $\theta_{\mathrm{vib}}$ factors,
- **bond strengths** enter through the dissociation energies $D_0$ in the exponential.

---

## Key takeaways

1. At constant $T,V$, equilibrium corresponds to minimizing Helmholtz free energy $A$, giving
   $\sum_i \nu_i\mu_i=0$.
2. For an ideal-gas mixture, the canonical partition function factorizes and
   $\mathcal{Q}_i = q_i^{N_i}/N_i!$.
3. The chemical potential can be obtained directly from $\mathcal{Q}$:
   $\mu_i = -RT\ln(q_i/N_i)$.
4. The equilibrium constant can be expressed in terms of single-particle partition functions
   (and therefore in terms of molecular masses and energy level spacings).
