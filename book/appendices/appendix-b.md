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

# Appendix B: Roots of Polynomials

Let

```{math}
p(z)=a_n z^n+a_{n-1}z^{n-1}+\cdots+a_1 z+a_0,\qquad a_n\neq 0,
```

with $a_k\in\mathbb{R}$ or $\mathbb{C}$. A **root** $r$ satisfies $p(r)=0$. If the coefficients are real, nonreal roots occur in complex-conjugate pairs.

This appendix summarizes (i) the *structure* of analytic formulas for $n\le 4$, (ii) why numerical methods are essential for $n\ge 5$ (and often preferable even for $n\le 4$), and (iii) the main numerical approaches and their practical caveats.

---

## A.1 Analytical solutions for degrees $\boldsymbol{n\le 4}$: structure and limitations

### Linear ($n=1$)

For $a_1 z + a_0 = 0$,

```{math}
z=-\frac{a_0}{a_1}.
```

### Quadratic ($n=2$)

For $a z^2+b z+c=0$ with $a\neq 0$,

```{math}
z=\frac{-b\pm\sqrt{b^2-4ac}}{2a}.
```

````{admonition} Numerical caveat (cancellation)
:class: dropdown
If $b^2\gg 4ac$, one root suffers catastrophic cancellation in $(-b\pm\sqrt{b^2-4ac})$. A stable variant is:

```{math}
q=-\tfrac12\left(b+\operatorname{sign}(b)\sqrt{b^2-4ac}\right),\qquad
z_1=\frac{q}{a},\quad z_2=\frac{c}{q}.
```

This uses the fact that $z_1 z_2=c/a$ and avoids subtracting nearly equal numbers.
````

### Cubic ($n=3$): Cardano via "depressed cubic"

For $a z^3+b z^2+c z+d=0$, shift to remove the quadratic term:

```{math}
z=y-\frac{b}{3a}\quad\Rightarrow\quad y^3+py+q=0,
```

where $p,q$ are explicit functions of $a,b,c,d$. Cardano's method seeks $y=u+v$ with $u^3$ and $v^3$ solving a quadratic ("resolvent") equation, giving a closed-form expression in radicals.

```{admonition} Key limitation (casus irreducibilis)
:class: dropdown
When the cubic has three real roots, Cardano's radical expression generally passes through complex numbers even though the final roots are real. This is algebraically correct, but can be numerically delicate in finite precision.
```

```{admonition} Practical note
:class: dropdown
Even though cubic formulas exist, robust numerical evaluation is nontrivial; in floating-point computations, it is common to prefer iterative methods or library routines designed for numerical stability.
```

### Quartic ($n=4$): Ferrari via "depressed quartic" and resolvent cubic

For $a z^4+b z^3+c z^2+d z+e=0$, shift

```{math}
z=y-\frac{b}{4a}
```

to obtain a depressed quartic

```{math}
y^4+p y^2+q y+r=0.
```

Ferrari's method introduces an auxiliary parameter to factor the quartic into two quadratics; determining this parameter reduces to solving a **resolvent cubic**, followed by two quadratic solves.

```{admonition} Key limitation
:class: dropdown
Quartic formulas are long and can be numerically unstable due to cancellation and sensitivity when roots cluster. For numerical work, "closed form" does not automatically imply "good floating-point behavior."
```

<!-- ---

## A.2 Why numerical methods become necessary for higher degrees

1. **No general radical formula for $n\ge 5$.** No formula in radicals solves all quintics (and higher) in terms of the coefficients. Specific families can be solvable, but a general-purpose method must be numerical.

2. **Conditioning can dominate.** Even when a root can be expressed analytically, computing it accurately may be hard because roots can be **highly sensitive** to small coefficient perturbations. A standard local sensitivity estimate for a *simple* root $r$ is

   ```{math}
   \delta r \approx -\frac{\delta p(r)}{p'(r)}.
   ```

   Thus, if $|p'(r)|$ is small (near a multiple root or tight cluster), small perturbations $\delta p$ can cause large $\delta r$.

3. **Physical-science inputs are rarely exact.** Coefficients often come from measured parameters, truncated series, or discretizations. In that setting, numerical methods paired with **backward-error** checks are typically more meaningful than symbolic expressions.

---

## A.3 Core numerical strategies

A recurring theme: methods differ in (i) what they require (brackets? derivatives?), (ii) their convergence rate when they work, and (iii) how they fail.

### Efficient evaluation: Horner's rule (recommended baseline)

Direct evaluation of (p(z)) from powers is inefficient and can amplify rounding error. Use Horner's form:

```{math}
p(z)=\Bigl(\bigl((a_n z+a_{n-1})z+a_{n-2}\bigr)z+\cdots+a_0\Bigr).
```

You can evaluate $p'(z)$ with a coupled recurrence (useful for Newton's method) without explicitly forming derivative coefficients:

```{math}
\begin{aligned}
b_n&=a_n,\quad &b_{k}&=a_k+z,\quad &b_{k+1},\quad &k=n-1,\dots,0,\\
c_n&=b_n,\quad &c_{k}&=b_k+z,\quad &c_{k+1},\quad &k=n-1,\dots,1,
\end{aligned}
```

then $p(z)=b_0$ and $p'(z)=c_1$.

---

### A.3.1 Bracketing methods for **real** roots (robust, slower)

Assume $p:\mathbb{R}\to\mathbb{R}$ is continuous. If you find an interval $[a,b]$ with a sign change,

```{math}
p(a)\,p(b)<0,
```

then at least one root lies in $[a,b]$ (Intermediate Value Theorem).

#### Bisection

Update $[a,b]$ by halving:

```{math}
m=\frac{a+b}{2},\quad \text{replace the endpoint that has the same sign as }p(m).
```

* **Convergence:** guaranteed, linear. After $k$ steps,

   ```{math}
   |b_k-a_k|=\frac{|b_0-a_0|}{2^k}.
   ```

* **Strength:** cannot diverge if a valid bracket exists.
* **Failure mode:** needs a sign change; **even-multiplicity roots** (touching the axis) do *not* change sign, so bisection cannot detect them from sign alone.

#### Regula falsi (false position) and Brent-type hybrids

Regula falsi uses the secant line on $[a,b]$ but maintains a bracket. Pure regula falsi can stagnate if one endpoint barely moves.

**Best practice (real roots):** use a safeguarded hybrid such as **Brent–Dekker** (bisection + secant + inverse quadratic interpolation). It combines reliability (bracketing) with faster local steps, and is often the default "industrial" choice for 1D real root finding.

---

### A.3.2 Newton–Raphson (fast locally; sensitive globally)

Newton iteration:

```{math}
z_{k+1}=z_k-\frac{p(z_k)}{p'(z_k)}.
```

#### Convergence behavior

If $r$ is a **simple root** ($p(r)=0$, $p'(r)\neq 0$) and $z_0$ is sufficiently close to $r$, then Newton converges **quadratically**:

```{math}
e_{k+1}\approx C e_k^2,\qquad e_k=z_k-r.
```

#### Failure modes and caveats

* **Bad initial guess:** Newton is not globally convergent; it can diverge or converge to an unintended root.
* **Derivative near zero:** if $p'(z_k)$ is small, the step can be huge and unstable.
* **Multiple roots:** if $r$ has multiplicity $m>1$, then $p'(r)=0$ and Newton loses quadratic convergence; typically it becomes **linear**:

   ```{math}
   e_{k+1}\approx \left(1-\frac1m\right)e_k.
   ```

   If $m$ is known (rare in noisy data), a modified step restores quadratic behavior:

   ```{math}
   z_{k+1}=z_k- m\,\frac{p(z_k)}{p'(z_k)}.
   ```

* **Complex plane dynamics:** for complex roots, Newton's basins of attraction can be intricate; iterates may wander before converging.

#### Practical guidance

* Prefer Newton when you have a **good initial guess**, or can obtain one from physics (asymptotics, perturbation limits, continuity in parameters).
* Use **damping/safeguards**: if a full Newton step increases $|p|$, reduce step length (line search) or fall back to a bracketing step (for real roots).

---

### A.3.3 Secant method (derivative-free, superlinear, less stable)

Secant iteration uses two previous points:

```{math}
z_{k+1} = z_k - p(z_k)\,\frac{z_k-z_{k-1}}{p(z_k)-p(z_{k-1})}.
```

* **Convergence:** superlinear with order $\varphi=\tfrac{1+\sqrt5}{2}\approx 1.618$ near a simple root.
* **Strength:** no derivative required; often competitive when $p'$ is expensive or noisy.
* **Failure modes:**

  * If $p(z_k)\approx p(z_{k-1})$, the denominator becomes small and the step can blow up.
  * It does **not** preserve a bracket (unless modified), so it can jump away from the root.
* **Best practice:** use secant inside a safeguarded routine (e.g., Brent), or monitor steps and fall back to bisection when needed.

---

### A.3.4 Companion-matrix/eigenvalue approach (all roots at once)

For polynomials, a powerful global approach is to convert root finding into an eigenvalue problem.

First scale to a **monic** polynomial:

```{math}
p(z)=a_n\left(z^n+c_{n-1}z^{n-1}+\cdots+c_0\right),\qquad c_k=\frac{a_k}{a_n}.
```

Define the $n\times n$ companion matrix

```{math}
C=
\begin{pmatrix}
0 & 0 & \cdots & 0 & -c_0\\
1 & 0 & \cdots & 0 & -c_1\\
0 & 1 & \cdots & 0 & -c_2\\
\vdots & & \ddots & \vdots & \vdots\\
0 & 0 & \cdots & 1 & -c_{n-1}
\end{pmatrix}.
```

Then the eigenvalues of $C$ are exactly the roots of $p$.

#### Pros

* Returns **all roots** (real and complex) simultaneously.
* Leverages mature eigenvalue algorithms (e.g., QR iterations) that are highly optimized and generally robust.

#### Cons/caveats

* **Conditioning:** the polynomial-to-root map can be ill-conditioned; eigenvalue routines cannot "fix" fundamental sensitivity. Clusters and near-multiple roots remain difficult.
* **Scaling matters:** large variations in coefficient magnitudes can lead to numerical issues. Rescale $z$ to reduce dynamic range before forming $C$.
* **Cost:** dense eigenvalue computation is typically $O(n^3)$, which can be expensive for very large $n$.

````{admonition} Practical guidance
:class: dropdown
* Use this approach when you need **the full root set** and $n$ is moderate.
* After computing eigenvalues $\hat r_i$, assess quality via a **backward error**:

  ```{math}
  \eta(\hat r)=\frac{|p(\hat r)|}{\sum_{k=0}^n |a_k|\,|\hat r|^k}.
  ```

  Small $\eta$ indicates $\hat r$ is the exact root of a nearby polynomial with small relative coefficient perturbations.
````

---

## A.4 Summary comparison (what to choose and why)

| Method            |                                 Needs |  Typical convergence (simple root) | Strength                      | Common failure mode                                                          |
| ----------------- | ------------------------------------: | ---------------------------------: | ----------------------------- | ---------------------------------------------------------------------------- |
| Bisection         | real bracket $[a,b]$ with sign change |                             linear | guaranteed, simple            | cannot start without sign change; misses even-multiplicity roots             |
| Brent-type hybrid |                               bracket |            superlinear in practice | "best of both": robust + fast | still requires a bracket; even-multiplicity roots need extra logic           |
| Newton–Raphson    |                  initial guess + $p'$ |                  quadratic (local) | very fast when it works       | diverges with poor guess; unstable if $p'\approx 0$; slow for multiple roots |
| Secant            |                   two initial guesses |               order $\approx1.618$ | no derivative; often fast     | unstable if denominator small; not globally convergent                       |
| Companion matrix  |                          coefficients | depends on eigenvalue conditioning | all roots at once             | sensitive to scaling/ill-conditioning; $O(n^3)$ cost                         |

---

## A.5 Key caveats that matter in practice

### 1. Multiple roots and near-multiple roots

* A root $r$ has multiplicity $m>1$ if $p(r)=p'(r)=\cdots=p^{(m-1)}(r)=0$.
* These are **ill-conditioned**: tiny perturbations can split a multiple root into a cluster.
* Newton slows dramatically unless modified; bracketing by sign change may fail.

### 2. Initial guesses and basins of attraction

* Iterative methods (Newton, secant) are **local**; convergence depends on the starting point.
* In parameter-dependent problems, continuation is effective: use the solution at parameter $\lambda$ as the initial guess for nearby $(\lambda+\Delta\lambda)$.

### 3. Deflation: useful but potentially dangerous

Once a root $\hat r$ is found, one can reduce the degree by dividing:

```{math}
p(z)=(z-\hat r)\,q(z)+\text{(remainder)}.
```

However, **deflation can amplify errors**, especially when $\hat r$ is inaccurate or roots are clustered. If you must deflate, refine $\hat r$ first and monitor backward error; otherwise, prefer methods that compute all roots simultaneously (companion matrix or simultaneous refinement methods).

### 4. Scaling and shifting often determine success

A simple rescaling $z=\alpha x$ changes coefficients and can greatly improve numerical behavior. Choose $\alpha$ so that "typical" root magnitudes are $O(1)$ when possible (using bounds; see below).

### 5. Root bounds (useful for bracketing and sanity checks)

A common bound (Cauchy) for any root $r$ of $p$ is

```{math}
|r|\le 1+\max_{0\le k\le n-1}\left|\frac{a_k}{a_n}\right|.
```

This gives a disk containing all complex roots; for real-root searches, it can guide interval selection.

---

## A.6 Practical workflow recommendations (physical-science setting)

1. **Preprocess**

   * Scale to monic form if feasible.
   * Rescale variables to reduce coefficient dynamic range.
   * Use Horner evaluation for $p$ (and $p'$ if needed).

2. **If you need one real root in a known interval**

   * Use a **bracketing + hybrid** method (e.g., Brent-type).
   * If a root may be even-multiplicity, consider also checking for local minima where $p\approx 0$ without a sign change (e.g., detect $|p|$ small and $p'\approx 0$).

3. **If you have a good initial guess**

   * Use **Newton** with damping/safeguards; monitor $|p|$ and step size.

4. **If you need all roots**

   * Use a **companion-matrix/eigenvalue** method (moderate $n$), followed by residual/backward-error checks.
   * Be cautious with clustered roots; interpret results through conditioning/backward error.

5. **Stopping criteria**
   Use both a step-based and a residual-based test, e.g.

   ```{math}
   |z_{k+1}-z_k|\le \varepsilon(1+|z_{k+1}|),
   \quad\text{and}\quad
   \eta(z_{k+1})\le \varepsilon,
   ```

   where $\eta$ is the backward error defined above. Relying on $|p(z)|$ alone can be misleading due to scaling and cancellation. -->
