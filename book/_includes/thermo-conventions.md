#

````{admonition} Course convention: heat & work signs
:class: important

We use the **chemistry sign convention** throughout CHEM 4020/5020:

- **Heat** $q>0$ when heat is **absorbed by the system**.
- **Work** $w>0$ when work is **done on the system**.

**First Law (this convention)**

```{math}
\Delta U = q + w
\qquad\text{and}\qquad
 dU = \delta q + \delta w.
```

**$PV$ work (expansion/compression against an external pressure $P_{\mathrm{ext}}$)**

```{math}
\delta w_{PV} = -P_{\mathrm{ext}}\,dV
\qquad\Rightarrow\qquad
w_{PV} = -\int P_{\mathrm{ext}}\,dV.
```

So:

- **Expansion** ($dV>0$) $\Rightarrow$ $w_{PV}<0$ (**system does work on surroundings**).
- **Compression** ($dV<0$) $\Rightarrow$ $w_{PV}>0$ (**surroundings do work on system**).

For a *reversible* $PV$ process, $P_{\mathrm{ext}}=P$ (the system pressure), so you will often see
$w_{PV,\,\mathrm{rev}}=-\int P\,dV$.

**Translating to the “work done *by* the system” convention**

Some texts define $W_{\text{by}}>0$ when work is **done by the system**. The two conventions are related by

```{math}
W_{\text{by}} = -w\;\;\;\text{(and likewise }\delta W_{\text{by}} = -\delta w\text{)}.
```

With $W_{\text{by}}$ the First Law is equivalently written as

```{math}
\Delta U = q - W_{\text{by}}.
```
````
