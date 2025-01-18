# SimpleDualAscent

> [!WARNING]
> This package is still under construction. Use with caution.

Inspired by [SimplePDHG.jl](https://github.com/Shuvomoy/SimplePDHG.jl), this package provides a minimal implementation of
the dual ascent algorithm for solving LPs of the form:

$$
\begin{align*}
\min_{x} \quad &  c^T x \\
\text{s.t.} \quad &  Ax = b \\
&  l \leq x \leq u
\end{align*}
$$

where $l<u$ and $l,u\in\mathbb{R}^n$ (consistent and finite).

> [!NOTE]
> Since not all LPs can be bridged to the above form, SimpleDualAscent may raise errors upon `optimize!`, since we only check for consistent/finite condition then.
> 
> This may require users to reformulate their problem such that JuMP bridges it correctly. In general, it is good practice to define variables with upper/lower bounds, e.g. using `@variable(model, l ≤ x ≤ u)`, and using only equality constraints. If an inequality constraint is required, manually define a slack variable (with double-sided bounds) and reformulate manually so JuMP sees an equality constraint.
> 
> Future versions of SimpleDualAscent may include some logic to deduce bounds for slack variables coming from bridges.

## Use with JuMP

```julia
using JuMP, SimpleDualAscent

model = Model(SimpleDualAscent.Optimizer)
...
```