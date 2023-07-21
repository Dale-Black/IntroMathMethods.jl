### A Pluto.jl notebook ###
# v0.19.26

using Markdown
using InteractiveUtils

# ╔═╡ ae0f795d-3237-4fd7-9f3f-a6b143c928f3
using DrWatson

# ╔═╡ fea70653-d6c7-438d-8a4f-940b8e7f21c3
# ╠═╡ show_logs = false
@quickactivate "graduate-courses"

# ╔═╡ a129cc0c-73b6-4c18-ac17-d609d0b1bda4
using PlutoUI, Symbolics, SymbolicNumericIntegration, SymbolicUtils

# ╔═╡ 4d46902e-7245-4987-9639-658f10005066
include(srcdir("helpers.jl"));

# ╔═╡ b37912b9-a9b0-4d69-952f-4d280f4ae382
TableOfContents()

# ╔═╡ e4ab3255-4fde-458a-a430-7dd6b3e1d521
md"""
# Calculus of Many Variables

**Homework Problems**
- 3.1.5
- 3.1.6
- 3.1.7
- 3.2.4
"""

# ╔═╡ ecb118f8-3a4f-44b9-adc3-421d94cb7514
md"""
## Problem 3.1.5

Find the shortest distance from the origin to any point on the line ``x + 2y = 4`` by using Lagrange multipliers. Check this by more elementary means: by first finding the equation for the line which is perpendicular to the given line and passing through the origin.

"""

# ╔═╡ d85c77a5-4b36-46ed-99f0-cb631c34baae
md"""
```math
\begin{align*}
f(x, y) &= x^2 + y^2 \\
g(x, y) &= x + 2y - 4 = 0 \\
L(x, y, λ) &= x^2 + y^2 - λ(x + 2y - 4) \\
\frac{\partial L}{\partial x} &= 2x - λ = 0 \\
\frac{\partial L}{\partial y} &= 2y - 2λ = 0\\
\frac{\partial L}{\partial λ} &= -(x + 2y - 4) = 0
\end{align*}
```

Solving the system gives:
```math
\begin{align*}
x &= 4/5 \\
y &= 8/5 \\
\lambda &= 8/5 \\
\end{align*}
```
"""

# ╔═╡ 722bda22-2044-4e6c-b17b-df1fa9801faf
@variables x y

# ╔═╡ 09bd6290-b690-46e7-8a17-6c0ed0e5627f
"""
    lagrangian(f, g)

Solve a constrained optimization problem using the Lagrangian method.

Given an objective function `f` and a constraint `g`, find the critical points 
that optimize `f` subject to `g` using the Lagrangian approach.

The Lagrangian `L` is formed as:

L(x, y, λ) = f(x, y) - λ(g(x, y))

The optimal points are found by taking derivatives of `L` w.r.t. x, y, λ 
and solving the resulting system of equations.

# Arguments
- `f`: Objective function to optimize
- `g`: Constraint equation 

# Returns
- Optimal values for x, y, λ that maximize/minimize `f` subject to `g`

# Examples
```julia
using Symbolics

@variables x y λ
f = x + y
g = x^2 + y^2 - 1
lagrangian(f, g)
```
"""
function lagrangian(f, g)
	@variables λ
	L = f - λ*g
	
	dLdx = derivative(L, x)
	dLdy = derivative(L, y)
	dLdλ = derivative(L, λ)

	sys = [dLdx ~ 0, dLdy ~ 0, dLdλ ~ 0]
	
	return solve_for(sys, [x, y, λ])
end

# ╔═╡ fab2caac-b604-4fad-97c0-ff7e9579d9fc
let
	f = x^2 + y^2
	g = x + 2y - 4	
	x, y, _ = lagrangian(f, g)
	x, y
end

# ╔═╡ 9d52232c-1f7a-40fe-82eb-4f2efa7014d6
md"""
## Problem 3.1.6
Show using the Lagrangian, that among all rectangles of a given perimeter, the square has the greatest area
"""

# ╔═╡ b20edb50-db8c-45d3-95b6-70e0402cab0a
let
	@variables P
	f = x*y
	g = 2x + 2y - P
	x, y, _ = lagrangian(f, g)
	x, y
end

# ╔═╡ 2d99d71f-0293-44e1-bc01-0416f33fa7ef
md"""
## Problem 3.1.7

Consider ``N`` particles in a box. According to quantum mechanics, the energies of the particles are quantized to some set of values or energy levels (``\epsilon_1, \epsilon_2, ...``). Let ``n_i`` be the number of particles in level ``i`` with energy ``\epsilon_i``. The multiplicity or number of distinct rearrangements of the particles consistent with any given distribution ``n_i``, is given by

```math
\begin{align*}
W(n_1, n_2, \ldots) &= \frac{N!}{n_1!n_2!\ldots}
\end{align*}
```

The question is this: which distribution of particles, subject to the constraint that the total number equal ``N`` and the total energy equal ``E``, gives the biggest ``W``? Use below

```math
\begin{align*}
S &= \ln W \qquad (\text{Note: } \ln n! \approx n \ln n - n)
\end{align*}
```

Then:
- (1) write the constraints on the ``n_i``'s due to the total number ``N`` and energy ``E``
- (2) treat all ``n_i``'s as continuous variables
- (3) introduce Lagrange multipliers ``\alpha`` and ``\beta`` for ``N`` and ``E``
- (4) maximize ``S``
"""

# ╔═╡ 192cf150-6baf-4246-a0de-15ed1a93f8ce
md"""
Total particles and total energy (constraint equations)
```math
\begin{align*}
N &= \sum_i n_i \\
g_1 &= \sum_i n_i - N \\
E &= \sum_i n_i \epsilon_i \\
g_2 &= \sum_i n_i \epsilon_i - E
\end{align*}
```

Entropy (maximization equation)
```math
\begin{align*}
S &= \ln W \\
&= \ln(\frac{N!}{n_1! n_2! ...}) \\
&= \ln (N!) - \ln (n_1! n_2! ...) \\
&= N \ln N - N - [\ln(n_1 !) + \ln(n_2 !) + ...] \\
&= N \ln N - N - [n_1 \ln n_1 - n_1 + n_2 \ln n_2 - n_2 + ...] \\
&= N \ln N - N - [\sum_i (n_i \ln n_i - n_i)]
\end{align*}
```

Lagrangian
```math
\begin{align*}
L &= S - \alpha(g_1) - \beta(g_2) \\
&= N \ln N - N - [\sum_i (n_i \ln n_i - n_i)] - [\alpha (\sum_i n_i - N)] - [\beta (\sum_i n_i \epsilon_i - E)]
\end{align*}
```
"""

# ╔═╡ 889ed804-81bd-4e87-97d5-d5b6f4986a1f
md"""
Answer
```math
\begin{align*}
\frac{\partial L}{\partial n_i} &= 0 \\
\frac{\partial }{\partial n_i}[-\int(n_i \ln n_i) + \int n_i - \alpha \int n_i - \beta \int (n_i \epsilon_i)] &= 0 \\
-\ln n_i - \alpha - \beta \epsilon &= 0 \\
n_i &= \boxed{e^{- \alpha - \beta \epsilon}}
\end{align*}
```
"""

# ╔═╡ 50b3d6e1-803a-497e-a8f0-01dbe2aec886
md"""
#### Symbolics.jl
"""

# ╔═╡ b1b895cc-0b3a-43cc-b15d-85383f646cea
@variables nᵢ ϵᵢ N E α β

# ╔═╡ e06e61b9-8830-4750-8560-89a088ad2b7a
@register_symbolic sum(..)

# ╔═╡ 551f7515-7f3e-4591-a237-35fa16788b1d
g₁ = sum(nᵢ) - N

# ╔═╡ c5bdcab2-7e11-42d8-be11-545482709ba6
g₂ = sum(nᵢ * ϵᵢ) - E

# ╔═╡ 63b3c4b4-72b3-4bc8-ac26-289d59c26bfd
S = N*log(N) - N - (sum(nᵢ * log(nᵢ)-nᵢ))

# ╔═╡ 3b359c98-3a5a-4268-92b5-3f5b8289e758
L = S - (α*g₁) - (β*g₂)

# ╔═╡ 4ae1e515-641a-4c8e-b56a-bca7185a8a29
function ∫(f, x)
	integrate(f, x)
end

# ╔═╡ af7a8f75-4ccd-42c7-9152-43473aad5910
function ∫(f)
	integrate(f)
end

# ╔═╡ 56369cef-eb80-4a0a-91b9-0b68d2363c5c
L_continous = substitute(L, sum => ∫)

# ╔═╡ 18d07a4d-7cd7-4a6e-8c24-57ecc4a4eab4
dLdnᵢ = derivative(L_continous, nᵢ) ~ 0

# ╔═╡ 46a4a70c-305c-4c50-ab93-992b3002ce7e
md"""
## Problem 3.2.4

Show that the volume of a sphere is ``V(r) = \frac{4}{3} \pi r^3`` by integrating ``f = 1`` over a sphere
"""

# ╔═╡ 952f53f4-555d-4074-904a-a0377e3e8b69
md"""
Given spherical coordinates:
```math
\begin{align*}
x &= r \sin \theta \cos \phi \\
y &= r \sin \theta \sin \phi \\
z &= r \cos \theta
\end{align*}
```

```math
\begin{align*}
r &= \sqrt{x^2 + y^2 + z^2} \\
\theta &= \cos^{-1}(z / r) \\ 
\phi &= \tan^{-1}(y / x) \\ 
\end{align*}
```

"""

# ╔═╡ 35d4162c-f6d1-4350-9337-836861fbcbad
md"""
```math
\begin{align*}
V &=\int_r  \int_{\phi}  \int_{\theta}r^2 \sin\theta d\theta d\phi dr \\
V &= \int_0^R r^2 dr \int_0^{2\pi} d\phi \int_0^\pi \sin\theta d\theta \\
V &= (\frac{1}{3} R^3) (2\pi) (2) \\
V &= \boxed{\frac{4}{3} \pi R^3}
\end{align*}
```
"""

# ╔═╡ 94ac2cc2-706b-4359-97df-36fddfd17367
@syms r θ ϕ π z

# ╔═╡ c5047d2f-4786-4dff-80b0-71cc48c6c366
f = r^2 + sin(θ)

# ╔═╡ ef23f392-fa34-4020-b699-efbbfd4e3990
function integrate_multivariate(f)
	args = arguments(f)
    for a in args
		@info a
    end
end

# ╔═╡ 3b4a3c28-08c1-4f22-9fa0-4b315437c94b
integrate_multivariate(f)

# ╔═╡ Cell order:
# ╠═ae0f795d-3237-4fd7-9f3f-a6b143c928f3
# ╠═fea70653-d6c7-438d-8a4f-940b8e7f21c3
# ╠═a129cc0c-73b6-4c18-ac17-d609d0b1bda4
# ╠═4d46902e-7245-4987-9639-658f10005066
# ╠═b37912b9-a9b0-4d69-952f-4d280f4ae382
# ╟─e4ab3255-4fde-458a-a430-7dd6b3e1d521
# ╟─ecb118f8-3a4f-44b9-adc3-421d94cb7514
# ╟─d85c77a5-4b36-46ed-99f0-cb631c34baae
# ╠═09bd6290-b690-46e7-8a17-6c0ed0e5627f
# ╠═722bda22-2044-4e6c-b17b-df1fa9801faf
# ╠═fab2caac-b604-4fad-97c0-ff7e9579d9fc
# ╟─9d52232c-1f7a-40fe-82eb-4f2efa7014d6
# ╠═b20edb50-db8c-45d3-95b6-70e0402cab0a
# ╟─2d99d71f-0293-44e1-bc01-0416f33fa7ef
# ╟─192cf150-6baf-4246-a0de-15ed1a93f8ce
# ╟─889ed804-81bd-4e87-97d5-d5b6f4986a1f
# ╟─50b3d6e1-803a-497e-a8f0-01dbe2aec886
# ╠═b1b895cc-0b3a-43cc-b15d-85383f646cea
# ╠═e06e61b9-8830-4750-8560-89a088ad2b7a
# ╠═551f7515-7f3e-4591-a237-35fa16788b1d
# ╠═c5bdcab2-7e11-42d8-be11-545482709ba6
# ╠═63b3c4b4-72b3-4bc8-ac26-289d59c26bfd
# ╠═3b359c98-3a5a-4268-92b5-3f5b8289e758
# ╠═4ae1e515-641a-4c8e-b56a-bca7185a8a29
# ╠═af7a8f75-4ccd-42c7-9152-43473aad5910
# ╠═56369cef-eb80-4a0a-91b9-0b68d2363c5c
# ╠═18d07a4d-7cd7-4a6e-8c24-57ecc4a4eab4
# ╟─46a4a70c-305c-4c50-ab93-992b3002ce7e
# ╟─952f53f4-555d-4074-904a-a0377e3e8b69
# ╟─35d4162c-f6d1-4350-9337-836861fbcbad
# ╠═94ac2cc2-706b-4359-97df-36fddfd17367
# ╠═c5047d2f-4786-4dff-80b0-71cc48c6c366
# ╠═ef23f392-fa34-4020-b699-efbbfd4e3990
# ╠═3b4a3c28-08c1-4f22-9fa0-4b315437c94b
