### A Pluto.jl notebook ###
# v0.19.26

using Markdown
using InteractiveUtils

# ╔═╡ b18ea998-4db1-4d19-9c83-bfc85301fc1b
using DrWatson

# ╔═╡ 1ffce464-c5e8-4da8-8e5f-ea357b34787f
# ╠═╡ show_logs = false
@quickactivate "graduate-courses"

# ╔═╡ cb753c10-df5f-4a5d-a93d-f49e6e22f3d0
using PlutoUI, Symbolics, SymbolicNumericIntegration

# ╔═╡ 599ed66b-49fc-460b-a1ed-280a827b3b2d
include(srcdir("helpers.jl"));

# ╔═╡ aefdfb86-c5df-44cc-ae56-00a7486ee14b
TableOfContents()

# ╔═╡ 6cbb915e-2299-11ee-0a63-3f883aafe2f2
md"""
# Integral Calculus

**Homework Problems**
- ✅ 2.1.3
- ✅ 2.2.1
- ✅ 2.2.8
- ✅ 2.2.9
"""

# ╔═╡ 0d949853-bf76-45ba-8fb6-3870d78ffc5e
md"""
## Problem 2.1.3

Given the function:
```math
\begin{align*}
F(n) = \int_0^{\infty} x^n e^{-x}dx
\end{align*}
```

Show using integration by parts:
```math
\begin{align*}
F(n) &= nF(n-1) \\
F(n) &= n!
\end{align*}
```

"""

# ╔═╡ 67ceccec-8229-47d2-9ce7-0c3ef67a82a2
@variables x n F(..); Dx = Differential(x);

# ╔═╡ 08146807-be78-4b5e-824b-9339ba74e5dd
f = x^n*exp(-x)

# ╔═╡ 2838d88a-5d04-4762-90a8-08f5412db9bd
F(n) ~ f

# ╔═╡ b34c15a1-502a-4197-aa04-5a7f04ed94fd
u = x^n

# ╔═╡ b2ec0ab8-c0e8-413d-9f44-caa064ecdc53
dv = exp(-x)

# ╔═╡ 847771d7-afdc-48f4-9a5e-1f7dfd7520a2
du = expand_derivatives(Dx(u))

# ╔═╡ fa9ac258-8d20-469b-b8ec-1473a5bd0d95
v, _ = integrate(dv, x)

# ╔═╡ e0c40bdb-01a3-47b1-97b5-d2b5c3fff3df
uv = substitute(u*v, x => 0) - substitute(u*v, x => Inf)

# ╔═╡ f61f1fac-cb1a-4a31-a5f2-d574b904e91e
begin
	vdu = -n * x^(n-1)exp(-x)
	vdu = -n * F(n-1)
end

# ╔═╡ 69599054-68dc-41c9-a179-2adb8da373b2
vdu

# ╔═╡ 977ccea8-4c15-48c6-808a-87c46ac6a178
Fn = uv - vdu

# ╔═╡ 674146c2-fc6c-4623-ac15-260a1ef4b487
begin
	F0, _ = integrate(substitute(f, n => 0), x)
	F0 = substitute(F0, x=>Inf) - substitute(F0, x=>0)
end

# ╔═╡ c57d145e-3e10-4582-a9d5-864f09356cb8
F1 = 1(F0)

# ╔═╡ 9ab18383-5872-4fbd-b69c-8c3acf8c25fe
F2 = 2(F1)

# ╔═╡ d0181d00-ca8f-4388-947f-db2090d6887f
F3 = 3(F2)

# ╔═╡ 6831bfdc-d66c-41ed-8da5-d1cbf915ef51
F4 = 4(F3)

# ╔═╡ afc9193c-77df-4aa4-bef1-8032e6552541
md"""
```math
\therefore F(n) = \boxed{n!}
```
"""

# ╔═╡ 1af99b41-85dc-463d-8a02-35600a5c82e3
md"""
## Problem 2.2.1
Evaluate ``\int_{x_1}^{x_2} \frac{dx}{\sqrt{a^2 - x^2}}`` by switching to ``\theta`` defined by ``x = a \sin \theta``, assume ``0 \leq x \leq \frac{\pi}{2}``
"""

# ╔═╡ e6643a9d-63c3-4684-a518-923c90cd785d
md"""
```math
\begin{align*}
& \int_{x_1}^{x_2} \frac{1}{\sqrt{a^2 - x^2}} dx \\
&\text{Let } x = a\sin\theta, \quad 0 \leq x \leq \frac{\pi}{2} \\
& \int_{\theta_1}^{\theta_2} \frac{a\cos\theta}{\sqrt{a^2\cos^2\theta}} d\theta \\
&= \int_{\theta_1}^{\theta_2} d\theta \\
&= \theta\big|_{\theta_1}^{\theta_2} \\
&= \theta_2 - \theta_1 \\
&= a - 0 = \boxed{a}
\end{align*}
```
"""

# ╔═╡ 60333e2a-c56e-4682-8e70-8bdaca13f6ab
md"""
## Problem 2.2.8
Given:

```math
\begin{align*}
I_1(a) = \int_{0}^{\infty} e^{-ax^2} x dx
\end{align*}
```

Show:
```math
\begin{align*}
I_1(a) = \frac{1}{2a}
\end{align*}
```

"""

# ╔═╡ 334a0ac1-f645-472d-ad83-93ea04c31e2d
md"""
```math
\begin{align*}
I_1(a) &= \int_{0}^{\infty} e^{-ax^2} x dx \\
&\text{Let } u = ax^2 \\
\therefore I_1(a) &= \int_{0}^{\infty} \frac{1}{2a} e^{-u} du \\
&= \left[-\frac{1}{2a}e^{-u}\right]_0^\infty \\
&= \frac{1}{2a}
\end{align*}
```
"""

# ╔═╡ 177eef94-2a63-4afe-8025-50d984b51306
# let
# 	@variables x a
# 	f = exp(-a * x^2) * x
# 	F = integrate(f, x)[1]

# 	# Fₓ₁ = substitute(F, x => Inf)
# 	# Fₓ₀ = substitute(F, x => 0)
# 	# answer = Fₓ₁ - Fₓ₀
# end

# ╔═╡ a8010ef2-53ba-4c59-9d52-b076ca69012a
md"""
## Problem 2.2.9
"""

# ╔═╡ 4e2756b5-9a0c-42d4-aadf-01e5fd5af718
md"""
Given:

```math
I_n(a) = \int_{0}^{\infty} (x^n) e^{-ax^2} dx
```

Evaluate:

```math
I_3(a) \ \text{and} \ I_4(a)
```
"""

# ╔═╡ 4123bd88-5de4-4233-ad93-edb2fc26f766
md"""
```math
\begin{align*}
I_3(a) &= \int_0^\infty x^3 e^{-ax^2} dx \\
&\text{Let } u=ax^2 , \frac{1}{2a}du = xdx \\
&= \frac{1}{2a} \int_0^\infty u^{3/2} e^{-u} du \\
&= \frac{1}{2a} \cdot \frac{1}{2} \Gamma\left(\frac{5}{2}\right) \\
&= \boxed{\frac{3}{4a}}
\end{align*}
```math

```math
\begin{align*}
I_4(a) &= \int_0^\infty x^4 e^{-ax^2} dx \\
&\text{Let } u=ax^2 \frac{1}{2a}du = xdx \\
\therefore \qquad I_4(a) &= \frac{1}{2a} \int_0^\infty u^2 e^{-u} du \\
&= \frac{1}{2a} \cdot 2\Gamma(3) \\
&= \boxed{\frac{2}{a}} \\
\end{align*}
```
"""

# ╔═╡ Cell order:
# ╠═b18ea998-4db1-4d19-9c83-bfc85301fc1b
# ╠═1ffce464-c5e8-4da8-8e5f-ea357b34787f
# ╠═cb753c10-df5f-4a5d-a93d-f49e6e22f3d0
# ╠═599ed66b-49fc-460b-a1ed-280a827b3b2d
# ╠═aefdfb86-c5df-44cc-ae56-00a7486ee14b
# ╟─6cbb915e-2299-11ee-0a63-3f883aafe2f2
# ╟─0d949853-bf76-45ba-8fb6-3870d78ffc5e
# ╠═67ceccec-8229-47d2-9ce7-0c3ef67a82a2
# ╠═08146807-be78-4b5e-824b-9339ba74e5dd
# ╠═2838d88a-5d04-4762-90a8-08f5412db9bd
# ╠═b34c15a1-502a-4197-aa04-5a7f04ed94fd
# ╠═b2ec0ab8-c0e8-413d-9f44-caa064ecdc53
# ╠═847771d7-afdc-48f4-9a5e-1f7dfd7520a2
# ╠═fa9ac258-8d20-469b-b8ec-1473a5bd0d95
# ╠═e0c40bdb-01a3-47b1-97b5-d2b5c3fff3df
# ╠═f61f1fac-cb1a-4a31-a5f2-d574b904e91e
# ╠═69599054-68dc-41c9-a179-2adb8da373b2
# ╠═977ccea8-4c15-48c6-808a-87c46ac6a178
# ╠═674146c2-fc6c-4623-ac15-260a1ef4b487
# ╠═c57d145e-3e10-4582-a9d5-864f09356cb8
# ╠═9ab18383-5872-4fbd-b69c-8c3acf8c25fe
# ╠═d0181d00-ca8f-4388-947f-db2090d6887f
# ╠═6831bfdc-d66c-41ed-8da5-d1cbf915ef51
# ╟─afc9193c-77df-4aa4-bef1-8032e6552541
# ╟─1af99b41-85dc-463d-8a02-35600a5c82e3
# ╟─e6643a9d-63c3-4684-a518-923c90cd785d
# ╟─60333e2a-c56e-4682-8e70-8bdaca13f6ab
# ╟─334a0ac1-f645-472d-ad83-93ea04c31e2d
# ╠═177eef94-2a63-4afe-8025-50d984b51306
# ╟─a8010ef2-53ba-4c59-9d52-b076ca69012a
# ╟─4e2756b5-9a0c-42d4-aadf-01e5fd5af718
# ╟─4123bd88-5de4-4233-ad93-edb2fc26f766
