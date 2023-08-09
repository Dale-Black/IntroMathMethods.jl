### A Pluto.jl notebook ###
# v0.19.26

using Markdown
using InteractiveUtils

# ╔═╡ fe301fc0-d32b-40e7-8ab8-17d81b6bcf8a
# ╠═╡ show_logs = false
using CondaPkg; CondaPkg.add("SymPy")

# ╔═╡ 94e6fcb1-d73f-467b-97fa-82d536f1703b
using PlutoUI, PythonCall, CairoMakie

# ╔═╡ 0b3fe3e4-ae8e-4a58-955f-8c396f2047a8
sp = pyimport("sympy");

# ╔═╡ 908a8f10-6375-427a-97a5-cb2e5caf5cfb
TableOfContents()

# ╔═╡ b61002b8-38ae-40a7-8344-a58ecb0ebffa
md"""
# Introduction

In this notebook, we explore the world of vector calculus, an essential mathematical tool for physics and engineering. We will tackle concepts like vector fields, gradients, curls, and divergences. 

We will be using the Python library SymPy and the Julia packages PlutoUI and PythonCall to work through homework assignments interactively. These assignments cover topics like:

- Time derivatives of vectors
- Line and surface integrals  
- The gradient and its properties
- The curl and its interpretation
- The divergence theorem

Vector calculus unlocks a powerful way of thinking about multivariate functions and vector fields. While the computations can be tricky, the payoff is immense: vector calculus provides the mathematical backbone for subjects like electromagnetism, fluid dynamics, and relativity.

My goal here is to develop an intuitive understanding of these topics. We will lean heavily on visualizations and examples to build up an understanding before diving into the math. With hard work and a little luck, we just might gain some insight into this beautiful subject.
"""

# ╔═╡ 13cdd42b-25c8-4f61-bc4e-1b5588aa9918
md"""
# Vector Calculus
**Homework Assignments**
- 7.1.3
- 7.2.4
- 7.4.3
- 7.5.1
- 7.5.2
- 7.5.3
- 7.5.7
- 7.6.1
- 7.6.6
- 7.6.12
- 7.7.3
"""

# ╔═╡ 5161556a-deae-48c9-bf6c-36a8aaf5882d
md"""
## Review of Vectors Analysis
"""

# ╔═╡ 9a39c2cb-d2e2-424d-a04d-ec52b47cd2dc
md"""
!!! info "Problem 7.1.3"

	Consider three vectors ``\vec{a}``, ``\vec{b}``, ``\vec{c}`` not in the same plane. Show that the box product, also called scalar triple product (``\vec{c} \cdot \vec{a} \times \vec{b}``) gives the volume of the parallelepiped with these vectors as three adjacent edges. At least verify for the case when the ``\vec{c}`` is perpendicular to the ``\vec{a} - \vec{b}`` plane. Show that ``\vec{a} \cdot \vec{b} \times \vec{c} = \vec{b} \cdot \vec{c} \times \vec{a} = \vec{c} \cdot \vec{a} \times \vec{b}`` either geometrically or algebraically. What happens to the box product when two of the vectors are parallel.
"""

# ╔═╡ 53fbf184-231b-488e-aaa9-0b0698be5a9e
md"""
!!! warning "By Hand"

	i)
	```math
	\begin{align*}
	\vec{b} \cdot \vec{c} \times \vec{a} &= \vec{a} \cdot \vec{b} \times \vec{c} \\ \\


	\vec{b} \cdot (\vec{c} \times \vec{a}) &= \begin{vmatrix}
    b_x & b_y & b_z \\
    c_x & c_y & c_z \\
    a_x & a_y & a_z
    \end{vmatrix} \\
    &= b_x \begin{vmatrix}
    c_y & c_z \\
    a_y & a_z
    \end{vmatrix} - b_y \begin{vmatrix}
    c_x & c_z \\
    a_x & a_z
    \end{vmatrix} + b_z \begin{vmatrix}
    c_x & c_y \\
    a_x & a_y
    \end{vmatrix} \\
    &= b_x (c_y a_z - c_z a_y) - b_y (c_x a_z - c_z a_x) + b_z (c_x a_y - c_y a_x) \\
    &= a_x (b_y c_z - b_z c_y) - a_y (b_x c_z - b_z c_x) + a_z (b_x c_y - b_y c_x) \\
    &= \vec{a} \cdot (\vec{b} \times \vec{c})
    \end{align*}
	```

	ii) 

	```math
    \begin{align*}
	\vec{c} \cdot \vec{a} \times \vec{b} &= \vec{b} \cdot \vec{c} \times \vec{a} \\ \\

    \vec{c} \cdot (\vec{a} \times \vec{b}) &= \begin{vmatrix}
    c_x & c_y & c_z \\
    a_x & a_y & a_z \\
    b_x & b_y & b_z
    \end{vmatrix} \\
    &= c_x \begin{vmatrix}
    a_y & a_z \\
    b_y & b_z
    \end{vmatrix} - c_y \begin{vmatrix}
    a_x & a_z \\
    b_x & b_z
    \end{vmatrix} + c_z \begin{vmatrix}
    a_x & a_y \\
    b_x & b_y
    \end{vmatrix} \\
    &= c_x (a_y b_z - a_z b_y) - c_y (a_x b_z - a_z b_x) + c_z (a_x b_y - a_y b_x) \\
    &= a_x (b_y c_z - b_z c_y) - a_y (b_x c_z - b_z c_x) + a_z (b_x c_y - b_y c_x) \\
    &= \vec{a} \cdot (\vec{b} \times \vec{c})
    \end{align*}
    ```

	iii)
    ```math
    \begin{align*}
    \vec{c} \cdot (\vec{a} \times \vec{b}) &= \vec{c} \cdot 0 \\
    &= 0
    \end{align*}
    ```
"""

# ╔═╡ 5d7e3199-dc42-4870-904c-1906a99cfa97
md"""
## Time Derivatives of Vectors
"""

# ╔═╡ bef5c3ec-8f59-4f4b-8b6a-ec22f5c9b8eb
md"""
!!! info "Problem 7.2.4"

	A particle has a position vector ``\vec{r} = \cos(\omega t) \hat{i} + \sin(\omega t) \hat{j}``. (1) Describe its motion in cartesian coordinates with a sketch and in words. (2) Compute its velocity and accelaration. (3) Find the magnitudes of both. (4) Switch to polar coordinates and find the ``r`` and ``\theta`` for this problem. (5) What is ``\vec{a}``?
"""

# ╔═╡ 4fadadf2-5b23-4038-92f2-8693918ca740
md"""
!!! warning "By Hand"
	i)
	```math
	\begin{align*}
	\vec{r} &= \cos(\omega t) \hat{i} + \sin(\omega t) \hat{j}
	\end{align*}
	```

	ii)
	```math
	\begin{align*}
	\vec{v} &= -\omega \sin(\omega t) \hat{i} + \omega \cos(\omega t) \hat{j} \\
	\vec{a} &= -\omega^2 \cos(\omega t) \hat{i} - \omega^2 \sin(\omega t) \hat{j}
	\end{align*}
	```

	iii)
	```math
	\begin{align*}
	|\vec{v}| &= |\omega| \\
	|\vec{a}| &= \omega^2
	\end{align*}
	```

	iv)
	```math
	\begin{align*}
	r &= 1 \\
	\theta &= \omega t
	\end{align*}
	```

	v)
	```math
	\begin{align*}
	\vec{a} &= -\omega^2 \cos(\omega t) \hat{i} - \omega^2 \sin(\omega t) \hat{j}
	\end{align*}
	```
"""

# ╔═╡ b839496c-bc56-4504-ab53-ae2f257bb4ba
md"""
#### With SymPy
"""

# ╔═╡ 4cacb360-625c-40e7-bdf9-437e700a1657
md"""
#### (i)
This is a circle with radius 1 and frequency ``\omega`` centered at the origin
"""

# ╔═╡ e1f7df8a-57c7-499e-b605-4be6604fae30
t0 = range(0, 2π, length = 1000)

# ╔═╡ 09a2bcab-1237-469e-ba0f-23cd9aad28d3
ω0 = 1

# ╔═╡ 72576547-e113-44aa-9b17-34bb5a575277
x0, y0  = cos.(ω0*t0), sin.(ω0*t0)

# ╔═╡ a3f36b49-526a-47cf-a55c-2648d4a367f1
let
	f = Figure()
	ax = CairoMakie.Axis(
		f[1, 1]
	)
	scatterlines!(x0, y0, markersize = 1)
	f
end

# ╔═╡ 68487535-e065-4778-94e7-db8344cc6237
md"""
#### (ii)
"""

# ╔═╡ 4bf65708-3747-4eba-b523-641223ff4d99
ω, t = sp.symbols("ω, t")

# ╔═╡ 47220365-3466-4fb3-b86e-f0428e0d9fa9
x₀, y₀ = sp.cos(ω * t), sp.sin(ω * t)

# ╔═╡ 9f3f770a-3ee5-4aec-9953-1b73641074c5
vᵢ, vⱼ = sp.diff(x₀, t), sp.diff(y₀, t)

# ╔═╡ 544eef1f-d204-425b-85ff-e25884721f5f
aᵢ, aⱼ = sp.diff(vᵢ, t), sp.diff(vⱼ, t)

# ╔═╡ b6c7ac3d-f2aa-4c4f-b53d-75a8d5886839
md"""
#### (iii)
"""

# ╔═╡ 24cd17dc-678b-4c7d-ab4a-bd5c0eeb8e76
v_magnitude = sp.sqrt(vᵢ^2 + vⱼ^2)

# ╔═╡ 0d5dceac-4230-4a2c-8783-1563e1738df4
sp.trigsimp(v_magnitude)

# ╔═╡ f8daa2ec-5e87-48a8-92ac-609996d302b4
a_magnitude = sp.sqrt(aᵢ^2 + aⱼ^2)

# ╔═╡ 44ffdc2f-3ce5-4883-b1d5-555dc8a0dddd
sp.trigsimp(a_magnitude)

# ╔═╡ 5487b403-64cb-46f5-b007-ba563f9fc218
md"""
#### (iv)
"""

# ╔═╡ 7afa8504-6f1e-4b5e-9d86-67227bf1bd7d
r = sp.sqrt(x₀^2 + y₀^2)

# ╔═╡ 5ef84561-4daa-4b5d-ba97-d26b521554f3
θ = sp.trigsimp(sp.atan(y₀ / x₀))

# ╔═╡ 75d86070-6ab0-488d-b53e-50f69cda76dc
md"""
#### (v)
"""

# ╔═╡ fa4d3c70-f60d-4639-a176-9a61897c7088
vᵣ, v₀ = sp.diff(r, t), sp.diff(θ, t)

# ╔═╡ 419354e8-4975-43dc-a119-d3fadc02004a
aᵣ, a₀ = sp.diff(vᵣ, t), sp.diff(v₀, t)

# ╔═╡ 08bbdbdf-3778-415d-813b-f81922f64215
md"""
## Scalar and Vector Fields
"""

# ╔═╡ 3b55ae9d-855f-4dca-a293-ad235ed7c8ff
md"""
## Line and Surface Integrals
"""

# ╔═╡ acd8c7a1-7ba1-4819-a91b-711cd9a3ac3b
md"""
!!! info "Problem 7.4.3"
	(i) Calculate the line integral of ``\vec{F} = 2xy^2\hat{i} + x^2\hat{j}`` between ``(0, 0)`` and ``(1, 1)``, but along ``y = x^2``. (ii) What is the line integral around a unit circle centered at the origin. First make a sensible choice for the parameter ``t`` that will move the coordinates along the circle. Note that even for a non-conservative field, the line integral along closed paths can vanish
"""

# ╔═╡ f159263e-72a9-4417-9cf3-cfd9e18214ac
md"""
!!! warning "By Hand"
	i)
	```math
	\begin{align*}
	y &= x^2 \\
	f_{i} &= 2 \cdot x \cdot y^2 \\
	f_{j} &= x^2 \\
	\frac{\partial y}{\partial x} &= 2x \\
	dy &= \frac{\partial y}{\partial x} \cdot dx \\
	F_{i} &= \int_{0}^{1} f_{i} dx \\
	F_{j} &= \int_{0}^{1} f_{j} dy \\
	F &= f_{i} + f_{j} = \frac{1}{3} + \frac{1}{2} = \frac{5}{6}
	\end{align*}
	```

	ii)
	```math
	\begin{align*}
	x(t) &= \cos(t) \\
	y(t) &= \sin(t) \\
	\frac{dx}{dt} &= -\sin(t) \\
	\frac{dy}{dt} &= \cos(t) \\
	f_{it} &= 2 x(t) y(t)^2 \\
	f_{jt} &= x(t)^2 \\
	F_{it} &= \int_{0}^{2\pi} f_{it} \cdot \frac{dx}{dt} \,dt \\
	F_{jt} &= \int_{0}^{2\pi} f_{jt} \cdot \frac{dy}{dt} \,dt \\
	F &= F_{it} + F_{jt} = 0
	\end{align*}
	```
"""

# ╔═╡ 8a9fa30c-49f8-4ae0-81a7-93ce3b3792fc
md"""
#### With SymPy
"""

# ╔═╡ 25f81dd0-27e2-4613-be93-fb8a9b4d0c57
md"""
**(i)**
"""

# ╔═╡ c139e863-0723-4cec-8131-3eda35a50e7a
x, dx = sp.symbols("x, dx")

# ╔═╡ 9e603a65-87d1-48ba-8052-37c9cd2b3a43
y1 = x^2

# ╔═╡ 14dfbbe8-5d0a-4572-8d22-3b73bbd22f93
fᵢ = 2 * x * y1^2

# ╔═╡ 9718a6c4-73ef-4fb7-a204-dd803fa5dc52
fⱼ = x^2

# ╔═╡ 24c4be9c-7c05-4925-a1c9-5a9aaaf6d878
∂y∂x = sp.diff(y1, x)

# ╔═╡ 117b734a-0d5a-4a6f-b239-d6067f399ef5
dy = ∂y∂x

# ╔═╡ 394f87f0-1046-4835-a9a6-3682fb8a1f03
Fᵢ = sp.integrate(fᵢ, (x, 0, 1))

# ╔═╡ 25877abf-5289-465e-8fab-7e147f1ca5b3
Fⱼ = sp.integrate(fⱼ * dy.subs(dx, 1), (x, 0, 1))

# ╔═╡ c7e9ba5f-5ce0-423c-86cd-015efca7137f
F = Fᵢ + Fⱼ

# ╔═╡ bc8c413f-8184-4d6d-a87b-a0ef3969c93a
md"""
**(ii)**
"""

# ╔═╡ 293e81fe-1744-40f6-bd49-c7476d281eb2
x2, y2 = sp.cos(t), sp.sin(t)

# ╔═╡ 7e47e8d2-65f3-4627-815c-d94212bac67b
dx2dt, dy2dt = sp.diff(x2, t), sp.diff(y2, t)

# ╔═╡ 451ae9c6-5ef6-4651-99bb-44132e844a89
fᵢₜ = fᵢ.subs(x, x2)

# ╔═╡ 5da3cbc8-6e1f-436a-a3b9-58c2992e2152
fⱼₜ = fⱼ.subs(x, x2)

# ╔═╡ da87a788-7538-43bd-abc6-510acb8999ff
Fᵢₜ = sp.integrate(fᵢₜ * dx2dt, (t, 0, 2 * sp.pi))

# ╔═╡ 61120a25-bd56-4d22-8431-dc0734117db7
Fⱼₜ = sp.integrate(fⱼₜ * dy2dt, (t, 0, 2 * sp.pi))

# ╔═╡ 1a722ba9-266f-4fc4-a521-c5038550b45e
F2 = Fᵢₜ + Fⱼₜ

# ╔═╡ 873551d7-aebb-41fc-8746-00c045d3a10f
md"""
## Scalar Field and Gradient
"""

# ╔═╡ 3ece13e7-18e9-4449-88d1-1d6e39b001ec
md"""
!!! info "Problem 7.5.1"
	Consider a sphere of radius ``R`` centered at the origin. Write a formula for ``h(x, y)`` the height of the hemisphere above the point ``x, y``. Calculate its gradient and compare the results against your expectations.
"""

# ╔═╡ 3ece9a12-8071-409e-a313-c68135e20503
md"""
!!! warning "By Hand"
	Formula for height of the hemisphere above the point ``(x, y)``:

	```math
	h(x, y) = \sqrt{R^2 - x^2 - y^2}
	```

	Gradient of ``h(x, y)``:

	```math
	\nabla h(x, y) = \left( -\frac{x}{\sqrt{R^2 - x^2 - y^2}}, -\frac{y}{\sqrt{R^2 - x^2 - y^2}} \right)
	```
"""


# ╔═╡ c7fe4490-fd11-4c51-89ac-770d3adaadc5
y, R = sp.symbols("y , R")

# ╔═╡ 3446cb4b-3dad-4659-8e53-2ffab246e151
h = sp.sqrt(R^2 - x^2 - y^2)

# ╔═╡ fc2f149c-3aa2-4f8d-b1f3-59a69d073cf8
∇h = [sp.diff(h, i) for i in (x, y)]

# ╔═╡ 946cb3e1-ba80-4b91-8afb-a9939ec6ef0b
md"""
!!! info "Problem 7.5.2"
	Show that
	```math
	\begin{align*}
	\vec{\nabla} \phi \chi = \phi \vec{\nabla \chi} + \chi \vec{\nabla} \phi
	\end{align*}
	```
"""

# ╔═╡ c13032f9-d9ac-4856-a9eb-537c8c102827
md"""
!!! warning "By Hand"
	```math
	\begin{align*}
	\vec{\nabla} (\phi \chi) &= \vec{\nabla} (\phi \chi_i + \phi \chi_j + \phi \chi_k) \\
	&= \vec{\nabla} (\phi \chi_i) + \vec{\nabla} (\phi \chi_j) + \vec{\nabla} (\phi \chi_k) \\
	&= \phi \vec{\nabla} \chi_i + \chi_i \vec{\nabla} \phi + \phi \vec{\nabla} \chi_j + \chi_j \vec{\nabla} \phi + \phi \vec{\nabla} \chi_k + \chi_k \vec{\nabla} \phi \\
	&= \phi \vec{\nabla} \chi + \chi \vec{\nabla} \phi
	\end{align*}
	```
"""

# ╔═╡ 77bfe4d8-c6ae-4832-8f00-a9e93bd1f540
md"""
#### With SymPy
"""

# ╔═╡ f1d8eaa1-9741-4988-8135-ac439c5f3005
z = sp.symbols("z")

# ╔═╡ 09c91686-5fab-4dfb-a906-3209abbbfa7d
χ = sp.Function("χ")(x, y, z)

# ╔═╡ 5a972777-a650-445f-b930-a05818559415
ϕ = sp.Function("ϕ")(x, y, z)

# ╔═╡ 3778e89a-b996-4fc2-87e1-b7a7af30cb56
∇ϕ = [sp.diff(ϕ, i) for i in (x, y, z)]

# ╔═╡ 2de36b1d-175f-496b-b5aa-6bc8dce4c064
∇χ = [sp.diff(χ, i) for i in (x, y, z)]

# ╔═╡ 37a2ffe4-49a3-46d6-87e5-dac31b9e82a8
∇ϕχ = [ϕ * ∇χ[i] + χ * ∇ϕ[i] for i in 1:3]

# ╔═╡ e5f68df7-144c-43b2-bef7-0df8905997ca
∇_product = [sp.diff(ϕ * χ, i) for i in (x, y, z)]

# ╔═╡ 3f5cabb9-546e-4cf3-8a79-2563d9d3ab61
∇ϕχ .== ∇_product

# ╔═╡ 27984f40-0d25-4216-95cb-8f15e1badd97
md"""
!!! info "Problem 7.5.3"
	Recall the equations for cylindrical coordinates ``(ρ, φ, z)`` and spherical coordinates ``(r, θ, φ)`` in terms of cartesian coordinates. 

	Verify that these are orthogonal coordinates, which means, for example, that the direction of purely increasing ``r`` (at fixed ``θ`` and ``φ``) is perpendicular to the direction of increasing ``θ`` and ``φ`` with the other two coordinates similarly held fixed. 

	Argue that the direction in which just one coordinate increases is perpendicular to the gradients of the other two coordinates (with all three spherical or cylindrical coordinates expressed as a function of ``x, y, z``). 

	Compute the three gradients for each coordinate system and verify orthogonality.
"""

# ╔═╡ d898d7fd-b6f8-4ec4-bcd4-6fd6611a8e70
md"""
!!! warning "By Hand"
	The cylindrical coordinates ``(ρ, φ, z)`` are defined as follows:

    ```math
    \begin{align*}
    \rho &= \sqrt{x^2 + y^2} \\
    \phi &= \tan^{-1}\left(\frac{y}{x}\right) \\
    z &= z
    \end{align*}
    ```

	We compute the gradient of each of these coordinates:

    ```math
    \begin{align*}
    \nabla \rho &= \left(\frac{\partial \rho}{\partial x}, \frac{\partial \rho}{\partial y}, \frac{\partial \rho}{\partial z}\right) \\
    \nabla \phi &= \left(\frac{\partial \phi}{\partial x}, \frac{\partial \phi}{\partial y}, \frac{\partial \phi}{\partial z}\right) \\
    \nabla z &= \left(\frac{\partial z}{\partial x}, \frac{\partial z}{\partial y}, \frac{\partial z}{\partial z}\right)
    \end{align*}
    ```

	We then verify orthogonality by computing the dot product of each pair of gradients and confirming that the result is zero:

    ```math
    \begin{align*}
    \nabla \rho \cdot \nabla \phi &= 0 \\
    \nabla \phi \cdot \nabla z &= 0 \\
    \nabla \rho \cdot \nabla z &= 0
    \end{align*}
    ```

	The spherical coordinates ``(r, θ, φ)`` are defined as follows:

    ```math
    \begin{align*}
    r &= \sqrt{x^2 + y^2 + z^2} \\
    \theta &= \cos^{-1}\left(\frac{z}{\sqrt{x^2 + y^2 + z^2}}\right) \\
    \phi &= \tan^{-1}\left(\frac{y}{x}\right)
    \end{align*}
    ```

	We compute the gradient of each of these coordinates:

    ```math
    \begin{align*}
    \nabla r &= \left(\frac{\partial r}{\partial x}, \frac{\partial r}{\partial y}, \frac{\partial r}{\partial z}\right) \\
    \nabla \theta &= \left(\frac{\partial \theta}{\partial x}, \frac{\partial \theta}{\partial y}, \frac{\partial \theta}{\partial z}\right) \\
    \nabla \phi &= \left(\frac{\partial \phi}{\partial x}, \frac{\partial \phi}{\partial y}, \frac{\partial \phi}{\partial z}\right)
    \end{align*}
    ```

	We then verify orthogonality by computing the dot product of each pair of gradients and confirming that the result is zero:

    ```math
    \begin{align*}
    \nabla r \cdot \nabla \theta &= 0 \\
    \nabla r \cdot \nabla \phi &= 0 \\
    \nabla \theta \cdot \nabla \phi &= 0
    \end{align*}
    ```
"""

# ╔═╡ 47894665-ac20-48f4-aa43-f6b27e854042
md"""
#### With SymPy
"""

# ╔═╡ 37f965aa-717d-43ec-adf5-8bc2790979b8
md"""
**Cylindrical**
"""

# ╔═╡ 185b0b11-7862-49c5-a28b-6b2d8c97e0b0
ρ_cyl = sp.sqrt(x^2 + y^2)

# ╔═╡ 0ebf5c1a-6287-4298-94f2-50361604bb0a
ϕ_cyl = sp.atan(y / x)

# ╔═╡ a45b5225-1c71-48a0-8cc4-0886dfe85cc3
z_cyl = z

# ╔═╡ 34bc2428-06d4-404c-b0b8-ac5d978119a5
∇ρ_cyl = [sp.diff(ρ_cyl, coord) for coord in (x, y, z)]

# ╔═╡ 16b9644c-aa6a-4665-a96d-ae4aa86d0b88
∇ϕ_cyl = [sp.diff(ϕ_cyl, coord) for coord in (x, y, z)]

# ╔═╡ 541fb80a-63a6-4b66-8a03-69ae1e0f4553
∇z_cyl = [sp.diff(z_cyl, coord) for coord in (x, y, z)]

# ╔═╡ cfd88052-e979-435b-bb41-2dc26552aa4a
function dot_product(v1, v2)
    return sum([v1[i] * v2[i] for i in 1:3])
end

# ╔═╡ 160e03a0-545e-4d86-bd8f-5fe0180580c5
dot_product(∇ρ_cyl, ∇ϕ_cyl)

# ╔═╡ 5bc43ec7-9720-411f-900e-cd5f3af68824
dot_product(∇ϕ_cyl, ∇z_cyl)

# ╔═╡ a3947c1e-3a5f-42fa-9ea3-7f0da19261c2
dot_product(∇ρ_cyl, ∇z_cyl)

# ╔═╡ 0c950c0c-a209-421f-8d25-b691a482334f
md"""
**Spherical**
"""

# ╔═╡ de382bf0-783d-487b-b3c0-83441bf7662a
r_sph = sp.sqrt(x^2 + y^2 + z^2)

# ╔═╡ c61be673-d114-47fc-9bfa-0457cf16aacf
θ_sph = sp.acos(z / sp.sqrt(x^2 + y^2 + z^2))

# ╔═╡ 8e2f5d10-a031-431a-a5a5-a475aa46533b
ϕ_sph = sp.atan(y / x)

# ╔═╡ d84d7137-cb72-426a-8866-a5156e5334de
∇r_sph = [sp.diff(r_sph, coord) for coord in (x, y, z)]

# ╔═╡ 2bc5ce49-f52a-48d6-be64-a040052927af
∇θ_sph = [sp.diff(θ_sph, coord) for coord in (x, y, z)]

# ╔═╡ 2e34d2ed-14c2-4746-8b3e-58510a5edf86
∇ϕ_sph = [sp.diff(ϕ_sph, coord) for coord in (x, y, z)]

# ╔═╡ 13deb593-c204-4b77-808f-62a7e333418e
sp.simplify(dot_product(∇r_sph, ∇θ_sph))

# ╔═╡ 7e81d4a3-2e75-4321-879a-df8259a93c1e
dot_product(∇r_sph, ∇ϕ_sph)

# ╔═╡ 5078c8e3-50ae-4312-ad4d-0e6416d8b360
dot_product(∇θ_sph, ∇ϕ_sph)

# ╔═╡ 32a3470d-2629-4e74-94fe-5383cb4af5b9
md"""
!!! info "Problem 7.5.7"
	You are on a hot volcanic mountain where the temperature is given by ``T(x, y) = x^2 + xy^3``. (i) If you are located at ``(x = 1, y = 1)``, in which direction will you run to beat the heat? (ii) If your steps are ``1/10`` units long, by how much will the temperature drop after the first step? (Work to first order.)
"""

# ╔═╡ bfb3a382-f362-4b20-afb4-ae82af191a57
md"""
!!! warning "By Hand"
	```math
	\begin{align*}
	T &= x^2 + xy^3 \\
	\nabla T &= \left(\frac{\partial T}{\partial x}, \frac{\partial T}{\partial y}\right) = (2x + y^3, 3xy^2) \\
	\nabla T|_{(1,1)} &= (2 + 1, 3) = (3, 3) \\
	-\nabla T|_{(1,1)} &= (-3, -3) \\
	\hat{d} &= \frac{-\nabla T|_{(1,1)}}{||\nabla T|_{(1,1)}||} = \left(-\frac{\sqrt{2}}{2}, -\frac{\sqrt{2}}{2}\right) \\
	\Delta T &= \text{step size} \times ||\nabla T|_{(1,1)}|| = \frac{1}{10} \times \sqrt{18} = 0.6
	\end{align*}
	```
"""

# ╔═╡ c00c7da9-b548-46c1-b53a-970539cc9780
md"""
#### With SymPy
"""

# ╔═╡ 088bd06b-5d6a-41c1-8a58-a4edbca9879e
md"""
**(i)**
"""

# ╔═╡ 933d3db0-470e-4e6b-8175-94290774af0c
T = x^2 + x*y^3

# ╔═╡ 710c5f83-871c-4eda-b10a-4f798c6b6fac
∇T = sp.Matrix([sp.diff(T, i) for i in (x, y)])

# ╔═╡ 94df1dd9-ee4c-4d23-bbb6-b5331f4196b4
∇T₀ = ∇T.subs([(x, 1), (y, 1)]) * -1

# ╔═╡ d85ad2a0-4cbc-45fb-a5dd-fbf0980181c0
direction = ∇T₀ / ∇T₀.norm()

# ╔═╡ 81e6ddce-0e8a-4d57-8e9c-f827dd613b9b
md"""
**(ii)**
"""

# ╔═╡ b7151a26-b9c9-4329-a317-2131c35f3dbc
step_size = 1/10

# ╔═╡ fc6e6bbd-8744-4369-bd7c-29e8f52c8048
ΔT = step_size * ∇T₀.norm()

# ╔═╡ 25112cba-34af-4ad7-9af2-27bf6b560374
md"""
## Curl of a Vector Field
"""

# ╔═╡ baa68d3f-8cd9-4f28-a8c6-a957ec5b3fb6
md"""
!!! info "Problem 7.6.1"
	Consider the complex integral ``\oint f(z) dz`` where ``f = u + iv`` and ``dz = dx + i dy`` over a contour that lies on the domain of analyticity of ``f``. Write the real and imaginary parts as circulations of two real vector fields and verify that the CRE ensure that both fields are conservative.
"""

# ╔═╡ 9efddcb3-6c01-494b-9512-e9d098868aae
md"""
!!! warning "By Hand"
	```math
	\begin{align*}
	&\oint f(z) dz = \oint (u + iv)(dx + idy) \\
	&\quad = \oint (u dx - v dy) + i \oint (v dx + u dy) \\
	&\quad = \oint F \cdot dr + i \oint G \cdot dr \\
	&\quad = \oint (u dx - v dy) + i \oint (v dx + u dy) \\
	&\text{where } F = (u, -v) \text{ and } G = (v, u) \\
	&\text{The curls of } F \text{ and } G \text{ are: } \\
	&\nabla \times \mathbf{F} = \frac{\partial (-v)}{\partial x} - \frac{\partial u}{\partial y} \\
	&\nabla \times \mathbf{G} = \frac{\partial u}{\partial x} - \frac{\partial v}{\partial y} \\
	&\text{Substituting the CRE into the expressions for the curls, we find: } \\
	&\nabla \times \mathbf{F} = - \frac{\partial u}{\partial y} - \frac{\partial u}{\partial y} = 0 \\
	&\nabla \times \mathbf{G} = \frac{\partial u}{\partial x} - \frac{\partial u}{\partial x} = 0 \\
	&\text{Therefore, both vector fields are conservative.}
	\end{align*}
	```
"""

# ╔═╡ a81a50ff-5d1a-4e77-9567-46ff2be1be80
md"""
#### With SymPy
"""

# ╔═╡ 528d1bb7-d56b-4f63-84a9-0a8035878c95
u, v = sp.Function("u")(x, y), sp.Function("v")(x, y)

# ╔═╡ f4010731-0793-4342-a65f-21dfe418d608
f = u + sp.I*v

# ╔═╡ cefd5ebb-041d-416d-8e05-390829e215f8
f_x, f_y = sp.diff(f, x), sp.diff(f, y)

# ╔═╡ a3f4fa24-314a-42d3-b932-ce17cbb766c8
F1, G1 = sp.Matrix([u, -v]), sp.Matrix([v, u])

# ╔═╡ 796ddfba-ef94-450c-9ce1-a87eeb9670a9
curl_F = sp.diff(F1[1], x) - sp.diff(F1[0], y)

# ╔═╡ 75b74bde-ab81-4414-af83-f1868b59657a
curl_G = sp.diff(G1[1], x) - sp.diff(G1[0], y)

# ╔═╡ a1d4a670-ccb6-4bc5-8d35-4011d6003812
curl_F_CRE = curl_F.subs([(sp.diff(u, x), sp.diff(v, y)), (sp.diff(u, y), sp.diff(-v, x))])

# ╔═╡ c6b8eb02-d3f8-419f-b3dd-0c67e810254a
curl_G_CRE = curl_G.subs([(sp.diff(u, x), sp.diff(v, y)), (sp.diff(u, y), sp.diff(-v, x))])

# ╔═╡ 512817e6-114d-4b80-b7a0-bc77c5c7c1ab
md"""
!!! info "Problem 7.6.6"
	If ``u + iv`` is an analytic function, show that the curves of constant ``u`` intersect the curves of constant ``v`` orthogonally at each point. (Bring in the gradient to solve this). For the case of ``f(z) = z^2``, sketch the curves of constant ``u`` and ``v`` and check this claim
"""

# ╔═╡ 2e283fd3-7cec-4926-96fa-b719eb288328
md"""
!!! warning "By Hand"
	```math
	\begin{align*}
	f(z) &= u + iv \\
	F &= \nabla u = \left(\frac{\partial u}{\partial x}, \frac{\partial u}{\partial y}\right) \\
	G &= \nabla v = \left(\frac{\partial v}{\partial x}, \frac{\partial v}{\partial y}\right). \\
	F \cdot G &= \frac{\partial u}{\partial x} \frac{\partial v}{\partial x} + \frac{\partial u}{\partial y} \frac{\partial v}{\partial y}. \\
	\frac{\partial u}{\partial x} &= \frac{\partial v}{\partial y} \text{ and } \frac{\partial u}{\partial y} = -\frac{\partial v}{\partial x}. \\
	F &\cdot G = \frac{\partial v}{\partial y} \frac{\partial v}{\partial x} - \frac{\partial v}{\partial x} \frac{\partial v}{\partial y} = 0. \\
	\end{align*}
	```
"""

# ╔═╡ 760e58cc-d859-4ce5-8d6f-1f105ded77fe
f2 = (u + sp.I*v)^2

# ╔═╡ 76087030-4287-41a8-ae10-805188dc7f65
begin
	u_x = sp.diff(u, x)
	u_y = sp.diff(u, y)
	v_x = sp.diff(v, x)
	v_y = sp.diff(v, y)
end

# ╔═╡ d5d3381c-bbcf-41e6-9cbd-b821126d4570
begin
	F3 = sp.Matrix([u_x, u_y])
	G3 = sp.Matrix([v_x, v_y])
end

# ╔═╡ b3ae3846-26a3-4da0-9ca1-23400a22e2fe
F_dot_G = F3.dot(G3)

# ╔═╡ effd8948-af6c-4040-ba06-3e31a784df8d
orthogonality_CRE = F_dot_G.subs([(u_x, v_y), (u_y, -v_x)])

# ╔═╡ b51c8800-e410-4481-b853-0763caf20baa
md"""
!!! info "Problem 7.6.12"
	Show that the area enclosed by a counter-clockwise curve ``C`` in the plane is given by ``\frac{1}{2} \oint_C (xdy - ydx)``
"""

# ╔═╡ a47e7c7a-bb5b-462a-87b3-69a6ddbe18b7
md"""
!!! warning "By Hand"
	```math
	\begin{align*}
	\oint_C \vec{F} \cdot \vec{dr} &= \oint P dx + Q dy \\
	&= \int \int (\frac{\partial Q}{\partial x} - \frac{\partial P}{\partial y}) dA \quad (\text{where } A = \frac{1}{2} \oint_C (xdy - ydx)) \\
	P &= -y \\ 
	Q &= x \\ \\
	\int \int (\frac{\partial Q}{\partial x} - \frac{\partial P}{\partial y}) dA &= \int \int (1 - (-1)) dA \\
	&= 2 \int \int dA \\
	&= 2A
	\end{align*}
	```
"""

# ╔═╡ 54468bc7-3115-4582-af0c-4d66f1cb0d3e
md"""
## The Divergence of a Vector Field
"""

# ╔═╡ b1633978-2110-4bef-a704-2284b34bcadc
md"""
!!! info "Problem 7.7.3"
	Find the surface integral of ``\vec{W} = x^3y \hat{i} + y^2 x \hat{j} + z \hat{k}`` over a unit cube centered at the origin with edges parallel to the axes by direct computation and by Gauss' theorem. Use symmetries to save time
"""

# ╔═╡ dae97f2e-a2a5-4265-8af6-d61c2bbd78ac
md"""
!!! warning "By Hand"
    ```math
    \begin{align*}
    \vec{W} &= x^3y \hat{i} + y^2 x \hat{j} + z \hat{k} \\
    \Phi_{\text{top}} &= \int_{-1/2}^{1/2} \int_{-1/2}^{1/2} \frac{1}{2} \, dx \, dy = \left[\frac{1}{2} x\right]_{-1/2}^{1/2} \bigg|_{-1/2}^{1/2} = \frac{1}{2} - \left(-\frac{1}{2}\right) = \frac{1}{2} \\
    \Phi_{\text{bottom}} &= \int_{-1/2}^{1/2} \int_{-1/2}^{1/2} -\left(-\frac{1}{2}\right) \, dx \, dy = \left[\frac{1}{2} x\right]_{-1/2}^{1/2} \bigg|_{-1/2}^{1/2} = \frac{1}{2} - \left(-\frac{1}{2}\right) = \frac{1}{2} \\
    \Phi_{\text{total}} &= \Phi_{\text{top}} + \Phi_{\text{bottom}} = \frac{1}{2} + \frac{1}{2} = 1 \\
    \nabla \cdot \vec{W} &= 3x^2y + 2y^2x + 1 \\
    \int_V (\nabla \cdot \vec{W}) \, dV &= \int_{-1/2}^{1/2} \int_{-1/2}^{1/2} \int_{-1/2}^{1/2} (3x^2y + 2y^2x + 1) \, dx \, dy \, dz = 1 
    \end{align*}
    ```
"""

# ╔═╡ 7b80df03-eeb4-4130-9b56-0a765a7d232b
md"""
#### With SymPy
"""

# ╔═╡ 37203e2b-8c2c-4fa0-997a-a8f59da022cc
W = sp.Matrix([x^3 * y, y^2 * x, z])

# ╔═╡ c967a4b7-7c4d-41c0-b0f7-2b326c86ebde
limits = ((x, -1/2, 1/2), (y, -1/2, 1/2), (z, -1/2, 1/2))

# ╔═╡ 49106849-ba23-45b4-90c4-54a13aa95b8d
div_W = W[0].diff(x) + W[1].diff(y) + W[2].diff(z)

# ╔═╡ fb10a998-87cb-476e-96d9-a232fa6fa5b8
volume_integral = sp.integrate(div_W, limits...)

# ╔═╡ c607bfcc-fa55-46eb-bc5f-b47fb5f2f34b
md"""
# Summary

As always, the book provides a great summary

---
Know the dot product, its definition in terms of componenets in an arthonormal basis and its properties such as linearity. Same for cross product

---
If ``\vec{r}(t)`` is a function of time then know

```math
\begin{align*}
\vec{r} &= \vec{i}x(t) + \vec{j}y(t) = \vec{e_r}r \\
\frac{d \vec{r}}{dt} &= \vec{i} \dot{x}(t) + \vec{j} \dot{y}(t) = \vec{e_r}\dot{r} + r \omega\vec{e_{\theta}}
\end{align*}
```

---

The line integral of a vector field ``\vec{F}(\vec{r})`` is defined as
```math
\begin{align*}
\int_1^2 \vec{F}(\vec{r}) \cdot \vec{dr} = \lim_{n \to \infty} \sum_{i=1}^{n} \vec{F}(\vec{r_i}) \cdot \vec{dr_i}
\end{align*}
```

- If the answer is path independent, the field is conservative

```math
\begin{align*}
\oint \vec{F}(\vec{r}) \cdot \vec{dr} = 0
\end{align*}
```

---

The surface integral of a vector field ``\vec{V}(\vec{r})`` over a surface ``S`` is defined as

```math
\begin{align*}
\int_S \vec{V} \cdot \vec{dS} = \lim_{i \to 100} \sum_i \vec{V}(\vec{r_i}) \cdot \vec{dS_i}
\end{align*}
```

---

The gradient enters as follows

```math
\begin{align*}
d \phi &= \frac{\partial \phi}{\partial x} dx + \frac{\partial \phi}{\partial y} dy + \frac{\partial \phi}{\partial z} dz \equiv \vec{\nabla \phi} \cdot \vec{dr} \text{ where} \\

\vec{\nabla \phi} &= \frac{\partial \phi}{\partial x} \hat{i} + \frac{\partial \phi}{\partial y} \hat{j} + \frac{\partial \phi}{\partial z} \hat{k} \text{ and} \\

\vec{dr} &= \hat{i} dx + \hat{j} dy + \hat{k} dz

\end{align*}
```

- since ``d \phi = \vec{\nabla } \cdot \vec{dr} = |\vec{\nabla }| |\vec{dr}| \cos \theta``, the gradient gives the direction of greatest rate of change and equals in magnitude that rate of change

---

The integral of the gradient is path independent

```math
\begin{align*}
\int_1^2 \vec{\nabla \phi} \cdot \vec{dr} = \phi(2) - \phi(1) = |\phi|_{\partial P}
\end{align*}
```

- where ``1`` and ``2`` are short for ``\vec{r_1}`` and ``\vec{r_2}`` respectively

---

In non cartesian coordinates

```math
\begin{align*}
\vec{\nabla \phi} = \sum_i \vec{e_i} \frac{1}{h_i} \frac{\partial \phi}{\partial u_i}
\end{align*}
```

---

Green's theorem says if ``\vec{W}`` and the loop ``C`` lie in a plane (say the x-y plane)

```math
\begin{align*}
\oint_{C = \partial S} \vec{W} \cdot \vec{dr} = \int \int (\frac{\partial W_y}{\partial x} - \frac{\partial W_x}{\partial y}) dx dy
\end{align*}
```

- ``\vec{W}`` is conservative ``\to (\frac{\partial W_y}{\partial x} - \frac{\partial W_x}{\partial y}) = 0``

---

The curl is given by
```math
\begin{align*}
\vec{\nabla} \times \vec{W} = \left(\frac{\partial W_z}{\partial y} - \frac{\partial W_y}{\partial z}\right) \hat{i} - \left(\frac{\partial W_x}{\partial z} - \frac{\partial W_z}{\partial x}\right) \hat{j} + \left(\frac{\partial W_y}{\partial x} - \frac{\partial W_x}{\partial y}\right) \hat{k}
\end{align*}
```

- If a field is conservative, its curl vanishes everywhere and vice versa. Look up the curl in general coordinates if needed

---

Stokes' theorem
```math
\begin{align*}
\oint_{C = \partial S} \vec{W} \cdot \vec{dr} = \int_S \vec{\nabla } \times \vec{W} \cdot \vec{dS}
\end{align*}
```

---

Gauss's law

```math
\begin{align*}
\oint_{S = \partial V} \vec{W} \cdot \vec{dS} = \int_V \vec{\nabla} \cdot \vec{W} dxdydz
\end{align*}
```

- Where ``\vec{\nabla} \cdot \vec{W}`` is the divergence of ``\vec{W}``

---

The laplacian ``\nabla^2`` is defined as follows

```math
\begin{align*}
\nabla \cdot \nabla \phi = \frac{\partial^2 \phi}{\partial x^2} + \frac{\partial^2 \phi}{\partial y^2} + \frac{\partial^2 \phi}{\partial z^2}
\end{align*}
```

- For general coordinates, look it up when needed

"""

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
CairoMakie = "13f3f980-e62b-5c42-98c6-ff1f3baf88f0"
CondaPkg = "992eb4ea-22a4-4c89-a5bb-47a3300528ab"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
PythonCall = "6099a3de-0909-46bc-b1f4-468b9a2dfc0d"

[compat]
CairoMakie = "~0.10.7"
CondaPkg = "~0.2.18"
PlutoUI = "~0.7.52"
PythonCall = "~0.9.13"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.9.2"
manifest_format = "2.0"
project_hash = "d23030f7311d94861c8bd7107c1347a09cb2fc7d"

[[deps.AbstractFFTs]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "cad4c758c0038eea30394b1b671526921ca85b21"
uuid = "621f4979-c628-5d54-868e-fcf4e3e8185c"
version = "1.4.0"
weakdeps = ["ChainRulesCore"]

    [deps.AbstractFFTs.extensions]
    AbstractFFTsChainRulesCoreExt = "ChainRulesCore"

[[deps.AbstractLattices]]
git-tree-sha1 = "f35684b7349da49fcc8a9e520e30e45dbb077166"
uuid = "398f06c4-4d28-53ec-89ca-5b2656b7603d"
version = "0.2.1"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "91bd53c39b9cbfb5ef4b015e8b582d344532bd0a"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.2.0"

[[deps.AbstractTrees]]
git-tree-sha1 = "faa260e4cb5aba097a73fab382dd4b5819d8ec8c"
uuid = "1520ce14-60c1-5f80-bbc7-55ef81b5835c"
version = "0.4.4"

[[deps.Adapt]]
deps = ["LinearAlgebra", "Requires"]
git-tree-sha1 = "76289dc51920fdc6e0013c872ba9551d54961c24"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "3.6.2"
weakdeps = ["StaticArrays"]

    [deps.Adapt.extensions]
    AdaptStaticArraysExt = "StaticArrays"

[[deps.Animations]]
deps = ["Colors"]
git-tree-sha1 = "e81c509d2c8e49592413bfb0bb3b08150056c79d"
uuid = "27a7e980-b3e6-11e9-2bcd-0b925532e340"
version = "0.4.1"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.1"

[[deps.ArrayInterface]]
deps = ["Adapt", "LinearAlgebra", "Requires", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "f83ec24f76d4c8f525099b2ac475fc098138ec31"
uuid = "4fba245c-0d91-5ea0-9b3e-6abc04ee57a9"
version = "7.4.11"

    [deps.ArrayInterface.extensions]
    ArrayInterfaceBandedMatricesExt = "BandedMatrices"
    ArrayInterfaceBlockBandedMatricesExt = "BlockBandedMatrices"
    ArrayInterfaceCUDAExt = "CUDA"
    ArrayInterfaceGPUArraysCoreExt = "GPUArraysCore"
    ArrayInterfaceStaticArraysCoreExt = "StaticArraysCore"
    ArrayInterfaceTrackerExt = "Tracker"

    [deps.ArrayInterface.weakdeps]
    BandedMatrices = "aae01518-5342-5314-be14-df237901396f"
    BlockBandedMatrices = "ffab5731-97b5-5995-9138-79e8c1846df0"
    CUDA = "052768ef-5323-5732-b1bb-66c8b64840ba"
    GPUArraysCore = "46192b85-c4d5-4398-a991-12ede77f4527"
    StaticArraysCore = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Automa]]
deps = ["ScanByte", "TranscodingStreams"]
git-tree-sha1 = "48e54446df62fdf9ef76959c32dc33f3cff659ee"
uuid = "67c07d97-cdcb-5c2c-af73-a7f9c32a568b"
version = "0.8.3"

[[deps.AxisAlgorithms]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "WoodburyMatrices"]
git-tree-sha1 = "66771c8d21c8ff5e3a93379480a2307ac36863f7"
uuid = "13072b0f-2c55-5437-9ae7-d433b7a33950"
version = "1.0.1"

[[deps.AxisArrays]]
deps = ["Dates", "IntervalSets", "IterTools", "RangeArrays"]
git-tree-sha1 = "16351be62963a67ac4083f748fdb3cca58bfd52f"
uuid = "39de3d68-74b9-583c-8d2d-e117c070f3a9"
version = "0.4.7"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "19a35467a82e236ff51bc17a3a44b69ef35185a2"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.8+0"

[[deps.CEnum]]
git-tree-sha1 = "eb4cb44a499229b3b8426dcfb5dd85333951ff90"
uuid = "fa961155-64e5-5f13-b03f-caf6b980ea82"
version = "0.4.2"

[[deps.CRC32c]]
uuid = "8bf52ea8-c179-5cab-976a-9e18b702a9bc"

[[deps.CRlibm]]
deps = ["CRlibm_jll"]
git-tree-sha1 = "32abd86e3c2025db5172aa182b982debed519834"
uuid = "96374032-68de-5a5b-8d9e-752f78720389"
version = "1.0.1"

[[deps.CRlibm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e329286945d0cfc04456972ea732551869af1cfc"
uuid = "4e9b3aee-d8a1-5a3d-ad8b-7d824db253f0"
version = "1.0.1+0"

[[deps.Cairo]]
deps = ["Cairo_jll", "Colors", "Glib_jll", "Graphics", "Libdl", "Pango_jll"]
git-tree-sha1 = "d0b3f8b4ad16cb0a2988c6788646a5e6a17b6b1b"
uuid = "159f3aea-2a34-519c-b102-8c37f9878175"
version = "1.0.5"

[[deps.CairoMakie]]
deps = ["Base64", "Cairo", "Colors", "FFTW", "FileIO", "FreeType", "GeometryBasics", "LinearAlgebra", "Makie", "PrecompileTools", "SHA"]
git-tree-sha1 = "e041782fed7614b1726fa250f2bf24fd5c789689"
uuid = "13f3f980-e62b-5c42-98c6-ff1f3baf88f0"
version = "0.10.7"

[[deps.Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "CompilerSupportLibraries_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "4b859a208b2397a7a623a03449e4636bdb17bcf2"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.16.1+1"

[[deps.Calculus]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "f641eb0a4f00c343bbc32346e1217b86f3ce9dad"
uuid = "49dc2e85-a5d0-5ad3-a950-438e2897f1b9"
version = "0.5.1"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "e30f2f4e20f7f186dc36529910beaedc60cfa644"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.16.0"

[[deps.ColorBrewer]]
deps = ["Colors", "JSON", "Test"]
git-tree-sha1 = "61c5334f33d91e570e1d0c3eb5465835242582c4"
uuid = "a2cac450-b92f-5266-8821-25eda20663c8"
version = "0.4.0"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "ColorVectorSpace", "Colors", "FixedPointNumbers", "PrecompileTools", "Random"]
git-tree-sha1 = "dd3000d954d483c1aad05fe1eb9e6a715c97013e"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.22.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "eb7f0f8307f71fac7c606984ea5fb2817275d6e4"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.4"

[[deps.ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "SpecialFunctions", "Statistics", "TensorCore"]
git-tree-sha1 = "600cc5508d66b78aae350f7accdb58763ac18589"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.9.10"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "fc08e5930ee9a4e03f84bfb5211cb54e7769758a"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.10"

[[deps.Combinatorics]]
git-tree-sha1 = "08c8b6831dc00bfea825826be0bc8336fc369860"
uuid = "861a8166-3701-5b0c-9a16-15d98fcdc6aa"
version = "1.0.2"

[[deps.CommonSubexpressions]]
deps = ["MacroTools", "Test"]
git-tree-sha1 = "7b8a93dba8af7e3b42fecabf646260105ac373f7"
uuid = "bbf7d656-a473-5ed7-a52c-81e309532950"
version = "0.3.0"

[[deps.Compat]]
deps = ["UUIDs"]
git-tree-sha1 = "5ce999a19f4ca23ea484e92a1774a61b8ca4cf8e"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.8.0"
weakdeps = ["Dates", "LinearAlgebra"]

    [deps.Compat.extensions]
    CompatLinearAlgebraExt = "LinearAlgebra"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.0.5+0"

[[deps.CondaPkg]]
deps = ["JSON3", "Markdown", "MicroMamba", "Pidfile", "Pkg", "TOML"]
git-tree-sha1 = "741146cf2ced5859faae76a84b541aa9af1a78bb"
uuid = "992eb4ea-22a4-4c89-a5bb-47a3300528ab"
version = "0.2.18"

[[deps.ConstructionBase]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "fe2838a593b5f776e1597e086dcd47560d94e816"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.5.3"
weakdeps = ["IntervalSets", "StaticArrays"]

    [deps.ConstructionBase.extensions]
    ConstructionBaseIntervalSetsExt = "IntervalSets"
    ConstructionBaseStaticArraysExt = "StaticArrays"

[[deps.Contour]]
git-tree-sha1 = "d05d9e7b7aedff4e5b51a029dced05cfb6125781"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.6.2"

[[deps.DataAPI]]
git-tree-sha1 = "8da84edb865b0b5b0100c0666a9bc9a0b71c553c"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.15.0"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "cf25ccb972fec4e4817764d01c82386ae94f77b4"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.14"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.DelaunayTriangulation]]
deps = ["DataStructures", "EnumX", "ExactPredicates", "Random", "SimpleGraphs"]
git-tree-sha1 = "a1d8532de83f8ce964235eff1edeff9581144d02"
uuid = "927a84f5-c5f4-47a5-9785-b46e178433df"
version = "0.7.2"
weakdeps = ["MakieCore"]

    [deps.DelaunayTriangulation.extensions]
    DelaunayTriangulationMakieCoreExt = "MakieCore"

[[deps.DiffResults]]
deps = ["StaticArraysCore"]
git-tree-sha1 = "782dd5f4561f5d267313f23853baaaa4c52ea621"
uuid = "163ba53b-c6d8-5494-b064-1a9d43ac40c5"
version = "1.1.0"

[[deps.DiffRules]]
deps = ["IrrationalConstants", "LogExpFunctions", "NaNMath", "Random", "SpecialFunctions"]
git-tree-sha1 = "23163d55f885173722d1e4cf0f6110cdbaf7e272"
uuid = "b552c78f-8df3-52c6-915a-8e097449b14b"
version = "1.15.1"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[deps.Distributions]]
deps = ["FillArrays", "LinearAlgebra", "PDMats", "Printf", "QuadGK", "Random", "SparseArrays", "SpecialFunctions", "Statistics", "StatsAPI", "StatsBase", "StatsFuns", "Test"]
git-tree-sha1 = "e76a3281de2719d7c81ed62c6ea7057380c87b1d"
uuid = "31c24e10-a181-5473-b8eb-7969acd0382f"
version = "0.25.98"

    [deps.Distributions.extensions]
    DistributionsChainRulesCoreExt = "ChainRulesCore"
    DistributionsDensityInterfaceExt = "DensityInterface"

    [deps.Distributions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    DensityInterface = "b429d917-457f-4dbc-8f4c-0cc954292b1d"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "2fb1e02f2b635d0845df5d7c167fec4dd739b00d"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.3"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.DualNumbers]]
deps = ["Calculus", "NaNMath", "SpecialFunctions"]
git-tree-sha1 = "5837a837389fccf076445fce071c8ddaea35a566"
uuid = "fa6b7ba4-c1ee-5f82-b5fc-ecf0adba8f74"
version = "0.6.8"

[[deps.EarCut_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e3290f2d49e661fbd94046d7e3726ffcb2d41053"
uuid = "5ae413db-bbd1-5e63-b57d-d24a61df00f5"
version = "2.2.4+0"

[[deps.EnumX]]
git-tree-sha1 = "bdb1942cd4c45e3c678fd11569d5cccd80976237"
uuid = "4e289a0a-7415-4d19-859d-a7e5c4648b56"
version = "1.0.4"

[[deps.ErrorfreeArithmetic]]
git-tree-sha1 = "d6863c556f1142a061532e79f611aa46be201686"
uuid = "90fa49ef-747e-5e6f-a989-263ba693cf1a"
version = "0.5.2"

[[deps.ExactPredicates]]
deps = ["IntervalArithmetic", "Random", "StaticArraysCore", "Test"]
git-tree-sha1 = "276e83bc8b21589b79303b9985c321024ffdf59c"
uuid = "429591f6-91af-11e9-00e2-59fbe8cec110"
version = "2.2.5"

[[deps.Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "4558ab818dcceaab612d1bb8c19cee87eda2b83c"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.5.0+0"

[[deps.Extents]]
git-tree-sha1 = "5e1e4c53fa39afe63a7d356e30452249365fba99"
uuid = "411431e0-e8b7-467b-b5e0-f676ba4f2910"
version = "0.1.1"

[[deps.FFMPEG]]
deps = ["FFMPEG_jll"]
git-tree-sha1 = "b57e3acbe22f8484b4b5ff66a7499717fe1a9cc8"
uuid = "c87230d0-a227-11e9-1b43-d7ebe4e7570a"
version = "0.4.1"

[[deps.FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "PCRE2_jll", "Pkg", "Zlib_jll", "libaom_jll", "libass_jll", "libfdk_aac_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "74faea50c1d007c85837327f6775bea60b5492dd"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "4.4.2+2"

[[deps.FFTW]]
deps = ["AbstractFFTs", "FFTW_jll", "LinearAlgebra", "MKL_jll", "Preferences", "Reexport"]
git-tree-sha1 = "b4fbdd20c889804969571cc589900803edda16b7"
uuid = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
version = "1.7.1"

[[deps.FFTW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c6033cc3892d0ef5bb9cd29b7f2f0331ea5184ea"
uuid = "f5851436-0d7a-5f13-b9de-f02708fd171a"
version = "3.3.10+0"

[[deps.FastRounding]]
deps = ["ErrorfreeArithmetic", "LinearAlgebra"]
git-tree-sha1 = "6344aa18f654196be82e62816935225b3b9abe44"
uuid = "fa42c844-2597-5d31-933b-ebd51ab2693f"
version = "0.3.1"

[[deps.FileIO]]
deps = ["Pkg", "Requires", "UUIDs"]
git-tree-sha1 = "299dc33549f68299137e51e6d49a13b5b1da9673"
uuid = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
version = "1.16.1"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.FillArrays]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "Statistics"]
git-tree-sha1 = "f0af9b12329a637e8fba7d6543f915fff6ba0090"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "1.4.2"

[[deps.FiniteDiff]]
deps = ["ArrayInterface", "LinearAlgebra", "Requires", "Setfield", "SparseArrays"]
git-tree-sha1 = "c6e4a1fbe73b31a3dea94b1da449503b8830c306"
uuid = "6a86dc24-6348-571c-b903-95158fe2bd41"
version = "2.21.1"

    [deps.FiniteDiff.extensions]
    FiniteDiffBandedMatricesExt = "BandedMatrices"
    FiniteDiffBlockBandedMatricesExt = "BlockBandedMatrices"
    FiniteDiffStaticArraysExt = "StaticArrays"

    [deps.FiniteDiff.weakdeps]
    BandedMatrices = "aae01518-5342-5314-be14-df237901396f"
    BlockBandedMatrices = "ffab5731-97b5-5995-9138-79e8c1846df0"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[deps.Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "21efd19106a55620a188615da6d3d06cd7f6ee03"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.13.93+0"

[[deps.Formatting]]
deps = ["Printf"]
git-tree-sha1 = "8339d61043228fdd3eb658d86c926cb282ae72a8"
uuid = "59287772-0a20-5a39-b81b-1366585eb4c0"
version = "0.4.2"

[[deps.ForwardDiff]]
deps = ["CommonSubexpressions", "DiffResults", "DiffRules", "LinearAlgebra", "LogExpFunctions", "NaNMath", "Preferences", "Printf", "Random", "SpecialFunctions"]
git-tree-sha1 = "00e252f4d706b3d55a8863432e742bf5717b498d"
uuid = "f6369f11-7733-5829-9624-2563aa707210"
version = "0.10.35"
weakdeps = ["StaticArrays"]

    [deps.ForwardDiff.extensions]
    ForwardDiffStaticArraysExt = "StaticArrays"

[[deps.FreeType]]
deps = ["CEnum", "FreeType2_jll"]
git-tree-sha1 = "cabd77ab6a6fdff49bfd24af2ebe76e6e018a2b4"
uuid = "b38be410-82b0-50bf-ab77-7b57e271db43"
version = "4.0.0"

[[deps.FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "87eb71354d8ec1a96d4a7636bd57a7347dde3ef9"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.10.4+0"

[[deps.FreeTypeAbstraction]]
deps = ["ColorVectorSpace", "Colors", "FreeType", "GeometryBasics"]
git-tree-sha1 = "38a92e40157100e796690421e34a11c107205c86"
uuid = "663a7486-cb36-511b-a19d-713bb74d65c9"
version = "0.10.0"

[[deps.FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "aa31987c2ba8704e23c6c8ba8a4f769d5d7e4f91"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.10+0"

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"

[[deps.GPUArraysCore]]
deps = ["Adapt"]
git-tree-sha1 = "2d6ca471a6c7b536127afccfa7564b5b39227fe0"
uuid = "46192b85-c4d5-4398-a991-12ede77f4527"
version = "0.1.5"

[[deps.GeoInterface]]
deps = ["Extents"]
git-tree-sha1 = "bb198ff907228523f3dee1070ceee63b9359b6ab"
uuid = "cf35fbd7-0cd7-5166-be24-54bfbe79505f"
version = "1.3.1"

[[deps.GeometryBasics]]
deps = ["EarCut_jll", "Extents", "GeoInterface", "IterTools", "LinearAlgebra", "StaticArrays", "StructArrays", "Tables"]
git-tree-sha1 = "424a5a6ce7c5d97cca7bcc4eac551b97294c54af"
uuid = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
version = "0.4.9"

[[deps.Gettext_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "9b02998aba7bf074d14de89f9d37ca24a1a0b046"
uuid = "78b55507-aeef-58d4-861c-77aaff3498b1"
version = "0.21.0+0"

[[deps.Glib_jll]]
deps = ["Artifacts", "Gettext_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE2_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "d3b3624125c1474292d0d8ed0f65554ac37ddb23"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.74.0+2"

[[deps.Graphics]]
deps = ["Colors", "LinearAlgebra", "NaNMath"]
git-tree-sha1 = "d61890399bc535850c4bf08e4e0d3a7ad0f21cbd"
uuid = "a2bd30eb-e257-5431-a919-1863eab51364"
version = "1.1.2"

[[deps.Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "344bf40dcab1073aca04aa0df4fb092f920e4011"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.14+0"

[[deps.GridLayoutBase]]
deps = ["GeometryBasics", "InteractiveUtils", "Observables"]
git-tree-sha1 = "678d136003ed5bceaab05cf64519e3f956ffa4ba"
uuid = "3955a311-db13-416c-9275-1d80ed98e5e9"
version = "0.9.1"

[[deps.Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[deps.HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg"]
git-tree-sha1 = "129acf094d168394e80ee1dc4bc06ec835e510a3"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "2.8.1+1"

[[deps.HypergeometricFunctions]]
deps = ["DualNumbers", "LinearAlgebra", "OpenLibm_jll", "SpecialFunctions"]
git-tree-sha1 = "83e95aaab9dc184a6dcd9c4c52aa0dc26cd14a1d"
uuid = "34004b35-14d8-5ef3-9330-4cdb6864b03a"
version = "0.3.21"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "8d511d5b81240fc8e6802386302675bdf47737b9"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.4"

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "c47c5fa4c5308f27ccaac35504858d8914e102f9"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.4"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "d75853a0bdbfb1ac815478bacd89cd27b550ace6"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.3"

[[deps.ImageAxes]]
deps = ["AxisArrays", "ImageBase", "ImageCore", "Reexport", "SimpleTraits"]
git-tree-sha1 = "2e4520d67b0cef90865b3ef727594d2a58e0e1f8"
uuid = "2803e5a7-5153-5ecf-9a86-9b4c37f5f5ac"
version = "0.6.11"

[[deps.ImageBase]]
deps = ["ImageCore", "Reexport"]
git-tree-sha1 = "b51bb8cae22c66d0f6357e3bcb6363145ef20835"
uuid = "c817782e-172a-44cc-b673-b171935fbb9e"
version = "0.1.5"

[[deps.ImageCore]]
deps = ["AbstractFFTs", "ColorVectorSpace", "Colors", "FixedPointNumbers", "Graphics", "MappedArrays", "MosaicViews", "OffsetArrays", "PaddedViews", "Reexport"]
git-tree-sha1 = "acf614720ef026d38400b3817614c45882d75500"
uuid = "a09fc81d-aa75-5fe9-8630-4744c3626534"
version = "0.9.4"

[[deps.ImageIO]]
deps = ["FileIO", "IndirectArrays", "JpegTurbo", "LazyModules", "Netpbm", "OpenEXR", "PNGFiles", "QOI", "Sixel", "TiffImages", "UUIDs"]
git-tree-sha1 = "bca20b2f5d00c4fbc192c3212da8fa79f4688009"
uuid = "82e4d734-157c-48bb-816b-45c225c6df19"
version = "0.6.7"

[[deps.ImageMetadata]]
deps = ["AxisArrays", "ImageAxes", "ImageBase", "ImageCore"]
git-tree-sha1 = "355e2b974f2e3212a75dfb60519de21361ad3cb7"
uuid = "bc367c6b-8a6b-528e-b4bd-a4b897500b49"
version = "0.9.9"

[[deps.Imath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "3d09a9f60edf77f8a4d99f9e015e8fbf9989605d"
uuid = "905a6f67-0a94-5f89-b386-d35d92009cd1"
version = "3.1.7+0"

[[deps.IndirectArrays]]
git-tree-sha1 = "012e604e1c7458645cb8b436f8fba789a51b257f"
uuid = "9b13fd28-a010-5f03-acff-a1bbcff69959"
version = "1.0.0"

[[deps.Inflate]]
git-tree-sha1 = "5cd07aab533df5170988219191dfad0519391428"
uuid = "d25df0c9-e2be-5dd7-82c8-3ad0b3e990b9"
version = "0.1.3"

[[deps.IntegerMathUtils]]
git-tree-sha1 = "b8ffb903da9f7b8cf695a8bead8e01814aa24b30"
uuid = "18e54dd8-cb9d-406c-a71d-865a43cbb235"
version = "0.1.2"

[[deps.IntelOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "0cb9352ef2e01574eeebdb102948a58740dcaf83"
uuid = "1d5cc7b8-4909-519e-a0f8-d0f5ad9712d0"
version = "2023.1.0+0"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.Interpolations]]
deps = ["Adapt", "AxisAlgorithms", "ChainRulesCore", "LinearAlgebra", "OffsetArrays", "Random", "Ratios", "Requires", "SharedArrays", "SparseArrays", "StaticArrays", "WoodburyMatrices"]
git-tree-sha1 = "721ec2cf720536ad005cb38f50dbba7b02419a15"
uuid = "a98d9a8b-a2ab-59e6-89dd-64a1c18fca59"
version = "0.14.7"

[[deps.IntervalArithmetic]]
deps = ["CRlibm", "FastRounding", "LinearAlgebra", "Markdown", "Random", "RecipesBase", "RoundingEmulator", "SetRounding", "StaticArrays"]
git-tree-sha1 = "5ab7744289be503d76a944784bac3f2df7b809af"
uuid = "d1acc4aa-44c8-5952-acd4-ba5d80a2a253"
version = "0.20.9"

[[deps.IntervalSets]]
deps = ["Dates", "Random", "Statistics"]
git-tree-sha1 = "16c0cc91853084cb5f58a78bd209513900206ce6"
uuid = "8197267c-284f-5f27-9208-e0e47529a953"
version = "0.7.4"

[[deps.IrrationalConstants]]
git-tree-sha1 = "630b497eafcc20001bba38a4651b327dcfc491d2"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.2.2"

[[deps.Isoband]]
deps = ["isoband_jll"]
git-tree-sha1 = "f9b6d97355599074dc867318950adaa6f9946137"
uuid = "f1662d9f-8043-43de-a69a-05efc1cc6ff4"
version = "0.1.1"

[[deps.IterTools]]
git-tree-sha1 = "4ced6667f9974fc5c5943fa5e2ef1ca43ea9e450"
uuid = "c8e1da08-722c-5040-9ed9-7db0dc04731e"
version = "1.8.0"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.JLLWrappers]]
deps = ["Preferences"]
git-tree-sha1 = "abc9885a7ca2052a736a600f7fa66209f96506e1"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.4.1"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "31e996f0a15c7b280ba9f76636b3ff9e2ae58c9a"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.4"

[[deps.JSON3]]
deps = ["Dates", "Mmap", "Parsers", "PrecompileTools", "StructTypes", "UUIDs"]
git-tree-sha1 = "5b62d93f2582b09e469b3099d839c2d2ebf5066d"
uuid = "0f8b85d8-7281-11e9-16c2-39a750bddbf1"
version = "1.13.1"

[[deps.JpegTurbo]]
deps = ["CEnum", "FileIO", "ImageCore", "JpegTurbo_jll", "TOML"]
git-tree-sha1 = "327713faef2a3e5c80f96bf38d1fa26f7a6ae29e"
uuid = "b835a17e-a41a-41e7-81f0-2f016b05efe0"
version = "0.1.3"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "6f2675ef130a300a112286de91973805fcc5ffbc"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "2.1.91+0"

[[deps.KernelDensity]]
deps = ["Distributions", "DocStringExtensions", "FFTW", "Interpolations", "StatsBase"]
git-tree-sha1 = "90442c50e202a5cdf21a7899c66b240fdef14035"
uuid = "5ab0869b-81aa-558d-bb23-cbf5423bbe9b"
version = "0.6.7"

[[deps.LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "f6250b16881adf048549549fba48b1161acdac8c"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.1+0"

[[deps.LLVMOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "f689897ccbe049adb19a065c495e75f372ecd42b"
uuid = "1d63c593-3942-5779-bab2-d838dc0a180e"
version = "15.0.4+0"

[[deps.LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e5b909bcf985c5e2605737d2ce278ed791b89be6"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.1+0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "f2355693d6778a178ade15952b7ac47a4ff97996"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.0"

[[deps.LazyArtifacts]]
deps = ["Artifacts", "Pkg"]
uuid = "4af54fe1-eca0-43a8-85a7-787d91b784e3"

[[deps.LazyModules]]
git-tree-sha1 = "a560dd966b386ac9ae60bdd3a3d3a326062d3c3e"
uuid = "8cdb02fc-e678-4876-92c5-9defec4f444e"
version = "0.3.1"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.3"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "7.84.0+0"

[[deps.LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.10.2+0"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "0b4a5d71f3e5200a7dff793393e09dfc2d874290"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.2.2+1"

[[deps.Libgcrypt_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgpg_error_jll", "Pkg"]
git-tree-sha1 = "64613c82a59c120435c067c2b809fc61cf5166ae"
uuid = "d4300ac3-e22c-5743-9152-c294e39db1e4"
version = "1.8.7+0"

[[deps.Libgpg_error_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c333716e46366857753e273ce6a69ee0945a6db9"
uuid = "7add5ba3-2f88-524e-9cd5-f83b8a55f7b8"
version = "1.42.0+0"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c7cb1f5d892775ba13767a87c7ada0b980ea0a71"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.16.1+2"

[[deps.Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9c30530bf0effd46e15e0fdcf2b8636e78cbbd73"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.35.0+0"

[[deps.Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "7f3efec06033682db852f8b3bc3c1d2b0a0ab066"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.36.0+0"

[[deps.LightXML]]
deps = ["Libdl", "XML2_jll"]
git-tree-sha1 = "e129d9391168c677cd4800f5c0abb1ed8cb3794f"
uuid = "9c8b4983-aa76-5018-a973-4c85ecc9e179"
version = "0.9.0"

[[deps.LineSearches]]
deps = ["LinearAlgebra", "NLSolversBase", "NaNMath", "Parameters", "Printf"]
git-tree-sha1 = "7bbea35cec17305fc70a0e5b4641477dc0789d9d"
uuid = "d3d80556-e9d4-5f37-9878-2ab0fcc64255"
version = "7.2.0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.LinearAlgebraX]]
deps = ["LinearAlgebra", "Mods", "Permutations", "Primes", "SimplePolynomials"]
git-tree-sha1 = "558a338f1eeabe933f9c2d4052aa7c2c707c3d52"
uuid = "9b3f67b0-2d00-526e-9884-9e4938f8fb88"
version = "0.1.12"

[[deps.LogExpFunctions]]
deps = ["DocStringExtensions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "c3ce8e7420b3a6e071e0fe4745f5d4300e37b13f"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.24"

    [deps.LogExpFunctions.extensions]
    LogExpFunctionsChainRulesCoreExt = "ChainRulesCore"
    LogExpFunctionsChangesOfVariablesExt = "ChangesOfVariables"
    LogExpFunctionsInverseFunctionsExt = "InverseFunctions"

    [deps.LogExpFunctions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    ChangesOfVariables = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.MIMEs]]
git-tree-sha1 = "65f28ad4b594aebe22157d6fac869786a255b7eb"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "0.1.4"

[[deps.MKL_jll]]
deps = ["Artifacts", "IntelOpenMP_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "Pkg"]
git-tree-sha1 = "154d7aaa82d24db6d8f7e4ffcfe596f40bff214b"
uuid = "856f044c-d86e-5d09-b602-aeab76dc8ba7"
version = "2023.1.0+0"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "42324d08725e200c23d4dfb549e0d5d89dede2d2"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.10"

[[deps.Makie]]
deps = ["Animations", "Base64", "ColorBrewer", "ColorSchemes", "ColorTypes", "Colors", "Contour", "DelaunayTriangulation", "Distributions", "DocStringExtensions", "Downloads", "FFMPEG", "FileIO", "FixedPointNumbers", "Formatting", "FreeType", "FreeTypeAbstraction", "GeometryBasics", "GridLayoutBase", "ImageIO", "InteractiveUtils", "IntervalSets", "Isoband", "KernelDensity", "LaTeXStrings", "LinearAlgebra", "MacroTools", "MakieCore", "Markdown", "Match", "MathTeXEngine", "Observables", "OffsetArrays", "Packing", "PlotUtils", "PolygonOps", "PrecompileTools", "Printf", "REPL", "Random", "RelocatableFolders", "Setfield", "ShaderAbstractions", "Showoff", "SignedDistanceFields", "SparseArrays", "StableHashTraits", "Statistics", "StatsBase", "StatsFuns", "StructArrays", "TriplotBase", "UnicodeFun"]
git-tree-sha1 = "729640354756782c89adba8857085a69e19be7ab"
uuid = "ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a"
version = "0.19.7"

[[deps.MakieCore]]
deps = ["Observables"]
git-tree-sha1 = "87a85ff81583bd392642869557cb633532989517"
uuid = "20f20a25-4f0e-4fdf-b5d1-57303727442b"
version = "0.6.4"

[[deps.MappedArrays]]
git-tree-sha1 = "2dab0221fe2b0f2cb6754eaa743cc266339f527e"
uuid = "dbb5928d-eab1-5f90-85c2-b9b0edb7c900"
version = "0.4.2"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.Match]]
git-tree-sha1 = "1d9bc5c1a6e7ee24effb93f175c9342f9154d97f"
uuid = "7eb4fadd-790c-5f42-8a69-bfa0b872bfbf"
version = "1.2.0"

[[deps.MathTeXEngine]]
deps = ["AbstractTrees", "Automa", "DataStructures", "FreeTypeAbstraction", "GeometryBasics", "LaTeXStrings", "REPL", "RelocatableFolders", "Test", "UnicodeFun"]
git-tree-sha1 = "8f52dbaa1351ce4cb847d95568cb29e62a307d93"
uuid = "0a4f8689-d25c-4efe-a92b-7142dfc1aa53"
version = "0.5.6"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.2+0"

[[deps.MicroMamba]]
deps = ["Pkg", "Scratch", "micromamba_jll"]
git-tree-sha1 = "6f0e43750a94574c18933e9456b18d4d94a4a671"
uuid = "0b3b1443-0f03-428d-bdfb-f27f9c1191ea"
version = "0.1.13"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "f66bdc5de519e8f8ae43bdc598782d35a25b1272"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.1.0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.Mods]]
git-tree-sha1 = "61be59e4daffff43a8cec04b5e0dc773cbb5db3a"
uuid = "7475f97c-0381-53b1-977b-4c60186c8d62"
version = "1.3.3"

[[deps.MosaicViews]]
deps = ["MappedArrays", "OffsetArrays", "PaddedViews", "StackViews"]
git-tree-sha1 = "7b86a5d4d70a9f5cdf2dacb3cbe6d251d1a61dbe"
uuid = "e94cdb99-869f-56ef-bcf0-1ae2bcbe0389"
version = "0.3.4"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2022.10.11"

[[deps.Multisets]]
git-tree-sha1 = "8d852646862c96e226367ad10c8af56099b4047e"
uuid = "3b2b4ff1-bcff-5658-a3ee-dbcf1ce5ac09"
version = "0.4.4"

[[deps.NLSolversBase]]
deps = ["DiffResults", "Distributed", "FiniteDiff", "ForwardDiff"]
git-tree-sha1 = "a0b464d183da839699f4c79e7606d9d186ec172c"
uuid = "d41bc354-129a-5804-8e4c-c37616107c6c"
version = "7.8.3"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "0877504529a3e5c3343c6f8b4c0381e57e4387e4"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.0.2"

[[deps.Netpbm]]
deps = ["FileIO", "ImageCore", "ImageMetadata"]
git-tree-sha1 = "d92b107dbb887293622df7697a2223f9f8176fcd"
uuid = "f09324ee-3d7c-5217-9330-fc30815ba969"
version = "1.1.1"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.Observables]]
git-tree-sha1 = "6862738f9796b3edc1c09d0890afce4eca9e7e93"
uuid = "510215fc-4207-5dde-b226-833fc4488ee2"
version = "0.5.4"

[[deps.OffsetArrays]]
deps = ["Adapt"]
git-tree-sha1 = "2ac17d29c523ce1cd38e27785a7d23024853a4bb"
uuid = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
version = "1.12.10"

[[deps.Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "887579a3eb005446d514ab7aeac5d1d027658b8f"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.5+1"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.21+4"

[[deps.OpenEXR]]
deps = ["Colors", "FileIO", "OpenEXR_jll"]
git-tree-sha1 = "327f53360fdb54df7ecd01e96ef1983536d1e633"
uuid = "52e1d378-f018-4a11-a4be-720524705ac7"
version = "0.3.2"

[[deps.OpenEXR_jll]]
deps = ["Artifacts", "Imath_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "a4ca623df1ae99d09bc9868b008262d0c0ac1e4f"
uuid = "18a262bb-aa17-5467-a713-aee519bc75cb"
version = "3.1.4+0"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.1+0"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1aa4b74f80b01c6bc2b89992b861b5f210e665b5"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "1.1.21+0"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[deps.Optim]]
deps = ["Compat", "FillArrays", "ForwardDiff", "LineSearches", "LinearAlgebra", "NLSolversBase", "NaNMath", "Parameters", "PositiveFactorizations", "Printf", "SparseArrays", "StatsBase"]
git-tree-sha1 = "e3a6546c1577bfd701771b477b794a52949e7594"
uuid = "429524aa-4258-5aef-a3af-852621145aeb"
version = "1.7.6"

[[deps.Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "51a08fb14ec28da2ec7a927c4337e4332c2a4720"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.3.2+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "2e73fe17cac3c62ad1aebe70d44c963c3cfdc3e3"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.6.2"

[[deps.PCRE2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "efcefdf7-47ab-520b-bdef-62a2eaa19f15"
version = "10.42.0+0"

[[deps.PDMats]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "67eae2738d63117a196f497d7db789821bce61d1"
uuid = "90014a1f-27ba-587c-ab20-58faa44d9150"
version = "0.11.17"

[[deps.PNGFiles]]
deps = ["Base64", "CEnum", "ImageCore", "IndirectArrays", "OffsetArrays", "libpng_jll"]
git-tree-sha1 = "9b02b27ac477cad98114584ff964e3052f656a0f"
uuid = "f57f5aa1-a3ce-4bc8-8ab9-96f992907883"
version = "0.4.0"

[[deps.Packing]]
deps = ["GeometryBasics"]
git-tree-sha1 = "ec3edfe723df33528e085e632414499f26650501"
uuid = "19eb6ba3-879d-56ad-ad62-d5c202156566"
version = "0.5.0"

[[deps.PaddedViews]]
deps = ["OffsetArrays"]
git-tree-sha1 = "0fac6313486baae819364c52b4f483450a9d793f"
uuid = "5432bcbf-9aad-5242-b902-cca2824c8663"
version = "0.5.12"

[[deps.Pango_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "FriBidi_jll", "Glib_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "84a314e3926ba9ec66ac097e3635e270986b0f10"
uuid = "36c8627f-9965-5494-a995-c6b170f724f3"
version = "1.50.9+0"

[[deps.Parameters]]
deps = ["OrderedCollections", "UnPack"]
git-tree-sha1 = "34c0e9ad262e5f7fc75b10a9952ca7692cfc5fbe"
uuid = "d96e819e-fc66-5662-9728-84c9c7592b0a"
version = "0.12.3"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "4b2e829ee66d4218e0cef22c0a64ee37cf258c29"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.7.1"

[[deps.Permutations]]
deps = ["Combinatorics", "LinearAlgebra", "Random"]
git-tree-sha1 = "6e6cab1c54ae2382bcc48866b91cf949cea703a1"
uuid = "2ae35dd2-176d-5d53-8349-f30d82d94d4f"
version = "0.4.16"

[[deps.Pidfile]]
deps = ["FileWatching", "Test"]
git-tree-sha1 = "2d8aaf8ee10df53d0dfb9b8ee44ae7c04ced2b03"
uuid = "fa939f87-e72e-5be4-a000-7fc836dbe307"
version = "1.3.0"

[[deps.Pixman_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "LLVMOpenMP_jll", "Libdl"]
git-tree-sha1 = "64779bc4c9784fee475689a1752ef4d5747c5e87"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.42.2+0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.9.2"

[[deps.PkgVersion]]
deps = ["Pkg"]
git-tree-sha1 = "f6cf8e7944e50901594838951729a1861e668cb8"
uuid = "eebad327-c553-4316-9ea0-9fa01ccd7688"
version = "0.3.2"

[[deps.PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "PrecompileTools", "Printf", "Random", "Reexport", "Statistics"]
git-tree-sha1 = "f92e1315dadf8c46561fb9396e525f7200cdc227"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.3.5"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "e47cd150dbe0443c3a3651bc5b9cbd5576ab75b7"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.52"

[[deps.PolygonOps]]
git-tree-sha1 = "77b3d3605fc1cd0b42d95eba87dfcd2bf67d5ff6"
uuid = "647866c9-e3ac-4575-94e7-e3d426903924"
version = "0.1.2"

[[deps.Polynomials]]
deps = ["LinearAlgebra", "RecipesBase"]
git-tree-sha1 = "3aa2bb4982e575acd7583f01531f241af077b163"
uuid = "f27b6e38-b328-58d1-80ce-0feddd5e7a45"
version = "3.2.13"

    [deps.Polynomials.extensions]
    PolynomialsChainRulesCoreExt = "ChainRulesCore"
    PolynomialsMakieCoreExt = "MakieCore"
    PolynomialsMutableArithmeticsExt = "MutableArithmetics"

    [deps.Polynomials.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    MakieCore = "20f20a25-4f0e-4fdf-b5d1-57303727442b"
    MutableArithmetics = "d8a4904e-b15c-11e9-3269-09a3773c0cb0"

[[deps.PositiveFactorizations]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "17275485f373e6673f7e7f97051f703ed5b15b20"
uuid = "85a6dd25-e78a-55b7-8502-1745935b8125"
version = "0.2.4"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "9673d39decc5feece56ef3940e5dafba15ba0f81"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.1.2"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "7eb1686b4f04b82f96ed7a4ea5890a4f0c7a09f1"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.4.0"

[[deps.Primes]]
deps = ["IntegerMathUtils"]
git-tree-sha1 = "4c9f306e5d6603ae203c2000dd460d81a5251489"
uuid = "27ebfcd6-29c5-5fa9-bf4b-fb8fc14df3ae"
version = "0.5.4"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.ProgressMeter]]
deps = ["Distributed", "Printf"]
git-tree-sha1 = "d7a7aef8f8f2d537104f170139553b14dfe39fe9"
uuid = "92933f4c-e287-5a05-a399-4b506db050ca"
version = "1.7.2"

[[deps.PythonCall]]
deps = ["CondaPkg", "Dates", "Libdl", "MacroTools", "Markdown", "Pkg", "REPL", "Requires", "Serialization", "Tables", "UnsafePointers"]
git-tree-sha1 = "0d15cb32f52654921169b4305dae8f66a0e345dc"
uuid = "6099a3de-0909-46bc-b1f4-468b9a2dfc0d"
version = "0.9.13"

[[deps.QOI]]
deps = ["ColorTypes", "FileIO", "FixedPointNumbers"]
git-tree-sha1 = "18e8f4d1426e965c7b532ddd260599e1510d26ce"
uuid = "4b34888f-f399-49d4-9bb3-47ed5cae4e65"
version = "1.0.0"

[[deps.QuadGK]]
deps = ["DataStructures", "LinearAlgebra"]
git-tree-sha1 = "6ec7ac8412e83d57e313393220879ede1740f9ee"
uuid = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
version = "2.8.2"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.RangeArrays]]
git-tree-sha1 = "b9039e93773ddcfc828f12aadf7115b4b4d225f5"
uuid = "b3c3ace0-ae52-54e7-9d0b-2c1406fd6b9d"
version = "0.3.2"

[[deps.Ratios]]
deps = ["Requires"]
git-tree-sha1 = "1342a47bf3260ee108163042310d26f2be5ec90b"
uuid = "c84ed2f1-dad5-54f0-aa8e-dbefe2724439"
version = "0.4.5"
weakdeps = ["FixedPointNumbers"]

    [deps.Ratios.extensions]
    RatiosFixedPointNumbersExt = "FixedPointNumbers"

[[deps.RecipesBase]]
deps = ["PrecompileTools"]
git-tree-sha1 = "5c3d09cc4f31f5fc6af001c250bf1278733100ff"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.3.4"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.RelocatableFolders]]
deps = ["SHA", "Scratch"]
git-tree-sha1 = "90bc7a7c96410424509e4263e277e43250c05691"
uuid = "05181044-ff0b-4ac5-8273-598c1e38db00"
version = "1.0.0"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.RingLists]]
deps = ["Random"]
git-tree-sha1 = "9712ebc42e91850f35272b48eb840e60c0270ec0"
uuid = "286e9d63-9694-5540-9e3c-4e6708fa07b2"
version = "0.2.7"

[[deps.Rmath]]
deps = ["Random", "Rmath_jll"]
git-tree-sha1 = "f65dcb5fa46aee0cf9ed6274ccbd597adc49aa7b"
uuid = "79098fc4-a85e-5d69-aa6a-4863f24498fa"
version = "0.7.1"

[[deps.Rmath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "6ed52fdd3382cf21947b15e8870ac0ddbff736da"
uuid = "f50d1b31-88e8-58de-be2c-1cc44531875f"
version = "0.4.0+0"

[[deps.RoundingEmulator]]
git-tree-sha1 = "40b9edad2e5287e05bd413a38f61a8ff55b9557b"
uuid = "5eaf0fd0-dfba-4ccb-bf02-d820a40db705"
version = "0.2.1"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.SIMD]]
deps = ["PrecompileTools"]
git-tree-sha1 = "0e270732477b9e551d884e6b07e23bb2ec947790"
uuid = "fdea26ae-647d-5447-a871-4b548cad5224"
version = "3.4.5"

[[deps.ScanByte]]
deps = ["Libdl", "SIMD"]
git-tree-sha1 = "d49e35f413186528f1d7cc675e67d0ed16fd7800"
uuid = "7b38b023-a4d7-4c5e-8d43-3f3097f304eb"
version = "0.4.0"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "30449ee12237627992a99d5e30ae63e4d78cd24a"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.2.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.SetRounding]]
git-tree-sha1 = "d7a25e439d07a17b7cdf97eecee504c50fedf5f6"
uuid = "3cc68bcd-71a2-5612-b932-767ffbe40ab0"
version = "0.2.1"

[[deps.Setfield]]
deps = ["ConstructionBase", "Future", "MacroTools", "StaticArraysCore"]
git-tree-sha1 = "e2cc6d8c88613c05e1defb55170bf5ff211fbeac"
uuid = "efcf1570-3423-57d1-acb7-fd33fddbac46"
version = "1.1.1"

[[deps.ShaderAbstractions]]
deps = ["ColorTypes", "FixedPointNumbers", "GeometryBasics", "LinearAlgebra", "Observables", "StaticArrays", "StructArrays", "Tables"]
git-tree-sha1 = "0d15c3e7b2003f4451714f08ffec2b77badc2dc4"
uuid = "65257c39-d410-5151-9873-9b3e5be5013e"
version = "0.3.0"

[[deps.SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"

[[deps.Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "91eddf657aca81df9ae6ceb20b959ae5653ad1de"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "1.0.3"

[[deps.SignedDistanceFields]]
deps = ["Random", "Statistics", "Test"]
git-tree-sha1 = "d263a08ec505853a5ff1c1ebde2070419e3f28e9"
uuid = "73760f76-fbc4-59ce-8f25-708e95d2df96"
version = "0.4.0"

[[deps.SimpleGraphs]]
deps = ["AbstractLattices", "Combinatorics", "DataStructures", "IterTools", "LightXML", "LinearAlgebra", "LinearAlgebraX", "Optim", "Primes", "Random", "RingLists", "SimplePartitions", "SimplePolynomials", "SimpleRandom", "SparseArrays", "Statistics"]
git-tree-sha1 = "b608903049d11cc557c45e03b3a53e9260579c19"
uuid = "55797a34-41de-5266-9ec1-32ac4eb504d3"
version = "0.8.4"

[[deps.SimplePartitions]]
deps = ["AbstractLattices", "DataStructures", "Permutations"]
git-tree-sha1 = "dcc02923a53f316ab97da8ef3136e80b4543dbf1"
uuid = "ec83eff0-a5b5-5643-ae32-5cbf6eedec9d"
version = "0.3.0"

[[deps.SimplePolynomials]]
deps = ["Mods", "Multisets", "Polynomials", "Primes"]
git-tree-sha1 = "d073c45302132b324ca653e1053966b4beacc2a5"
uuid = "cc47b68c-3164-5771-a705-2bc0097375a0"
version = "0.2.11"

[[deps.SimpleRandom]]
deps = ["Distributions", "LinearAlgebra", "Random"]
git-tree-sha1 = "3a6fb395e37afab81aeea85bae48a4db5cd7244a"
uuid = "a6525b86-64cd-54fa-8f65-62fc48bdc0e8"
version = "0.3.1"

[[deps.SimpleTraits]]
deps = ["InteractiveUtils", "MacroTools"]
git-tree-sha1 = "5d7e3f4e11935503d3ecaf7186eac40602e7d231"
uuid = "699a6c99-e7fa-54fc-8d76-47d257e15c1d"
version = "0.9.4"

[[deps.Sixel]]
deps = ["Dates", "FileIO", "ImageCore", "IndirectArrays", "OffsetArrays", "REPL", "libsixel_jll"]
git-tree-sha1 = "2da10356e31327c7096832eb9cd86307a50b1eb6"
uuid = "45858cf5-a6b0-47a3-bbea-62219f50df47"
version = "0.1.3"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "c60ec5c62180f27efea3ba2908480f8055e17cee"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.1.1"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.SpecialFunctions]]
deps = ["IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "7beb031cf8145577fbccacd94b8a8f4ce78428d3"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.3.0"
weakdeps = ["ChainRulesCore"]

    [deps.SpecialFunctions.extensions]
    SpecialFunctionsChainRulesCoreExt = "ChainRulesCore"

[[deps.StableHashTraits]]
deps = ["CRC32c", "Compat", "Dates", "SHA", "Tables", "TupleTools", "UUIDs"]
git-tree-sha1 = "0b8b801b8f03a329a4e86b44c5e8a7d7f4fe10a3"
uuid = "c5dd0088-6c3f-4803-b00e-f31a60c170fa"
version = "0.3.1"

[[deps.StackViews]]
deps = ["OffsetArrays"]
git-tree-sha1 = "46e589465204cd0c08b4bd97385e4fa79a0c770c"
uuid = "cae243ae-269e-4f55-b966-ac2d0dc13c15"
version = "0.1.1"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "Random", "StaticArraysCore"]
git-tree-sha1 = "9cabadf6e7cd2349b6cf49f1915ad2028d65e881"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.6.2"
weakdeps = ["Statistics"]

    [deps.StaticArrays.extensions]
    StaticArraysStatisticsExt = "Statistics"

[[deps.StaticArraysCore]]
git-tree-sha1 = "36b3d696ce6366023a0ea192b4cd442268995a0d"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.4.2"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.9.0"

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "45a7769a04a3cf80da1c1c7c60caf932e6f4c9f7"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.6.0"

[[deps.StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "75ebe04c5bed70b91614d684259b661c9e6274a4"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.34.0"

[[deps.StatsFuns]]
deps = ["HypergeometricFunctions", "IrrationalConstants", "LogExpFunctions", "Reexport", "Rmath", "SpecialFunctions"]
git-tree-sha1 = "f625d686d5a88bcd2b15cd81f18f98186fdc0c9a"
uuid = "4c63d2b9-4356-54db-8cca-17b64c39e42c"
version = "1.3.0"

    [deps.StatsFuns.extensions]
    StatsFunsChainRulesCoreExt = "ChainRulesCore"
    StatsFunsInverseFunctionsExt = "InverseFunctions"

    [deps.StatsFuns.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.StructArrays]]
deps = ["Adapt", "DataAPI", "GPUArraysCore", "StaticArraysCore", "Tables"]
git-tree-sha1 = "521a0e828e98bb69042fec1809c1b5a680eb7389"
uuid = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
version = "0.6.15"

[[deps.StructTypes]]
deps = ["Dates", "UUIDs"]
git-tree-sha1 = "ca4bccb03acf9faaf4137a9abc1881ed1841aa70"
uuid = "856f2bd8-1eba-4b0a-8007-ebc267875bd4"
version = "1.10.0"

[[deps.SuiteSparse]]
deps = ["Libdl", "LinearAlgebra", "Serialization", "SparseArrays"]
uuid = "4607b0f0-06f3-5cda-b6b1-a6196a1729e9"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "Pkg", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "5.10.1+6"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[deps.Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "LinearAlgebra", "OrderedCollections", "TableTraits", "Test"]
git-tree-sha1 = "1544b926975372da01227b382066ab70e574a3ec"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.10.1"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.TensorCore]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1feb45f88d133a655e001435632f019a9a1bcdb6"
uuid = "62fd8b95-f654-4bbd-a8a5-9c27f68ccd50"
version = "0.1.1"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.TiffImages]]
deps = ["ColorTypes", "DataStructures", "DocStringExtensions", "FileIO", "FixedPointNumbers", "IndirectArrays", "Inflate", "Mmap", "OffsetArrays", "PkgVersion", "ProgressMeter", "UUIDs"]
git-tree-sha1 = "8621f5c499a8aa4aa970b1ae381aae0ef1576966"
uuid = "731e570b-9d59-4bfa-96dc-6df516fadf69"
version = "0.6.4"

[[deps.TranscodingStreams]]
deps = ["Random", "Test"]
git-tree-sha1 = "9a6ae7ed916312b41236fcef7e0af564ef934769"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.9.13"

[[deps.Tricks]]
git-tree-sha1 = "aadb748be58b492045b4f56166b5188aa63ce549"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.7"

[[deps.TriplotBase]]
git-tree-sha1 = "4d4ed7f294cda19382ff7de4c137d24d16adc89b"
uuid = "981d1d27-644d-49a2-9326-4793e63143c3"
version = "0.1.0"

[[deps.TupleTools]]
git-tree-sha1 = "3c712976c47707ff893cf6ba4354aa14db1d8938"
uuid = "9d95972d-f1c8-5527-a6e0-b4b365fa01f6"
version = "1.3.0"

[[deps.URIs]]
git-tree-sha1 = "074f993b0ca030848b897beff716d93aca60f06a"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.4.2"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.UnPack]]
git-tree-sha1 = "387c1f73762231e86e0c9c5443ce3b4a0a9a0c2b"
uuid = "3a884ed6-31ef-47d7-9d2a-63182c4928ed"
version = "1.0.2"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.UnicodeFun]]
deps = ["REPL"]
git-tree-sha1 = "53915e50200959667e78a92a418594b428dffddf"
uuid = "1cfade01-22cf-5700-b092-accc4b62d6e1"
version = "0.4.1"

[[deps.UnsafePointers]]
git-tree-sha1 = "c81331b3b2e60a982be57c046ec91f599ede674a"
uuid = "e17b2a0c-0bdf-430a-bd0c-3a23cae4ff39"
version = "1.0.0"

[[deps.WoodburyMatrices]]
deps = ["LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "de67fa59e33ad156a590055375a30b23c40299d3"
uuid = "efce3f68-66dc-5838-9240-27a6d6f5f9b6"
version = "0.5.5"

[[deps.XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "93c41695bc1c08c46c5899f4fe06d6ead504bb73"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.10.3+0"

[[deps.XSLT_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgcrypt_jll", "Libgpg_error_jll", "Libiconv_jll", "Pkg", "XML2_jll", "Zlib_jll"]
git-tree-sha1 = "91844873c4085240b95e795f692c4cec4d805f8a"
uuid = "aed1982a-8fda-507f-9586-7b0439959a61"
version = "1.1.34+0"

[[deps.Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "afead5aba5aa507ad5a3bf01f58f82c8d1403495"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.8.6+0"

[[deps.Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "6035850dcc70518ca32f012e46015b9beeda49d8"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.11+0"

[[deps.Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "34d526d318358a859d7de23da945578e8e8727b7"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.4+0"

[[deps.Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "b7c0aa8c376b31e4852b360222848637f481f8c3"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.4+4"

[[deps.Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "19560f30fd49f4d4efbe7002a1037f8c43d43b96"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.10+4"

[[deps.Xorg_libpthread_stubs_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8fdda4c692503d44d04a0603d9ac0982054635f9"
uuid = "14d82f49-176c-5ed1-bb49-ad3f5cbd8c74"
version = "0.1.1+0"

[[deps.Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "XSLT_jll", "Xorg_libXau_jll", "Xorg_libXdmcp_jll", "Xorg_libpthread_stubs_jll"]
git-tree-sha1 = "b4bfde5d5b652e22b9c790ad00af08b6d042b97d"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.15.0+0"

[[deps.Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e92a1a012a10506618f10b7047e478403a046c77"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.5.0+0"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.13+0"

[[deps.isoband_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "51b5eeb3f98367157a7a12a1fb0aa5328946c03c"
uuid = "9a68df92-36a6-505f-a73e-abb412b6bfb4"
version = "0.2.3+0"

[[deps.libaom_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "3a2ea60308f0996d26f1e5354e10c24e9ef905d4"
uuid = "a4ae2306-e953-59d6-aa16-d00cac43593b"
version = "3.4.0+0"

[[deps.libass_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "5982a94fcba20f02f42ace44b9894ee2b140fe47"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.15.1+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.8.0+0"

[[deps.libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "daacc84a041563f965be61859a36e17c4e4fcd55"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.2+0"

[[deps.libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "94d180a6d2b5e55e447e2d27a29ed04fe79eb30c"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.38+0"

[[deps.libsixel_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Pkg", "libpng_jll"]
git-tree-sha1 = "d4f63314c8aa1e48cd22aa0c17ed76cd1ae48c3c"
uuid = "075b6546-f08a-558a-be8f-8157d0f608a5"
version = "1.10.3+0"

[[deps.libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll", "Pkg"]
git-tree-sha1 = "b910cb81ef3fe6e78bf6acee440bda86fd6ae00c"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.7+1"

[[deps.micromamba_jll]]
deps = ["Artifacts", "JLLWrappers", "LazyArtifacts", "Libdl"]
git-tree-sha1 = "087555b0405ed6adf526cef22b6931606b5af8ac"
uuid = "f8abcde7-e9b7-5caa-b8af-a437887ae8e4"
version = "1.4.1+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.48.0+0"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+0"

[[deps.x264_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fea590b89e6ec504593146bf8b988b2c00922b2"
uuid = "1270edf5-f2f9-52d2-97e9-ab00b5d0237a"
version = "2021.5.5+0"

[[deps.x265_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "ee567a171cce03570d77ad3a43e90218e38937a9"
uuid = "dfaa095f-4041-5dcd-9319-2fabd8486b76"
version = "3.5.0+0"
"""

# ╔═╡ Cell order:
# ╠═fe301fc0-d32b-40e7-8ab8-17d81b6bcf8a
# ╠═94e6fcb1-d73f-467b-97fa-82d536f1703b
# ╠═0b3fe3e4-ae8e-4a58-955f-8c396f2047a8
# ╠═908a8f10-6375-427a-97a5-cb2e5caf5cfb
# ╟─b61002b8-38ae-40a7-8344-a58ecb0ebffa
# ╟─13cdd42b-25c8-4f61-bc4e-1b5588aa9918
# ╟─5161556a-deae-48c9-bf6c-36a8aaf5882d
# ╟─9a39c2cb-d2e2-424d-a04d-ec52b47cd2dc
# ╟─53fbf184-231b-488e-aaa9-0b0698be5a9e
# ╟─5d7e3199-dc42-4870-904c-1906a99cfa97
# ╟─bef5c3ec-8f59-4f4b-8b6a-ec22f5c9b8eb
# ╟─4fadadf2-5b23-4038-92f2-8693918ca740
# ╟─b839496c-bc56-4504-ab53-ae2f257bb4ba
# ╟─4cacb360-625c-40e7-bdf9-437e700a1657
# ╠═e1f7df8a-57c7-499e-b605-4be6604fae30
# ╠═09a2bcab-1237-469e-ba0f-23cd9aad28d3
# ╠═72576547-e113-44aa-9b17-34bb5a575277
# ╟─a3f36b49-526a-47cf-a55c-2648d4a367f1
# ╟─68487535-e065-4778-94e7-db8344cc6237
# ╠═4bf65708-3747-4eba-b523-641223ff4d99
# ╠═47220365-3466-4fb3-b86e-f0428e0d9fa9
# ╠═9f3f770a-3ee5-4aec-9953-1b73641074c5
# ╠═544eef1f-d204-425b-85ff-e25884721f5f
# ╟─b6c7ac3d-f2aa-4c4f-b53d-75a8d5886839
# ╠═24cd17dc-678b-4c7d-ab4a-bd5c0eeb8e76
# ╠═0d5dceac-4230-4a2c-8783-1563e1738df4
# ╠═f8daa2ec-5e87-48a8-92ac-609996d302b4
# ╠═44ffdc2f-3ce5-4883-b1d5-555dc8a0dddd
# ╟─5487b403-64cb-46f5-b007-ba563f9fc218
# ╠═7afa8504-6f1e-4b5e-9d86-67227bf1bd7d
# ╠═5ef84561-4daa-4b5d-ba97-d26b521554f3
# ╟─75d86070-6ab0-488d-b53e-50f69cda76dc
# ╠═fa4d3c70-f60d-4639-a176-9a61897c7088
# ╠═419354e8-4975-43dc-a119-d3fadc02004a
# ╟─08bbdbdf-3778-415d-813b-f81922f64215
# ╟─3b55ae9d-855f-4dca-a293-ad235ed7c8ff
# ╟─acd8c7a1-7ba1-4819-a91b-711cd9a3ac3b
# ╟─f159263e-72a9-4417-9cf3-cfd9e18214ac
# ╟─8a9fa30c-49f8-4ae0-81a7-93ce3b3792fc
# ╟─25f81dd0-27e2-4613-be93-fb8a9b4d0c57
# ╟─c139e863-0723-4cec-8131-3eda35a50e7a
# ╠═9e603a65-87d1-48ba-8052-37c9cd2b3a43
# ╠═14dfbbe8-5d0a-4572-8d22-3b73bbd22f93
# ╠═9718a6c4-73ef-4fb7-a204-dd803fa5dc52
# ╠═24c4be9c-7c05-4925-a1c9-5a9aaaf6d878
# ╠═117b734a-0d5a-4a6f-b239-d6067f399ef5
# ╠═394f87f0-1046-4835-a9a6-3682fb8a1f03
# ╠═25877abf-5289-465e-8fab-7e147f1ca5b3
# ╠═c7e9ba5f-5ce0-423c-86cd-015efca7137f
# ╟─bc8c413f-8184-4d6d-a87b-a0ef3969c93a
# ╠═293e81fe-1744-40f6-bd49-c7476d281eb2
# ╠═7e47e8d2-65f3-4627-815c-d94212bac67b
# ╠═451ae9c6-5ef6-4651-99bb-44132e844a89
# ╠═5da3cbc8-6e1f-436a-a3b9-58c2992e2152
# ╠═da87a788-7538-43bd-abc6-510acb8999ff
# ╠═61120a25-bd56-4d22-8431-dc0734117db7
# ╠═1a722ba9-266f-4fc4-a521-c5038550b45e
# ╟─873551d7-aebb-41fc-8746-00c045d3a10f
# ╟─3ece13e7-18e9-4449-88d1-1d6e39b001ec
# ╟─3ece9a12-8071-409e-a313-c68135e20503
# ╟─c7fe4490-fd11-4c51-89ac-770d3adaadc5
# ╠═3446cb4b-3dad-4659-8e53-2ffab246e151
# ╠═fc2f149c-3aa2-4f8d-b1f3-59a69d073cf8
# ╟─946cb3e1-ba80-4b91-8afb-a9939ec6ef0b
# ╟─c13032f9-d9ac-4856-a9eb-537c8c102827
# ╟─77bfe4d8-c6ae-4832-8f00-a9e93bd1f540
# ╠═f1d8eaa1-9741-4988-8135-ac439c5f3005
# ╠═09c91686-5fab-4dfb-a906-3209abbbfa7d
# ╠═5a972777-a650-445f-b930-a05818559415
# ╠═3778e89a-b996-4fc2-87e1-b7a7af30cb56
# ╠═2de36b1d-175f-496b-b5aa-6bc8dce4c064
# ╠═37a2ffe4-49a3-46d6-87e5-dac31b9e82a8
# ╠═e5f68df7-144c-43b2-bef7-0df8905997ca
# ╠═3f5cabb9-546e-4cf3-8a79-2563d9d3ab61
# ╟─27984f40-0d25-4216-95cb-8f15e1badd97
# ╟─d898d7fd-b6f8-4ec4-bcd4-6fd6611a8e70
# ╟─47894665-ac20-48f4-aa43-f6b27e854042
# ╟─37f965aa-717d-43ec-adf5-8bc2790979b8
# ╠═185b0b11-7862-49c5-a28b-6b2d8c97e0b0
# ╠═0ebf5c1a-6287-4298-94f2-50361604bb0a
# ╠═a45b5225-1c71-48a0-8cc4-0886dfe85cc3
# ╠═34bc2428-06d4-404c-b0b8-ac5d978119a5
# ╠═16b9644c-aa6a-4665-a96d-ae4aa86d0b88
# ╠═541fb80a-63a6-4b66-8a03-69ae1e0f4553
# ╠═cfd88052-e979-435b-bb41-2dc26552aa4a
# ╠═160e03a0-545e-4d86-bd8f-5fe0180580c5
# ╠═5bc43ec7-9720-411f-900e-cd5f3af68824
# ╠═a3947c1e-3a5f-42fa-9ea3-7f0da19261c2
# ╟─0c950c0c-a209-421f-8d25-b691a482334f
# ╠═de382bf0-783d-487b-b3c0-83441bf7662a
# ╠═c61be673-d114-47fc-9bfa-0457cf16aacf
# ╠═8e2f5d10-a031-431a-a5a5-a475aa46533b
# ╠═d84d7137-cb72-426a-8866-a5156e5334de
# ╠═2bc5ce49-f52a-48d6-be64-a040052927af
# ╠═2e34d2ed-14c2-4746-8b3e-58510a5edf86
# ╠═13deb593-c204-4b77-808f-62a7e333418e
# ╠═7e81d4a3-2e75-4321-879a-df8259a93c1e
# ╠═5078c8e3-50ae-4312-ad4d-0e6416d8b360
# ╟─32a3470d-2629-4e74-94fe-5383cb4af5b9
# ╟─bfb3a382-f362-4b20-afb4-ae82af191a57
# ╟─c00c7da9-b548-46c1-b53a-970539cc9780
# ╟─088bd06b-5d6a-41c1-8a58-a4edbca9879e
# ╠═933d3db0-470e-4e6b-8175-94290774af0c
# ╠═710c5f83-871c-4eda-b10a-4f798c6b6fac
# ╠═94df1dd9-ee4c-4d23-bbb6-b5331f4196b4
# ╠═d85ad2a0-4cbc-45fb-a5dd-fbf0980181c0
# ╟─81e6ddce-0e8a-4d57-8e9c-f827dd613b9b
# ╠═b7151a26-b9c9-4329-a317-2131c35f3dbc
# ╠═fc6e6bbd-8744-4369-bd7c-29e8f52c8048
# ╟─25112cba-34af-4ad7-9af2-27bf6b560374
# ╟─baa68d3f-8cd9-4f28-a8c6-a957ec5b3fb6
# ╟─9efddcb3-6c01-494b-9512-e9d098868aae
# ╟─a81a50ff-5d1a-4e77-9567-46ff2be1be80
# ╠═528d1bb7-d56b-4f63-84a9-0a8035878c95
# ╠═f4010731-0793-4342-a65f-21dfe418d608
# ╠═cefd5ebb-041d-416d-8e05-390829e215f8
# ╠═a3f4fa24-314a-42d3-b932-ce17cbb766c8
# ╠═796ddfba-ef94-450c-9ce1-a87eeb9670a9
# ╠═75b74bde-ab81-4414-af83-f1868b59657a
# ╠═a1d4a670-ccb6-4bc5-8d35-4011d6003812
# ╠═c6b8eb02-d3f8-419f-b3dd-0c67e810254a
# ╟─512817e6-114d-4b80-b7a0-bc77c5c7c1ab
# ╟─2e283fd3-7cec-4926-96fa-b719eb288328
# ╠═760e58cc-d859-4ce5-8d6f-1f105ded77fe
# ╠═76087030-4287-41a8-ae10-805188dc7f65
# ╠═d5d3381c-bbcf-41e6-9cbd-b821126d4570
# ╠═b3ae3846-26a3-4da0-9ca1-23400a22e2fe
# ╠═effd8948-af6c-4040-ba06-3e31a784df8d
# ╟─b51c8800-e410-4481-b853-0763caf20baa
# ╟─a47e7c7a-bb5b-462a-87b3-69a6ddbe18b7
# ╟─54468bc7-3115-4582-af0c-4d66f1cb0d3e
# ╟─b1633978-2110-4bef-a704-2284b34bcadc
# ╟─dae97f2e-a2a5-4265-8af6-d61c2bbd78ac
# ╟─7b80df03-eeb4-4130-9b56-0a765a7d232b
# ╠═37203e2b-8c2c-4fa0-997a-a8f59da022cc
# ╠═c967a4b7-7c4d-41c0-b0f7-2b326c86ebde
# ╠═49106849-ba23-45b4-90c4-54a13aa95b8d
# ╠═fb10a998-87cb-476e-96d9-a232fa6fa5b8
# ╟─c607bfcc-fa55-46eb-bc5f-b47fb5f2f34b
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002