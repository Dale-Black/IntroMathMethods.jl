"""
    taylor_series(f, x, n)

Compute the Taylor series expansion of a symbolic function `f` up to the `n`-th term, about the origin (0).

#### Arguments
- `f`: A symbolic expression representing a function. For example, if you have defined a symbolic variable `x` with `@syms x`, then `f` could be `sin(x)`.
- `x`: A symbolic variable representing the point about which the Taylor series is to be expanded. This variable is usually defined using the `@syms` macro.
- `n`: A non-negative integer representing the degree of the Taylor series expansion. The function will compute the Taylor series up to the `x^n` term.

#### Returns
- The function returns a symbolic expression representing the Taylor series expansion of `f` up to the `x^n` term.

#### Examples
```julia
@syms x
f = sin(x)
taylor_series(f, x, 5)
```
"""
function taylor_series(f, x, n::Int)
	Dx = Differential(x)
	
    f_ts = substitute(f, x => 0)
    for i in 1:n
		Dxi = Dx^i
		df = expand_derivatives(Dxi(f))
		df_0 = substitute(df, x => 0)
        f_ts += (df_0 * x^i) / factorial(i)
    end
    return f_ts
end


"""
    taylor_series_exp(x::Num, n::Int)

Compute the Taylor series expansion of `exp(x)` up to the `n`th term about the origin (0).

#### Arguments
- `x::Num` : The point around which to expand the function.
- `n::Int` : The number of terms in the Taylor series expansion.

#### Returns
- `::Number` : The calculated value of the Taylor series expansion of `exp(x)` up to the `n`th term.

#### Example
```julia
using Symbolics

@variables x
taylor_series_exp(x, 5)  

# Returns
julia > 1 + x + x^2/2 + x^3/6 + x^4/24 + x^5/120
```
"""
function taylor_series_exp(x::Num, n::Int)
    sum_expr = 0
    for i in 0:n
        sum_expr += x^i // factorial(i)
    end
    return sum_expr
end

export taylor_series, taylor_series_exp

