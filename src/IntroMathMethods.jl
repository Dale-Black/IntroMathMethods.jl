module IntroMathMethods

using Symbolics

export taylor_series_exp


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

end
