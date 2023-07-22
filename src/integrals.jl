"""
    integrate_definite(f, x, x0, x1)

Simple wrapper around `integrate()` that allows for definite integrals to be computed.
"""
function integrate_definite(f, x, x0, x1)
	F = SymbolicNumericIntegration.integrate(f, x)[1]
	F = substitute(F, x => x1) - substitute(F, x => x0)
end

export integrate_definite