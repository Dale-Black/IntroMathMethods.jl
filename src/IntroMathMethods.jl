module IntroMathMethods

using PythonCall

const sp = PythonCall.pynew()

function __init__()
    PythonCall.pycopy!(sp, pyimport("sympy"))
end

export sp

include("lagrangian.jl")


end