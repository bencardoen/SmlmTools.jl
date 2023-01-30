# using Pkg;
# using Logging;
# @info "Initiating build"
# ## We want the Conda local Python env, anything else is out of control
# ENV["PYTHON"] = ""
# # Conda and PyCall are dependencies, but we need to make sure they get prebuilt first.
# # We're in our own env, so explicitly adding them now does not harm.
# Pkg.add("Conda")
# Pkg.add("PyCall")
# ## --> Initiates an PyConda env local to us
# Pkg.build("PyCall")
# # Precompile
# using PyCall
# using Conda
# ## Add the two packages we need
# Conda.add("gcc=12.1.0"; channel="conda-forge")
# Conda.add("kneed"; channel="conda-forge")
# Conda.add("scikit-image")
# Conda.add("scipy=1.8.0")
# PyCall.pyimport("kneed");
# PyCall.pyimport("skimage");
# @info "Success!"

using Pkg;
using Logging;
@info "Initiating build"
## We want the Conda local Python env, anything else is out of control
install_p = false
if !haskey(ENV, "PYTHON")
    install_p = true
    ENV["PYTHON"] = ""
end
# ENV["R_HOME"] = "*"
# ENV["LD_PRELOAD"] = joinpath(ENV["HOME"], ".julia/conda/3/lib/libstdc++.so.6.0.30")
# Conda and PyCall are dependencies, but we need to make sure they get prebuilt first.
# We're in our own env, so explicitly adding them now does not harm.
Pkg.add("Conda")
Pkg.add("PyCall")
## --> Initiates an PyConda env local to us
Pkg.build("PyCall")
# Precompile
using PyCall
using Logging
using Conda
## Add the two packages we need
# Conda.pip_interop(true)
# Conda.add("gcc=12.1.0"; channel="conda-forge")
# Pin this version, to avoid clashes with libgcc.34
# Conda.add("scipy=1.8.0"))
try
    PyCall.pyimport("smlmvis");
    PyCall.pyimport("meshio");
catch e
    @warn "Failed import $e -- installing"
    println("Failed import $e -- installing")
    Conda.pip_interop(true)
    Conda.add("gcc=12.1.0"; channel="conda-forge")
    #Pin this version, to avoid clashes with libgcc.34
    Conda.add("scipy=1.8.0")
    Conda.pip("install", "smlmvis")
end
PyCall.pyimport("smlmvis");
@info "Success!"
