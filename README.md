# Norma
A Julia prototype for testing algorithms and ideas for coupling and multiphysics mainly in solid mechanics and heat conduction.

To install within the package manager (press `]` in the Julia REPL):

    pkg> activate .
    pkg> registry update
    pkg> update
    pkg> instantiate
 
Then press `delete` to exit the package manager.

On MacOS, it is necessary to ignore package hashes for the dependence on Exodus.jl:

    ENV["JULIA_PKG_IGNORE_HASHES"] = 1

To run the code, assuming that Julia is in the executable path:

    julia --project=@. /path/to/src/norma.jl input.yaml
