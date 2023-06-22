# Norma
A Julia prototype for testing algorithms and ideas for coupling and multiphysics mainly in solid mechanics and heat conduction.

To install within the package manager (press `]` in the Julia REPL):

    ```
    pkg> add https://github.com/lxmota/norma
    ```

On MacOS, it is necessary to ignore package hashes for the dependence on Exodus.jl:

    ```
    ENV["JULIA_PKG_IGNORE_HASHES"] = 1
    ```
