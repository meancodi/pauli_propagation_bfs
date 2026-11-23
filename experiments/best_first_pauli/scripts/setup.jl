using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
try
    if !isfile(joinpath(@__DIR__, "..", "Project.toml"))
        Pkg.add(["PauliPropagation", "CSV", "DataFrames", "JSON3"])
    else
        Pkg.instantiate()
    end
    Pkg.precompile()
    include(joinpath(@__DIR__,"verify_environment.jl"))
catch err
    @error "Setup failed. Check network access to Julia registries/package servers and writable JULIA_DEPOT_PATH; rerun this script." exception=(err, catch_backtrace())
    exit(1)
end
