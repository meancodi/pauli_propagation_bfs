module BestFirstPauli
using Random, LinearAlgebra
import PauliPropagation as PP
export PauliWord, word, symbols, commutes, product, Rotation, Clifford, Circuit,
       bestfirst, depthfirst, make_circuit, observable, upstream, dense_expectation
include("encoding.jl")
include("gates.jl")
include("frontier.jl")
include("search.jl")
include("dfs.jl")
include("circuits.jl")
include("baselines.jl")
end
