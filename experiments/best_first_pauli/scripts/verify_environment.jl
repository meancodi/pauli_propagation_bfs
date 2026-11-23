using Pkg, TOML, JSON3, Dates
import PauliPropagation
root = normpath(joinpath(@__DIR__,".."))
manifest = TOML.parsefile(joinpath(root,"Manifest.toml"))
entry = only(manifest["deps"]["PauliPropagation"])
source = dirname(dirname(pathof(PauliPropagation)))
actual = bytes2hex(Pkg.GitTools.tree_hash(source))
expected = entry["git-tree-sha1"]
actual == expected || error("Installed upstream source differs from the registered manifest tree: $actual != $expected")
mkpath(joinpath(root,"results","raw"))
write(joinpath(root,"results","raw","upstream_verification.json"),JSON3.write((;
    package="PauliPropagation",version=entry["version"],source_tree=actual,manifest_tree=expected,
    unmodified=true,source_path=source,official_repository="https://github.com/SparqleSim/PauliPropagation.jl",
    timestamp=string(now(UTC))*"Z")))
println("Verified unmodified PauliPropagation $(entry["version"]) tree $actual")
