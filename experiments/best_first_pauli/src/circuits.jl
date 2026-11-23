"Open-chain circuits, gates stored in Schrodinger order; seeded parameters."
function make_circuit(n::Int,family::Symbol,layers::Int,seed::Int)
    n >= 2 && layers >= 0 || throw(ArgumentError("invalid circuit size/depth"))
    rng = MersenneTwister(seed)
    gates = SupportedGate[]
    for l in 1:layers
        if family == :ising
            # First-order tilted Ising step: ZZ bonds, then X field, then Z field.
            for q in 1:n-1
                push!(gates,Rotation(word([:Z,:Z],[q,q+1]),0.31+0.04rand(rng)))
            end
            for q in 1:n
                push!(gates,Rotation(word([:X],[q]),0.47+0.04rand(rng)))
                push!(gates,Rotation(word([:Z],[q]),0.23+0.04rand(rng)))
            end
        elseif family == :random
            for q in 1:n
                push!(gates,Clifford(rand(rng,(:H,:S)),q))
                push!(gates,Rotation(word([rand(rng,(:X,:Y,:Z))],[q]),rand(rng)*1.2-0.6))
            end
            for q in (isodd(l) ? 1 : 2):2:n-1
                push!(gates,Clifford(:CNOT,q,q+1))
                push!(gates,Rotation(word([rand(rng,(:X,:Y,:Z)),rand(rng,(:X,:Y,:Z))],[q,q+1]),rand(rng)*1.2-0.6))
            end
        else
            throw(ArgumentError("unknown circuit family $family"))
        end
    end
    Circuit(n,gates)
end
observable(n,kind::Symbol) = kind == :Z_center ? word([:Z],[cld(n,2)]) :
    kind == :ZZ_center ? word([:Z,:Z],[cld(n,2),cld(n,2)+1]) : throw(ArgumentError("unknown observable"))
