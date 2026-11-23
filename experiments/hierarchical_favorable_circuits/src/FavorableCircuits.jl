module FavorableCircuits

using Random, LinearAlgebra
const B=Main.BestFirstPauli

export echo_circuit, repeated_block_circuit, clifford_sandwich_circuit,
       trotter_circuit, periodic_circuit, inverse_gate, observables,
       supported_local_block, family_case

function inverse_gate(g)
    g isa B.Rotation && return B.Rotation(g.generator,-g.theta)
    g.kind==:H && return B.Clifford(:H,g.a)
    g.kind==:CNOT && return B.Clifford(:CNOT,g.a,g.b)
    # S† = S^3, represented using the supported Clifford gate set.
    [B.Clifford(:S,g.a),B.Clifford(:S,g.a),B.Clifford(:S,g.a)]
end
function inverse_circuit(gates)
    out=B.SupportedGate[]
    for g in reverse(gates)
        x=inverse_gate(g); x isa AbstractVector ? append!(out,x) : push!(out,x)
    end
    out
end
function supported_local_block(k::Int,ngates::Int,seed::Int;kind=:mixed)
    k>=2 || throw(ArgumentError("local block needs at least two qubits")); rng=MersenneTwister(seed); gs=B.SupportedGate[]
    for i in 1:ngates
        q=mod1(i,k)
        if kind==:clifford
            push!(gs,B.Clifford(rand(rng,(:H,:S)),q))
            iseven(i) && push!(gs,B.Clifford(:CNOT,q,mod1(q+1,k)))
        elseif kind==:nonclifford
            push!(gs,B.Rotation(B.word([rand(rng,(:X,:Y,:Z))],[q]),.08+.5rand(rng)))
            iseven(i) && push!(gs,B.Rotation(B.word([:Z,:Z],[q,mod1(q+1,k)]),.05+.3rand(rng)))
        else
            push!(gs,B.Clifford(rand(rng,(:H,:S)),q))
            push!(gs,B.Rotation(B.word([rand(rng,(:X,:Y,:Z))],[q]),rand(rng)*.8-.4))
            iseven(i) && push!(gs,B.Clifford(:CNOT,q,mod1(q+1,k)))
        end
    end
    gs
end
function echo_circuit(n::Int,vdepth::Int,seed::Int; perturb=:identity,theta=.1)
    v=supported_local_block(n,vdepth,seed;kind=:mixed); w=B.SupportedGate[]
    perturb==:single && push!(w,B.Rotation(B.word([:X],[1]),theta))
    perturb==:two && (push!(w,B.Rotation(B.word([:X],[1]),theta));push!(w,B.Rotation(B.word([:Z,:Z],[1,min(2,n)]),theta*.7)))
    perturb==:noncommuting && (push!(w,B.Rotation(B.word([:X],[1]),theta));push!(w,B.Rotation(B.word([:Y],[min(2,n)]),theta*.8)))
    B.Circuit(n,vcat(v,w,inverse_circuit(v)))
end
function repeated_block_circuit(k::Int,ngates::Int,reps::Int,seed::Int;kind=:mixed,perturbed=false)
    b=supported_local_block(k,ngates,seed;kind=kind); rng=MersenneTwister(seed+1000); gs=B.SupportedGate[]
    for _ in 1:reps
        if perturbed
            append!(gs,[g isa B.Rotation ? B.Rotation(g.generator,g.theta+.001rand(rng)) : g for g in b])
        else; append!(gs,b) end
    end
    B.Circuit(k,gs)
end
function clifford_block(n,depth,seed)
    rng=MersenneTwister(seed); gs=B.SupportedGate[]
    for i in 1:depth
        q=mod1(i,n); push!(gs,B.Clifford(rand(rng,(:H,:S)),q));
        iseven(i) && push!(gs,B.Clifford(:CNOT,q,mod1(q+1,n)))
    end
    gs
end
function clifford_sandwich_circuit(n::Int,cdepth::Int,nc::Int,seed::Int;theta=.2)
    c=clifford_block(n,cdepth,seed); r=B.SupportedGate[]
    for i in 1:nc
        q=mod1(seed+i,n); push!(r,B.Rotation(B.word([isodd(i) ? :X : :Z],[q]),theta/(1+.2i)))
    end
    B.Circuit(n,vcat(c,r,inverse_circuit(c)))
end
function ising_block(n,dt,perturb)
    gs=B.SupportedGate[]
    for q in 1:n-1; push!(gs,B.Rotation(B.word([:Z,:Z],[q,q+1]),.31dt)) end
    for q in 1:n
        push!(gs,B.Rotation(B.word([:X],[q]),.47dt)); push!(gs,B.Rotation(B.word([:Z],[q]),.23dt))
    end
    perturb!=0 && push!(gs,B.Rotation(B.word([:Y],[cld(n,2)]),perturb))
    gs
end
function trotter_circuit(n::Int,reps::Int,dt::Real,perturb::Real;vary=false,seed=11)
    gs=B.SupportedGate[]; rng=MersenneTwister(seed)
    for i in 1:reps
        b=ising_block(n,dt,perturb)
        vary ? append!(gs,[g isa B.Rotation ? B.Rotation(g.generator,g.theta*(1+.02rand(rng))) : g for g in b]) : append!(gs,b)
    end
    B.Circuit(n,gs)
end
function periodic_circuit(n::Int,reps::Int,perturb::Real;order=4,seed=11)
    gs=B.SupportedGate[]; q=1
    for _ in 1:reps
        push!(gs,B.Rotation(B.word([:X],[q]),2π/order))
        n>=2 && push!(gs,B.Rotation(B.word([:Z,:Z],[1,2]),perturb))
    end
    B.Circuit(n,gs)
end
observables(n)=(X=B.word([:X],[cld(n,2)]),Y=B.word([:Y],[cld(n,2)]),Z=B.word([:Z],[cld(n,2)]),ZZ=B.word([:Z,:Z],[cld(n,2),cld(n,2)+1]),mixed=n>=3 ? B.word([:X,:Y,:Z],[1,2,3]) : B.word([:X,:Y],[1,2]))
function family_case(family::Symbol; n=4, depth=4, reps=4, seed=11, theta=.1, perturb=:identity, kind=:mixed, vary=false)
    if family==:echo; c=echo_circuit(n,depth,seed;perturb=perturb,theta=theta)
    elseif family==:repeated; c=repeated_block_circuit(n,depth,reps,seed;kind=kind,perturbed=vary)
    elseif family==:clifford_sandwich; c=clifford_sandwich_circuit(n,depth,reps,seed;theta=theta)
    elseif family==:trotter; c=trotter_circuit(n,reps,theta,perturb;vary=vary,seed=seed)
    elseif family==:periodic; c=periodic_circuit(n,reps,perturb;seed=seed)
    else; throw(ArgumentError("unknown family $family")) end
    c
end
end
