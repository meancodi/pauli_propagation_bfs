abstract type LocalGate end
struct Rotation <: LocalGate
    generator::PauliWord
    theta::Float64
    sine::Float64
    cosine::Float64
    function Rotation(p::PauliWord, theta::Real)
        t = Float64(theta)
        isfinite(t) || throw(ArgumentError("angle must be finite"))
        s,c = sincos(t)
        new(p,t,s,c)
    end
end
struct Clifford <: LocalGate
    kind::Symbol
    a::Int
    b::Int
    function Clifford(kind::Symbol,a::Int,b::Int=0)
        kind in (:H,:S,:CNOT) || throw(ArgumentError("unsupported Clifford $kind"))
        1 <= a <= 64 || throw(ArgumentError("invalid site"))
        (kind == :CNOT ? (1 <= b <= 64 && a != b) : b == 0) || throw(ArgumentError("invalid Clifford sites"))
        new(kind,a,b)
    end
end
const SupportedGate = Union{Rotation,Clifford}
struct Circuit
    nqubits::Int
    gates::Vector{SupportedGate}
    function Circuit(n::Int,gates)
        1 <= n <= 64 || throw(ArgumentError("nqubits must be in 1:64"))
        for g in gates
            g isa SupportedGate || throw(ArgumentError("unsupported gate $(typeof(g))"))
            g isa Rotation ? validword(g.generator,n) : (max(g.a,g.b) <= n || throw(ArgumentError("gate exceeds circuit width")))
        end
        new(n,SupportedGate[g for g in gates])
    end
end
function expand(g::Rotation,p::PauliWord,c::Float64)
    commutes(g.generator,p) && return ((p,c),(p,0.0))
    q, phase = product(g.generator,p)
    # U' P U = cos(theta) P + i sin(theta) QP.
    ((p,c*g.cosine),(q,c*g.sine*real(im*phase)))
end
function expand(g::Clifford,p::PauliWord,c::Float64)
    x,z = p.x,p.z
    xa,za = bit(x,g.a),bit(z,g.a)
    ba = UInt64(1) << (g.a-1)
    if g.kind == :H
        xa && za && (c = -c)
        xa != za && (x ⊻= ba; z ⊻= ba)
    elseif g.kind == :S
        xa && !za && (c = -c)
        xa && (z ⊻= ba)
    else
        xb,zb = bit(x,g.b),bit(z,g.b)
        xa && zb && !(xb ⊻ za) && (c = -c)
        xa && (x ⊻= UInt64(1) << (g.b-1))
        zb && (z ⊻= ba)
    end
    ((PauliWord(x,z),c),(p,0.0))
end
