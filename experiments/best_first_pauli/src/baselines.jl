function upstream_inputs(c::Circuit,p::PauliWord)
    gates = PP.Gate[]
    theta = Float64[]
    for g in c.gates
        if g isa Rotation
            syms = symbols(g.generator,c.nqubits)
            inds = findall(!=(:I),syms)
            # Identity rotations are a global phase and may be omitted.
            isempty(inds) && continue
            push!(gates,PP.PauliRotation(syms[inds],inds))
            push!(theta,g.theta)
        else
            push!(gates,PP.CliffordGate(g.kind,g.kind == :CNOT ? [g.a,g.b] : [g.a]))
        end
    end
    obs = PP.PauliString(c.nqubits,symbols(p,c.nqubits),collect(1:c.nqubits),1.0)
    gates,theta,obs
end
function upstream(c::Circuit,p::PauliWord,backend::Symbol=:PauliSum)
    gates,theta,obs = upstream_inputs(c,p)
    initial = backend == :PauliSum ? PP.PauliSum(obs) : backend == :VectorPauliSum ? PP.VectorPauliSum(obs) : throw(ArgumentError("unknown backend"))
    out = PP.propagate(gates,initial,theta; min_abs_coeff=0.0,max_weight=Inf,thread=false)
    (;expectation=PP.overlapwithzero(out),unique_terms=length(out))
end

# Independent reference: explicit 2x2 matrices and state-vector gate actions.
# No expansion kernel or Pauli product/sign formula is used here.
const PAULI_MATRICES = Dict(:I=>ComplexF64[1 0;0 1],:X=>ComplexF64[0 1;1 0],
    :Y=>ComplexF64[0 -im;im 0],:Z=>ComplexF64[1 0;0 -1])
function dense_pauli(p::PauliWord,n)
    reduce(kron,[PAULI_MATRICES[s] for s in reverse(symbols(p,n))])
end
function dense_gate(g::Rotation,n)
    cos(g.theta/2)*Matrix{ComplexF64}(I,2^n,2^n)-im*sin(g.theta/2)*dense_pauli(g.generator,n)
end
function dense_gate(g::Clifford,n)
    if g.kind != :CNOT
        localmat = g.kind == :H ? ComplexF64[1 1;1 -1]/sqrt(2) : ComplexF64[1 0;0 im]
        return reduce(kron,[q == g.a ? localmat : PAULI_MATRICES[:I] for q in n:-1:1])
    end
    u = zeros(ComplexF64,2^n,2^n)
    for j in 0:2^n-1
        k = (j >> (g.a-1)) & 1 == 1 ? j ⊻ (1 << (g.b-1)) : j
        u[k+1,j+1] = 1
    end
    u
end
function dense_expectation(c::Circuit,p::PauliWord)
    c.nqubits <= 8 || throw(ArgumentError("dense reference limited to 8 qubits"))
    psi = zeros(ComplexF64,2^c.nqubits); psi[1] = 1
    for g in c.gates
        psi = dense_gate(g,c.nqubits)*psi
    end
    real(dot(psi,dense_pauli(p,c.nqubits)*psi))
end
