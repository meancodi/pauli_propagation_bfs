"Hermitian Pauli: i^count_ones(x & z) X^x Z^z, qubit 1 is the low bit."
struct PauliWord
    x::UInt64
    z::UInt64
end
function word(syms, sites=eachindex(syms))
    length(syms) == length(sites) || throw(ArgumentError("symbols/sites length mismatch"))
    allunique(sites) || throw(ArgumentError("duplicate sites"))
    x = z = UInt64(0)
    for (s, q) in zip(syms, sites)
        1 <= q <= 64 || throw(ArgumentError("sites must be in 1:64"))
        s in (:I, :X, :Y, :Z) || throw(ArgumentError("invalid Pauli $s"))
        b = UInt64(1) << (q-1)
        s in (:X, :Y) && (x |= b)
        s in (:Z, :Y) && (z |= b)
    end
    PauliWord(x,z)
end
bit(a, q) = !iszero(a & (UInt64(1) << (q-1)))
symbols(p::PauliWord, n) = [bit(p.x,q) ? (bit(p.z,q) ? :Y : :X) : (bit(p.z,q) ? :Z : :I) for q in 1:n]
commutes(a::PauliWord,b::PauliWord) = iseven(count_ones(a.x & b.z) + count_ones(a.z & b.x))
function product(a::PauliWord,b::PauliWord)
    p = PauliWord(a.x ⊻ b.x, a.z ⊻ b.z)
    e = mod(count_ones(a.x & a.z) + count_ones(b.x & b.z) - count_ones(p.x & p.z) + 2count_ones(a.z & b.x), 4)
    p, (1.0+0im, 0.0+1im, -1.0+0im, 0.0-1im)[e+1]
end
function validword(p, n)
    1 <= n <= 64 || throw(ArgumentError("nqubits must be in 1:64"))
    n == 64 || iszero((p.x | p.z) >> n) || throw(ArgumentError("Pauli exceeds circuit width"))
end
