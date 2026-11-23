using LinearAlgebra, Dates
include(joinpath(@__DIR__,"run_v2.jl"))

function local_index(p,sites)
    code=0; stride=1
    for q in sites
        x=B.bit(p.x,q); z=B.bit(p.z,q)
        code += ((x ? 1 : 0)+(z ? 2 : 0))*stride
        stride *= 4
    end
    code+1
end
function schmidt_rank(d,n)
    left=1:div(n,2); right=div(n,2)+1:n
    mat=zeros(Float64,4^length(left),4^length(right))
    for (p,v) in d
        mat[local_index(p,left),local_index(p,right)]+=v
    end
    s=svdvals(mat)
    isempty(s) || s[1]==0 ? 0 : count(x->x>1e-10*s[1],s)
end
function main_rank()
    root=joinpath(V2ROOT,"raw"); mkpath(root)
    path=joinpath(root,"v2_rank_"*Dates.format(now(),"yyyymmdd_HHMMSS")*".csv")
    open(path,"w") do io
        writecsv(io,("family","n","layers","observable","seed","fingerprint","cut","schmidt_rank","s2_peak_components","s3_root_sum_arity","s3_reuse_events"))
        for (family,n,layers,kind,seedidx) in ((:pairs,4,1,:factorized,1),(:pairs,8,1,:factorized,1),(:cross,8,1,:factorized,1),(:cross,8,2,:zz,1))
            seed=100_000+FAMILY_CODE[family]*10_000+n*100+layers*10+seedidx
            c=v2circuit(n,family,layers,seed,layers); terms=v2observable(n,kind)
            ref=sequential(c,terms;limits=REFLIMITS)
            s2=spatial(c,terms;leaf_size=2,limits=V2LIMITS)
            s3=recursive_spatial(c,terms;leaf_size=2,limits=V2LIMITS)
            rootarity=s3.expr isa SpatialPauli.RSum ? length(s3.expr.terms) : s3.expr===nothing ? 0 : 1
            writecsv(io,(family,n,layers,kind,seed,fingerprint(c,terms),div(n,2),schmidt_rank(ref,n),s2.peak_components,rootarity,s3.reuse_events))
        end
    end
    println(path)
end
main_rank()
