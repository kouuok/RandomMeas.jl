# 実務では真の密度は使えない。UHF が自分で予言する密度プロファイルでも
# 「どこで自分が破綻するか」を予測できるかを確かめる。
using DelimitedFiles, Statistics, Printf, LinearAlgebra, Random
const W=4
include("crm_2d_site.jl_uhf_only")
for (f, U, nel, lab) in (("crm_2d_site_U8p0_h4.tsv", 8.0, 28, "U=8 δ=1/8"),
                         ("crm_2d_site_U4p0_h4.tsv", 4.0, 28, "U=4 δ=1/8"))
    uhf = solve_uhf_best(LX, W, 1.0, U; Nup=(nel+1)÷2, Ndn=nel÷2)
    ndens = uhf.nup .+ uhf.ndn
    raw,hdr = readdlm(f,'\t';header=true)
    c=Dict(String(h)=>j for (j,h) in enumerate(vec(hdr)))
    x=Int.(raw[:,c["x"]]); y=Int.(raw[:,c["y"]]); ob=String.(raw[:,c["observable"]])
    pr=String.(raw[:,c["prior"]]); re=Float64.(raw[:,c["relerr"]]); dn=Float64.(raw[:,c["dens"]])
    println("\n=== $lab ===")
    @printf("%-14s %18s %18s %14s\n","observable","corr(lnε, |n_UHF-1|)","corr(lnε, |n_true-1|)","corr(n_UHF,n_true)")
    for o in ["ZZ onsite","DoubleOcc","SzSz nb"]
        idx=findall(i-> ob[i]==o && pr[i]=="UHF" && !isnan(re[i]), 1:length(x))
        length(idx)<8 && continue
        nu = [abs(ndens[(x[i]-1)*W + y[i]] - 1) for i in idx]
        nt = abs.(dn[idx] .- 1); e = log.(re[idx])
        @printf("%-14s %18.3f %18.3f %14.3f\n", o, cor(e,nu), cor(e,nt), cor(nu,nt))
    end
end
