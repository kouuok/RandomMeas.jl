# crm_ed_pbc.jl の結果を §2h 用に集計する
#   julia --project=Hubbard_MPS_Env_v2 crm_ed_pbc_analyze.jl
using DelimitedFiles, Statistics, Printf

raw, hdr = readdlm(joinpath(@__DIR__, "crm_ed_pbc_results.tsv"), '\t'; header=true)
c = Dict(String(h) => j for (j, h) in enumerate(vec(hdr)))
L   = Int.(raw[:, c["L"]]);        bc  = String.(raw[:, c["bc"]])
rho = String.(raw[:, c["rho"]]);   site= Int.(raw[:, c["site"]])
pr  = String.(raw[:, c["prior"]]); ob  = String.(raw[:, c["observable"]])
Δ   = Float64.(raw[:, c["Delta"]]);G   = Float64.(raw[:, c["G"]])
Gmx = Float64.(raw[:, c["G_max"]])
sd  = Int.(raw[:, c["seed"]]);      Ghi = Float64.(raw[:, c["G_hi"]])
const S0 = minimum(sd)             # 表は代表シードで作る
n = length(L)
Ls = sort(unique(L)); bcs = ["OBC", "PBC"]

meta, mh = readdlm(joinpath(@__DIR__, "crm_ed_pbc_meta.tsv"), '\t'; header=true)
mc = Dict(String(h) => j for (j, h) in enumerate(vec(mh)))

println("="^92)
println("(1) DMRG は PBC でもこの規模で厳密か")
@printf("%-4s %-5s %18s %12s %12s %12s %10s\n", "L","bc","E0(ED)","|ΔE|","1-F","局所量max差","χ_max")
for r in 1:size(meta,1)
    f = meta[r, mc["fid"]]
    @printf("%-4d %-5s %18.10f %12.2e %12s %12.2e %10d\n",
            meta[r,mc["L"]], meta[r,mc["bc"]], meta[r,mc["E_ed"]], meta[r,mc["dE"]],
            f isa Number && !isnan(f) ? @sprintf("%.2e", 1-f) : "—",
            meta[r,mc["max_obs_diff"]], meta[r,mc["maxlinkdim"]])
end

println("\n", "="^92)
println("(2) ED を参照にしたときと DMRG を参照にしたときで利得は変わるか")
@printf("%-4s %-5s %4s %8s %14s %12s %10s %12s\n",
        "L","bc","seed","点数","max|ΔG|","max|ΔG|/G","ずれた点数","G の1σ幅")
key = Dict{Tuple,Int}()
for i in 1:n; key[(L[i],bc[i],sd[i],rho[i],site[i],pr[i],ob[i])] = i; end
for l in Ls, b in bcs, s0 in sort(unique(sd))
    gd = Float64[]; rel = Float64[]; err = Float64[]
    for i in 1:n
        (L[i]==l && bc[i]==b && sd[i]==s0 && rho[i]=="ED") || continue
        j = get(key, (l,b,s0,"DMRG",site[i],pr[i],ob[i]), 0); j == 0 && continue
        push!(gd, abs(G[i]-G[j])); push!(rel, abs(G[i]-G[j])/max(G[i],1e-12))
        push!(err, (Ghi[i]-G[i]))
    end
    isempty(gd) && continue
    @printf("%-4d %-5s %4d %8d %14.3e %12.2e %10d %12.3f\n", l, b, s0, length(gd),
            maximum(gd), maximum(rel), count(>(0.0), gd), median(err))
end

println("\n", "="^92)
println("(3) 境界条件で利得はどう変わるか (ρ=ED、全サイト中央値)")
obs_order = ["n","Sz","DoubleOcc","ZZ onsite","ZZ up-up r=1","SzSz r=1","hop up r=1"]
for l in Ls
    @printf("\n-- L=%d --\n", l)
    @printf("%-14s %-5s %8s %8s %8s %8s %8s %8s | %8s\n",
            "observable","bc","χ=2","χ=4","χ=8","χ=16","χ=32","UHFsym","天井")
    for o in obs_order
        for b in bcs
            f(lb) = (h=[G[i] for i in 1:n if L[i]==l&&bc[i]==b&&rho[i]=="ED"&&sd[i]==S0&&pr[i]==lb&&ob[i]==o];
                     isempty(h) ? NaN : median(h))
            gm = (h=[Gmx[i] for i in 1:n if L[i]==l&&bc[i]==b&&rho[i]=="ED"&&sd[i]==S0&&ob[i]==o];
                  isempty(h) ? NaN : median(h))
            isnan(gm) && continue
            @printf("%-14s %-5s %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f | %8.2f\n",
                    b=="OBC" ? o : "", b, f("chi2"),f("chi4"),f("chi8"),
                    f("chi16"),f("chi32"),f("UHFsym"), gm)
        end
    end
end

println("\n", "="^92)
println("(4) 損をする点の数 (ブートストラップCI上端<0.9、ρ=ED、全観測量×全サイト)")
plist = ["chi2","chi4","chi8","chi16","chi32","UHF","UHFsym"]
@printf("%-4s %-5s %8s", "L","bc","総点数")
for p in plist; @printf("%10s", p); end; println()
for l in Ls, b in bcs
    tot = count(i -> L[i]==l&&bc[i]==b&&rho[i]=="ED"&&sd[i]==S0&&pr[i]=="chi2", 1:n)
    tot == 0 && continue
    @printf("%-4d %-5s %8d", l, b, tot)
    for p in plist
        nb = count(i -> L[i]==l&&bc[i]==b&&rho[i]=="ED"&&sd[i]==S0&&pr[i]==p&&Ghi[i]<0.9, 1:n)
        @printf("%10s", "$nb/$tot")
    end
    println()
end

println("\n", "="^92)
println("(5) PBC でのサイト依存性 (並進対称なら消えるはず。ρ=ED、χ_p=8)")
@printf("%-4s %-5s %-14s %10s %10s %10s\n", "L","bc","observable","min G","max G","max/min")
for l in Ls, b in bcs, o in obs_order
    h = [G[i] for i in 1:n if L[i]==l&&bc[i]==b&&rho[i]=="ED"&&sd[i]==S0&&pr[i]=="chi8"&&ob[i]==o]
    (isempty(h) || median(h) < 1.01) && continue
    @printf("%-4d %-5s %-14s %10.2f %10.2f %10.2f\n", l, b, o,
            minimum(h), maximum(h), maximum(h)/minimum(h))
end

println("\n", "="^92)
println("(6) 同じ χ_p で PBC は OBC よりどれだけ利得を失うか (ρ=ED、中央値の比)")
@printf("%-4s %-14s %10s %10s %10s %10s\n", "L","observable","χ=4","χ=8","χ=16","χ=32")
for l in Ls, o in obs_order
    g(b,lb) = (h=[G[i] for i in 1:n if L[i]==l&&bc[i]==b&&rho[i]=="ED"&&sd[i]==S0&&pr[i]==lb&&ob[i]==o];
               isempty(h) ? NaN : median(h))
    g("OBC","chi8") < 1.01 && continue
    @printf("%-4d %-14s %10.2f %10.2f %10.2f %10.2f\n", l, o,
            g("PBC","chi4")/g("OBC","chi4"), g("PBC","chi8")/g("OBC","chi8"),
            g("PBC","chi16")/g("OBC","chi16"), g("PBC","chi32")/g("OBC","chi32"))
end
