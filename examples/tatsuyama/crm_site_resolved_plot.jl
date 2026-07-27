# サイト分解CRM利得の解析・作図 (crm_site_resolved.jl の出力を読む)
#   julia --project=Hubbard_MPS_Env_v2 crm_site_resolved_plot.jl
using DelimitedFiles, Statistics, Printf, Plots

const FILE = joinpath(@__DIR__, "crm_site_resolved_results.tsv")
raw, hdr = readdlm(FILE, '\t'; header=true)
col = Dict(String(h) => j for (j, h) in enumerate(vec(hdr)))
getc(c) = raw[:, col[c]]

L    = Int.(getc("L"));         δ    = Float64.(getc("doping"))
site = Int.(getc("site"));      chi  = Int.(getc("chi_prior"))
obs  = String.(getc("observable"))
Δ    = Float64.(getc("Delta")); Ptrue = Float64.(getc("true"))
G    = Float64.(getc("G_emp")); Glo  = Float64.(getc("G_lo")); Ghi = Float64.(getc("G_hi"))
Gmax = Float64.(getc("G_max")); Sb   = Float64.(getc("S_bond"))
fid  = Float64.(getc("prior_fid"))

sel(f) = findall(f, 1:length(L))
Ls = sort(unique(L)); chis = sort(unique(chi)); δs = sort(unique(δ)); obss = unique(obs)

# ============================================================
# A. 主張の検定: 「利得はグローバル忠実度が落ちても、鎖のどのサイトでも生き残るか」
# ============================================================
println("="^92)
println("A. 全サイト最小利得 (これが 1 を大きく上回れば「どこでも効く」と言える)")
@printf("%6s %6s %5s %11s %9s %9s %9s %9s\n",
        "L", "doping", "chi_p", "prior_fid", "min G", "median G", "max G", "min G_lo")
for Lv in Ls, dv in δs, c in chis
    idx = sel(i -> L[i]==Lv && δ[i]==dv && chi[i]==c)
    isempty(idx) && continue
    @printf("%6d %6.3f %5d %11.2e %9.2f %9.2f %9.2f %9.2f\n",
            Lv, dv, c, fid[idx[1]], minimum(G[idx]), median(G[idx]),
            maximum(G[idx]), minimum(Glo[idx]))
end

# ============================================================
# B. サイト依存性: 端 vs 中央、およびばらつき
# ============================================================
println("\n", "="^92)
println("B. サイト依存性 (chi_p=8): 端(i≤2 or i≥L-1)と中央(|i-L/2|≤1)の利得比")
@printf("%-14s %6s %6s %9s %9s %7s %10s %12s\n",
        "observable", "L", "doping", "G_edge", "G_center", "比", "G_max中央", "median G/Gmax")
for o in obss, Lv in Ls, dv in δs
    idx = sel(i -> obs[i]==o && L[i]==Lv && δ[i]==dv && chi[i]==8)
    length(idx) < 4 && continue
    g = G[idx]; s = site[idx]
    edge = mean(g[(s .<= 2) .| (s .>= Lv-1)])
    ctr  = mean(g[abs.(s .- Lv/2) .<= 1])
    @printf("%-14s %6d %6.3f %9.2f %9.2f %7.2f %10.2f %12.3f\n",
            o, Lv, dv, edge, ctr, edge/ctr, median(Gmax[idx]), median(g ./ Gmax[idx]))
end

# ============================================================
# C. サイト依存性は何で決まるか: Δ とボンドエントロピーの相関
# ============================================================
println("\n", "="^92)
println("C. 局所prior誤差 |Δ(i)| と ボンドエントロピー S(i) の相関 (サイト分解)")
@printf("%6s %6s %5s %10s %10s %14s\n", "L", "doping", "chi_p", "n", "median|Δ|", "corr(|Δ|,S)")
for Lv in Ls, dv in δs, c in chis
    idx = sel(i -> L[i]==Lv && δ[i]==dv && chi[i]==c && abs(Δ[i]) > 1e-12)
    length(idx) < 10 && continue
    @printf("%6d %6.3f %5d %10d %10.2e %+14.3f\n",
            Lv, dv, c, length(idx), median(abs.(Δ[idx])), cor(abs.(Δ[idx]), Sb[idx]))
end

# ============================================================
# D. 2レジーム診断のサイト分解
# ============================================================
println("\n", "="^92)
println("D. 2レジーム診断: G/G_max が 1 に近い = ショットノイズ律速, ≪1 = Δ律速")
@printf("%6s %6s %5s %14s %14s %14s\n", "L", "doping", "chi_p",
        "median G/Gmax", "min G/Gmax", "Δ律速の割合(<0.5)")
for Lv in Ls, dv in δs, c in chis
    idx = sel(i -> L[i]==Lv && δ[i]==dv && chi[i]==c)
    isempty(idx) && continue
    r = G[idx] ./ Gmax[idx]
    @printf("%6d %6.3f %5d %14.3f %14.3f %14.3f\n",
            Lv, dv, c, median(r), minimum(r), mean(r .< 0.5))
end

# ============================================================
# 作図
# ============================================================
Lmax = maximum(Ls)

# 図1: G vs サイト位置 (半充填, chi_p=4)
p1 = plot(xlabel="site position  i / L", ylabel="G = Var_std / Var_CRM", yscale=:log10,
          legend=:bottomleft, title="Site-resolved gain, ZZ onsite (δ=0, χ_p=4)")
for Lv in Ls
    idx = sel(i -> obs[i]=="ZZ onsite" && L[i]==Lv && δ[i]==0.0 && chi[i]==4)
    isempty(idx) && continue
    o = sortperm(site[idx]); idx = idx[o]
    plot!(p1, site[idx]./Lv, G[idx], ribbon=(G[idx].-Glo[idx], Ghi[idx].-G[idx]),
          marker=:circle, ms=2.5, lw=1.8, label="L=$Lv (F=$(@sprintf("%.1e", fid[idx[1]])))")
end
hline!(p1, [1.0], color=:black, ls=:dot, label="no gain")

# 図2: ドープ系 — 悪いpriorではサイトによって「損」になる
p2 = plot(xlabel="site i", ylabel="G", yscale=:log10, legend=:bottomright,
          title="Doped chain δ=1/8, L=$Lmax: ZZ onsite")
hspan!(p2, [1e-2, 1.0], color=:red, alpha=0.08, label="CRM loses (G<1)")
for c in [2, 4, 8]
    idx = sel(i -> obs[i]=="ZZ onsite" && L[i]==Lmax && δ[i]==0.125 && chi[i]==c)
    isempty(idx) && continue
    s = sortperm(site[idx]); idx = idx[s]
    plot!(p2, site[idx], G[idx], marker=:circle, ms=2, lw=1.4, label="χ_p=$c")
end
hline!(p2, [1.0], color=:black, ls=:dot, label="")

# 図3: 局所prior誤差はサイトにも L にもほとんど依存しない (＝局所性の実体)
p3 = plot(xlabel="site position  i / L", ylabel="|Δ(i)|  (ZZ onsite)", yscale=:log10,
          legend=:bottomright, title="Local prior error is uniform and L-independent")
for c in [2, 4], Lv in Ls
    idx = sel(i -> obs[i]=="ZZ onsite" && L[i]==Lv && δ[i]==0.0 && chi[i]==c)
    isempty(idx) && continue
    s = sortperm(site[idx]); idx = idx[s]
    plot!(p3, site[idx]./Lv, abs.(Δ[idx]), lw=1.4, alpha=0.85,
          ls = c==2 ? :solid : :dash, label="χ_p=$c, L=$Lv")
end

# 図4: 2レジーム診断のサイト分解
p4 = plot(xlabel="site position  i / L", ylabel="G / G_max", ylims=(0, 1.15),
          legend=:bottomleft, title="Regime diagnostic, L=$Lmax (δ=0)")
for c in chis
    idx = sel(i -> obs[i]=="ZZ onsite" && L[i]==Lmax && δ[i]==0.0 && chi[i]==c)
    isempty(idx) && continue
    s = sortperm(site[idx]); idx = idx[s]
    plot!(p4, site[idx]./Lmax, G[idx]./Gmax[idx], marker=:circle, ms=2.5, lw=1.6,
          label="χ_p=$c")
end
hline!(p4, [1.0], color=:black, ls=:dot, label="shot-noise limit")

fig = plot(p1, p2, p3, p4, layout=(2,2), size=(1400, 900), margin=6Plots.mm)
out = joinpath(@__DIR__, "crm_site_resolved.png")
savefig(fig, out)
println("\nfigure saved: $out")
