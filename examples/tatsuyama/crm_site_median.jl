# 各観測量・各サイトごとの利得 G の中央値
#   julia --project=Hubbard_MPS_Env_v2 crm_site_median.jl
#
# 中央値は χ_p ∈ {2,4,8,16,32} について取る（そのサイト・その観測量の
# 「prior品質に対して代表的な利得」）。サイト添字 i は L ごとに意味が
# 変わるので、L と δ は集約せず層別のまま残す。
using DelimitedFiles, Statistics, Printf, Plots

raw, hdr = readdlm(joinpath(@__DIR__, "crm_site_resolved_results.tsv"), '\t'; header=true)
col = Dict(String(h) => j for (j, h) in enumerate(vec(hdr)))
gc(c) = raw[:, col[c]]

L = Int.(gc("L")); δ = Float64.(gc("doping")); site = Int.(gc("site"))
chi = Int.(gc("chi_prior")); obs = String.(gc("observable"))
G = Float64.(gc("G_emp")); Gmax = Float64.(gc("G_max")); Δ = Float64.(gc("Delta"))

# (L, δ, observable, site) ごとに χ_p で集約
key = Dict{Tuple{Int,Float64,String,Int},Vector{Int}}()
for i in 1:length(L)
    push!(get!(key, (L[i], δ[i], obs[i], site[i]), Int[]), i)
end

out = joinpath(@__DIR__, "crm_site_median_results.tsv")
open(out, "w") do io
    println(io, "L\tdoping\tobservable\tsite\tG_median\tG_min\tG_max_over_chi\tGmax_ceiling\tratio_median")
    for k in sort(collect(keys(key)); by = x -> (x[1], x[2], x[3], x[4]))
        idx = key[k]
        g = G[idx]
        @printf(io, "%d\t%.3f\t%s\t%d\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n",
                k[1], k[2], k[3], k[4], median(g), minimum(g), maximum(g),
                median(Gmax[idx]), median(g ./ Gmax[idx]))
    end
end
println("saved: $out\n")

obss = ["n", "Sz", "DoubleOcc", "ZZ onsite", "ZZ up-up r=1", "SzSz r=1", "hop up r=1"]

# --- 表示1: L=16 の全サイト（画面に収まる大きさ）------------------------
for dv in (0.0, 0.125)
    @printf("\n=== 各観測量・各サイトの G 中央値 (χ_p で集約)  L=16, δ=%.3f ===\n", dv)
    @printf("%-14s", "observable")
    for s in 1:16; @printf("%7d", s); end
    println()
    for o in obss
        @printf("%-14s", o)
        for s in 1:16
            k = (16, dv, o, s)
            if haskey(key, k)
                @printf("%7.1f", median(G[key[k]]))
            else
                @printf("%7s", "-")
            end
        end
        println()
    end
end

# --- 表示2: L=128 は端/バルクに要約 -------------------------------------
println("\n\n=== L=128: 端(i≤4) と バルク(5≤i≤L-4) の G中央値 ===")
@printf("%-14s %7s | %9s %9s %9s | %9s\n",
        "observable", "doping", "端 中央値", "バルク中央値", "バルク最小", "端/バルク")
for o in obss, dv in (0.0, 0.125)
    edge = Float64[]; bulk = Float64[]
    for s in 1:128
        k = (128, dv, o, s)
        haskey(key, k) || continue
        m = median(G[key[k]])
        (s <= 4 || s >= 125) ? push!(edge, m) : push!(bulk, m)
    end
    isempty(bulk) && continue
    @printf("%-14s %7.3f | %9.2f %9.2f %9.2f | %9.2f\n",
            o, dv, median(edge), median(bulk), minimum(bulk), median(edge)/median(bulk))
end

# --- 表示3: サイト方向のばらつきの大きさ --------------------------------
println("\n\n=== サイト方向のばらつき (L=128, G中央値の site 間 変動係数) ===")
@printf("%-14s %8s %10s %10s %10s\n", "observable", "doping", "site平均", "site標準偏差", "変動係数")
for o in obss, dv in (0.0, 0.125)
    v = Float64[]
    for s in 1:128
        k = (128, dv, o, s); haskey(key, k) && push!(v, median(G[key[k]]))
    end
    isempty(v) && continue
    @printf("%-14s %8.3f %10.2f %10.2f %10.3f\n", o, dv, mean(v), std(v), std(v)/mean(v))
end

# --- 図: 観測量 × サイト の G中央値ヒートマップ (log スケール) -----------
function heat(Lv, dv)
    M = fill(NaN, length(obss), Lv)
    for (r, o) in enumerate(obss), s in 1:Lv
        k = (Lv, dv, o, s)
        haskey(key, k) && (M[r, s] = median(G[key[k]]))
    end
    heatmap(1:Lv, 1:length(obss), log10.(M), c=:viridis,
            yticks=(1:length(obss), obss), xlabel="site i",
            colorbar_title="log₁₀ G (median over χ_p)",
            title="L=$Lv, δ=$(dv==0 ? "0" : "1/8")")
end

fig = plot(heat(128, 0.0), heat(128, 0.125), layout=(2,1), size=(1250, 620),
           left_margin=14Plots.mm, bottom_margin=5Plots.mm)
outp = joinpath(@__DIR__, "crm_site_median.png")
savefig(fig, outp)
println("\nfigure saved: $outp")
