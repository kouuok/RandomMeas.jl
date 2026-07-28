# 2Dサイト分解の解析・作図: 安価な prior はどこで破綻するか
#   julia --project=Hubbard_MPS_Env_v2 crm_2d_site_plot.jl
using DelimitedFiles, Statistics, Printf, Plots

const FILES = [("crm_2d_site_U8p0_h0.tsv",  "U=8, half"),
               ("crm_2d_site_U8p0_h4.tsv",  "U=8, δ=1/8"),
               ("crm_2d_site_U4p0_h4.tsv",  "U=4, δ=1/8")]
const W = 4

function load(f)
    raw, hdr = readdlm(joinpath(@__DIR__, f), '\t'; header=true)
    c = Dict(String(h) => j for (j, h) in enumerate(vec(hdr)))
    (x = Int.(raw[:, c["x"]]), y = Int.(raw[:, c["y"]]),
     ob = String.(raw[:, c["observable"]]), pr = String.(raw[:, c["prior"]]),
     re = Float64.(raw[:, c["relerr"]]), G = Float64.(raw[:, c["G"]]),
     dens = Float64.(raw[:, c["dens"]]), tru = Float64.(raw[:, c["true"]]))
end

# ---- 定量サマリ -----------------------------------------------------------
println("="^100)
println("prior ごとの相対誤差の空間変動 (max/min over sites)")
@printf("%-14s %-14s %10s %10s %10s %12s\n", "config", "observable",
        "chi8", "chi32", "UHF", "UHF/chi32")
for (f, lab) in FILES
    d = load(f)
    for o in ["ZZ onsite", "DoubleOcc", "SzSz nb"]
        vals = Float64[]
        for p in ["chi8", "chi32", "UHF"]
            idx = findall(i -> d.ob[i]==o && d.pr[i]==p && !isnan(d.re[i]), 1:length(d.x))
            push!(vals, isempty(idx) ? NaN : maximum(d.re[idx])/minimum(d.re[idx]))
        end
        @printf("%-14s %-14s %10.1f %10.1f %10.1f %12.1f\n", lab, o, vals..., vals[3]/vals[2])
    end
end

println("\n", "="^100)
println("UHF prior の相対誤差と局所電子密度の相関 (ドープ系のみ意味を持つ)")
@printf("%-14s %-14s %14s %14s\n", "config", "observable", "corr(ε, n)", "corr(ε, |n-1|)")
for (f, lab) in FILES
    d = load(f)
    for o in ["ZZ onsite", "DoubleOcc"]
        idx = findall(i -> d.ob[i]==o && d.pr[i]=="UHF" && !isnan(d.re[i]), 1:length(d.x))
        length(idx) < 8 && continue
        n = d.dens[idx]; e = log.(d.re[idx])
        std(n) < 1e-6 && (@printf("%-14s %-14s %14s %14s\n", lab, o, "(一様)", "(一様)"); continue)
        @printf("%-14s %-14s %+14.3f %+14.3f\n", lab, o, cor(e, n), cor(e, abs.(n .- 1)))
    end
end

# ---- 図: 相対誤差の空間地図 ----------------------------------------------
function heat(d, o, p, ttl; Lx)
    M = fill(NaN, W, Lx)
    for i in 1:length(d.x)
        (d.ob[i]==o && d.pr[i]==p) || continue
        M[d.y[i], d.x[i]] = d.re[i]
    end
    heatmap(1:Lx, 1:W, log10.(M), c=:magma, xlabel="x", ylabel="y",
            colorbar_title="log₁₀ ε", title=ttl, titlefontsize=9)
end

d0 = load(FILES[1][1]); d8 = load(FILES[2][1]); d4 = load(FILES[3][1])
Lx = maximum(d8.x)
plots = Any[]
for (d, lab) in ((d0, "U=8 half"), (d8, "U=8 δ=1/8"))
    for p in ("chi8", "UHF")
        push!(plots, heat(d, "ZZ onsite", p, "$lab — $p (ZZ onsite)"; Lx))
    end
end

# 密度プロファイルと UHF 誤差の対応
pd = plot(xlabel="column x", ylabel="⟨n⟩", legend=:bottomright,
          title="Charge profile vs UHF error", titlefontsize=9)
for (d, lab) in ((d8, "U=8 δ=1/8"), (d4, "U=4 δ=1/8"))
    ys = [mean(d.dens[findall(i -> d.x[i]==x, 1:length(d.x))]) for x in 1:Lx]
    plot!(pd, 1:Lx, ys, marker=:circle, lw=2, label="⟨n⟩ $lab")
end
pe = twinx(pd)
for (d, lab) in ((d8, "U=8 δ=1/8"), (d4, "U=4 δ=1/8"))
    ys = [median(filter(!isnan, d.re[findall(i -> d.x[i]==x && d.pr[i]=="UHF" &&
                                             d.ob[i]=="ZZ onsite", 1:length(d.x))])) for x in 1:Lx]
    plot!(pe, 1:Lx, ys, marker=:square, ls=:dash, lw=2, color=:red,
          ylabel="UHF ε (ZZ onsite)", yscale=:log10, label="")
end
push!(plots, pd)

# prior 間の空間変動の比較
pv = plot(xlabel="column x", ylabel="median ε (ZZ onsite)", yscale=:log10,
          legend=:topright, title="U=8 δ=1/8: prior comparison", titlefontsize=9)
for p in ("chi8", "chi32", "UHF")
    ys = [median(filter(!isnan, d8.re[findall(i -> d8.x[i]==x && d8.pr[i]==p &&
                                              d8.ob[i]=="ZZ onsite", 1:length(d8.x))])) for x in 1:Lx]
    plot!(pv, 1:Lx, ys, marker=:circle, lw=2, label=p)
end
push!(plots, pv)

fig = plot(plots..., layout=(3,2), size=(1300, 1100), margin=5Plots.mm)
out = joinpath(@__DIR__, "crm_2d_site.png")
savefig(fig, out)
println("\nfigure saved: $out")
