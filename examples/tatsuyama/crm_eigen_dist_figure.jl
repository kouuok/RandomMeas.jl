# 距離 r の Z↑Z↑ 相関で「固有状態への近さ」と利得を見る図(crm_eigen_corr_analysis.py の出力を読む)。
#   (a) 1次元鎖: 距離とともに真の相関は減衰し、UHF は固有状態の近くに留まる → r=2 から全対で損
#   (b) 2次元: 相関の減衰が遅いので、群平均した UHF は損をしない
#   (c) 全サイト対の (|<P>ρ|, |<P>σ|) 平面: UHF・UHF-sym・MPS はそれぞれ水平の帯・水平の帯・対角線付近に並ぶ
#   (d) 損得は y/x = 2 で厳密に切り替わる
#   julia --project=. examples/tatsuyama/crm_eigen_dist_figure.jl
using DelimitedFiles, Printf, Statistics, Plots
gr(fontfamily="Helvetica", legendfontsize=7, guidefontsize=9, titlefontsize=10, tickfontsize=8)
const DIR = @__DIR__
A, h = readdlm(joinpath(DIR, "crm_eigcorr_pairs.tsv"), '\t'; header=true)
h = vec(h); col(n) = A[:, findfirst(==(n), h)]
Wc, geo, LXc, Uc, rc = Int.(col("W")), string.(col("geometry")), Int.(col("LX")), Float64.(col("U")), Int.(col("r"))
pri, xv, yv, Gv = string.(col("prior")), Float64.(col("x")), Float64.(col("y")), Float64.(col("G"))

hexc(s) = RGB(parse(Int, s[2:3]; base=16)/255, parse(Int, s[4:5]; base=16)/255, parse(Int, s[6:7]; base=16)/255)
const INK = hexc("#52514e")
const SER = [("UHF", hexc("#2a78d6"), "UHF (as is)"),
             ("UHF-sym", hexc("#eb6834"), "UHF averaged over SU(2)"),
             ("chi8", hexc("#1baf7a"), "truncated MPS, χ = 8")]
const LOGT = [1e-2, 1e-1, 1, 10, 100]

function vsr!(p, W, g, L, U)
    m0 = (Wc .== W) .& (geo .== g) .& (LXc .== L) .& (Uc .== U)
    any(m0) || return false
    for (q, c, lab) in SER
        m = m0 .& (pri .== q); rs = sort(unique(rc[m]))
        med = [median(Gv[m .& (rc .== r)]) for r in rs]
        lo  = [minimum(Gv[m .& (rc .== r)]) for r in rs]; hi = [maximum(Gv[m .& (rc .== r)]) for r in rs]
        plot!(p, rs, med; ribbon=(med .- lo, hi .- med), fillalpha=0.12, color=c, lw=2, marker=:circle, ms=6, msw=0, label=lab)
    end
    true
end

# (a) 1次元: あれば L=64、なければ L=32
La = any((Wc .== 1) .& (LXc .== 64)) ? 64 : 32
pa = plot(yscale=:log10, yticks=LOGT, ylims=(5e-3, 200), legend=:topright,
          xlabel="distance r   (both sites in the central half)", ylabel="gain G   (median; band = min–max over pairs)",
          title=@sprintf("(a) 1D chain L=%d, U=8: raw UHF loses from r = 2 on", La))
hline!(pa, [1.0], color=:black, ls=:dot, lw=1.2, label="")
vsr!(pa, 1, "cylinder", La, 8.0)
let m = (Wc .== 1) .& (LXc .== La) .& (Uc .== 8.0) .& (pri .== "UHF")
    x1, x2 = median(abs.(xv[m .& (rc .== 1)])), median(abs.(xv[m .& (rc .== 2)]))
    y1, y2 = median(abs.(yv[m .& (rc .== 1)])), median(abs.(yv[m .& (rc .== 2)]))
    annotate!(pa, 2.4, 0.012, text(@sprintf("r = 1 → 2:  truth |⟨P⟩ρ| %.2f → %.2f,\nraw UHF |⟨P⟩σ| %.2f → %.2f  (stays near the eigenstate)", x1, x2, y1, y2), 7, :left, INK))
end

# (b) 2次元: あれば 2本脚梯子 L=8、なければ幅4シリンダー L=3
has_w2 = any((Wc .== 2) .& (LXc .== 8))
(Wb, Lb) = has_w2 ? (2, 8) : (4, 3)
pb = plot(yscale=:log10, yticks=LOGT, ylims=(5e-3, 200), legend=:topright,
          xlabel="distance r   (Manhattan, periodic in y)", ylabel="gain G",
          title=@sprintf("(b) %s, U=8: correlations decay slowly", has_w2 ? "2-leg ladder L=8" : "width-4 cylinder L=3"))
hline!(pb, [1.0], color=:black, ls=:dot, lw=1.2, label="")
vsr!(pb, Wb, "cylinder", Lb, 8.0)

# (c) 平面
pc = plot(xlims=(-0.02, 1.02), ylims=(-0.02, 1.02), aspect_ratio=:equal, legend=:bottomright,
          xlabel="|⟨P⟩ρ|   (truth)", ylabel="|⟨P⟩σ|   (prior)",
          title="(c) every pair, every system and U: three bands")
plot!(pc, [0, 0.5, 0], [0, 1, 1]; seriestype=:shape, fillcolor=:gray, fillalpha=0.13, lw=0, label="CRM loses (G < 1)")
plot!(pc, [0, 1], [0, 1]; color=:black, lw=1.2, label="σ = ρ  (G = Gmax)")
plot!(pc, [0, 0.5], [0, 1]; color=:black, lw=1.2, ls=:dash, label="break-even |⟨P⟩σ| = 2|⟨P⟩ρ|")
for (q, c, lab) in SER
    m = (pri .== q) .& (xv .* yv .> 0)
    scatter!(pc, abs.(xv[m]), abs.(yv[m]); color=c, ms=3.2, msw=0, alpha=0.45, label=lab)
end

# (d) y/x で切り替わる
pd = plot(xscale=:log10, yscale=:log10, yticks=LOGT, xlims=(0.03, 40), ylims=(5e-3, 200), legend=:bottomleft,
          xlabel="⟨P⟩σ / ⟨P⟩ρ", ylabel="gain G", title="(d) the switch is exactly at ⟨P⟩σ/⟨P⟩ρ = 2")
hline!(pd, [1.0], color=:black, ls=:dot, lw=1.2, label="")
vline!(pd, [1.0], color=:gray, lw=1, label="")
vline!(pd, [2.0], color=:black, ls=:dash, lw=1.2, label="")
nskip = 0
for (q, c, lab) in SER
    m = (pri .== q) .& (yv ./ xv .> 0)
    global nskip += sum((pri .== q) .& .!(yv ./ xv .> 0))
    scatter!(pd, yv[m] ./ xv[m], Gv[m]; color=c, ms=3.2, msw=0, alpha=0.45, label=lab)
end
annotate!(pd, 0.9, 0.3, text("prior less extreme\nthan the truth: G ≥ 1", 7, :right, INK))
annotate!(pd, 2.3, 3.0, text("more than twice\nthe truth: G < 1", 7, :left, INK))
annotate!(pd, 38, 120, text(@sprintf("%d pairs (%d with ⟨P⟩σ/⟨P⟩ρ ≤ 0 not shown)", length(Gv), nskip), 7, :right, INK))

fig = plot(pa, pb, pc, pd; layout=(2, 2), size=(1500, 1180), margin=6Plots.mm)
out = joinpath(DIR, "crm_new_fig_eigen_dist.png"); savefig(fig, out); println("書き出し: ", out)
