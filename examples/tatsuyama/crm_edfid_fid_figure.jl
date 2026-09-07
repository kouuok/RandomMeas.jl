# 「何が利得を決めているのか」を全3912点で比べる。
#   (a) 大域忠実度      (b) 台の忠実度      (c) 相対誤差 ε
#   julia --project=. examples/tatsuyama/crm_edfid_fid_figure.jl
using DelimitedFiles, Printf, Statistics, Plots
gr(fontfamily="Helvetica", legendfontsize=6, guidefontsize=9,
   titlefontsize=10, tickfontsize=8)
const DIR = @__DIR__
A, h = readdlm(joinpath(DIR,"crm_edfid_all.tsv"), '\t'; header=true)
h = vec(h); c(n) = findfirst(==(n), h)
S(n) = String.(A[:, c(n)]); F(n) = Float64.(A[:, c(n)])
obs, pri, exa = S("observable"), S("prior"), S("exact")
gf, fs, eps, G, Gm = F("global_fid"), F("F_supp"), F("eps"), F("G"), F("G_max")
ok = (exa .== "yes") .& isfinite.(Gm) .& isfinite.(eps) .& (eps .> 0)
@printf("点数 %d\n", sum(ok))

function spearman(x, y)
    n = length(x); rx = zeros(n); ry = zeros(n)
    for (r,i) in enumerate(sortperm(x)); rx[i] = r; end
    for (r,i) in enumerate(sortperm(y)); ry[i] = r; end
    rx .-= mean(rx); ry .-= mean(ry)
    (rx' * ry) / sqrt((rx' * rx) * (ry' * ry))
end

OBS = ["ZZ onsite","ZZ up-up nb","SzSz nb","SxSx nb","DoubleOcc","Sz"]
COL = [:darkorange, :steelblue, :mediumpurple, :seagreen, :goldenrod, :gray55]
MK  = [:circle, :utriangle, :diamond, :star5, :rect, :xcross]

function panel(xv, yv, xlab, ylab, ttl; xl, yl, leg=false)
    p = plot(xscale=:log10, yscale=:log10, xlabel=xlab, ylabel=ylab, title=ttl,
             xlims=xl, ylims=yl, legend=leg ? :bottomleft : false)
    for (j,o) in enumerate(OBS)
        m = ok .& (obs .== o)
        any(m) || continue
        scatter!(p, clamp.(xv[m], xl[1]*1.02, xl[2]*0.98), clamp.(yv[m], yl[1]*1.02, yl[2]*0.98),
                 ms=2.6, mc=COL[j], msw=0, alpha=0.42, marker=MK[j], label=o)
    end
    p
end

ρ1 = spearman(gf[ok], G[ok]); ρ2 = spearman(fs[ok], G[ok])
ratio = G ./ Gm
ρ3 = spearman(eps[ok], ratio[ok])

pa = panel(gf, G, "global fidelity  F(ρ, σ)", "gain G",
           "(a) global fidelity"; xl=(3e-5,2.0), yl=(0.3,3e3), leg=true)
hline!(pa, [1.0], color=:black, ls=:dot, lw=1.3, label="")
annotate!(pa, 4e-5, 1.6e3, text(@sprintf("Spearman ρ = %+.3f", ρ1), 8, :left, :black))

pb = panel(fs, G, "fidelity on the support  F_A", "gain G",
           "(b) fidelity on the observable's own support"; xl=(3e-3,1.6), yl=(0.3,3e3))
hline!(pb, [1.0], color=:black, ls=:dot, lw=1.3, label="")
annotate!(pb, 3.4e-3, 1.6e3, text(@sprintf("Spearman ρ = %+.3f", ρ2), 8, :left, :black))
annotate!(pb, 3.4e-3, 5.5e2, text("no better than (a)", 8, :left, :firebrick))

pc = panel(eps, ratio, "relative error  ε = |Δ| / |⟨P⟩|", "G / G_max",
           "(c) the relative error, against the ceiling"; xl=(1e-6,20.0), yl=(2e-4,2.0))
annotate!(pc, 1.4e-6, 8e-4, text(@sprintf("Spearman ρ = %+.3f", ρ3), 8, :left, :black))

fig = plot(pa, pb, pc, layout=(1,3), size=(1680,470), margin=7Plots.mm, bottom_margin=9Plots.mm)
savefig(fig, joinpath(DIR,"crm_new_fig_edfid_what_predicts.png"))
@printf("保存 (ρ: 大域 %.3f / 台 %.3f / ε %.3f)\n", ρ1, ρ2, ρ3)
