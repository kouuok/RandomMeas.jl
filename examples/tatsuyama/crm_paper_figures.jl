# ============================================================
# 論文向け清書図の生成
#
# これまでの数値実験のTSVから3枚の図を作る:
#  Fig 1 (crm_fig1_gainlaw.png):  利得法則のマスタープロット
#         全実験 (1D L掃引 / nm掃引 / 2D U=4,8) の純Pauli列について
#         経験利得 G_emp vs 理論 G_theo (5桁にわたり y=x)
#  Fig 2 (crm_fig2_locality.png): (a) 局所利得は L に依存せず、
#         グローバル忠実度の崩壊を生き延びる  (b) ショット配分 nm 依存と理論曲線
#  Fig 3 (crm_fig3_priors2d.png): 2Dシリンダーでの prior 比較
#         (UHF / 対称性回復UHF / 切断MPS chi=8..64)
#
# スタイル: Okabe-Ito 色覚多様性対応パレット + 系列ごとに固有マーカー
# 実行: JULIA_LOAD_PATH="@:@v#.#:@stdlib" julia --project=Hubbard_MPS_Env_v2 crm_paper_figures.jl
# ============================================================

using DelimitedFiles
using Printf
using Plots

# ---- スタイル ----
const OI = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#56B4E9", "#E69F00"]  # Okabe-Ito
const MK = [:circle, :square, :diamond, :utriangle, :dtriangle, :pentagon]
default(fontfamily="Helvetica", framestyle=:box, grid=true, gridalpha=0.12,
        gridlinewidth=0.5, linewidth=1.8, markersize=5.5, markerstrokewidth=0.6,
        markerstrokecolor=:white, guidefontsize=11, tickfontsize=9,
        legendfontsize=8, titlefontsize=11, dpi=200)

dir = @__DIR__
read_tsv(f) = readdlm(joinpath(dir, f), '\t'; header=true)

# ============================================================
# Fig 1: 利得法則マスタープロット
# ============================================================
sets = [("crm_mps_scaling_results.tsv", "1D chain (L = 8–32)", 8, 9),
        ("crm_sweep_nm.tsv",            "1D shot sweep (nm = 1–1000)", 7, 8),
        ("crm_2d_symuhf_results.tsv",   "2D cylinder (W = 4, U = 4, 8)", 8, 9)]

p1 = plot(xscale=:log10, yscale=:log10, legend=:topleft, size=(560, 520),
          xlabel="theoretical gain  G_theory", ylabel="empirical gain  G",
          title="Variance gain: theory vs experiment")
allv = Float64[]
for (i, (file, lab, gcol, tcol)) in enumerate(sets)
    d, hdr = read_tsv(file)
    pure_col = findfirst(==("pure"), vec(hdr))
    gx = Float64[]; gy = Float64[]
    for r in 1:size(d, 1)
        string(d[r, pure_col]) == "true" || continue
        gt = d[r, tcol]; ge = d[r, gcol]
        (gt isa Number && ge isa Number && isfinite(gt) && isfinite(ge)) || continue
        push!(gx, gt); push!(gy, ge)
    end
    append!(allv, gx); append!(allv, gy)
    scatter!(p1, gx, gy, color=OI[i], marker=MK[i], alpha=0.75, label=lab)
    @printf("Fig1: %-32s %3d points\n", lab, length(gx))
end
lims = (0.6 * minimum(allv), 1.8 * maximum(allv))
plot!(p1, [lims...], [lims...], color=:black, ls=:dash, lw=1.2, label="G = G_theory")
xlims!(p1, lims); ylims!(p1, lims)
savefig(p1, joinpath(dir, "crm_fig1_gainlaw.png"))

# ============================================================
# Fig 2: 局所性 (G vs L + 忠実度崩壊) と ショット配分 (G vs nm)
# ============================================================
d, hdr = read_tsv("crm_mps_scaling_results.tsv")
h = Dict(string(k) => i for (i, k) in enumerate(vec(hdr)))
Ls = [8, 16, 32]

pa = plot(xlabel="chain length  L", ylabel="G,  fidelity", yscale=:log10,
          legend=:left, title="(a) Local gains survive fidelity collapse",
          xticks=Ls)
series_a = [("ZZ up-up r=1",  2,  "ZZ r = 1  (chi_p = 2)"),
            ("SzSz r=1",      2,  "SzSz r = 1  (chi_p = 2)"),
            ("ZZ onsite(i0)", 16, "on-site ZZ  (chi_p = 16)")]
for (i, (name, chi, lab)) in enumerate(series_a)
    g = [only([d[r, h["G_emp"]] for r in 1:size(d,1)
               if d[r, h["L"]] == L && d[r, h["chi_prior"]] == chi &&
                  string(d[r, h["observable"]]) == name]) for L in Ls]
    plot!(pa, Ls, g, color=OI[i], marker=MK[i], label=lab)
end
for (j, chi) in enumerate((2, 16))
    f = [only(unique([d[r, h["prior_fid"]] for r in 1:size(d,1)
                      if d[r, h["L"]] == L && d[r, h["chi_prior"]] == chi])) for L in Ls]
    plot!(pa, Ls, f, color=:gray40, ls=(j == 1 ? :dash : :dot), marker=:x,
          markerstrokecolor=:gray40, lw=1.4, label="global fidelity (chi_p = $chi)")
end
hline!(pa, [1.0], color=:black, ls=:dot, lw=1, label="")
annotate!(pa, 30.5, 0.011, text("F = 0.007", 8, :gray30, :center))

dn, hdrn = read_tsv("crm_sweep_nm.tsv")
hn = Dict(string(k) => i for (i, k) in enumerate(vec(hdrn)))
nms = [10, 30, 100, 300, 1000]
theory_G(absA, P, Δ, nm) = ((3.0^absA - 1)*P^2 + 3.0^absA*(1-P^2)/nm) /
                           ((3.0^absA - 1)*Δ^2 + 3.0^absA*(1-P^2)/nm)

pb = plot(xlabel="shots per setting  nm", ylabel="G", xscale=:log10, yscale=:log10,
          legend=:topleft, title="(b) Gain vs shot allocation (L = 8, U = 4)")
series_b = [("ZZ onsite(i0)", 32, "on-site ZZ  (chi_p = 32)"),
            ("ZZ up-up r=1",  32, "ZZ r = 1  (chi_p = 32)"),
            ("ZZ onsite(i0)", 4,  "on-site ZZ  (chi_p = 4)")]
for (i, (name, chi, lab)) in enumerate(series_b)
    sel = [r for r in 1:size(dn,1)
           if dn[r, hn["chi_prior"]] == chi && string(dn[r, hn["observable"]]) == name]
    g = [only([dn[r, hn["G_emp"]] for r in sel if dn[r, hn["nm"]] == nm]) for nm in nms]
    plot!(pb, nms, g, color=OI[i], marker=MK[i], label=lab)
    r0 = sel[1]
    P = dn[r0, hn["true"]]; Δ = dn[r0, hn["Delta"]]
    nmc = 10 .^ range(0.9, 3.1; length=60)
    plot!(pb, nmc, theory_G.(2, P, Δ, nmc), color=OI[i], ls=:dot, lw=1.1, label="")
end
plot!(pb, [NaN], [NaN], color=:gray, ls=:dot, lw=1.1, label="theory")  # 凡例用ダミー
hline!(pb, [1.0], color=:black, ls=:dot, lw=1, label="")

p2 = plot(pa, pb, layout=(1, 2), size=(1120, 460), margin=7Plots.mm)
savefig(p2, joinpath(dir, "crm_fig2_locality.png"))

# ============================================================
# Fig 3: 2D prior 比較
# ============================================================
d2, hdr2 = read_tsv("crm_2d_symuhf_results.tsv")
h2 = Dict(string(k) => i for (i, k) in enumerate(vec(hdr2)))
prior_order = ["UHF", "UHF-sym", "chi=8", "chi=16", "chi=32", "chi=64"]
obs_sel = [("ZZ onsite",   "on-site ZZ"),
           ("DoubleOcc",   "double occupancy"),
           ("hop y-bond",  "kinetic bond"),
           ("ZZ up-up x1", "ZZ  x-neighbor"),
           ("SzSz x-bond", "SzSz  x-neighbor"),
           ("ZZ up-up y2", "ZZ  y-distance 2")]

panels = []
for (iu, U) in enumerate((4.0, 8.0))
    p = plot(xticks=(1:6, prior_order), xrotation=25, yscale=:log10,
             ylabel=(iu == 1 ? "G" : ""), legend=(iu == 1 ? :none : :outerright),
             title=@sprintf("W = 4 cylinder,  U = %.0f", U))
    for (i, (name, lab)) in enumerate(obs_sel)
        g = [only([d2[r, h2["G_emp"]] for r in 1:size(d2,1)
                   if d2[r, h2["U"]] == U && string(d2[r, h2["prior"]]) == pr &&
                      string(d2[r, h2["observable"]]) == name]) for pr in prior_order]
        plot!(p, 1:6, g, color=OI[i], marker=MK[i], label=lab)
    end
    hline!(p, [1.0], color=:black, ls=:dot, lw=1, label="")
    vline!(p, [2.5], color=:gray70, ls=:dash, lw=1, label="")
    push!(panels, p)
end
p3 = plot(panels..., layout=(1, 2), size=(1240, 460), margin=7Plots.mm)
savefig(p3, joinpath(dir, "crm_fig3_priors2d.png"))

println("figures saved: crm_fig1_gainlaw.png, crm_fig2_locality.png, crm_fig3_priors2d.png")
