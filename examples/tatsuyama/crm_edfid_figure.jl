# §10 の主要結果を3枚で示す。
#   (a) 系を伸ばすと大域忠実度だけが崩れ、台の忠実度と利得は動かない
#   (b) U を上げると電荷量とスピン量が逆を向く
#   (c) 全3920点: 大域忠実度は利得を決めていない
#   julia --project=. examples/tatsuyama/crm_edfid_figure.jl
using DelimitedFiles, Printf, Statistics, Plots
gr(fontfamily="Helvetica", legendfontsize=7, guidefontsize=9,
   titlefontsize=10, tickfontsize=8)
const DIR = @__DIR__
A, h = readdlm(joinpath(DIR,"crm_edfid_all.tsv"), '\t'; header=true)
h = vec(h); c(n) = findfirst(==(n), h)
col(n) = A[:, c(n)]
S(n) = String.(col(n)); F(n) = Float64.(col(n))
dim, geo, pri, obs, exa = S("dim"), S("geometry"), S("prior"), S("observable"), S("exact")
W, LX, ns, nA = Int.(col("W")), Int.(col("LX")), Int.(col("nsites")), Int.(col("nA"))
U, gf, fs, ds, eps, G, Gm = F("U"), F("global_fid"), F("F_supp"), F("D_supp"), F("eps"), F("G"), F("G_max")
ok = (exa .== "yes") .& isfinite.(Gm)        # 厳密かつ利得が定義できる点
@printf("使用点数 %d / %d\n", sum(ok), length(ok))

# ---------- (a) サイズ走査 ----------
Ls = [8,10,12,14,16]
pick(L, o, q) = begin
    m = (W.==1) .& (LX.==L) .& (geo.=="cylinder") .& (U.==8.0) .& (pri.=="UHF") .& (obs.==o)
    q[findfirst(m)]
end
gfv = [pick(L,"ZZ onsite",gf) for L in Ls]
fsv = [pick(L,"ZZ onsite",fs) for L in Ls]
gv  = [pick(L,"ZZ onsite",G)  for L in Ls]
pa = plot(xlabel="chain length L   (2L qubits)", ylabel="value relative to L = 8",
          title="(a) only the global fidelity notices the system size",
          legend=:left, yscale=:log10, ylims=(0.15, 2.2), xticks=Ls,
          yticks=([0.25,0.5,1.0,2.0], ["0.25","0.5","1.0","2.0"]))
plot!(pa, Ls, gfv./gfv[1], marker=:circle,  ms=6, lw=2.4, color=:firebrick,
      label="global fidelity  F(ρ,σ)")
plot!(pa, Ls, fsv./fsv[1], marker=:rect,    ms=6, lw=2.4, color=:steelblue,
      label="fidelity on the support")
plot!(pa, Ls, gv./gv[1],   marker=:diamond, ms=6, lw=2.4, ls=:dash, color=:seagreen,
      label="gain G  (onsite ZZ)")
hline!(pa, [1.0], color=:black, ls=:dot, lw=1.2, label="")
annotate!(pa, 14.6, 0.30, text(@sprintf("%.3f → %.3f\n(÷%.1f)", gfv[1], gfv[end], gfv[1]/gfv[end]),
                               7, :right, :firebrick))
annotate!(pa, 11.6, 1.55, text(@sprintf("F support: %.4f at every L\ngain: %.0f → %.0f",
                                        fsv[1], gv[1], gv[end]), 7, :left, :steelblue))

# ---------- (b) U 走査: 電荷 vs スピン ----------
Us = [2.0,4.0,8.0,12.0]
pu(o,q) = [ (m = (W.==1) .& (LX.==12) .& (geo.=="cylinder") .& (U.==u) .& (pri.=="UHF") .& (obs.==o);
             q[findfirst(m)]) for u in Us ]
pb = plot(xscale=:log2, yscale=:log10, xlabel="interaction U / t",
          ylabel="relative error  ε = |Δ| / |⟨P⟩|",
          title="(b) charge and spin observables move in opposite directions",
          legend=:left, xticks=(Us, ["2","4","8","12"]), ylims=(8e-3, 3.0))
for (o,lab,cl,mk) in (("ZZ onsite","ZZ onsite = double occ. (charge)",:darkorange,:circle),)
    plot!(pb, Us, pu(o,eps), marker=mk, ms=6, lw=2.4, color=cl, label=lab)
end
for (o,lab,cl,mk) in (("ZZ up-up nb","ZZ up-up (spin)",:steelblue,:utriangle),
                      ("SzSz nb","SzSz (spin)",:mediumpurple,:diamond),
                      ("SxSx nb","SxSx (spin)",:seagreen,:star5))
    plot!(pb, Us, pu(o,eps), marker=mk, ms=6, lw=2.4, ls=:dash, color=cl, label=lab)
end
hline!(pb, [1.0], color=:black, ls=:dot, lw=1.5, label="ε = 1 (break-even)")
annotate!(pb, 8.4, 1.4e-2, text("÷34", 8, :left, :darkorange))
annotate!(pb, 3.0, 1.35, text("crosses into loss", 7, :left, :steelblue))

# ---------- (c) 全点: 大域忠実度と利得 ----------
pc = plot(xscale=:log10, yscale=:log10, xlabel="global fidelity  F(ρ, σ)", ylabel="gain G",
          title=@sprintf("(c) all %d exact single-Pauli points: F does not determine G", sum(ok)),
          legend=:bottomleft, xlims=(3e-5, 2.0), ylims=(0.3, 3e3))
mps = ok .& (pri .!= "UHF"); uhf = ok .& (pri .== "UHF")
scatter!(pc, max.(gf[mps],3e-5), max.(G[mps],0.3), ms=2.6, mc=:steelblue, msw=0,
         alpha=0.38, label="MPS priors (χ_p = 2–256)")
scatter!(pc, max.(gf[uhf],3e-5), max.(G[uhf],0.3), ms=4.2, mc=:firebrick, msw=0.3,
         alpha=0.8, marker=:star5, label="UHF (mean field)")
hline!(pc, [1.0], color=:black, ls=:dot, lw=1.5, label="G = 1")
annotate!(pc, 2e-4, 1.1e3, text("same G, F spread over\nfour decades →", 7, :left, :black))

fig = plot(pa, pb, pc, layout=(1,3), size=(1680,470), margin=7Plots.mm, bottom_margin=9Plots.mm)
savefig(fig, joinpath(DIR,"crm_new_fig_edfid.png"))
println("保存: crm_new_fig_edfid.png")
