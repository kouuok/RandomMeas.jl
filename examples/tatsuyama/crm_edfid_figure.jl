# §10 の主要結果を3枚で示す。
#   (a1) 系を伸ばすと大域忠実度だけが崩れ、台の忠実度は動かない(実数値)
#   (a2) 同じ prior の利得も動かない(実数値)
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
# 比(L=8 基準)ではなく実数値で描く。忠実度(≲1)と利得(≈190)は単位も桁も違うので
# 上下2枚に分ける。二重軸は2つの目盛りの合わせ方しだいで見かけの相関を作るので使わない。
Ls = [8,10,12,14,16]
pick(L, o, q) = begin
    m = (W.==1) .& (LX.==L) .& (geo.=="cylinder") .& (U.==8.0) .& (pri.=="UHF") .& (obs.==o)
    q[findfirst(m)]
end
gfv = [pick(L,"ZZ onsite",gf) for L in Ls]
fsv = [pick(L,"ZZ onsite",fs) for L in Ls]
gv  = [pick(L,"ZZ onsite",G)  for L in Ls]
@printf("(a) 大域F %s\n    台F  %s  (幅 %.1e)\n    G    %s\n",
        join((@sprintf("%.4f",x) for x in gfv), " "), join((@sprintf("%.5f",x) for x in fsv), " "),
        maximum(fsv)-minimum(fsv), join((@sprintf("%.2f",x) for x in gv), " "))
const BLUE = "#1F77B4"   # steelblue は彩度が足りずパレット検証に落ちる
const INK  = :gray20     # 数値ラベルは系列色でなく文字色
const XL   = (6.3, 17.7) # 両端に数値ラベルを置く余白。上下で x 軸をそろえる
spread(v) = replace(@sprintf("%.0e", maximum(v)-minimum(v)), "e-0" => "e-")
pa1 = plot(ylabel="fidelity", yscale=:log10, ylims=(0.01, 1.0), xlims=XL,
           yticks=([0.01,0.02,0.05,0.2,0.5,1.0], ["0.01","0.02","0.05","0.2","0.5","1"]),  # "0.1" は GR で点が潰れて "01" に見える
           xticks=(Ls, fill("", length(Ls))), legend=:bottomleft,
           title="(a1) fidelity: only the global one falls with L")
plot!(pa1, Ls, fsv, marker=:rect,   ms=5, lw=2, color=BLUE,       label="fidelity on the support  F_A")
plot!(pa1, Ls, gfv, marker=:circle, ms=5, lw=2, color=:firebrick, label="global fidelity  F(ρ,σ)")
sig3(x) = x >= 0.1 ? @sprintf("%.3f", x) : @sprintf("%.4f", x)   # 有効数字3桁
for (v, f) in ((gfv, sig3), (fsv, x -> @sprintf("%.4f", x)))
    annotate!(pa1, 7.75,  v[1],   text(f(v[1]),   7, :right, INK))
    annotate!(pa1, 16.25, v[end], text(f(v[end]), 7, :left,  INK))
end
annotate!(pa1, 12.0, 0.30, text("F_A moves by only " * spread(fsv) * " across L", 7, :center, INK))
pa2 = plot(xlabel="chain length L   (2L qubits)", ylabel="gain G  (onsite ZZ)",
           ylims=(0, 250), xlims=XL, xticks=Ls, legend=false,
           title="(a2) gain for the same prior: flat")
plot!(pa2, Ls, gv, marker=:diamond, ms=6, lw=2, color=:seagreen)
annotate!(pa2, 7.75,  gv[1],   text(@sprintf("%.1f", gv[1]),   7, :right, INK))
annotate!(pa2, 16.25, gv[end], text(@sprintf("%.1f", gv[end]), 7, :left,  INK))

# ---------- (b) U 走査: 電荷 vs スピン ----------
Us = [2.0,4.0,8.0,12.0]
pu(o,q) = [ (m = (W.==1) .& (LX.==12) .& (geo.=="cylinder") .& (U.==u) .& (pri.=="UHF") .& (obs.==o);
             q[findfirst(m)]) for u in Us ]
pb = plot(xscale=:log2, yscale=:log10, xlabel="interaction U / t",
          ylabel="relative error  ε = |Δ| / |⟨P⟩|",
          title="(b) charge and spin observables move in opposite directions",
          legend=:right, xticks=(Us, ["2","4","8","12"]), ylims=(8e-3, 3.0))
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

fig = plot(pa1, pa2, pb, pc, layout=@layout([[a1; a2] b c]), size=(1680,680),
           margin=6Plots.mm, bottom_margin=7Plots.mm)
savefig(fig, joinpath(DIR,"crm_new_fig_edfid.png"))
println("保存: crm_new_fig_edfid.png")
