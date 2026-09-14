# §10 の主要結果の PBC(周期境界/トーラス)版。OBC 版 crm_edfid_figure.jl と同じ3枚を
# そのまま並べられるように作り、境界条件で何が変わるかの4枚目を足す。
#   (a1) 系を伸ばすと大域忠実度だけが崩れ、台の忠実度は動かない(実数値)← PBC でも成立
#   (a2) 同じ prior の利得も動かない(実数値)                         ← PBC でも成立
#   (b) U を上げると電荷量とスピン量が逆を向く                        ← PBC では損に入らない
#   (c) PBC 全点: 大域忠実度は利得を決めていない                      ← PBC でも成立
#   (d) (a)(b) の差の出どころ: 開放鎖の中央ボンドの偶奇による交替
#   julia --project=. examples/tatsuyama/crm_edfid_pbc_figure.jl
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

pbc = (exa .== "yes") .& (geo .== "torus")
obc = (exa .== "yes") .& (geo .== "cylinder")
okP = pbc .& isfinite.(Gm)
@printf("PBC 厳密点 %d / うち単一パウリ %d\n", sum(pbc), sum(okP))

# 1D の (幾何, L, U, prior, 観測量) を1点引く
pick(g, L, u, o, q; p="UHF") = begin
    m = (W.==1) .& (LX.==L) .& (geo.==g) .& (U.==u) .& (pri.==p) .& (obs.==o)
    any(m) ? q[findfirst(m)] : NaN
end

# ---------- (a) サイズ走査(PBC) ----------
# 比(L=8 基準)ではなく実数値で描く。忠実度(≲1)と利得(≈190)は単位も桁も違うので
# 上下2枚に分ける。二重軸は2つの目盛りの合わせ方しだいで見かけの相関を作るので使わない。
Ls = [8,10,12,14,16]
gfv = [pick("torus",L,8.0,"ZZ onsite",gf) for L in Ls]
fsv = [pick("torus",L,8.0,"ZZ onsite",fs) for L in Ls]
gv  = [pick("torus",L,8.0,"ZZ onsite",G)  for L in Ls]
gfo = [pick("cylinder",L,8.0,"ZZ onsite",gf) for L in Ls]
@printf("(a) 大域F %s\n    台F  %s  (幅 %.1e)\n    G    %s\n",
        join((@sprintf("%.4f",x) for x in gfv), " "), join((@sprintf("%.5f",x) for x in fsv), " "),
        maximum(fsv)-minimum(fsv), join((@sprintf("%.2f",x) for x in gv), " "))
const BLUE = "#1F77B4"   # steelblue は彩度が足りずパレット検証に落ちる
const INK  = :gray20     # 数値ラベルは系列色でなく文字色
const XL   = (6.3, 17.7) # 両端に数値ラベルを置く余白。上下で x 軸をそろえる
spread(v) = replace(@sprintf("%.0e", maximum(v)-minimum(v)), "e-0" => "e-")
sig3(x) = x >= 0.1 ? @sprintf("%.3f", x) : @sprintf("%.4f", x)   # 有効数字3桁
pa1 = plot(ylabel="fidelity", yscale=:log10, ylims=(0.01, 1.0), xlims=XL,
           yticks=([0.01,0.02,0.05,0.2,0.5,1.0], ["0.01","0.02","0.05","0.2","0.5","1"]),  # "0.1" は GR で点が潰れる
           xticks=(Ls, fill("", length(Ls))), legend=:bottomleft,
           title="(a1) PBC fidelity: only the global one falls with L")
plot!(pa1, Ls, gfo, lw=1.6, ls=:dashdot, color=:gray55, label="open chain, global F (reference)")
plot!(pa1, Ls, fsv, marker=:rect,   ms=5, lw=2, color=BLUE,       label="fidelity on the support  F_A")
plot!(pa1, Ls, gfv, marker=:circle, ms=5, lw=2, color=:firebrick, label="global fidelity  F(ρ,σ)")
for (v, f) in ((gfv, sig3), (fsv, x -> @sprintf("%.4f", x)))
    annotate!(pa1, 7.75,  v[1],   text(f(v[1]),   7, :right, INK))
    annotate!(pa1, 16.25, v[end], text(f(v[end]), 7, :left,  INK))
end
annotate!(pa1, 12.0, 0.30, text("F_A moves by only " * spread(fsv) * " across L", 7, :center, INK))
pa2 = plot(xlabel="ring length L   (2L qubits)", ylabel="gain G  (onsite ZZ)",
           ylims=(0, 250), xlims=XL, xticks=Ls, legend=false,
           title="(a2) PBC gain for the same prior: flat")
plot!(pa2, Ls, gv, marker=:diamond, ms=6, lw=2, color=:seagreen)
annotate!(pa2, 7.75,  gv[1],   text(@sprintf("%.1f", gv[1]),   7, :right, INK))
annotate!(pa2, 16.25, gv[end], text(@sprintf("%.1f", gv[end]), 7, :left,  INK))

# ---------- (b) U 走査(PBC): 電荷 vs スピン ----------
Us = [2.0,4.0,8.0,12.0]
pb = plot(xscale=:log2, yscale=:log10, xlabel="interaction U / t",
          ylabel="relative error  ε = |Δ| / |⟨P⟩ − c_0|",
          title="(b) PBC: the spin observables never cross into loss",
          legend=(0.56,0.60), xticks=(Us, ["2","4","8","12"]), ylims=(6e-3, 3.5))
for (i,(o,cl)) in enumerate((("ZZ up-up nb",:steelblue), ("SzSz nb",:mediumpurple)))
    plot!(pb, Us, [pick("cylinder",12,u,o,eps) for u in Us], lw=1.8, ls=:dashdot,
          color=cl, alpha=0.75, label=(i==1 ? "same observable, open chain" : ""))
end
for (o,lab,cl,mk) in (("ZZ onsite","ZZ onsite = double occ. (charge)",:darkorange,:circle),)
    plot!(pb, Us, [pick("torus",12,u,o,eps) for u in Us], marker=mk, ms=6, lw=2.4, color=cl, label=lab)
end
for (o,lab,cl,mk) in (("ZZ up-up nb","ZZ up-up (spin)",:steelblue,:utriangle),
                      ("SzSz nb","SzSz (spin)",:mediumpurple,:diamond),
                      ("SxSx nb","SxSx (spin)",:seagreen,:star5))
    plot!(pb, Us, [pick("torus",12,u,o,eps) for u in Us], marker=mk, ms=6, lw=2.4, ls=:dash,
          color=cl, label=lab)
end
hline!(pb, [1.0], color=:black, ls=:dot, lw=1.5, label="ε = 1 (break-even)")
annotate!(pb, 2.1, 1.55, text("open chain crosses here →", 7, :left, :gray35))
annotate!(pb, 8.6, 0.45, text("ring stays below 1\nat every U", 7, :left, :steelblue))

# ---------- (c) PBC 全点: 大域忠実度と利得 ----------
pc = plot(xscale=:log10, yscale=:log10, xlabel="global fidelity  F(ρ, σ)", ylabel="gain G",
          title=@sprintf("(c) all %d exact PBC single-Pauli points: F does not determine G", sum(okP)),
          legend=:bottomleft, xlims=(1e-4, 3.0), ylims=(0.3, 3e3))
mps = okP .& (pri .!= "UHF"); uhf = okP .& (pri .== "UHF")
scatter!(pc, max.(gf[mps],1.3e-4), max.(G[mps],0.3), ms=2.6, mc=:steelblue, msw=0,
         alpha=0.38, label="MPS priors (χ_p = 2–256)")
scatter!(pc, max.(gf[uhf],1.3e-4), max.(G[uhf],0.3), ms=4.2, mc=:firebrick, msw=0.3,
         alpha=0.8, marker=:star5, label="UHF (mean field)")
hline!(pc, [1.0], color=:black, ls=:dot, lw=1.5, label="G = 1")
# 同じ G の帯のなかで F がどれだけ広がるかを、外れ値に頼らず 5–95 パーセンタイルで測る
let b = Dict{Int,Vector{Float64}}()
    for i in eachindex(okP)
        okP[i] && G[i] > 0 || continue
        push!(get!(b, round(Int, log10(G[i])*4), Float64[]), gf[i])
    end
    q(v,p) = (w = sort(v); k = (length(w)-1)*p; lo = floor(Int,k)+1;
              hi = min(lo+1,length(w)); w[lo] + (w[hi]-w[lo])*(k-floor(k)))
    best = (0.0, 0.0, 0, 0.0, 0.0)
    for (k,v) in b
        length(v) >= 15 || continue
        d = log10(q(v,0.95)/q(v,0.05))
        d > best[1] && (best = (d, 10.0^(k/4), length(v), q(v,0.05), q(v,0.95)))
    end
    @printf("(c) 最大の同G帯: G≈%.0f (%d点) で F の 5–95%%tile = [%.2e, %.3f] → %.1f 桁\n",
            best[2], best[3], best[4], best[5], best[1])
    annotate!(pc, 1.5e-4, 1.1e3,
        text(@sprintf("at G ≈ %.0f (%d pts), F spans %.1f decades\n(5–95%%: %.0e to %.2f)",
                      best[2], best[3], best[1], best[4], best[5]), 7, :left, :black))
end
annotate!(pc, 1.5e-4, 3.2e2, text("2 UHF pts clamped to the left edge:\nF = 4e−15 (4×3 torus, U = 4 — anomalous,\nsee note in README)", 6, :left, :gray40))

# ---------- (d) 差の出どころ: 開放鎖の中央ボンドの偶奇 ----------
Lo = [8,10,12,14,16,24,32,48,64]
pd = plot(xscale=:log2, xlabel="chain / ring length L", ylabel="relative error  ε  (U = 8, UHF)",
          title="(d) why (b) differs: the open chain's central bond alternates",
          legend=:topright, xticks=(Lo, string.(Lo)), ylims=(0.0, 1.80))
for (o,cl,mk) in (("SzSz nb",:mediumpurple,:diamond), ("ZZ up-up nb",:steelblue,:utriangle))
    eo = [pick("cylinder",L,8.0,o,eps) for L in Lo]
    et = [pick("torus",L,8.0,o,eps)    for L in Lo]
    plot!(pd, Lo, eo, marker=mk, ms=5, lw=2.2, color=cl, label="open chain, $o")
    plot!(pd, Lo, et, marker=mk, ms=5, lw=2.2, ls=:dash, color=cl, alpha=0.55,
          label="ring, $o")
end
hline!(pd, [1.0], color=:black, ls=:dot, lw=1.5, label="ε = 1 (break-even)")
annotate!(pd, 10.0, 0.18, text("L ≡ 2 (mod 4):\ncentral bond is odd", 7, :center, :gray25))
annotate!(pd, 9.6, 1.68, text("L ≡ 0 (mod 4): even", 7, :left, :gray25))
annotate!(pd, 21.0, 0.38, text("ring: every bond equivalent\n→ smooth, never crosses 1", 7, :left, :gray25))

# 上段は (a1)(a2) を縦に積むので高さを多めに配分する
top = plot(pa1, pa2, pb, layout=@layout([[a1; a2] b]))
bot = plot(pc, pd, layout=(1,2))
fig = plot(top, bot, layout=@layout([t{0.58h}; u]), size=(1300,1250),
           margin=6Plots.mm, bottom_margin=7Plots.mm)
savefig(fig, joinpath(DIR,"crm_new_fig_edfid_pbc.png"))
println("保存: crm_new_fig_edfid_pbc.png")
