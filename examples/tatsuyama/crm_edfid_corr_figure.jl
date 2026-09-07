# 観測量ごとに「忠実度は利得を予測するか」を調べる。
#   (a) 観測量ごと: 忠実度 vs prior の相対誤差 ε  — よく効く
#   (b) 同じ観測量でも prior 群で符号が変わる
#   (c) UHF では台の忠実度が利得を逆向きに予測する
#   julia --project=. examples/tatsuyama/crm_edfid_corr_figure.jl
using DelimitedFiles, Printf, Statistics, Plots
gr(fontfamily="Helvetica", legendfontsize=7, guidefontsize=9,
   titlefontsize=10, tickfontsize=8)
const DIR = @__DIR__
A, hh = readdlm(joinpath(DIR,"crm_edfid_all.tsv"), '\t'; header=true)
hh = vec(hh); c(n) = findfirst(==(n), hh)
S(n) = String.(A[:, c(n)]); F(n) = Float64.(A[:, c(n)])
obs, pri, exa = S("observable"), S("prior"), S("exact")
gf, fs, U = F("global_fid"), F("F_supp"), F("U")
epsr = A[:, c("eps")]; Gr = A[:, c("G")]
num(v) = v isa Number ? Float64(v) : NaN
eps = num.(epsr); G = num.(Gr)
base = exa .== "yes"

function spear(x, y)
    n = length(x); n < 5 && return NaN
    rx = zeros(n); ry = zeros(n)
    for (r,i) in enumerate(sortperm(x)); rx[i] = r; end
    for (r,i) in enumerate(sortperm(y)); ry[i] = r; end
    rx .-= mean(rx); ry .-= mean(ry)
    d = sqrt((rx'*rx)*(ry'*ry)); d == 0 ? NaN : (rx'*ry)/d
end

OB   = ["ZZ onsite","ZZ up-up nb","SzSz nb","SxSx nb","DoubleOcc"]
LAB  = ["ZZ onsite","ZZ up-up","SzSz","SxSx","double occ."]

# ---------- (a) 観測量ごと: 忠実度 vs -log ε ----------
ca_g = Float64[]; ca_s = Float64[]
for o in OB
    m = base .& (obs .== o) .& isfinite.(eps) .& (eps .> 0)
    t = -log.(eps[m])
    push!(ca_g, spear(gf[m], t)); push!(ca_s, spear(fs[m], t))
end
pa = plot(xticks=(1:length(OB), LAB), xrotation=18, ylabel="Spearman ρ",
          title="(a) within one observable, fidelity does track ε", ylims=(0,1.05), legend=:bottomright)
w = 0.8/2
bar!(pa, (1:length(OB)) .- w/2, ca_g, bar_width=w*0.9, color=:firebrick, lw=0, label="global fidelity")
bar!(pa, (1:length(OB)) .+ w/2, ca_s, bar_width=w*0.9, color=:steelblue, lw=0, label="fidelity on the support")
for (i,v) in enumerate(ca_g); annotate!(pa, i-w/2, v+0.045, text(@sprintf("%.2f",v), 6, :center)); end
for (i,v) in enumerate(ca_s); annotate!(pa, i+w/2, v+0.045, text(@sprintf("%.2f",v), 6, :center)); end

# ---------- (b) 見かけの逆相関は天井の効果 ----------
Us = [2.0,4.0,8.0,12.0]
raw = Float64[]; norm = Float64[]
for u in Us
    m = base .& (obs .== "ZZ onsite") .& (pri .== "UHF") .& (U .== u) .& isfinite.(G)
    push!(raw,  spear(fs[m], G[m]))
    push!(norm, spear(fs[m], G[m] ./ F("G_max")[m]))
end
pb = plot(xticks=(1:4, ["2","4","8","12"]), xlabel="interaction U / t", ylabel="Spearman ρ",
          title="(b) the apparent inversion is the ceiling, not the prior",
          ylims=(-1.15,1.0), legend=:bottomright)
w = 0.8/2
bar!(pb, (1:4) .- w/2, raw,  bar_width=w*0.9, color=:firebrick, lw=0, label="F_A  vs  G")
bar!(pb, (1:4) .+ w/2, norm, bar_width=w*0.9, color=:steelblue, lw=0, label="F_A  vs  G / G_max")
hline!(pb, [0.0], color=:black, lw=1.2, label="")
for (i,v) in enumerate(raw);  annotate!(pb, i-w/2, v-0.09, text(@sprintf("%+.2f",v), 6, :center)); end
for (i,v) in enumerate(norm); annotate!(pb, i+w/2, v+(v<0 ? -0.09 : 0.06), text(@sprintf("%+.2f",v), 6, :center)); end
annotate!(pb, 2.5, 0.72, text("dividing out the ceiling\nflips the sign at U = 8, 12", 7, :center, :black))

# ---------- (c) 1D では F_A がほとんど動いていない ----------
pc = plot(xlabel="fidelity on the support  F_A", ylabel="gain G  (onsite ZZ)",
          title="(c) UHF at U = 8: the 1D points span 1e-4 in F_A",
          legend=:bottomleft, ylims=(120, 245), xlims=(0.554, 0.605),
          yticks=(120:20:240))
for (d,cl,mk,lab) in (("1D",:steelblue,:circle,"1D chains (10 systems)"),
                      ("2D",:firebrick,:rect,"2D lattices (6 systems)"))
    m = base .& (obs .== "ZZ onsite") .& (pri .== "UHF") .& (U .== 8.0) .& (S("dim") .== d)
    scatter!(pc, fs[m], G[m], ms=7, mc=cl, marker=mk, msw=0.3, alpha=0.85, label=lab)
end
annotate!(pc, 0.5588, 205, text("all ten 1D systems\nsit in this column:\nF_A = 0.55768 … 0.55780", 7, :left, :steelblue))
annotate!(pc, 0.5735, 238, text("2D: higher F_A (smaller UHF moment)\nbut lower G (|⟨ZZ⟩| sets a lower ceiling)\n— two unrelated trends, not a relation", 7, :left, :black))

fig = plot(pa, pb, pc, layout=(1,3), size=(1680,470), margin=7Plots.mm, bottom_margin=10Plots.mm)
savefig(fig, joinpath(DIR,"crm_new_fig_edfid_corr.png"))
println("保存: crm_new_fig_edfid_corr.png")
