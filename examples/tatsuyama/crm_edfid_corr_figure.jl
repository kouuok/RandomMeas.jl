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

# ---------- (b) prior 群で符号が変わる ----------
OB2 = ["ZZ onsite","ZZ up-up nb"]; LAB2 = ["ZZ onsite","ZZ up-up"]
grp = [("MPS priors", x -> x != "UHF"), ("UHF", x -> x == "UHF")]
vals = Dict{String,Vector{Float64}}()
for (gname, f2) in grp
    v = Float64[]
    for o in OB2
        m = base .& (obs .== o) .& isfinite.(G) .& [f2(p) for p in pri]
        push!(v, spear(fs[m], G[m]))
    end
    vals[gname] = v
end
pb = plot(xticks=(1:2, LAB2), ylabel="Spearman ρ  (support fidelity vs G)",
          title="(b) the sign flips with the kind of prior", ylims=(-1.15,1.0), legend=:bottomleft)
bar!(pb, (1:2) .- w/2, vals["MPS priors"], bar_width=w*0.9, color=:steelblue, lw=0, label="MPS priors")
bar!(pb, (1:2) .+ w/2, vals["UHF"],        bar_width=w*0.9, color=:firebrick, lw=0, label="UHF (mean field)")
hline!(pb, [0.0], color=:black, lw=1.2, label="")
for (i,v) in enumerate(vals["MPS priors"]); annotate!(pb, i-w/2, v+0.07, text(@sprintf("%+.2f",v), 7, :center)); end
for (i,v) in enumerate(vals["UHF"]);        annotate!(pb, i+w/2, v-0.09, text(@sprintf("%+.2f",v), 7, :center)); end

# ---------- (c) UHF: 忠実度が高いほど利得が低い ----------
pc = plot(yscale=:log10, xlabel="fidelity on the support  F_A", ylabel="gain G  (onsite ZZ)",
          title="(c) UHF on the on-site ZZ: the relation runs backwards",
          legend=:topright, ylims=(0.7, 900))
UC = [(2.0,:steelblue,:circle),(4.0,:seagreen,:rect),(8.0,:goldenrod,:utriangle),(12.0,:firebrick,:diamond)]
for (u,cl,mk) in UC
    m = base .& (obs .== "ZZ onsite") .& (pri .== "UHF") .& (U .== u) .& isfinite.(G)
    ρ = spear(fs[m], G[m])
    scatter!(pc, fs[m], G[m], ms=6, mc=cl, marker=mk, msw=0.3, alpha=0.85,
             label=@sprintf("U = %d   ρ = %+.3f", Int(u), ρ))
end
annotate!(pc, 0.60, 1.6, text("higher fidelity, lower gain — at every U.\n(ZZ up-up stays positive: see panel b)", 7, :left, :black))

fig = plot(pa, pb, pc, layout=(1,3), size=(1680,470), margin=7Plots.mm, bottom_margin=10Plots.mm)
savefig(fig, joinpath(DIR,"crm_new_fig_edfid_corr.png"))
println("保存: crm_new_fig_edfid_corr.png")
