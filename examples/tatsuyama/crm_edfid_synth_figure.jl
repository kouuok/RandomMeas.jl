# 検証で確定した2点を図にする。
#   図A: 忠実度が prior の質を測れるのはどういうときか(系内では測れる/系をまたぐと測れない)
#   図B: 平均場が MPS に勝つのはどういうときか(MPS が張り詰めているとき)
using DelimitedFiles, Printf, Statistics, Plots
gr(fontfamily="Helvetica", legendfontsize=7, guidefontsize=9, titlefontsize=10, tickfontsize=8)
const DIR = @__DIR__
A, hh = readdlm(joinpath(DIR,"crm_edfid_all.tsv"), '\t'; header=true)
hh = vec(hh); c(n) = findfirst(==(n), hh)
S(n) = String.(A[:, c(n)]); F(n) = Float64.(A[:, c(n)])
num(v) = v isa Number ? Float64(v) : NaN
dim, obs, pri, exa, geo = S("dim"), S("observable"), S("prior"), S("exact"), S("geometry")
gf, fs, U, W, LX = F("global_fid"), F("F_supp"), F("U"), F("W"), F("LX")
eps = num.(A[:, c("eps")]); G = num.(A[:, c("G")])
ok = (exa .== "yes") .& (obs .== "ZZ onsite") .& isfinite.(eps) .& (eps .> 0)

function spear(x,y)
    n=length(x); n<5 && return NaN
    rx=zeros(n); ry=zeros(n)
    for (r,i) in enumerate(sortperm(x)); rx[i]=r; end
    for (r,i) in enumerate(sortperm(y)); ry[i]=r; end
    rx.-=mean(rx); ry.-=mean(ry)
    d=sqrt((rx'*rx)*(ry'*ry)); d==0 ? NaN : (rx'*ry)/d
end
key = [(W[i],LX[i],geo[i],U[i]) for i in eachindex(W)]
ukeys = unique(key[ok])

# ---- (a) 系を固定して prior を変える ----
pa = plot(xscale=:log10, xlabel="global fidelity  F(ρ, σ)", ylabel="relative error  ε",
          yscale=:log10, title="(a) fix the system, vary the prior", legend=:topright,
          xlims=(2e-2,1.6), ylims=(1e-7,3))
shown = 0
for k in ukeys
    m = ok .& [key[i]==k for i in eachindex(key)]
    sum(m) < 6 && continue
    global shown += 1
    lab = shown <= 3 ? @sprintf("%dx%d %s, U=%d", Int(k[2]), Int(k[1]), k[3]=="torus" ? "PBC" : "OBC", Int(k[4])) : ""
    plot!(pa, gf[m][sortperm(gf[m])], eps[m][sortperm(gf[m])], marker=:circle, ms=3,
          lw=1.2, alpha=shown<=3 ? 0.95 : 0.16, color=shown<=3 ? shown : :gray55, label=lab)
end
cs = [spear(gf[ok .& [key[i]==k for i in eachindex(key)]], -log.(eps[ok .& [key[i]==k for i in eachindex(key)]])) for k in ukeys]
annotate!(pa, 0.30, 3e-6, text(@sprintf("each system: ρ = %+.2f (median over %d)", median(filter(isfinite,cs)), length(ukeys)), 7, :left))

# ---- (b) prior を固定して系を変える ----
pb = plot(xscale=:log10, yscale=:log10, xlabel="global fidelity  F(ρ, σ)", ylabel="relative error  ε",
          title="(b) fix the prior AND U, vary only the system", legend=:topright,
          xlims=(2e-2,1.6), ylims=(1e-7,3))
for (p,cl,mk,lab) in (("UHF",:firebrick,:star5,"UHF"),("chi8",:steelblue,:circle,"χ_p = 8"),
                      ("chi32",:seagreen,:rect,"χ_p = 32"))
    m = ok .& (pri .== p) .& (U .== 8.0)
    ρ = spear(gf[m], -log.(eps[m]))
    scatter!(pb, gf[m], eps[m], ms=4, mc=cl, marker=mk, msw=0.2, alpha=0.75,
             label=@sprintf("%s   ρ = %+.2f", lab, ρ))
end
annotate!(pb, 2.6e-2, 3e-7, text("at U = 8, across 16 systems:\nUHF\u0027s F moves 6x, its ε barely moves", 7, :left, :firebrick))

# ---- (c) 平均場が勝つのは MPS が張り詰めているとき ----
rat = Float64[]; strain = Float64[]; is2d = Bool[]
for k in ukeys
    m = ok .& [key[i]==k for i in eachindex(key)]
    iu = findfirst(m .& (pri .== "UHF")); ic = findfirst(m .& (pri .== "chi32"))
    (iu === nothing || ic === nothing) && continue
    push!(rat, G[iu]/G[ic]); push!(strain, gf[ic]); push!(is2d, dim[iu]=="2D")
end
pc = plot(xlabel="global fidelity of the χ_p = 32 MPS prior", ylabel="G(UHF) / G(χ_p = 32)",
          title="(c) mean field wins where the MPS is strained", legend=:topright,
          yscale=:log10, ylims=(0.3, 2.2), xlims=(0.25,1.03),
          yticks=([0.3,0.5,0.7,1.0,1.5,2.0], ["0.3","0.5","0.7","1.0","1.5","2.0"]))
scatter!(pc, strain[.!is2d], rat[.!is2d], ms=6, mc=:steelblue, marker=:circle, msw=0.2,
         alpha=0.8, label="1D chains")
scatter!(pc, strain[is2d], rat[is2d], ms=7, mc=:firebrick, marker=:rect, msw=0.2,
         alpha=0.85, label="2D lattices")
hline!(pc, [1.0], color=:black, ls=:dash, lw=1.3, label="UHF = χ_p=32")
annotate!(pc, 0.28, 1.75, text(@sprintf("UHF wins in %d of %d systems,\n7 of them 2D", sum(rat.>1), length(rat)), 7, :left, :black))

fig = plot(pa, pb, pc, layout=(1,3), size=(1680,470), margin=7Plots.mm, bottom_margin=9Plots.mm)
savefig(fig, joinpath(DIR,"crm_new_fig_edfid_synth.png"))
@printf("保存 (系内 ρ 中央値 %.3f, UHF 勝ち %d/%d)\n", median(filter(isfinite,cs)), sum(rat.>1), length(rat))
