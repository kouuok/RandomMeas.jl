# 128量子ビットまで厳密な ρ を使って、系サイズ依存を測り直す。
#   (a) 大域忠実度は最大2万倍崩れ、台の忠実度と利得は動かない
#   (b) 崩壊が激しい U ほど局所量は安定する
#   (c) 1次元では平均場は MPS に追いつかない — L を8倍にしても差が一定
using DelimitedFiles, Printf, Statistics, Plots
gr(fontfamily="Helvetica", legendfontsize=7, guidefontsize=9, titlefontsize=10, tickfontsize=8)
const DIR = @__DIR__
A, hh = readdlm(joinpath(DIR,"crm_edfid_all.tsv"), '\t'; header=true)
hh = vec(hh); c(n) = findfirst(==(n), hh)
S(n) = String.(A[:, c(n)]); F(n) = Float64.(A[:, c(n)])
num(v) = v isa Number ? Float64(v) : NaN
dim, geo, pri, obs, exa = S("dim"), S("geometry"), S("prior"), S("observable"), S("exact")
W, LX, U = F("W"), F("LX"), F("U")
gf, fs = F("global_fid"), F("F_supp")
G = num.(A[:, c("G")])
base = (exa .== "yes") .& (dim .== "1D") .& (geo .== "cylinder") .& (obs .== "ZZ onsite")
Ls = [8,10,12,14,16,24,32,48,64]
US = [2.0,4.0,8.0,12.0]
COL = Dict(2.0=>:steelblue, 4.0=>:seagreen, 8.0=>:goldenrod, 12.0=>:firebrick)
MK  = Dict(2.0=>:circle, 4.0=>:rect, 8.0=>:utriangle, 12.0=>:diamond)
pick(u,p,L,q) = (m = base .& (U.==u) .& (pri.==p) .& (LX.==L); any(m) ? q[findfirst(m)] : NaN)

# ---- (a) 大域忠実度 vs 台の忠実度 ----
pa = plot(yscale=:log10, xscale=:log2, xlabel="chain length L   (2L qubits)",
          ylabel="fidelity", title="(a) exact ρ up to 128 qubits: only F collapses",
          legend=:bottomleft, xticks=(Ls, string.(Ls)), ylims=(3e-6, 3))
for u in US
    g = [pick(u,"UHF",L,gf) for L in Ls]; s = [pick(u,"UHF",L,fs) for L in Ls]
    plot!(pa, Ls, g, marker=MK[u], ms=5, lw=2, color=COL[u], label=@sprintf("global F,  U=%d", Int(u)))
    plot!(pa, Ls, s, marker=MK[u], ms=4, lw=1.6, ls=:dash, color=COL[u], alpha=0.75,
          label=@sprintf("F on support,  U=%d", Int(u)))
end
annotate!(pa, 26, 1.1e-5, text("solid: global fidelity (falls 2e4x)\ndashed: support fidelity (flat)", 7, :left))

# ---- (b) 崩壊率 と 局所量の安定性 ----
drop = [pick(u,"UHF",Ls[1],gf)/pick(u,"UHF",Ls[end],gf) for u in US]
fvar = [begin s=[pick(u,"UHF",L,fs) for L in Ls]; maximum(s)-minimum(s) end for u in US]
pb = plot(xscale=:log10, yscale=:log10, xlabel="collapse of the global fidelity, L = 8 → 64",
          ylabel="drift of the support fidelity", legend=:topright,
          title="(b) the harder F falls, the stiller the local quantity")
for (i,u) in enumerate(US)
    scatter!(pb, [drop[i]], [fvar[i]], ms=9, mc=COL[u], marker=MK[u], msw=0.4,
             label=@sprintf("U = %d", Int(u)))
end
for (i,u) in enumerate(US)
    annotate!(pb, drop[i]*1.35, fvar[i], text(@sprintf("%.0fx / %.0e", drop[i], fvar[i]), 6, :left))
end
let x = log.(drop), y = log.(fvar)
    b = sum((x .- mean(x)).*(y .- mean(y)))/sum((x .- mean(x)).^2)
    a = mean(y) - b*mean(x)
    xs = exp.(range(log(minimum(drop))*0.95, log(maximum(drop))*1.02, length=50))
    plot!(pb, xs, exp.(a) .* xs.^b, lw=1.6, ls=:dash, color=:gray45,
          label=@sprintf("slope %.2f", b))
end

# ---- (c) 1次元では平均場は追いつかない ----
pc = plot(xscale=:log2, xlabel="chain length L", ylabel="G(UHF) / G(χ_p = 32)",
          title="(c) in 1D the gap does not close, even at 8x the length",
          legend=:right, xticks=(Ls, string.(Ls)), ylims=(0.15, 1.18))
for u in US
    r = [pick(u,"UHF",L,G)/pick(u,"chi32",L,G) for L in Ls]
    plot!(pc, Ls, r, marker=MK[u], ms=5, lw=2, color=COL[u], label=@sprintf("U = %d", Int(u)))
end
hline!(pc, [1.0], color=:black, ls=:dash, lw=1.3, label="mean field = χ_p=32")
annotate!(pc, 9.0, 0.36, text("2D is different: there the mean field\nwins (W=6: 157 vs 96)", 7, :left, :black))

fig = plot(pa, pb, pc, layout=(1,3), size=(1680,470), margin=7Plots.mm, bottom_margin=9Plots.mm)
savefig(fig, joinpath(DIR,"crm_new_fig_edfid_size.png"))
@printf("保存: 崩壊率 %s\n", join([@sprintf("U=%d:%.0fx", Int(u), d) for (u,d) in zip(US,drop)], " "))
