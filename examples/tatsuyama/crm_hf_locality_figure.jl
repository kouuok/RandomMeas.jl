# 「大域忠実度が低くても、局所観測量なら平均場(UHF)が prior として働く」を3枚で示す。
#   (a) 同じ系で大域忠実度と利得を並べる — UHF だけ傾向から外れる
#   (b) その UHF を観測量ごとに分解する — どこで効き、どこで損をするか
#   (c) 系を大きくすると、タダの UHF が χ=32 の MPS を追い越す
#   julia --project=. examples/tatsuyama/crm_hf_locality_figure.jl
using DelimitedFiles, Statistics, Printf, Plots
gr(fontfamily="Helvetica", legendfontsize=7, guidefontsize=9,
   titlefontsize=10, tickfontsize=8)
const DIR = @__DIR__
rd(f) = readdlm(joinpath(DIR,f), '\t'; header=true)

# ---------- (a) 1D L=8 (16量子ビット), U=8 半充填, ED相当のDMRG参照 ----------
A, ha = rd("crm_fid_vs_gain_L8_U8.0.tsv")
ha = vec(ha); col(n) = findfirst(==(n), ha)
pri = String.(A[:, col("prior")]); obs = String.(A[:, col("observable")])
fid = Float64.(A[:, col("global_fid")]); G = Float64.(A[:, col("G")])

obs_u = unique(obs)
mk = [:circle, :rect, :utriangle, :diamond, :star5, :hexagon]
pa = plot(xscale=:log10, yscale=:log10,
          xlabel="global fidelity  F(ρ, σ)", ylabel="gain G",
          title="(a) same system: global F vs gain   (1D L=8, U=8, half filling)",
          legend=:bottomright, legendfontsize=6, xlims=(0.09, 3.2), ylims=(2e-2, 6e2),
          xticks=([0.1,0.2,0.5,1.0], ["0.1","0.2","0.5","1.0"]))
for (j,o) in enumerate(obs_u)
    m = (obs .== o) .& (pri .!= "UHF")
    scatter!(pa, fid[m], max.(G[m],3e-2), label=o,
             marker=mk[mod1(j,length(mk))], markersize=5, color=j, markerstrokewidth=0.4)
end
mu = pri .== "UHF"
scatter!(pa, fid[mu], max.(G[mu],3e-2), label="UHF (mean field)",
         marker=:star6, markersize=10, color=:black, markerstrokewidth=1.2)
hline!(pa, [1.0], color=:gray, linestyle=:dash, label="G = 1 (no gain)")
vline!(pa, [fid[findfirst(mu)]], color=:black, linestyle=:dot, label="")
annotate!(pa, fid[findfirst(mu)]*1.1, 3.5e2,
          text(@sprintf("UHF\nF = %.3f", fid[findfirst(mu)]), 7, :left, :black))
annotate!(pa, 0.30, 2.2e2, text("same prior,\n4 decades of G", 7, :left, :black))

# ---------- (b) §1 (8量子ビット ED, プラケット): UHF を観測量ごとに ----------
B, hb = rd("crm_gain_verification_results.tsv")
hb = vec(hb); cb(n) = findfirst(==(n), hb)
m = (String.(B[:, cb("config")]) .== "uniform") .& (Float64.(B[:, cb("U")]) .== 8.0)
ob = String.(B[m, cb("observable")]); gb = Float64.(B[m, cb("G_emp")])
nmv = Int.(B[m, cb("nm")]); Fhf = Float64.(B[findfirst(m), cb("HF_fid")])
ord = ["ZZ(1u,1d) |A|=2","ZZZZ(1,3) |A|=4","DoubleOcc site1","Fidelity","Sz1*Sz3","Z(1u)  |A|=1"]
short = ["ZZ onsite","ZZZZ","double occ.","fidelity","SzSz","single Z"]
nms = sort(unique(nmv))
mat = [maximum(gb[(ob .== o) .& (nmv .== n)]) for o in ord, n in nms]
pb = plot(ylabel="log10 gain G", legend=:topright,
          title=@sprintf("(b) one prior, six observables   (UHF, global F = %.3f)", Fhf),
          xticks=(1:length(ord), short), xrotation=20, ylims=(-3.2, 3.6))
w = 0.8/length(nms)
for (j,n) in enumerate(nms)
    bar!(pb, (1:length(ord)) .+ (j-(length(nms)+1)/2)*w, log10.(max.(mat[:,j],1e-3)),
         bar_width=w*0.9, label="n_m = $n", color=j, linewidth=0)
end
hline!(pb, [0.0], color=:black, linestyle=:dash, label="G = 1")
annotate!(pb, 1.0, 3.2, text("CRM wins", 8, :left, :darkgreen))
annotate!(pb, 4.6, -2.6, text("CRM loses", 8, :left, :darkred))

# ---------- (c) 系を大きくすると UHF が χ=32 を追い越す ----------
function zz(file, U)
    D, h = rd(file); h = vec(h); c(n) = findfirst(==(n), h)
    m = (Float64.(D[:, c("U")]) .== U) .& (String.(D[:, c("observable")]) .== "ZZ onsite")
    Dict(String.(D[m, c("prior")]) .=> Float64.(D[m, c("G_emp")])),
    Dict(String.(D[m, c("prior")]) .=> D[m, c("prior_fid")])
end
g4, f4 = zz("crm_2d_results.tsv", 8.0)
g6, f6 = zz("crm_2d_w6_results_U8p0.tsv", 8.0)
labs = ["W=4 cylinder\n(64 qubits)", "W=6 cylinder\n(96 qubits)"]
series = ["chi=8","chi=16","chi=32","chi=64","chi=128","UHF"]
mat2 = [get(g4,s,NaN) for s in series]'
mat2 = vcat(mat2, [get(g6,s,NaN) for s in series]')
pc = plot(ylabel="gain G  (onsite ZZ)", legend=:topright,
          title="(c) scaling up: the free mean-field prior overtakes MPS",
          legendfontsize=6,
          xticks=(1:2, labs), ylims=(0, 290))
w = 0.8/length(series)
cols = [:steelblue, :seagreen, :orange, :firebrick, :purple, :black]
for j in 1:length(series)
    vals = [isnan(mat2[i,j]) ? 0.0 : mat2[i,j] for i in 1:2]
    bar!(pc, (1:2) .+ (j-(length(series)+1)/2)*w, vals, bar_width=w*0.9,
         label=series[j], color=cols[j], linewidth=0)
end
for (i,fd) in enumerate((f4,f6)), (j,s) in enumerate(series)
    v = get(fd, s, nothing); v === nothing && continue
    fv = v isa AbstractString ? NaN : Float64(v)
    (isnan(fv) || isnan(mat2[i,j]) || mat2[i,j] == 0) && continue
    annotate!(pc, i + (j-(length(series)+1)/2)*w, mat2[i,j] + 11,
              text(@sprintf("F=%.2f", fv), 5, :center))
end

fig = plot(pa, pb, pc, layout=(1,3), size=(1650, 480), margin=6Plots.mm,
           bottom_margin=10Plots.mm)
savefig(fig, joinpath(DIR, "crm_new_fig_hf_locality.png"))
println("保存: crm_new_fig_hf_locality.png")
