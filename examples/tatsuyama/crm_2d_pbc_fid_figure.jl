# 2次元(4×4 シリンダー / トーラス、32量子ビット)で、
# 「大域忠実度が低くても局所観測量なら平均場 prior が働く」を3枚で示す。
#   julia --project=. examples/tatsuyama/crm_2d_pbc_fid_figure.jl
using DelimitedFiles, Printf, Plots
gr(fontfamily="Helvetica", legendfontsize=6, guidefontsize=9,
   titlefontsize=10, tickfontsize=8)
const DIR = @__DIR__

function load(geo)
    A, h = readdlm(joinpath(DIR, "crm_2d_pbc_fid_W4L4_$(geo)_U8.0.tsv"), '\t'; header=true)
    h = vec(h); c(n) = findfirst(==(n), h)
    (prior=String.(A[:, c("prior")]), obs=String.(A[:, c("observable")]),
     fid=Float64.(A[:, c("global_fid")]), G=Float64.(A[:, c("G")]))
end
cyl = load("cylinder"); tor = load("torus")
chis = ["chi2","chi4","chi8","chi16","chi32","chi64","chi128"]
pick(d, p, o) = begin i = findfirst(k -> d.prior[k]==p && d.obs[k]==o, eachindex(d.prior)); i===nothing ? NaN : d.G[i] end
pf(d, p) = begin i = findfirst(==(p), d.prior); d.fid[i] end

# ---------- (a) 大域忠実度 vs 利得(onsite ZZ) ----------
pa = plot(xscale=:log10, xlabel="global fidelity  F(ρ, σ)", ylabel="gain G  (onsite ZZ)",
          title="(a) 4x4, 32 qubits: global F vs local gain", legend=:bottomright,
          xlims=(0.02, 1.6), ylims=(0, 175))
for (d, nm, col, mk) in ((cyl,"cylinder (MPS)",:steelblue,:circle), (tor,"torus (MPS)",:seagreen,:rect))
    plot!(pa, [pf(d,c) for c in chis], [pick(d,c,"ZZ onsite") for c in chis],
          marker=mk, markersize=5, color=col, label=nm, linewidth=1.5)
end
scatter!(pa, [pf(cyl,"UHF")], [pick(cyl,"UHF","ZZ onsite")], marker=:star6, markersize=13,
         color=:black, label="UHF, cylinder", markerstrokewidth=1)
scatter!(pa, [pf(tor,"UHF")], [pick(tor,"UHF","ZZ onsite")], marker=:star5, markersize=13,
         color=:firebrick, label="UHF, torus", markerstrokewidth=1)
annotate!(pa, 0.14, 163, text("mean field: F ≈ 0.11\nbut the best gain\nof every prior tried", 7, :left, :black))

# ---------- (b) UHF を観測量ごとに ----------
obs = ["ZZ onsite","ZZ up-up nb","SzSz nb","DoubleOcc","Sz"]
short = ["ZZ onsite","ZZ up-up","SzSz","double occ.","Sz"]
mat = hcat([pick(cyl,"UHF",o) for o in obs], [pick(tor,"UHF",o) for o in obs])
pb = plot(ylabel="log10 gain G", title="(b) the same mean-field prior, five observables",
          xticks=(1:length(obs), short), xrotation=18, legend=:topright, ylims=(-1.6, 2.6))
w = 0.8/2
for (j, nm, col) in ((1,"cylinder",:steelblue), (2,"torus",:firebrick))
    bar!(pb, (1:length(obs)) .+ (j-1.5)*w, log10.(max.(mat[:,j],1e-2)),
         bar_width=w*0.9, label=nm, color=col, linewidth=0)
end
hline!(pb, [0.0], color=:black, linestyle=:dash, label="G = 1")
annotate!(pb, 2.6, 2.3, text("CRM wins", 8, :left, :darkgreen))
annotate!(pb, 2.6, -1.35, text("CRM loses", 8, :left, :darkred))

# ---------- (c) MPS は PBC の代金を払うが、平均場は払わない ----------
xs = [2,4,8,16,32,64,128]
pc = plot(xscale=:log2, xlabel="MPS bond dimension χ_p", ylabel="global fidelity  F",
          title="(c) periodic boundaries cost the MPS prior, not the mean field",
          legend=:bottomright, ylims=(0, 1.05), xticks=(xs, string.(xs)))
plot!(pc, xs, [pf(cyl,c) for c in chis], marker=:circle, color=:steelblue,
      label="cylinder (MPS)", linewidth=2)
plot!(pc, xs, [pf(tor,c) for c in chis], marker=:rect, color=:seagreen,
      label="torus (MPS)", linewidth=2)
hline!(pc, [pf(cyl,"UHF")], color=:black, linestyle=:dash, label="UHF, cylinder")
hline!(pc, [pf(tor,"UHF")], color=:firebrick, linestyle=:dot, label="UHF, torus")
annotate!(pc, 3.0, 0.20, text(@sprintf("UHF: %.3f → %.3f\n(χ-free, boundary-insensitive)",
                                       pf(cyl,"UHF"), pf(tor,"UHF")), 7, :left, :black))

fig = plot(pa, pb, pc, layout=(1,3), size=(1650, 470), margin=6Plots.mm, bottom_margin=10Plots.mm)
savefig(fig, joinpath(DIR, "crm_new_fig_hf_locality_2d.png"))
println("保存: crm_new_fig_hf_locality_2d.png")
