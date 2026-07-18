# W=6 シリンダー結果 (crm_2d_w6_results_U*.tsv) の清書図
# 実行: JULIA_LOAD_PATH="@:@v#.#:@stdlib" julia --project=Hubbard_MPS_Env_v2 crm_2d_w6_plot.jl
using DelimitedFiles, Printf, Plots
const OI = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#56B4E9", "#E69F00"]
const MK = [:circle, :square, :diamond, :utriangle, :dtriangle, :pentagon]
default(fontfamily="Helvetica", framestyle=:box, grid=true, gridalpha=0.12,
        linewidth=1.8, markersize=6, markerstrokewidth=0.6, markerstrokecolor=:white,
        guidefontsize=11, tickfontsize=9, legendfontsize=8, titlefontsize=11, dpi=200)

dir = @__DIR__
rows = []
for f in filter(f -> startswith(f, "crm_2d_w6_results_U"), readdir(dir))
    d, hdr = readdlm(joinpath(dir, f), '\t'; header=true)
    for r in 1:size(d, 1)
        push!(rows, Tuple(d[r, :]))
    end
end
isempty(rows) && error("no crm_2d_w6_results_U*.tsv found")
# 列: U, prior, observable, pure, true, Delta, prior_fid, G_emp, G_theo
prior_order = ["UHF", "UHF-sym", "chi=8", "chi=16", "chi=32", "chi=64", "chi=128"]
obs_sel = [("ZZ onsite", "on-site ZZ"), ("DoubleOcc", "double occupancy"),
           ("hop y-bond", "kinetic bond"), ("ZZ up-up x1", "ZZ  x-neighbor"),
           ("SzSz x-bond", "SzSz  x-neighbor"), ("ZZ up-up y2", "ZZ  y-distance 2")]

Us = sort(unique(r[1] for r in rows))
panels = []
for (iu, U) in enumerate(Us)
    prs = [p for p in prior_order if any(r[1]==U && string(r[2])==p for r in rows)]
    p = plot(xticks=(1:length(prs), prs), xrotation=30, yscale=:log10,
             ylabel=(iu == 1 ? "G" : ""), legend=(iu == length(Us) ? :outerright : :none),
             title=@sprintf("W=6 x 8 cylinder (96 qubits),  U=%.0f", U))
    for (i, (name, lab)) in enumerate(obs_sel)
        g = [only([r[8] for r in rows if r[1]==U && string(r[2])==pr && string(r[3])==name])
             for pr in prs]
        plot!(p, 1:length(prs), g, color=OI[i], marker=MK[i], label=lab)
    end
    hline!(p, [1.0], color=:black, ls=:dot, lw=1, label="")
    vline!(p, [2.5], color=:gray70, ls=:dash, lw=1, label="")
    push!(panels, p)
end
final = plot(panels..., layout=(1, length(Us)), size=(660*length(Us)+240, 460), margin=7Plots.mm)
savefig(final, joinpath(dir, "crm_2d_w6.png"))
println("figure saved: crm_2d_w6.png")

# W=4 との比較表 (共通observable/priorのG比較を標準出力に)
d4, hdr4 = readdlm(joinpath(dir, "crm_2d_symuhf_results.tsv"), '\t'; header=true)
println("\n=== G comparison W=4 vs W=6 (same observable/prior) ===")
@printf("%-6s %-12s %-9s %10s %10s\n", "U", "observable", "prior", "G(W=4)", "G(W=6)")
for U in Us, (name, _) in obs_sel, pr in ["UHF", "UHF-sym", "chi=32"]
    g6 = [r[8] for r in rows if r[1]==U && string(r[2])==pr && string(r[3])==name]
    g4 = [d4[r, 8] for r in 1:size(d4,1)
          if d4[r,1]==U && string(d4[r,2])==pr && string(d4[r,3])==name]
    (isempty(g6) || isempty(g4)) && continue
    @printf("%-6.0f %-12s %-9s %10.2f %10.2f\n", U, name, pr, g4[1], g6[1])
end
