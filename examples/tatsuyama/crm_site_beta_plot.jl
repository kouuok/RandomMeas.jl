# 実験A(最適係数CRM)の解析・作図
#   julia --project=Hubbard_MPS_Env_v2 crm_site_beta_plot.jl
using DelimitedFiles, Statistics, Printf, Plots

raw, hdr = readdlm(joinpath(@__DIR__, "crm_site_beta_results.tsv"), '\t'; header=true)
col = Dict(String(h) => j for (j, h) in enumerate(vec(hdr)))
gc(x) = raw[:, col[x]]

L = Int.(gc("L")); δ = Float64.(gc("doping")); site = Int.(gc("site"))
cp = Int.(gc("chi_prior")); nu = Int.(gc("nu")); ob = String.(gc("observable"))
Gmax = Float64.(gc("G_max"))

const EST = ["crm", "beta", "shrunk", "split", "exactvar", "splitev"]
const LAB = ["CRM (β=1)", "plug-in β", "shrunk β", "split-sample", "exact Var(Y)", "split × exactVar"]
G  = [Float64.(gc("G_" * e)) for e in EST]
HI = [Float64.(gc(e == "crm" ? "G_crm_hi" : "G_" * e * "_hi")) for e in EST]
M  = [Float64.(gc("M_" * e)) for e in EST]
B  = [abs.(Float64.(gc("bias_" * e))) for e in EST]

sel(f) = findall(f, 1:length(L))

# ---- 図1: サイト分解 — 損がどこで起き、どう直るか ------------------------
p1 = plot(xlabel="site i", ylabel="G = Var_std / Var_est", yscale=:log10,
          legend=:bottomright, title="Loss regime and its repair\n(L=64, δ=1/8, χ_p=2, n_u=200)")
hspan!(p1, [1e-3, 1.0], color=:red, alpha=0.08, label="CRM loses")
for (e, cl) in ((1, :firebrick), (2, :steelblue), (5, :seagreen))
    idx = sel(i -> L[i]==64 && δ[i]==0.125 && cp[i]==2 && nu[i]==200 && ob[i]=="ZZ onsite")
    isempty(idx) && continue
    o = sortperm(site[idx]); idx = idx[o]
    plot!(p1, site[idx], G[e][idx], lw=1.5, color=cl, marker=:circle, ms=2, label=LAB[e])
end
hline!(p1, [1.0], color=:black, ls=:dot, label="")

# ---- 図2: MSE比 vs prior品質（分散比ではなくMSEで比較する）---------------
p2 = plot(xlabel="χ_prior", ylabel="median MSE ratio", xscale=:log10, yscale=:log10,
          xticks=([2,4,8,16,32], ["2","4","8","16","32"]),
          yticks=([1,2,5,10,20], ["1","2","5","10","20"]),
          legend=:topleft, title="Ranking by MSE, not variance\n(L=64, δ=1/8, n_u=200)")
for e in 1:6
    ys = [ (idx = sel(i -> L[i]==64 && δ[i]==0.125 && cp[i]==x && nu[i]==200);
            isempty(idx) ? NaN : median(M[e][idx])) for x in [2,4,8,16,32] ]
    plot!(p2, [2,4,8,16,32], ys, marker=:circle, ms=4, lw=2, label=LAB[e])
end
hline!(p2, [1.0], color=:black, ls=:dot, label="")

# ---- 図3: バイアスと「損なサイト数」のトレードオフ -----------------------
p3 = plot(xlabel="max |bias| / sd", ylabel="losing sites at n_u=50  (of 445)",
          legend=:topleft, title="Bias vs worst-case variance\n(L=64, δ=1/8, χ_p=2)")
for e in 1:6
    idx50 = sel(i -> L[i]==64 && δ[i]==0.125 && cp[i]==2 && nu[i]==50)
    isempty(idx50) && continue
    scatter!(p3, [maximum(B[e][idx50])], [count(i -> HI[e][i] < 0.9, idx50)],
             ms=8, label=LAB[e], markerstrokewidth=0)
end

# ---- 図4: priorフリー性 — 単一Pauli列 vs 多項観測量 ----------------------
p4 = plot(xlabel="χ_prior", ylabel="G (plug-in β)", xscale=:log10,
          xticks=([2,4,8,16,32], ["2","4","8","16","32"]),
          legend=:bottomright, title="Prior-free for single Pauli strings\n(L=64, δ=1/8, n_u=200)")
for (o, st) in (("ZZ onsite", :solid), ("ZZ up-up r=1", :solid),
                ("DoubleOcc", :dash), ("SzSz r=1", :dash))
    ys = [ (idx = sel(i -> L[i]==64 && δ[i]==0.125 && cp[i]==x && nu[i]==200 && ob[i]==o);
            isempty(idx) ? NaN : median(G[2][idx])) for x in [2,4,8,16,32] ]
    plot!(p4, [2,4,8,16,32], ys, marker=:circle, ms=4, lw=2, ls=st,
          label=o * (st == :solid ? " (1 Pauli string)" : " (multi-term)"))
end

fig = plot(p1, p2, p3, p4, layout=(2,2), size=(1400, 900), margin=6Plots.mm)
out = joinpath(@__DIR__, "crm_site_beta.png")
savefig(fig, out)
println("figure saved: $out")

# ---- 表示: 実務推奨の根拠となる要約 --------------------------------------
println("\n", "="^96)
println("推定量の選択指針 (L=64, δ=1/8, χ_p=2)")
@printf("%-18s %12s %12s %12s %12s\n", "estimator", "損@n_u=50", "損@n_u=200",
        "MSE比(χ_p=4)", "max|bias|/sd")
for e in 1:6
    i50 = sel(i -> L[i]==64 && δ[i]==0.125 && cp[i]==2 && nu[i]==50)
    i200 = sel(i -> L[i]==64 && δ[i]==0.125 && cp[i]==2 && nu[i]==200)
    i4 = sel(i -> L[i]==64 && δ[i]==0.125 && cp[i]==4 && nu[i]==200)
    (isempty(i50) || isempty(i200)) && continue
    @printf("%-18s %12s %12s %12.2f %12.2f\n", LAB[e],
            "$(count(i -> HI[e][i] < 0.9, i50))/$(length(i50))",
            "$(count(i -> HI[e][i] < 0.9, i200))/$(length(i200))",
            isempty(i4) ? NaN : median(M[e][i4]), maximum(B[e][i200]))
end
