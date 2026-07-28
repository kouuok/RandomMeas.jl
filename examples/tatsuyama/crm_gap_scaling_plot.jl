# 実験B(Δ の L 依存と gapped/gapless)の作図
#   julia --project=Hubbard_MPS_Env_v2 crm_gap_scaling_plot.jl
using DelimitedFiles, Statistics, Printf, Plots

raw, hdr = readdlm(joinpath(@__DIR__, "crm_gap_scaling_results.tsv"), '\t'; header=true)
col = Dict(String(h) => j for (j, h) in enumerate(vec(hdr))); gc(x) = raw[:, col[x]]
U = Float64.(gc("U")); δ = Float64.(gc("doping")); L = Int.(gc("L"))
cp = Int.(gc("chi_prior")); ob = String.(gc("observable"))
re = Float64.(gc("median_relerr"))
sr, _ = readdlm(joinpath(@__DIR__, "crm_gap_states.tsv"), '\t'; header=true)

Ls = [8, 16, 32, 64, 128]
cfg = [(0.0,0.0,"U=0 half"), (1.0,0.0,"U=1 half"), (2.0,0.0,"U=2 half"),
       (4.0,0.0,"U=4 half"), (8.0,0.0,"U=8 half"),
       (4.0,0.125,"U=4 δ=1/8"), (8.0,0.125,"U=8 δ=1/8")]

eps_of(u,d,o,x) = [ (i = findfirst(k -> U[k]==u && δ[k]==d && L[k]==l && cp[k]==x && ob[k]==o,
                                   1:length(L)); i === nothing ? NaN : re[i]) for l in Ls ]
S_of(u,d) = [ (i = findfirst(k -> sr[k,1]==u && sr[k,2]==d && sr[k,3]==l, 1:size(sr,1));
               i === nothing ? NaN : sr[i,7]) for l in Ls ]
function fit_c(u,d)
    y = S_of(u,d)[3:5]; x = log.(Ls[3:5])
    6*((mean(x.*y) - mean(x)*mean(y)) / (mean(x.^2) - mean(x)^2))
end

# --- 図1: 中央ボンドエントロピー — gapless セクターの同定 -----------------
p1 = plot(xlabel="L", ylabel="S_center", xscale=:log10,
          xticks=(Ls, string.(Ls)), legend=:topleft,
          title="Entanglement identifies the gapless sectors")
for (u,d,lab) in cfg
    plot!(p1, Ls, S_of(u,d), marker=:circle, ms=4, lw=2,
          ls = d==0 ? :solid : :dash,
          label=@sprintf("%s  (c=%.2f)", lab, fit_c(u,d)))
end

# --- 図2: 相対prior誤差の L 依存 -----------------------------------------
p2 = plot(xlabel="L", ylabel="median relative error  ε = |Δ/⟨P⟩|",
          xscale=:log10, yscale=:log10, xticks=(Ls, string.(Ls)),
          legend=:bottomright, title="Prior error vs L  (ZZ onsite, χ_p=8)")
for (u,d,lab) in cfg
    y = eps_of(u,d,"ZZ onsite",8)
    all(isnan, y) && continue
    plot!(p2, Ls, y, marker=:circle, ms=4, lw=2, ls = d==0 ? :solid : :dash, label=lab)
end

# --- 図3: エントロピーと prior 誤差の対応 --------------------------------
p3 = plot(xlabel="S_center(L)", ylabel="ε(L)", yscale=:log10, legend=:bottomright,
          title="Does entanglement predict the prior error?")
for (u,d,lab) in cfg
    y = eps_of(u,d,"ZZ onsite",8); s = S_of(u,d)
    all(isnan, y) && continue
    plot!(p3, s, y, marker=:circle, ms=4, lw=1.6, ls = d==0 ? :solid : :dash, label=lab)
end

# --- 図4: χ_p 依存 (L=128) — どこまで prior を良くすれば足りるか ---------
p4 = plot(xlabel="χ_prior", ylabel="ε (L=128)", xscale=:log10, yscale=:log10,
          xticks=([4,8,16,32], ["4","8","16","32"]), legend=:bottomleft,
          title="How good must the prior be?  (L=128)")
for (u,d,lab) in cfg
    y = [ (i = findfirst(k -> U[k]==u && δ[k]==d && L[k]==128 && cp[k]==x && ob[k]=="ZZ onsite",
                         1:length(L)); i === nothing ? NaN : re[i]) for x in [4,8,16,32] ]
    all(isnan, y) && continue
    plot!(p4, [4,8,16,32], y, marker=:circle, ms=4, lw=2, ls = d==0 ? :solid : :dash, label=lab)
end

fig = plot(p1, p2, p3, p4, layout=(2,2), size=(1400, 900), margin=6Plots.mm)
out = joinpath(@__DIR__, "crm_gap_scaling.png")
savefig(fig, out)
println("figure saved: $out")

# --- 要約表 ---------------------------------------------------------------
println("\n", "="^92)
@printf("%-14s %10s %12s %12s %12s\n", "config", "c (L≥32)", "ε(L=8)", "ε(L=128)", "L128/L32")
for (u,d,lab) in cfg
    y = eps_of(u,d,"ZZ onsite",8)
    @printf("%-14s %10.2f %12.3e %12.3e %12.2f\n", lab, fit_c(u,d), y[1], y[5], y[5]/y[3])
end
