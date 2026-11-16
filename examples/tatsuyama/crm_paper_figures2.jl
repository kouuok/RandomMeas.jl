# ============================================================
# 論文向けの追加図(既存データのみを使う。新しい計算はしない)
#   julia --project=Hubbard_MPS_Env_v2 crm_paper_figures2.jl
#
#   collapse    : 利得法則の普遍曲線  (G-1)/(G_max-1) = x/(1+x),  x = n_m/n_m*
#   designmap   : 設計図  G(ε, n_m) と、切替コストがあるときの最適 n_m
#   bc_prior    : 開放端と周期境界で prior の序列がどう変わるか(L=14, ED参照)
#   chistar     : 必要ボンド次元 χ*(L) と中心電荷
#   derand      : 貪欲 derandomization との比較
#   nu          : n_u を変えても利得は偏らない(8シード)
#   symbreak    : 4x4 トーラスで DMRG が SU(2) を自発的に破る
#
# ラベルは英語(crm_new_figures.jl と揃える)。
# 対数量の棒は使わず、線・点で描く(対数軸の棒は基準が -∞ になるため)。
# ============================================================
using DelimitedFiles, Statistics, Printf, Plots
gr(fontfamily="Helvetica", legendfontsize=7, guidefontsize=9,
   titlefontsize=10, tickfontsize=8)

const DIR = @__DIR__
function load(f)
    raw, hdr = readdlm(joinpath(DIR, f), '\t'; header=true)
    c = Dict(String(h) => j for (j, h) in enumerate(vec(hdr)))
    (raw = raw, n = size(raw, 1), col = name -> raw[:, c[name]])
end
num(v) = Float64.(v); str(v) = String.(strip.(string.(v)))
out(p, name) = (f = joinpath(DIR, name); savefig(p, f); println("wrote ", name))

# 単一 Pauli 列に対する利得法則の部品
gmax(P, Δ) = P^2 / Δ^2
nmstar(P, Δ, A) = 3.0^A * (1 - P^2) / ((3.0^A - 1) * Δ^2)
gain(P, Δ, A, nm) = ((3.0^A - 1) * P^2 + 3.0^A * (1 - P^2) / nm) /
                    ((3.0^A - 1) * Δ^2 + 3.0^A * (1 - P^2) / nm)

# ------------------------------------------------------------
# Fig 1: 普遍曲線への崩壊
#   G - 1 = (G_max - 1) x/(1+x),  x = n_m/n_m*  (単一 Pauli 列で厳密)
# ------------------------------------------------------------
d = load("crm_sweep_nm.tsv")
nm = num(d.col("nm")); chi = Int.(d.col("chi_prior")); ob = str(d.col("observable"))
P = num(d.col("true")); Δ = num(d.col("Delta")); Ge = num(d.col("G_emp"))
single = ["ZZ onsite(i0)", "ZZ up-up r=1", "ZZ up-up r=2", "ZZ up-up r=4"]   # |A| = 2
chis = sort(unique(chi))
ccol = Dict(c => cgrad(:viridis)[(k - 1) / (length(chis) - 1) * 0.9] for (k, c) in enumerate(chis))
mk = Dict(single[1] => :circle, single[2] => :square, single[3] => :diamond, single[4] => :utriangle)

p1a = plot(xscale=:log10, yscale=:log10, xlabel="shots per setting  n_m", ylabel="gain  G",
           title="(a) raw: 20 curves (4 Pauli strings x 5 priors)", legend=:topleft)
for o in single, c in chis
    k = findall(i -> ob[i] == o && chi[i] == c, 1:d.n)
    k = k[sortperm(nm[k])]
    plot!(p1a, nm[k], Ge[k], color=ccol[c], marker=mk[o], ms=3.5, lw=1, msw=0,
          label=(o == single[1] ? "chi_p = $c" : ""))
end
hline!(p1a, [1.0], color=:black, ls=:dot, label="")

xs = Float64[]; ys = Float64[]; cs = Int[]; os = String[]
for i in 1:d.n
    ob[i] in single || continue
    gm = gmax(P[i], Δ[i])
    abs(gm - 1) > 0.05 || continue          # G_max ≈ 1 では比が不定になるので除く
    push!(xs, nm[i] / nmstar(P[i], Δ[i], 2)); push!(ys, (Ge[i] - 1) / (gm - 1))
    push!(cs, chi[i]); push!(os, ob[i])
end
xl = 10 .^ range(log10(minimum(xs)) - 0.3, log10(maximum(xs)) + 0.3, length=200)
p1b = plot(xscale=:log10, yscale=:log10, xlabel="n_m / n_m*",
           ylabel="(G - 1) / (G_max - 1)", legend=:bottomright,
           title="(b) collapsed onto x/(1+x)")
plot!(p1b, xl, xl ./ (1 .+ xl), color=:black, lw=1.5, label="x / (1 + x)   (no fit parameter)")
for o in single
    k = findall(i -> os[i] == o && ys[i] > 0, 1:length(xs))
    scatter!(p1b, xs[k], ys[k], marker=mk[o], ms=4, msw=0,
             color=[ccol[c] for c in cs[k]], label="")
    scatter!(p1b, [NaN], [NaN], marker=mk[o], ms=4, msw=0, color=:gray40,
             label=replace(o, "(i0)" => "") * "  (L = 8)")
end
# 独立な実験(固定予算の配分, L=32)も同じ曲線に乗る
s = load("crm_shot_split_results.tsv")
ks = findall(i -> str(s.col("observable"))[i] == "ZZ onsite" && Int(s.col("B")[i]) == 5000 &&
                  str(s.col("prior"))[i] in ("chi2", "chi4", "chi8"), 1:s.n)
# 注意: このファイルの G_max 列は「その n_m での厳密 prior の利得」なので使わず、(P, Δ) から作り直す
sP = num(s.col("true"))[ks]; sΔ = num(s.col("Delta"))[ks]
sx = num(s.col("n_m"))[ks] ./ nmstar.(sP, sΔ, 2)
sy = (num(s.col("G_emp"))[ks] .- 1) ./ (gmax.(sP, sΔ) .- 1)
kk = findall(>(0), sy)
scatter!(p1b, sx[kk], sy[kk], marker=:xcross, ms=4, color=:firebrick,
         label="ZZ onsite, fixed-budget run  (L = 32)")
nneg = count(<=(0), ys) + count(<=(0), sy)
println("collapse: points = ", length(xs) + length(sx), "  (non-positive, not shown on log axis: ", nneg, ")")
resid = vcat([log10(ys[i] / (xs[i] / (1 + xs[i]))) for i in 1:length(xs) if ys[i] > 0],
             [log10(sy[i] / (sx[i] / (1 + sx[i]))) for i in kk])
@printf("collapse: x spans %.1f decades\n", log10(maximum(vcat(xs, sx)) / minimum(vcat(xs, sx))))
@printf("collapse: median |log10(data/curve)| = %.3f  (factor %.2f)\n",
        median(abs.(resid)), 10^median(abs.(resid)))
out(plot(p1a, p1b, layout=(1, 2), size=(1050, 420), margin=5Plots.mm), "crm_new_fig_collapse.png")

# ------------------------------------------------------------
# Fig 2: 設計図
# ------------------------------------------------------------
P0 = 0.6; A0 = 2                      # onsite ZZ, U = 4 の典型値 |<P>| ≈ 0.6
εs = 10 .^ range(-2.5, 0.5, length=240); nms = 10 .^ range(0, 4, length=200)
Z = [log10(gain(P0, ε * P0, A0, n)) for n in nms, ε in εs]
p2a = heatmap(εs, nms, Z, xscale=:log10, yscale=:log10, c=:RdBu, clims=(-3.2, 3.2),
              xlabel="relative prior error  |Delta / <P>|", ylabel="shots per setting  n_m",
              colorbar_title="log10 G", title="(a) gain map  (|A| = 2, |<P>| = 0.6)")
contour!(p2a, εs, nms, Z, levels=[0.5, 1.0, 2.0, 3.0], color=:black, lw=0.6, colorbar_entry=false)
plot!(p2a, εs, [nmstar(P0, ε * P0, A0) for ε in εs], color=:black, lw=2, ls=:dash, label="n_m*")
vline!(p2a, [1.0], color=:black, lw=1.5, label="G = 1")
ylims!(p2a, (1, 1e4)); xlims!(p2a, (εs[1], εs[end]))
# 実際の prior がどこにいるか(L=8, U=4, onsite ZZ)
for c in chis
    i = findfirst(j -> ob[j] == "ZZ onsite(i0)" && chi[j] == c, 1:d.n)
    ε = abs(Δ[i] / P[i])
    ε < εs[1] && continue
    scatter!(p2a, [ε], [100.0], marker=:star5, ms=7, color=:gold, msc=:black, label="")
    annotate!(p2a, ε, 190, text("chi_p=$c", 7, :black))
end
annotate!(p2a, 0.008, 2.5, text("n_m < n_m*:\nmore shots per setting help", 7, :black))
annotate!(p2a, 0.3, 2500, text("n_m > n_m*:\nsaturated,\nG ~ G_max", 7, :black))
annotate!(p2a, 2.2, 30, text("loss", 8, :black))

# 切替コスト c(設定を1回変えるコスト、ショット換算)があるときの固定予算分散
#   B Var = (n_m + c) [ (3^A-1) Δ² + 3^A(1-P²)/n_m ]   → 最小は n_m = sqrt(c n_m*)
i8 = findfirst(j -> ob[j] == "ZZ onsite(i0)" && chi[j] == 4, 1:d.n)
Pv, Δv = P[i8], Δ[i8]; ns = nmstar(Pv, Δv, 2)
varB(n, c, Δx) = (n + c) * ((3.0^2 - 1) * Δx^2 + 3.0^2 * (1 - Pv^2) / n)
p2b = plot(xscale=:log10, yscale=:log10, xlabel="shots per setting  n_m",
           ylabel="B x Var   (fixed total cost B)", legend=:topleft,
           title=@sprintf("(b) optimal n_m with switching cost c  (chi_p = 4, n_m* = %.0f)", ns))
for (c, col) in zip((0.0, 10.0, 100.0, 1000.0), (:gray60, :seagreen, :darkorange, :firebrick))
    plot!(p2b, nms, varB.(nms, c, Δv), color=col, lw=2, label="CRM,  c = $(Int(c))")
    plot!(p2b, nms, varB.(nms, c, Pv), color=col, lw=1, ls=:dash, label="")
    if c > 0
        no = sqrt(c * ns)
        scatter!(p2b, [no], [varB(no, c, Δv)], marker=:circle, ms=6, color=col, msc=:black, label="")
    end
end
plot!(p2b, [NaN], [NaN], color=:black, lw=1, ls=:dash, label="standard shadow (same c)")
scatter!(p2b, [NaN], [NaN], marker=:circle, ms=6, color=:white, msc=:black, label="n_m = sqrt(c n_m*)")
out(plot(p2a, p2b, layout=(1, 2), size=(1150, 440), margin=5Plots.mm), "crm_new_fig_designmap.png")

# ------------------------------------------------------------
# Fig 3: 境界条件と prior の序列(L=14, ED 参照, 3シード × 全サイト)
# ------------------------------------------------------------
e = load("crm_ed_pbc_results.tsv")
eL = Int.(e.col("L")); ebc = str(e.col("bc")); erho = str(e.col("rho"))
epr = str(e.col("prior")); eob = str(e.col("observable")); eG = num(e.col("G"))
ladder = ["chi2", "chi4", "chi8", "chi16", "chi32", "exact"]
ticks = ["2", "4", "8", "16", "32", "exact"]
function gstat(o, bc, p)
    v = eG[findall(i -> eL[i] == 14 && erho[i] == "ED" && ebc[i] == bc && epr[i] == p && eob[i] == o, 1:e.n)]
    (median(v), quantile(v, 0.1), quantile(v, 0.9))
end
panels = []
for (j, o) in enumerate(["ZZ onsite", "SzSz r=1", "hop up r=1"])
    p = plot(yscale=:log10, xticks=(1:8, vcat(ticks, ["UHF", "UHF-sym"])), xrotation=30,
             xlabel="prior:  MPS bond dimension chi_p  |  mean field", ylabel=(j == 1 ? "gain G" : ""),
             title="($("abc"[j])) $o", legend=(j == 1 ? :right : :none), ylims=(0.3, 120),
             xlims=(0.5, 8.5), bottom_margin=12Plots.mm)
    hline!(p, [1.0], color=:black, ls=:dot, label="")
    vline!(p, [6.5], color=:gray70, label="")
    for (bc, col, dx) in (("OBC", :steelblue, -0.08), ("PBC", :firebrick, 0.08))
        st = [gstat(o, bc, q) for q in ladder]
        m = first.(st); lo = getindex.(st, 2); hi = getindex.(st, 3)
        plot!(p, (1:6) .+ dx, m, yerror=(m .- lo, hi .- m), color=col, marker=:circle, ms=4, lw=2,
              msc=col, label="$bc  (median, 10-90% over sites and 3 seeds)")
        st2 = [gstat(o, bc, q) for q in ("UHF", "UHFsym")]
        m2 = first.(st2)
        scatter!(p, [7, 8] .+ dx, m2, yerror=(m2 .- getindex.(st2, 2), getindex.(st2, 3) .- m2),
                 color=col, marker=:diamond, ms=5, msc=col, label="")
    end
    push!(panels, p)
end
annotate!(panels[1], 3.5, 0.42, text("L = 14 chain, U = 4, exact (ED) reference", 7, :black))
out(plot(panels..., layout=(1, 3), size=(1300, 450), margin=6Plots.mm), "crm_new_fig_bc_prior.png")

# ------------------------------------------------------------
# Fig 4: 必要ボンド次元 χ*(L) と中心電荷
# ------------------------------------------------------------
c4 = load("crm_chistar_fixed.tsv")
cU = num(c4.col("U")); cd = num(c4.col("doping")); cL = Int.(c4.col("L"))
ct = num(c4.col("target")); cx = num(c4.col("chi_star"))
g4 = load("crm_gap_states.tsv")
Sc = Dict((Float64(g4.col("U")[i]), Float64(g4.col("doping")[i]), Int(g4.col("L")[i])) =>
          Float64(g4.col("S_center")[i]) for i in 1:g4.n)
cfgs = [(4.0, 0.0, "U=4 half filling", :steelblue, :solid),
        (8.0, 0.0, "U=8 half filling", :navy, :solid),
        (2.0, 0.0, "U=2 half filling", :seagreen, :solid),
        (1.0, 0.0, "U=1 half filling", :darkorange, :solid),
        (4.0, 0.125, "U=4 doped 1/8", :firebrick, :dash),
        (8.0, 0.125, "U=8 doped 1/8", :purple, :dash)]
Ls = [16, 32, 64, 128]
geo(v) = exp(mean(log.(v)))
p4a = plot(xscale=:log10, yscale=:log10, xticks=(Ls, string.(Ls)), yticks=([8, 12, 16, 24, 32, 48], string.([8, 12, 16, 24, 32, 48])),
           xlabel="chain length  L", ylabel="required bond dimension  chi*  (rel. error < 1%)",
           legend=:topleft, title="(a) prior quality needed vs system size", ylims=(7.5, 48))
ps = Float64[]; ccs = Float64[]
for (U, dp, lab, col, ls) in cfgs
    y = [geo(cx[findall(i -> cU[i] == U && cd[i] == dp && cL[i] == L && ct[i] == 0.01, 1:c4.n)]) for L in Ls]
    plot!(p4a, Ls, y, color=col, ls=ls, lw=2, marker=:circle, ms=4, msw=0, label=lab)
    push!(ps, log(y[end] / y[1]) / log(Ls[end] / Ls[1]))
    push!(ccs, 6 * (Sc[(U, dp, 128)] - Sc[(U, dp, 32)]) / log(4))
end
p4b = plot(xlabel="central charge c  (from entanglement entropy, L = 32 -> 128)",
           ylabel="growth exponent p   (chi* ~ L^p)", legend=:topleft,
           title=@sprintf("(b) exponent tracks c   (corr = %.2f)", cor(ccs, ps)),
           xlims=(0.8, 1.8), ylims=(0, 0.55),
           yticks=(0:0.1:0.5, ["0", "0.1", "0.2", "0.3", "0.4", "0.5"]))
cl = range(0.8, 1.8, length=10)
plot!(p4b, cl, cl ./ 6, color=:black, ls=:dash, label="p = c/6  (entropy alone)")
for (k, (U, dp, lab, col, ls)) in enumerate(cfgs)
    scatter!(p4b, [ccs[k]], [ps[k]], color=col, ms=7, msw=0,
             marker=(dp == 0 ? :circle : :diamond), label=lab)
end
for k in 1:length(cfgs)
    @printf("chistar: %-18s p = %.3f  c = %.3f\n", cfgs[k][3], ps[k], ccs[k])
end
out(plot(p4a, p4b, layout=(1, 2), size=(1100, 430), margin=5Plots.mm), "crm_new_fig_chistar.png")

# ------------------------------------------------------------
# Fig 5: 貪欲 derandomization との比較
# ------------------------------------------------------------
r = load("crm_derand_results.tsv")
rs = str(r.col("set")); es = num(r.col("err_shadow")); ec = num(r.col("err_crm")); ed = num(r.col("err_derand"))
cols3 = (:steelblue, :firebrick, :seagreen); labs3 = ("random shadow", "random shadow + CRM", "greedy derandomization")
function cdfplot(set, ttl; leg=:bottomright)
    k = findall(==(set), rs)
    p = plot(xscale=:log10, xlabel="error of the estimate", ylabel="fraction of observables below",
             legend=leg, title=ttl, ylims=(0, 1.02))
    for (v, col, lab) in zip((es, ec, ed), cols3, labs3)
        x = sort(v[k]); n = length(x)
        plot!(p, x, (1:n) ./ n, color=col, lw=2, lt=:steppost, label=@sprintf("%s  (median %.3f)", lab, median(x)))
    end
    @printf("derand %-18s n=%d  median shadow %.4f  CRM %.4f  derand %.4f\n", set, length(k),
            median(es[k]), median(ec[k]), median(ed[k]))
    p
end
p5a = cdfplot("SetA(structured)", "(a) structured set (23 physical observables)"; leg=:bottomright)
p5b = cdfplot("SetB(random)", "(b) random set (150 Pauli strings)"; leg=:bottomright)
out(plot(p5a, p5b, layout=(1, 2), size=(1100, 420), margin=5Plots.mm), "crm_new_fig_derand.png")

# ------------------------------------------------------------
# Fig 6: n_u を変えても利得は偏らない
# ------------------------------------------------------------
u = load("crm_nu_seed_results.tsv")
unu = Int.(u.col("n_u")); uob = str(u.col("observable")); upr = str(u.col("prior"))
uG = num(u.col("G")); uA = Int.(u.col("absA")); nrep = Int(u.col("n_repeat")[1])
obsu = [o for o in unique(uob) if !startswith(o, "Sz ")]       # Sz は G = 1 なので除く
nus = sort(unique(unu))
ocol = palette(:tab10)
p6a = plot(xscale=:log10, xticks=(nus, string.(nus)), xlabel="number of random settings  n_u",
           ylabel="G / (mean over all n_u and seeds)", legend=:outerright, ylims=(0.4, 1.6),
           title="(a) no drift with n_u  (chi_p = 8, 8 seeds each)", xlims=(30, 5500))
hline!(p6a, [1.0], color=:black, ls=:dot, label="")
p6b = plot(xscale=:log10, xticks=(nus, string.(nus)), xlabel="number of random settings  n_u",
           ylabel="relative seed-to-seed spread of G", legend=:topright, ylims=(0, 0.3),
           title="(b) spread vs prediction  (8 seeds: each point +-27%)", xlims=(30, 5500))
for (j, o) in enumerate(obsu)
    ka = findall(i -> uob[i] == o && upr[i] == "chi8", 1:u.n)
    gbar_ = mean(uG[ka]); A = uA[ka[1]]
    ms_ = Float64[]; sd_ = Float64[]
    for (q, n) in enumerate(nus)
        v = uG[ka[unu[ka] .== n]] ./ gbar_
        jit = n .* 10 .^ ((j - (length(obsu) + 1) / 2) * 0.035)
        scatter!(p6a, fill(jit, length(v)), v, color=ocol[j], ms=2.5, msw=0, alpha=0.6, label="")
        scatter!(p6a, [jit], [mean(v)], yerror=[std(v) / sqrt(length(v))], color=ocol[j], ms=5,
                 marker=:diamond, msc=:black, label=(q == 1 ? o : ""))
        push!(sd_, std(v))
    end
    plot!(p6b, nus, sd_, color=ocol[j], marker=:circle, ms=4, msw=0, lw=0.7, alpha=0.7, label="")
end
for (A, ls) in ((2, :dash), (3, :dashdot), (4, :solid))
    plot!(p6b, nus, [sqrt((2 + 3.0^A / n) / nrep) for n in nus], color=:black, ls=ls, lw=1.5,
          label="sqrt((2 + 3^|A|/n_u)/n_repeat),  |A| = $A")
end
out(plot(p6a, p6b, layout=grid(1, 2, widths=[0.58, 0.42]), size=(1200, 420), margin=5Plots.mm),
    "crm_new_fig_nu.png")

# ------------------------------------------------------------
# Fig 7: 4x4 トーラスで DMRG(χ=1024)が SU(2) を自発的に破る
# ------------------------------------------------------------
t = load("crm_2d_ed_W4L4_U8p0.tsv")
tg = str(t.col("geometry")); trho = str(t.col("rho")); tob = str(t.col("observable"))
tpr = str(t.col("prior")); ts = Int.(t.col("site")); tv = num(t.col("true"))
W = 4
val(geo_, rho, o) = (k = findall(i -> tg[i] == geo_ && trho[i] == rho && tob[i] == o && tpr[i] == "exact", 1:t.n);
                     (ts[k], tv[k]))
function szmap(geo_, rho, ttl)
    s, v = val(geo_, rho, "Sz")
    x = (s .- 1) .÷ W .+ 1; y = (s .- 1) .% W .+ 1
    p = plot(xlims=(0.4, 4.6), ylims=(0.4, 4.6), aspect_ratio=1, xticks=1:4, yticks=1:4,
             xlabel="x", ylabel="y", title=ttl, legend=false, grid=false)
    for a in 1:4, b in 1:4     # 格子の辺
        a < 4 && plot!(p, [a, a + 1], [b, b], color=:gray80)
        b < 4 && plot!(p, [a, a], [b, b + 1], color=:gray80)
    end
    scatter!(p, 1:4, fill(4, 4), marker=:circle, ms=10, color=:white, msc=:gray70)   # 未測定(y=4)
    scatter!(p, x, y, marker=:circle, ms=16, zcolor=v .* 100, c=:RdBu, clims=(-2, 2),
             msc=:black, colorbar=true, colorbar_title="<Sz_i>  x 100")
    @printf("symbreak %-8s %-4s max|Sz| = %.3e\n", geo_, rho, maximum(abs.(v)))
    p
end
p7a = szmap("torus", "ED", "4x4 torus (32 qubits), U = 8, half filling\n(a) exact (ED):  <Sz_i> = 0")
p7b = szmap("torus", "DMRG", "\n(b) DMRG, chi = 1024:  staggered <Sz_i>")
annotate!(p7a, 2.5, 4.42, text("y = 4: not in the data", 6, :gray40))
sE, zzE = val("torus", "ED", "SzSz nb"); _, xxE = val("torus", "ED", "SxSx nb")
sD, zzD = val("torus", "DMRG", "SzSz nb"); _, xxD = val("torus", "DMRG", "SxSx nb")
p7c = plot(xlabel="bond index (site of the window)", ylabel="nearest-neighbour spin correlation",
           title="\n(c) isotropy is lost:  <SzSz> != <SxSx>", legend=:right,
           ylims=(-0.1025, -0.0935))
hline!(p7c, [mean(zzE)], color=:black, lw=2, label="ED:  SzSz = SxSx")
scatter!(p7c, sD, zzD, color=:firebrick, marker=:circle, ms=5, msw=0, label="DMRG  SzSz")
scatter!(p7c, sD, xxD, color=:steelblue, marker=:square, ms=5, msw=0, label="DMRG  SxSx")
@printf("symbreak: ED SzSz %.4f SxSx %.4f | DMRG SzSz %.4f SxSx %.4f\n",
        mean(zzE), mean(xxE), mean(zzD), mean(xxD))
out(plot(p7a, p7b, p7c, layout=grid(1, 3, widths=[0.3, 0.3, 0.4]), size=(1400, 420), margin=6Plots.mm), "crm_new_fig_symbreak.png")
