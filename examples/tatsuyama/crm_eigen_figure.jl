# 仮説「HF 基底状態が観測量の固有状態に近いほど利得が大きい」を検証する図。
#   (a) 仮説そのもの: U を上げて UHF を固有状態に近づけると、利得は観測量によって逆を向く
#   (b) UHF は実際どこにいるか: (|<P>_ρ|, |<P>_σ|) 平面。対角線が理想、|y| = 2|x| が損益分岐
#   (c) ρ を固定して σ を動かす: 利得は <P>_σ = <P>_ρ で最大になり、固有状態では最大にならない
#   (d) σ を固定して ρ を動かす: 固有状態の prior が得をするのは |<P>_ρ| > 1/2 のときだけ
#   julia --project=. examples/tatsuyama/crm_eigen_figure.jl
using DelimitedFiles, Printf, Plots
gr(fontfamily="Helvetica", legendfontsize=7, guidefontsize=9, titlefontsize=10, tickfontsize=8)
const DIR = @__DIR__
A, h = readdlm(joinpath(DIR, "crm_edfid_all.tsv"), '\t'; header=true)
h = vec(h); col(n) = A[:, findfirst(==(n), h)]

# 参照パレットの先頭3色(全ペアで検証済みの枠)。観測量ごとに固定し、全パネルで同じ色にする
hexc(s) = RGB(parse(Int, s[2:3]; base=16)/255, parse(Int, s[4:5]; base=16)/255, parse(Int, s[6:7]; base=16)/255)
const C_ON = hexc("#2a78d6")    # ZZ onsite       (電荷)
const C_NB = hexc("#eb6834")    # ZZ up-up nb     (スピン, 最近接)
const C_Z1 = hexc("#1baf7a")    # Z_up 単一サイト (スピン, 局所磁化)
const INK  = hexc("#52514e")
const LOGT1 = [1e-2, 1e-1, 1, 10, 100, 1000]   # 全ての桁に目盛り(GR が 10^n で描く)
const LOGT2 = [0.1, 1, 10, 100]

gain(x, y; nA=2, nm=100) = begin
    K = 3.0^nA - 1; vs = 3.0^nA*(1 - x^2)/nm
    (K*x^2 + vs) / (K*(x - y)^2 + vs)
end

# (W, geometry, LX, U, prior) => observable => (ρ の値, σ の値)
V = Dict{Tuple{Int,String,Int,Float64,String}, Dict{String,Tuple{Float64,Float64}}}()
for r in 1:size(A, 1)
    col("exact")[r] == "yes" || continue
    k = (Int(col("W")[r]), string(col("geometry")[r]), Int(col("LX")[r]), Float64(col("U")[r]), string(col("prior")[r]))
    get!(V, k, Dict{String,Tuple{Float64,Float64}}())[string(col("observable")[r])] =
        (Float64(col("true")[r]), Float64(col("prior_val")[r]))
end
# 単一サイトの Z_up = 1 - n - 2 S^z(台 1 量子ビット)
for d in values(V)
    haskey(d, "n") && haskey(d, "Sz") || continue
    (nr, ns), (sr, ss) = d["n"], d["Sz"]
    d["Zup"] = (1 - nr - 2sr, 1 - ns - 2ss)
end
const SERIES = [("ZZ onsite", 2, C_ON, "ZZ onsite  (charge)"),
                ("ZZ up-up nb", 2, C_NB, "Z↑Z↑ nearest neighbour  (spin)"),
                ("Zup", 1, C_Z1, "Z↑ single site  (local moment)")]
lattices = sort(unique((k[1], k[2], k[3]) for k in keys(V)))
Us = [2.0, 4.0, 8.0, 12.0]

# ---------------------------------------------------------------- (a)
pa = plot(yscale=:log10, xlims=(-0.02, 1.02), ylims=(4e-3, 1.5e3), yticks=LOGT1, legend=:topleft,
          xlabel="|⟨P⟩σ| of UHF      (1 = eigenstate of P)", ylabel="gain G",
          title="(a) push UHF toward the eigenstate (raise U): the gain goes opposite ways")
hline!(pa, [1.0], color=:black, ls=:dot, lw=1.2, label="")
for (obs, nA, c, lab) in SERIES
    for lat in lattices
        pts = [(abs(V[(lat..., U, "UHF")][obs][2]), gain(V[(lat..., U, "UHF")][obs]...; nA))
               for U in Us if haskey(V, (lat..., U, "UHF"))]
        length(pts) == 4 || continue
        hero = lat == (1, "cylinder", 64)
        plot!(pa, first.(pts), last.(pts); color=c, lw=hero ? 2.2 : 0.8, alpha=hero ? 1.0 : 0.22,
              marker=hero ? :circle : :none, ms=7, msw=0, label=hero ? lab : "")
        if hero
            for (U, (px, py)) in zip(Us, pts)   # U=8 と U=12 は近いので上下に振り分ける
                below = U == 8.0 && obs != "ZZ onsite"
                annotate!(pa, px + (below ? -0.02 : 0.0), below ? py/2.1 : py*2.1,
                          text(@sprintf("U=%d", U), 7, INK, below ? :right : :center))
            end
        end
    end
end
annotate!(pa, 0.02, 1.5e-2, text("bold: 1D chain L=64 (128 qubits)\nfaint: the other 19 lattices", 7, :left, INK))

# ---------------------------------------------------------------- (b)
pb = plot(xlims=(-0.02, 1.02), ylims=(-0.02, 1.02), aspect_ratio=:equal, legend=:bottomright,
          xlabel="|⟨P⟩ρ|   (truth's closeness to the eigenstate)",
          ylabel="|⟨P⟩σ|   (UHF's closeness to the eigenstate)",
          title="(b) where UHF actually sits  (all exact systems, U = 2–12)")
plot!(pb, [0, 0.5, 0], [0, 1, 1]; seriestype=:shape, fillcolor=:gray, fillalpha=0.13, lw=0, label="")
plot!(pb, [0, 1], [0, 1]; color=:black, lw=1.2, label="")
plot!(pb, [0, 0.5], [0, 1]; color=:black, lw=1.2, ls=:dash, label="")
for (obs, nA, c, lab) in SERIES
    pts = [(abs(d[obs][1]), abs(d[obs][2])) for (k, d) in V if k[5] == "UHF" && haskey(d, obs) && d[obs][1]*d[obs][2] >= 0]
    scatter!(pb, first.(pts), last.(pts); color=c, ms=5.5, msw=0, alpha=0.75, label=lab)
end
annotate!(pb, 0.05, 0.80, text("CRM loses\n(G < 1)", 8, :left, INK))
annotate!(pb, 0.10, 0.30, text("break-even  |⟨P⟩σ| = 2|⟨P⟩ρ|", 7, :center, INK, rotation=63.4))
annotate!(pb, 0.72, 0.67, text("σ = ρ   (G = Gmax)", 7, :center, INK, rotation=45))
annotate!(pb, 0.63, 0.30, text("prior less extreme than truth:\nnever loses", 7, :left, INK))

# ---------------------------------------------------------------- (c)
lat_c, U_c = (1, "cylinder", 64), 8.0
pc = plot(yscale=:log10, xlims=(-1.04, 1.04), ylims=(0.1, 450), yticks=LOGT2, legend=:topright,
          xlabel="⟨P⟩σ   (the prior's value)", ylabel="gain G",
          title=@sprintf("(c) truth fixed (1D L=64, U=%d): the peak is at σ = ρ, not at the eigenstate", U_c))
hline!(pc, [1.0], color=:black, ls=:dot, lw=1.2, label="")
vline!(pc, [-1.0], color=:gray, lw=5, alpha=0.25, label="")
annotate!(pc, -0.97, 0.16, text("eigenstate\n⟨P⟩σ = -1", 7, :left, INK))
for (obs, nA, c, lab) in SERIES[1:2]
    x = V[(lat_c..., U_c, "UHF")][obs][1]
    ys = range(-1, 1; length=2001)
    plot!(pc, ys, [gain(x, y; nA) for y in ys]; color=c, lw=1.8, label=@sprintf("%s:  ⟨P⟩ρ = %.3f", obs, x))
    vline!(pc, [x]; color=c, ls=:dash, lw=1, label="")
    for p in ("UHF", "chi2", "chi4", "chi8", "chi16", "chi32")
        haskey(V, (lat_c..., U_c, p)) || continue
        y = V[(lat_c..., U_c, p)][obs][2]
        scatter!(pc, [y], [gain(x, y; nA)]; color=c, ms=p == "UHF" ? 12 : 7.5,
                 marker=p == "UHF" ? :star5 : :circle, msw=1.2, msc=:white, label="")
    end
end
# 点の位置は表から既知(χ=2: -1.000/0.000, χ=4: -0.919/-0.491, UHF: -0.881/-0.940, χ≥8: ≈<P>ρ)
annotate!(pc, -0.975, 19.0, text("χ=2 (eigenstate)", 7, :left, INK))
annotate!(pc, -0.905, 66.0, text("χ=4", 7, :left, INK))
annotate!(pc, -0.895, 320.0, text("UHF", 7, :right, INK))
annotate!(pc, -0.835, 330.0, text("χ≥8", 7, :left, INK))
annotate!(pc, -0.915, 2.3, text("UHF", 7, :left, INK))
annotate!(pc, -0.47, 50.0, text("χ≥4", 7, :left, INK))
annotate!(pc, 0.03, 0.62, text("χ=2", 7, :left, INK))
scatter!(pc, [NaN], [NaN]; color=:gray, marker=:star5, ms=10, msw=0, label="UHF prior")
scatter!(pc, [NaN], [NaN]; color=:gray, marker=:circle, ms=7, msw=0, label="truncated MPS prior (χ)")

# ---------------------------------------------------------------- (d)
pd = plot(yscale=:log10, xlims=(0, 1), ylims=(8e-3, 1.5e3), yticks=LOGT1, legend=:topleft,
          xlabel="|⟨P⟩ρ|   (truth's closeness to the eigenstate)", ylabel="gain G",
          title="(d) prior fixed near the eigenstate: it helps only if |⟨P⟩ρ| > |⟨P⟩σ|/2")
hline!(pd, [1.0], color=:black, ls=:dot, lw=1.2, label="")
vline!(pd, [0.5], color=:gray, lw=1, ls=:dash, label="")
xs = range(0.005, 0.995; length=800)
plot!(pd, xs, [gain(x, 1.0) for x in xs]; color=C_ON, lw=1.6, alpha=0.7, label="theory, ⟨P⟩σ = ±1 exactly")
p2 = [(abs(d["ZZ onsite"][1]), gain(d["ZZ onsite"]...)) for (k, d) in V if k[5] == "chi2"]
scatter!(pd, first.(p2), last.(p2); color=C_ON, ms=5.5, msw=0, alpha=0.8,
         label=@sprintf("χ=2 MPS, ZZ onsite: exactly an eigenstate (%d systems, all U)", length(p2)))
yu = abs(V[(1, "cylinder", 64, 12.0, "UHF")]["ZZ up-up nb"][2])
plot!(pd, xs, [gain(x, yu) for x in xs]; color=C_NB, lw=1.6, alpha=0.7, label=@sprintf("theory, |⟨P⟩σ| = %.4f", yu))
pu = [(abs(d["ZZ up-up nb"][1]), gain(d["ZZ up-up nb"]...)) for (k, d) in V if k[5] == "UHF" && k[4] == 12.0 && k[1] == 1]
scatter!(pd, first.(pu), last.(pu); color=C_NB, ms=11, marker=:star5, msw=0.8, msc=:white,
         label=@sprintf("UHF, Z↑Z↑ nb, U=12, 1D (%d systems: strong & weak bonds)", length(pu)))
annotate!(pd, 0.52, 4e-2, text("|⟨P⟩ρ| = 1/2:\nthe truth must put > 75% weight\non the prior's eigenspace", 7, :left, INK))

fig = plot(pa, pb, pc, pd; layout=(2, 2), size=(1560, 1240), margin=6Plots.mm)
out = joinpath(DIR, "crm_new_fig_eigen.png")
savefig(fig, out); println("書き出し: ", out)
