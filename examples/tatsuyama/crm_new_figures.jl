# ============================================================
# 今回のセッションで出した結果の図(§0.5b, §2i, §2j, §2k, §2L, §2m)
#   julia --project=Hubbard_MPS_Env_v2 crm_new_figures.jl
# 既存の図(§1〜§9)は crm_paper_figures.jl 側にある。
# ラベルは英語(日本語フォントが環境依存のため。既存図と揃える)。
# ============================================================
using DelimitedFiles, Statistics, Printf, Plots
gr(fontfamily="Helvetica", legendfontsize=8, guidefontsize=9,
   titlefontsize=10, tickfontsize=8)


"""StatsPlots に依存せずに横並びの棒グラフを描く。"""
function gbar(labels, mat; series, colors, kw...)
    n, k = size(mat)
    w = 0.8/k
    xs = 1:n
    plt = plot(; xticks=(collect(xs), labels), kw...)
    for j in 1:k
        bar!(plt, xs .+ (j - (k+1)/2)*w, mat[:, j], bar_width=w*0.92,
             label=series[j], color=colors[j], linewidth=0)
    end
    plt
end

const DIR = @__DIR__
function load(f)
    raw, hdr = readdlm(joinpath(DIR, f), '\t'; header=true)
    c = Dict(String(h) => j for (j, h) in enumerate(vec(hdr)))
    (raw = raw, col = name -> raw[:, c[name]])
end
num(v) = Float64.(v); str(v) = String.(v)

# ------------------------------------------------------------
# Fig A: 総ショット数を固定した配分 (§0.5b)
# ------------------------------------------------------------
d = load("crm_shot_split_results.tsv")
B = Int.(d.col("B")); nm = Int.(d.col("n_m")); ob = str(d.col("observable"))
pr = str(d.col("prior")); sds = num(d.col("sd_std")); sdc = num(d.col("sd_crm"))
sdst = num(d.col("sd_std_theory")); sdct = num(d.col("sd_crm_theory")); G = num(d.col("G_emp"))

sel(f) = findall(f, 1:length(B))
nms = sort(unique(nm[sel(i -> B[i] == 5000)]))
getv(v, p) = [ (k = sel(i -> B[i]==5000 && nm[i]==m && ob[i]=="ZZ onsite" && pr[i]==p);
                isempty(k) ? NaN : v[k[1]]) for m in nms ]

pA1 = plot(xscale=:log10, yscale=:log10, xlabel="n_m  (n_u = B/n_m)",
           ylabel="standard error of the estimate", legend=:topleft,
           title="Fixed total budget B = 5000\n(ZZ onsite, L=32, U=4)")
plot!(pA1, nms, getv(sds, "chi8"), marker=:circle, lw=2, color=:steelblue,
      label="standard shadow")
plot!(pA1, nms, getv(sdst, "chi8"), ls=:dash, lw=1, color=:steelblue, label="theory")
plot!(pA1, nms, getv(sdc, "chi2"), marker=:diamond, lw=2, color=:darkorange,
      label="CRM (chi_p = 2, coarse prior)")
plot!(pA1, nms, getv(sdc, "chi8"), marker=:square, lw=2, color=:firebrick,
      label="CRM (chi_p = 8, good prior)")
plot!(pA1, nms, getv(sdct, "chi8"), ls=:dash, lw=1, color=:firebrick, label="")
annotate!(pA1, 300, 0.62, text("standard: 20x worse", 7, :right, :steelblue))
annotate!(pA1, 300, 0.043, text("CRM (chi_p=8): flat", 7, :right, :firebrick))

pA2 = plot(xscale=:log10, yscale=:log10, xlabel="n_m", ylabel="gain G",
           legend=:topleft, title="G grows only because the\nbaseline is handicapped")
plot!(pA2, nms, getv(G, "chi8"), marker=:square, lw=2, color=:firebrick, label="chi_p = 8")
plot!(pA2, nms, getv(G, "chi2"), marker=:diamond, lw=2, color=:darkorange, label="chi_p = 2")
hline!(pA2, [1.559], ls=:dot, color=:black, label="bound 1/(1-<P>^2) at n_m=1")

figA = plot(pA1, pA2, layout=(1,2), size=(1000, 400), margin=5Plots.mm)
savefig(figA, joinpath(DIR, "crm_new_fig_shotsplit.png"))

# ------------------------------------------------------------
# Fig B: 幾何 — シリンダー vs トーラス (§2k)
# ------------------------------------------------------------
d = load("crm_2d_ed_W4L4_U8p0.tsv")
geo = str(d.col("geometry")); rho = str(d.col("rho")); ob = str(d.col("observable"))
pr = str(d.col("prior")); G = num(d.col("G"))
med(g, o, p) = (k = findall(i -> geo[i]==g && rho[i]=="ED" && ob[i]==o && pr[i]==p,
                            1:length(geo)); isempty(k) ? NaN : median(G[k]))
obs_b = ["DoubleOcc", "ZZ onsite", "ZZ up-up nb", "SzSz nb", "SxSx nb"]
lab_b = ["double occ", "ZZ onsite", "ZZ nb", "SzSz nb", "SxSx nb"]

r_cyl = [med("cylinder", o, "UHFsym")/med("cylinder", o, "chi32") for o in obs_b]
r_tor = [med("torus",    o, "UHFsym")/med("torus",    o, "chi32") for o in obs_b]
pB1 = gbar(lab_b, hcat(r_cyl, r_tor); series=["cylinder","torus"],
           colors=[:steelblue,:firebrick],
           ylabel="G(mean-field, sym) / G(chi=32 MPS)", legend=:topleft,
           title="4x4 Hubbard, 32 qubits, exact reference\nmean-field overtakes MPS on the torus",
           xrotation=20)
hline!(pB1, [1.0], ls=:dash, color=:black, label="equal")

# DMRG のエネルギー誤差(同じ chi、同じサイト数、境界条件だけ違う)
sys_lab = ["4x3\nU=8", "2x6\nU=8", "2x6\nU=4", "4x4\nU=8"]
dE_cyl  = [9.53e-7, 2.30e-10, 1.29e-10, 5.86e-4]
dE_tor  = [4.31e-4, 8.27e-7,  9.90e-6,  8.71e-2]
# 対数軸で棒を描くと基準が -inf になって壊れるので log10 を線形軸に描く
pB2 = gbar(sys_lab, hcat(log10.(dE_cyl), log10.(dE_tor));
           series=["cylinder / open","torus / periodic"],
           colors=[:steelblue,:firebrick],
           ylabel="log10 |E_DMRG - E_exact|  at fixed chi",
           legend=:bottomright, ylims=(-11.5, 1.9),
           title="Cost of periodicity to the MPS side\n(UHF does not move: m = 0.4546 -> 0.4539)")
for (i, r) in enumerate(dE_tor ./ dE_cyl)   # 悪化率は棒の上、パネル上端に揃えて置く
    annotate!(pB2, i, 0.75, text(@sprintf("x%.0f", r), 7, :center, :firebrick))
end

figB = plot(pB1, pB2, layout=(1,2), size=(1050, 420), margin=5Plots.mm)
savefig(figB, joinpath(DIR, "crm_new_fig_geometry.png"))

# ------------------------------------------------------------
# Fig C: 対称性の回復 (§2i, §2j)
# ------------------------------------------------------------
d = load("crm_transverse_results.tsv")
ob = str(d.col("observable")); pr = str(d.col("prior"))
tru = num(d.col("true")); pv = num(d.col("prior_val")); Ls = Int.(d.col("L"))
pick(o, p) = (k = findall(i -> ob[i]==o && pr[i]==p && Ls[i]==8, 1:length(ob));
              isempty(k) ? NaN : median(abs.(pv[k])))
tv(o) = (k = findall(i -> ob[i]==o && pr[i]=="exact" && Ls[i]==8, 1:length(ob));
         median(abs.(tru[k])))
vals = hcat([tv("SzSz r=1"), pick("SzSz r=1","UHF"), pick("SzSz r=1","UHFsym")],
            [tv("SxSx r=1"), pick("SxSx r=1","UHF"), pick("SxSx r=1","UHFsym")])
pC1 = gbar(["true state","UHF","UHF-sym"], vals;
           series=["<SzSz> (longitudinal)","<SxSx> (transverse)"],
           colors=[:steelblue,:darkorange], ylabel="|correlation|",
           title="SU(2) forces these equal.\nUHF gets them wrong by 7.1x", legend=:topleft)

d = load("crm_symmetry_prior_results.tsv")
ob = str(d.col("observable")); pr = str(d.col("prior")); G = num(d.col("G")); Ls = Int.(d.col("L"))
mg(o, p) = (k = findall(i -> ob[i]==o && pr[i]==p && Ls[i]==16, 1:length(ob));
            isempty(k) ? NaN : median(G[k]))
obs_c = ["ZZ up-up r=1", "SzSz r=1", "hop up r=1"]
plain = [mg(o, "chi4") for o in obs_c]; tavg = [mg(o, "chi4_tavg") for o in obs_c]
exact = [mg(o, "exact") for o in obs_c]
pC2 = gbar(["ZZ nb","SzSz nb","hop"], hcat(plain, tavg, exact);
           series=["chi_p = 4","chi_p = 4 + translation avg","exact prior"],
           colors=[:gray70,:seagreen,:black], ylabel="gain G", legend=:outerbottom,
           legendcolumns=3, ylims=(0, 40),
           title="Translation averaging is free\nand reaches the exact prior (L=16, PBC)")

figC = plot(pC1, pC2, layout=(1,2), size=(1050, 420), margin=5Plots.mm)
savefig(figC, joinpath(DIR, "crm_new_fig_symmetry.png"))

# ------------------------------------------------------------
# Fig D: 部分系の忠実度 (§2m)
# ------------------------------------------------------------
d = load("crm_support_fidelity_results.tsv")
ob = str(d.col("observable")); pr = str(d.col("prior"))
aD = num(d.col("absDelta")); Ds = num(d.col("D_supp")); Dw = num(d.col("D_trace"))
kk(p, o) = findall(i -> pr[i]==p && ob[i]==o && aD[i] > 0 && Ds[i] > 0, 1:length(ob))

pD1 = plot(xscale=:log10, yscale=:log10, xlabel="trace distance on the support  D",
           ylabel="|Delta|", legend=:topleft,
           title="Support-restricted: single Pauli strings\nsit exactly on |Delta| = 2D")
xs = 10 .^ range(-9, 0; length=50)
plot!(pD1, xs, 2 .* xs, ls=:dash, color=:black, lw=1.5, label="|Delta| = 2D (bound)")
for (p, c, m) in (("chi8", :steelblue, :circle), ("chi32", :seagreen, :utriangle),
                  ("UHF", :firebrick, :diamond), ("UHFsym", :darkorange, :square))
    k = kk(p, "ZZ onsite")
    scatter!(pD1, Ds[k], aD[k], ms=3, msw=0, color=c, marker=m, alpha=0.6, label=p)
end

# 上界にどれだけ近いか(1.0 なら飽和)を prior ごとに並べる
prs = ["chi2","chi4","chi8","chi16","chi32","UHF","UHFsym"]
pD2 = plot(xlabel="", ylabel="|Delta| / (2D)   (1.0 = bound saturated)",
           legend=false, xticks=(1:length(prs), prs), ylims=(-0.05, 1.15),
           title="Single Pauli string (ZZ onsite), half filling\n" *
                 "MPS priors saturate the bound; UHF does not")
for (j, p) in enumerate(prs)
    k = kk(p, "ZZ onsite")
    r = aD[k] ./ (2 .* Ds[k])
    scatter!(pD2, fill(j, length(r)) .+ 0.16 .* (rand(length(r)) .- 0.5), r,
             ms=2.5, msw=0, alpha=0.45,
             color = (p == "UHF" ? :firebrick : (p == "UHFsym" ? :darkorange : :steelblue)))
    scatter!(pD2, [j], [median(r)], ms=7, marker=:hline, color=:black)
end
hline!(pD2, [1.0], ls=:dash, color=:black)
annotate!(pD2, 6, 0.35, text("distance points along\nthe broken symmetry", 7, :center, :firebrick))

figD = plot(pD1, pD2, layout=(1,2), size=(1050, 430), margin=5Plots.mm)
savefig(figD, joinpath(DIR, "crm_new_fig_supportfid.png"))

# ------------------------------------------------------------
# Fig E: 1つのシャドウを全観測量で共有 (§2L)
# ------------------------------------------------------------
d = load("crm_global_shadow_results.tsv")
ob = str(d.col("observable")); pr = str(d.col("prior"))
nu = Int.(d.col("n_u")); G = num(d.col("G")); nt = Int.(d.col("nterms"))
gg(o, p) = (k = findall(i -> ob[i]==o && pr[i]==p && nu[i]==50, 1:length(ob));
            isempty(k) ? NaN : G[k[1]])
obs_e = ["ZZ onsite", "DoubleOcc", "SzSz r=1", "hop r=1", "構造因子 S(π)", "全エネルギー"]
lab_e = ["ZZ onsite\n(1 term)", "double occ\n(4)", "SzSz\n(4)", "hop\n(2)",
         "S(pi)\n(248)", "total energy\n(84)"]
g2 = [gg(o, "chi2") for o in obs_e]; g8 = [gg(o, "chi8") for o in obs_e]
# log10(G) を 0 基準で描く。負の棒がそのまま「損」になる
pE = gbar(lab_e, hcat(log10.(g2), log10.(g8));
          series=["chi_p = 2 (coarse)","chi_p = 8"],
          colors=[:darkorange,:steelblue],
          ylabel="log10 G     (below 0 = CRM loses)", legend=:topleft,
          title="One shadow, all observables (L=8, real protocol)\n" *
                "at chi_p=2 every part gains, but the 84-term total energy loses")
hline!(pE, [0.0], ls=:dash, color=:black, label="")
annotate!(pE, 5.72, 0.14, text(@sprintf("G = %.2f\n(loses)", g2[6]), 7, :center, :firebrick))
annotate!(pE, 4.0, 0.06, text("hop: G = 1.00 exactly\n(prior predicts 0)", 6, :center, :gray40))
figE = plot(pE, size=(760, 430), margin=5Plots.mm)
savefig(figE, joinpath(DIR, "crm_new_fig_globalshadow.png"))

println("figures saved:")
for f in ("crm_new_fig_shotsplit.png","crm_new_fig_geometry.png","crm_new_fig_symmetry.png",
          "crm_new_fig_supportfid.png","crm_new_fig_globalshadow.png")
    println("  ", f)
end
