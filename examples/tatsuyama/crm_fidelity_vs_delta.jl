# ============================================================
# 部分系の「忠実度」で同じ解析をする(§2m のトレース距離版との比較)
#
# §2m はトレース距離 D_tr で解析し、単一Pauli列では上界 |Δ| ≤ 2 D_tr が
# 厳密に飽和することを見た。忠実度でも同じことが言えるかを調べる。
#
# 理論的な予想:
#   Fuchs-van de Graaf より D_tr ≤ √(1-F) なので、忠実度版の上界
#       |Δ| ≤ 2 D_tr ≤ 2√(1-F)
#   は必ずトレース距離版より緩い。さらに半充填のサイト内2量子ビットでは
#   RDM が二重占有 D ひとつで決まり、差が e=(δ,-δ,-δ,δ) なので
#       1-F ≃ (1/4) Σ e_i²/p_i = (δ²/2)(1/D + 1/(1/2-D))
#   となり
#       |Δ|/(2√(1-F)) = 2 / √( (1/2)(1/D + 1/(1/2-D)) )
#   と **D に依存する係数**が付く。トレース距離版は D によらず厳密に 1 だった。
#   この違いが実測に出るかを確かめる。
#
# 入力: crm_support_fidelity_results.tsv(crm_support_fidelity.jl の出力)
# 実行: julia --project=Hubbard_MPS_Env_v2 crm_fidelity_vs_delta.jl
# ============================================================
using DelimitedFiles, Statistics, Printf

raw, hdr = readdlm(joinpath(@__DIR__, "crm_support_fidelity_results.tsv"), '\t'; header=true)
c = Dict(String(h) => j for (j, h) in enumerate(vec(hdr)))
col(n) = raw[:, c[n]]
L   = Int.(col("L"));            dop = Float64.(col("doping"))
site= Int.(col("site"));         ob  = String.(col("observable"))
pr  = String.(col("prior"));     tru = Float64.(col("true"))
Δ   = Float64.(col("Delta"));    aΔ  = Float64.(col("absDelta"))
Fs  = Float64.(col("F_supp"));   iFs = Float64.(col("infid_supp"))
Ds  = Float64.(col("D_supp"));   ns  = Int.(col("nsupp"))
Dw  = Float64.(col("D_trace"));  N = length(L)

PRIORS = ["chi2","chi4","chi8","chi16","chi32","UHF","UHFsym"]
OBS    = ["ZZ onsite","ZZ up-up r=1","DoubleOcc","SzSz r=1","hop up r=1","Sz","n"]

function corr(x, y)
    length(x) < 5 && return NaN
    mx, my = mean(x), mean(y)
    sx = sqrt(sum((a-mx)^2 for a in x)); sy = sqrt(sum((b-my)^2 for b in y))
    (sx*sy) > 0 ? sum((a-mx)*(b-my) for (a,b) in zip(x,y))/(sx*sy) : NaN
end

println("="^100)
println("(1) 上界のタイトさ: トレース距離版 |Δ|/(2 D_台) vs 忠実度版 |Δ|/(2√(1-F_台))")
println("    半充填のみ。単一Pauli列ではトレース距離版が厳密に飽和するはず")
@printf("%-16s %5s %20s %22s %16s\n", "観測量", "台", "|Δ|/(2 D_台)", "|Δ|/(2√(1-F_台))", "飽和(>0.99)")
for o in OBS
    idx = findall(i -> ob[i]==o && dop[i]==0.0 && aΔ[i]>0 && Ds[i]>0 && iFs[i]>0, 1:N)
    isempty(idx) && continue
    rt = [aΔ[i]/(2Ds[i]) for i in idx]
    rf = [aΔ[i]/(2sqrt(iFs[i])) for i in idx]
    @printf("%-16s %5d %20.4f %22.4f %15.0f%%\n", o, ns[idx[1]], median(rt), median(rf),
            100*count(x -> x>0.99, rt)/length(rt))
end

println("\n", "="^100)
println("(2) 忠実度版の係数は二重占有 D に依存するか(理論予想: 2/√((1/2)(1/D+1/(1/2-D))))")
println("    ZZ onsite、半充填。サイトごとの D は <ZZ>=4D-1 から復元する")
@printf("%-9s %12s %14s %16s %16s\n", "prior", "D 中央値", "予想の係数", "実測(忠実度版)", "実測(距離版)")
for p in ("chi2","chi4","chi8","chi32","UHFsym")
    idx = findall(i -> ob[i]=="ZZ onsite" && pr[i]==p && dop[i]==0.0 && aΔ[i]>0 && iFs[i]>0, 1:N)
    isempty(idx) && continue
    Dv = [(tru[i]+1)/4 for i in idx]                     # <ZZ> = 4D - 1
    pred = [2/sqrt(0.5*(1/d + 1/(0.5-d))) for d in Dv]
    rf = [aΔ[i]/(2sqrt(iFs[i])) for i in idx]
    rt = [aΔ[i]/(2Ds[i]) for i in idx]
    @printf("%-9s %12.5f %14.4f %16.4f %16.4f\n", p, median(Dv), median(pred), median(rf), median(rt))
end

println("\n", "="^100)
println("(3) prior を固定した相関: log|Δ| 対 log(1-F_台) と log D_台")
@printf("%-16s", "観測量")
for p in PRIORS; @printf("%9s", p); end
println("   ← 忠実度版")
for o in OBS
    @printf("%-16s", o)
    for p in PRIORS
        idx = findall(i -> ob[i]==o && pr[i]==p && aΔ[i]>0 && iFs[i]>0, 1:N)
        cv = corr([log(aΔ[i]) for i in idx], [log(iFs[i]) for i in idx])
        @printf("%9s", isnan(cv) ? "—" : @sprintf("%.3f", cv))
    end
    println()
end
@printf("\n%-16s", "観測量")
for p in PRIORS; @printf("%9s", p); end
println("   ← 距離版(比較)")
for o in OBS
    @printf("%-16s", o)
    for p in PRIORS
        idx = findall(i -> ob[i]==o && pr[i]==p && aΔ[i]>0 && Ds[i]>0, 1:N)
        cv = corr([log(aΔ[i]) for i in idx], [log(Ds[i]) for i in idx])
        @printf("%9s", isnan(cv) ? "—" : @sprintf("%.3f", cv))
    end
    println()
end

println("\n", "="^100)
println("(4) 1-F_台 と D_台 の関係(対数回帰の傾き。2 なら 1-F ∝ D²)")
idx = findall(i -> Ds[i]>0 && iFs[i]>0, 1:N)
seen = Set{Tuple{Int,Float64,Int,String}}(); u = Int[]
for i in idx
    k = (L[i], dop[i], site[i], pr[i])
    k in seen || (push!(seen, k); push!(u, i))
end
x = [log10(Ds[i]) for i in u]; y = [log10(iFs[i]) for i in u]
mx, my = mean(x), mean(y)
a = sum((xi-mx)*(yi-my) for (xi,yi) in zip(x,y)) / sum((xi-mx)^2 for xi in x)
res = [yi - (a*xi + (my-a*mx)) for (xi,yi) in zip(x,y)]
@printf("    傾き = %.4f   相関 = %.5f   残差sd = %.4f 桁   (%d 組)\n",
        a, corr(x,y), std(res), length(u))

println("\n", "="^100)
println("(5) prior ごとの部分系忠実度(半充填、ZZ onsite の台=2量子ビット)")
@printf("%-9s %14s %14s %14s %14s\n", "prior", "F_台", "1-F_台", "D_台", "|Δ|(ZZ)")
for p in PRIORS
    idx = findall(i -> ob[i]=="ZZ onsite" && pr[i]==p && dop[i]==0.0, 1:N)
    isempty(idx) && continue
    @printf("%-9s %14.6f %14.3e %14.3e %14.3e\n", p,
            median(Fs[idx]), median(iFs[idx]), median(Ds[idx]), median(aΔ[idx]))
end
