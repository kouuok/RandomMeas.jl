# ============================================================
# n_u と n_m をいくつにすべきか — 実測データからの設計
#
# 利得法則から、n_u と n_m の役割はまったく違う:
#   Var[標準] = (1/n_u)[(3^|A|-1)<P>^2 + v_s],  v_s = 3^|A|(1-<P>^2)/n_m
#   Var[CRM]  = (1/n_u)[(3^|A|-1)Δ^2   + v_s]
# なので n_u は比 G に現れず(推定量自身の精度だけを決める)、n_m だけが G を動かす。
#
# 本スクリプトは §2b の 17,210 点の実測 (<P>, Δ) から
#   (1) 観測量ごとの n_m* = 3^|A|(1-<P>^2)/[(3^|A|-1)Δ^2]
#   (2) n_m=100 での「prior の質に対する感度」|dlogG/dlogΔ|
#   (3) 総ショット数を固定したときの最適な n_m(基底切替コスト c 付き)
#   (4) 測定される G のばらつきが n_u と n_repeat にどう依存するか
# を出して、設計の根拠を数値で示す。
#
# 実行: julia --project=Hubbard_MPS_Env_v2 crm_budget_design.jl
# ============================================================
using DelimitedFiles, Statistics, Printf

const ABSA = Dict("n"=>1, "Sz"=>1, "DoubleOcc"=>2, "ZZ onsite"=>2,
                  "ZZ up-up r=1"=>2, "SzSz r=1"=>2, "hop up r=1"=>3)
# 利得法則が閉形式で使えるのは単一Pauli列だけ(多項観測量は項ごとのΔに依存)
const PURE = ("ZZ onsite", "ZZ up-up r=1")

raw, hdr = readdlm(joinpath(@__DIR__, "crm_site_resolved_results.tsv"), '\t'; header=true)
c = Dict(String(h) => j for (j, h) in enumerate(vec(hdr)))
col(name) = raw[:, c[name]]
L = Int.(col("L")); dop = Float64.(col("doping")); cp = Int.(col("chi_prior"))
ob = String.(col("observable")); tru = Float64.(col("true")); Δ = abs.(Float64.(col("Delta")))
Gem = Float64.(col("G_emp")); Glo = Float64.(col("G_lo")); Ghi = Float64.(col("G_hi"))
Gmx = Float64.(col("G_max"))

nmstar(a, P, d) = 3.0^a * (1 - P^2) / ((3.0^a - 1) * d^2)
gain(a, P, d, nm) = ((3.0^a-1)*P^2*nm + 3.0^a*(1-P^2)) / ((3.0^a-1)*d^2*nm + 3.0^a*(1-P^2))

println("="^96)
println("(1) 必要ショット数 n_m* — 「prior の誤差がショットノイズと釣り合う」点")
println("    n_m ≪ n_m* : ショットノイズ律速。G は天井に貼り付き、prior の質を見ていない")
println("    n_m = n_m* : G = G_max/2")
println("    n_m ≫ n_m* : prior 律速。G → (<P>/Δ)^2")
@printf("%-16s %10s", "観測量(単一Pauli列)", "χ_p")
for x in (2,4,8,16,32); @printf("%13d", x); end; println()
for o in PURE
    @printf("%-16s %10s", o, "n_m* 中央値")
    for x in (2,4,8,16,32)
        idx = findall(i -> ob[i]==o && cp[i]==x && L[i]==32 && dop[i]==0.0 && Δ[i]>0, 1:length(L))
        v = [nmstar(ABSA[o], tru[i], Δ[i]) for i in idx]
        @printf("%13.3g", isempty(v) ? NaN : median(v))
    end
    println()
end

println("\n", "="^96)
println("(2) n_m=100 のとき、測定は prior の質をどれだけ見分けているか")
println("    感度 = |dlogG/dlogΔ| = 2(3^|A|-1)Δ²/[(3^|A|-1)Δ²+v_s]  (0 なら prior に無反応)")
@printf("%6s %10s %22s %20s\n", "χ_p", "点数", "G/G_max>0.9 の割合", "感度の中央値")
for x in (2,4,8,16,32)
    idx = findall(i -> haskey(ABSA, ob[i]) && cp[i]==x && Δ[i]>0 && Gmx[i]>1.05, 1:length(L))
    isempty(idx) && continue
    s = Float64[]; f = Float64[]
    for i in idx
        a = ABSA[ob[i]]; vs = 3.0^a*(1-tru[i]^2)/100
        push!(s, 2*(3.0^a-1)*Δ[i]^2 / ((3.0^a-1)*Δ[i]^2 + vs))
        push!(f, Gem[i]/Gmx[i] > 0.9 ? 1.0 : 0.0)
    end
    @printf("%6d %10d %21.0f%% %20.3f\n", x, length(idx), 100*mean(f), median(s))
end

println("\n", "="^96)
println("(3) 総コストを固定したときの最適 n_m")
println("    基底切替に「ショット c 個ぶん」の手間がかかるとすると、総コスト T = n_u(c+n_m) のもとで")
println("    Var を最小にするのは n_m_opt = sqrt(c * n_m*)  (c=0 なら n_m=1、すなわち標準的なシャドウ)")
@printf("%-16s %6s %14s %12s %12s %12s\n",
        "観測量", "χ_p", "n_m*", "n_m_opt(c=1)", "n_m_opt(c=100)", "n_m=100 の損")
for o in PURE, x in (2, 4, 8, 32)
    idx = findall(i -> ob[i]==o && cp[i]==x && L[i]==32 && dop[i]==0.0 && Δ[i]>0, 1:length(L))
    isempty(idx) && continue
    a = ABSA[o]; P = median(tru[idx]); d = median(Δ[idx])
    ns = nmstar(a, P, d)
    A = (3.0^a-1)*d^2; B = 3.0^a*(1-P^2)
    f(nm, cc) = (cc+nm)*(A + B/nm)
    opt1 = sqrt(1*ns); opt100 = sqrt(100*ns)
    @printf("%-16s %6d %14.3g %12.1f %12.1f %12.1f倍\n",
            o, x, ns, opt1, opt100, f(100,1)/f(opt1,1))
end

println("\n", "="^96)
println("(4) n_u は G の値には影響しない。効くのは「測定された G のばらつき」だけ")
println("    相対誤差 ≈ sqrt((2 + 3^|A|/n_u)/n_repeat)")
println("    総設定数 N = n_repeat*n_u を固定すると sqrt((2n_u + 3^|A|)/N) なので **n_u は小さいほど良い**")
println("    ただし n_u ≳ 3^|A| でないと基底が一致する設定がほとんど出ない")
@printf("%-16s %5s %10s %14s %16s %16s\n",
        "観測量", "|A|", "3^|A|", "一致設定数", "CI相対幅(実測)", "予測値")
for o in sort(collect(keys(ABSA)))
    idx = findall(i -> ob[i]==o && Gem[i]>1.2, 1:length(L))
    isempty(idx) && continue
    a = ABSA[o]
    w = median([(Ghi[i]-Glo[i])/2/Gem[i] for i in idx])
    @printf("%-16s %5d %10d %14.2f %16.3f %16.3f\n",
            o, a, 3^a, 50/3^a, w, sqrt((2 + 3.0^a/50)/1000))
end

println("\n", "="^96)
println("(5) 同じ総ショット数での CRM の実力 — n_m を最適に割ったときの上限")
println("    n_m=1 では G ≤ 1/(1-<P>^2) で、prior の質によらず観測量だけで決まる")
@printf("%-16s %10s %12s %14s %14s %14s\n",
        "観測量", "<P>", "1/(1-<P>²)", "G(n_m=1)", "G(n_m=100)", "G(n_m=1000)")
for o in PURE
    idx = findall(i -> ob[i]==o && cp[i]==8 && L[i]==32 && dop[i]==0.0 && Δ[i]>0, 1:length(L))
    isempty(idx) && continue
    a = ABSA[o]; P = median(tru[idx]); d = median(Δ[idx])
    @printf("%-16s %10.4f %12.2f %14.2f %14.1f %14.1f\n",
            o, P, 1/(1-P^2), gain(a,P,d,1), gain(a,P,d,100), gain(a,P,d,1000))
end
println("\n→ n_m を大きくすると G は伸びるが、それは標準シャドウ側が不利な配分を強いられているため。")
println("  総ショット数を固定して両者を最適配分すると差は 1/(1-<P>²) 程度に縮む。")
println("  逆に基底切替コスト c>0 があれば n_m_opt=sqrt(c·n_m*) は 1 より十分大きくなり、")
println("  良い prior(n_m* が大きい)ほど大きな n_m が正当化される。")
