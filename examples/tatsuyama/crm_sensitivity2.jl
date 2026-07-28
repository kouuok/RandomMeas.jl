# ============================================================
# 観測量ごとの「切断への感度」は定義できるか
#
# ROADMAP §3.1 で、切断が捨てる重み w のうち観測量誤差 |Δ| に効く割合が
# 観測量と U によって 1% から 41% まで変わることが分かった。
# ならば χ を掃引したとき |Δ(χ)| と w(χ) の関係は観測量ごとに一定の
# 「感度」 s = |Δ|/w で書けるのか、それとも指数が違うのか?
#
#   |Δ(χ)| ≈ s · w(χ)^α
#
# α=1 で s が観測量ごとの定数なら、「この観測量にはχいくつ必要か」を
# 捨てる重みだけから見積もれることになる(prior側だけで完結する指標)。
#
# Monte Carlo 不要。DMRG 1本 + 窓RDM のみ。
# 実行: julia --project=Hubbard_MPS_Env_v2 crm_sensitivity.jl
# ============================================================
include(joinpath(@__DIR__, "crm_chain_common.jl"))

const L = parse(Int, get(ENV, "LSYS", "64"))
const CHIS = [2, 4, 8, 16, 32, 64, 128]

"""中央ボンドで χ を残したときに捨てる Schmidt 重み。"""
function discarded_weights(ψ::MPS, chis)
    N = length(ψ); ψo = copy(ψ); orthogonalize!(ψo, N ÷ 2)
    _, Sv, _ = svd(ψo[N÷2], (linkinds(ψo)[N÷2-1], siteinds(ψo)[N÷2]))
    p = sort(diag(Array(Sv, inds(Sv)...)).^2; rev=true)
    return [1 - sum(p[1:min(c, length(p))]) for c in chis]
end

function main()
    t = 1.0
    rows = []
    for U in (2.0, 4.0, 8.0), δ in (0.0, 0.125)
        E, ψ, _, _, nel = ground_state(L, t, U, U/2; chi_max=400, nsweeps=30, δ)
        w = discarded_weights(ψ, CHIS)
        @printf("\n=== U=%.1f, δ=%.3f : E0=%.6f, N_el=%d, chi=%d ===\n",
                U, δ, E, nel, maxlinkdim(ψ)); flush(stdout)
        @printf("捨てる重み w(χ): %s\n",
                join([@sprintf("χ=%d:%.2e", c, x) for (c,x) in zip(CHIS,w)], "  "))

        i0 = L ÷ 2
        obs = site_observables(i0; bond=true)
        q1, q2 = obs_support(obs); Nw = q2 - q1 + 1
        ow = shift_obs(obs, q1 - 1)
        Mats = [[term_window_matrix(tm, Nw) for tm in o.terms] for o in ow]
        rρ = window_rdms(ψ, [(q1,q2)])[1]
        vρ = [sum(tm.coeff * expect_rdm(rρ, Mats[k][ti]) for (ti,tm) in enumerate(ow[k].terms))
              for k in 1:length(obs)]

        @printf("%-14s %10s", "observable", "<P>")
        for c in CHIS; @printf("%11s", "χ=$c"); end
        @printf("%12s %10s\n", "感度 s", "指数 α")
        for (k, o) in enumerate(obs)
            ds = Float64[]
            for c in CHIS
                σ = truncate(ψ; maxdim=c); normalize!(σ)
                rσ = window_rdms(σ, [(q1,q2)])[1]
                vσ = sum(tm.coeff * expect_rdm(rσ, Mats[k][ti]) for (ti,tm) in enumerate(ow[k].terms))
                push!(ds, abs(vρ[k] - vσ))
            end
            # log|Δ| = log s + α log w の最小二乗
            ok = findall(i -> ds[i] > 1e-13 && w[i] > 1e-13, 1:length(CHIS))
            α = NaN; s = NaN
            if length(ok) >= 3
                X = log.(w[ok]); Y = log.(ds[ok])
                α = (mean(X.*Y) - mean(X)*mean(Y)) / (mean(X.^2) - mean(X)^2)
                s = exp(mean(Y) - α*mean(X))
            end
            @printf("%-14s %10.4f", o.name, vρ[k])
            for x in ds; @printf("%11.2e", x); end
            @printf("%12.3f %10.3f\n", s, α)
            for (c, x, ww) in zip(CHIS, ds, w)
                push!(rows, (U, δ, o.name, c, vρ[k], x, ww, s, α))
            end
        end
        flush(stdout)
    end
    out = joinpath(@__DIR__, "crm_sensitivity_results.tsv")
    open(out, "w") do io
        println(io, "U\tdoping\tobservable\tchi_prior\ttrue\tabsDelta\tdiscarded_w\tsens_s\talpha")
        for r in rows; println(io, join(r, "\t")); end
    end
    println("\nresults saved: $out")
end

main()
