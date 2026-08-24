# ============================================================
# 総ショット数を固定して、設定数 n_u とショット数 n_m の配分を掃引する
#
# 動機:
#   既存の掃引(crm_param_sweep.jl 実験A)は n_u=50 を固定して n_m を振ったので、
#   総ショット数 B = n_u*n_m 自体が 100 倍変わっていた。それでは「配分をどうすべきか」
#   という設計問題に答えられない。ここでは B を固定して配分だけを変える。
#
# 検証したい予言(README §0.5b):
#   Var = (1/B)[(3^|A|-1)Δ² n_m + 3^|A|(1-<P>²)]      → n_m に単調増加、n_m=1 が最良
#   G/G_max = 1/(1 + n_m/n_m*)                          → 天井までの距離が配分の無駄と一致
#   n_m* = 3^|A|(1-<P>²)/[(3^|A|-1)Δ²]
#   基底切替コスト c(ショット何個ぶんか)を入れると n_m_opt = sqrt(c·n_m*)
#
#   特に n_m=1(標準的な古典シャドウの設定)は本研究で一度も試していなかった。
#
# 実行: julia --project=Hubbard_MPS_Env_v2 crm_shot_split.jl
# 環境変数: L(既定 32), U(既定 4.0), CHI_EXP(既定 256), N_REPEAT(既定 400)
# ============================================================
include(joinpath(@__DIR__, "crm_chain_common.jl"))

"""台が単一Pauli列の観測量について、理論の分散(標準・CRM)を返す。"""
function theory_pair(absA, P, Δ, nu, nm)
    vs = 3.0^absA * (1 - P^2) / nm
    ((3.0^absA - 1)*P^2 + vs)/nu, ((3.0^absA - 1)*Δ^2 + vs)/nu
end

nm_star(absA, P, Δ) = 3.0^absA*(1-P^2) / ((3.0^absA - 1)*Δ^2)

function main()
    t  = 1.0
    U  = parse(Float64, get(ENV, "U", "4.0"))
    L  = parse(Int, get(ENV, "L", "32"))
    chi_exp  = parse(Int, get(ENV, "CHI_EXP", "256"))
    n_repeat = parse(Int, get(ENV, "N_REPEAT", "400"))
    chi_priors = [2, 4, 8, 32]

    println("=== validation ==="); flush(stdout)
    dense_chain_check()

    E, ψ, _, _, _ = ground_state(L, t, U, U/2; chi_max=chi_exp, nsweeps=20)
    @printf("\nL=%d U=%.1f  E0=%.6f  chi=%d\n", L, U, E, maxlinkdim(ψ)); flush(stdout)

    i0  = L ÷ 2
    obs = site_observables(i0; bond=true)
    q1, q2 = obs_support(obs); Nw = q2 - q1 + 1
    ow  = shift_obs(obs, q1 - 1)
    Mats = [[term_window_matrix(tm, Nw) for tm in o.terms] for o in ow]
    ρw  = window_rdms(ψ, [(q1, q2)])[1]
    Otrue = [sum(tm.coeff*expect_rdm(ρw, Mats[k][ti]) for (ti,tm) in enumerate(ow[k].terms))
             for k in 1:length(obs)]

    priors = MPS[]
    for cp in chi_priors
        σ = truncate(ψ; maxdim=cp); normalize!(σ); push!(priors, σ)
    end
    push!(priors, ψ)                       # 完全な prior = 天井 G_max の実測用
    labels = vcat(["chi$c" for c in chi_priors], ["exact"])
    Pσ_all = Vector{Vector{Vector{Float64}}}(); trOσ_all = Vector{Vector{Float64}}()
    for σ in priors
        rσ = window_rdms(σ, [(q1,q2)])[1]
        Pσ = [[expect_rdm(rσ, M) for M in Ms] for Ms in Mats]
        push!(Pσ_all, Pσ)
        push!(trOσ_all, [sum(tm.coeff*Pσ[k][ti] for (ti,tm) in enumerate(obs[k].terms))
                         for k in 1:length(obs)])
    end

    # 単一Pauli列の |A|(理論式が閉形式で使えるもの)
    absA = [o.pure ? length(o.terms[1].sup) : 0 for o in obs]

    budgets = [(5000, [1,2,5,10,25,50,100,250,500,1000], n_repeat),
               (50000, [1,10,100,1000], max(100, n_repeat ÷ 3))]

    rows = []
    for (B, nm_list, nrep) in budgets
        @printf("\n%s\n総ショット数 B=%d  (n_repeat=%d)\n%s\n", "="^92, B, nrep, "="^92)
        for nm in nm_list
            nu = B ÷ nm
            nu == 0 && continue
            t0 = time()
            es, ec = run_locals(ρw, Nw, ow, Pσ_all, trOσ_all;
                                nu, nm, n_repeat=nrep, seed = 31337 + nm)
            @printf("\n-- n_u=%-6d n_m=%-5d (%.0fs) --\n", nu, nm, time()-t0)
            @printf("  %-14s %9s %9s", "observable", "<P>", "sd(標準)")
            for lb in labels; @printf("%12s", "sd($lb)"); end
            @printf("%10s%10s%12s\n", "G(chi8)", "G_max", "n_m*(chi8)")
            for (k, o) in enumerate(obs)
                sds = std(es[:,k]); v = var(es[:,k])
                sdc = [std(ec[:,k,p]) for p in 1:length(priors)]
                g8  = v / var(ec[:,k,3]); gmx = v / var(ec[:,k,end])
                Δ8  = Otrue[k] - trOσ_all[3][k]
                ns  = (o.pure && abs(Δ8) > 0) ? nm_star(absA[k], Otrue[k], Δ8) : NaN
                @printf("  %-14s %9.4f %9.5f", o.name, Otrue[k], sds)
                for s in sdc; @printf("%12.5f", s); end
                @printf("%10.2f%10.2f%12.3g\n", g8, gmx, ns)
                for (p, lb) in enumerate(labels)
                    Δp = Otrue[k] - trOσ_all[p][k]
                    vth_s, vth_c = o.pure ? theory_pair(absA[k], Otrue[k], Δp, nu, nm) : (NaN, NaN)
                    push!(rows, (B, nu, nm, o.name, o.pure, absA[k], lb, Otrue[k], Δp,
                                 sds, sdc[p], v/var(ec[:,k,p]), gmx,
                                 sqrt(vth_s), sqrt(vth_c),
                                 o.pure && abs(Δp)>0 ? nm_star(absA[k], Otrue[k], Δp) : NaN))
                end
            end
            flush(stdout)
        end
    end

    out = joinpath(@__DIR__, "crm_shot_split_results.tsv")
    open(out, "w") do io
        println(io, "B\tn_u\tn_m\tobservable\tpure\tabsA\tprior\ttrue\tDelta\t" *
                    "sd_std\tsd_crm\tG_emp\tG_max\tsd_std_theory\tsd_crm_theory\tnm_star")
        for r in rows; println(io, join(r, "\t")); end
    end
    println("\nresults saved: $out")
end
main()
