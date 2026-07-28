# ============================================================
# 実験B改: 「目標誤差に到達するのに必要な χ」で局所性の限界を測る
#
# 動機:
#   実験B(crm_gap_scaling.jl)は固定 χ_p での相対誤差 ε の L 依存を見たが、
#   |Δ|(χ) が偶然の打ち消しで単調にならないため、ε(L) も その比も χ_p に
#   対して頑健でなかった(README §2e)。
#
# 修正:
#   固定 χ での ε ではなく、逆に
#       χ*(ε_target, L) ≡ min{ χ : ε(χ') < ε_target が全ての χ'≥χ で成立 }
#   を指標にする。「最後に閾値を上から横切った点」を取るので、非単調性に
#   影響されず単調で well-defined になる。
#
# 予言:
#   ギャップのある系  : χ*(L) は L とともに飽和
#   gapless な系      : χ*(L) は L とともに増加（面積則の対数補正を反映）
#
# Monte Carlo 不要。DMRG と窓RDMのみ。
# 実行: julia --project=Hubbard_MPS_Env_v2 crm_chistar.jl
# 環境変数: L_LIST(既定 "16,32,64,128"), CHI_EXP(既定 400), NSWEEPS(既定 30)
# ============================================================
include(joinpath(@__DIR__, "crm_chain_common.jl"))

const CHIS = [2, 4, 6, 8, 12, 16, 24, 32, 48, 64, 96, 128, 192, 256]
const TARGETS = [1e-1, 1e-2, 1e-3]

"""ε(χ) の列から χ*(target) を求める。非単調性に強い「最後の上向き交差」定義:
   χ' ≥ χ の全てで ε(χ') < target となる最小の χ。存在しなければ NaN。"""
function chi_star(chis, eps, target)
    n = length(chis)
    ok = true
    best = NaN
    for k in n:-1:1
        if isnan(eps[k]) || eps[k] >= target
            ok = false
        end
        ok && (best = chis[k])
    end
    return best
end

function main()
    t = 1.0
    L_list = parse.(Int, split(get(ENV, "L_LIST", "16,32,64,128"), ","))
    chi_exp = parse(Int, get(ENV, "CHI_EXP", "400"))
    nsweeps = parse(Int, get(ENV, "NSWEEPS", "30"))
    configs = [(1.0,0.0),(2.0,0.0),(4.0,0.0),(8.0,0.0),(4.0,0.125),(8.0,0.125)]
    obsnames = ["DoubleOcc","ZZ onsite","SzSz r=1","ZZ up-up r=1","hop up r=1"]

    println("=== validation (L=4) ==="); flush(stdout)
    dense_chain_check()

    rows = []
    for (U, δ) in configs, L in L_list
        tstart = time()
        E, ψ, _, _, nel = ground_state(L, t, U, U/2; chi_max=chi_exp, nsweeps, δ)
        @printf("\nU=%.1f δ=%.3f L=%3d : E0=%12.6f N_el=%3d chi=%3d (%.0fs)\n",
                U, δ, L, E, nel, maxlinkdim(ψ), time()-tstart); flush(stdout)

        # 全サイトの窓と真値
        obs_by_site = [site_observables(i; bond = i < L) for i in 1:L]
        windows = [obs_support(o) for o in obs_by_site]
        rρ = window_rdms(ψ, windows)
        Mats = [[[term_window_matrix(tm, windows[s][2]-windows[s][1]+1) for tm in o.terms]
                 for o in shift_obs(obs_by_site[s], windows[s][1]-1)] for s in 1:L]
        vρ = [[sum(tm.coeff * expect_rdm(rρ[s], Mats[s][k][ti])
                   for (ti,tm) in enumerate(obs_by_site[s][k].terms))
               for k in 1:length(obs_by_site[s])] for s in 1:L]

        # χ 掃引: 各 χ で全サイトの ε の中央値
        epsmed = Dict(o => Float64[] for o in obsnames)
        for cp in CHIS
            cp > maxlinkdim(ψ) && (for o in obsnames; push!(epsmed[o], NaN); end; continue)
            σ = truncate(ψ; maxdim=cp); normalize!(σ)
            rσ = window_rdms(σ, windows)
            acc = Dict(o => Float64[] for o in obsnames)
            for s in 1:L, (k, o) in enumerate(obs_by_site[s])
                o.name in obsnames || continue
                vσ = sum(tm.coeff * expect_rdm(rσ[s], Mats[s][k][ti])
                         for (ti,tm) in enumerate(shift_obs(obs_by_site[s], windows[s][1]-1)[k].terms))
                abs(vρ[s][k]) > 1e-3 && push!(acc[o.name], abs((vρ[s][k]-vσ)/vρ[s][k]))
            end
            for o in obsnames
                push!(epsmed[o], isempty(acc[o]) ? NaN : median(acc[o]))
            end
        end

        @printf("%-14s", "observable"); for c in CHIS; @printf("%9d", c); end
        for tg in TARGETS; @printf("%10s", "χ*($tg)"); end; println()
        for o in obsnames
            @printf("%-14s", o)
            for x in epsmed[o]; @printf("%9.2e", x); end
            for tg in TARGETS
                cs = chi_star(CHIS, epsmed[o], tg)
                @printf("%10s", isnan(cs) ? ">256" : string(Int(cs)))
                push!(rows, (U, δ, L, o, tg, cs))
            end
            println()
        end
        flush(stdout)
    end

    out = joinpath(@__DIR__, "crm_chistar_results.tsv")
    open(out, "w") do io
        println(io, "U\tdoping\tL\tobservable\ttarget\tchi_star")
        for r in rows; println(io, join(r, "\t")); end
    end
    println("\nresults saved: $out")

    println("\n", "="^92)
    println("χ*(ε<1e-2) の L 依存 — gapped なら飽和、gapless なら増加")
    @printf("%-16s %-14s", "config", "observable")
    for L in L_list; @printf("%9s", "L=$L"); end; println()
    for (U, δ) in configs, o in obsnames
        lab = δ == 0 ? "U=$(Int(U)) half" : "U=$(Int(U)) δ=1/8"
        @printf("%-16s %-14s", lab, o)
        for L in L_list
            i = findfirst(r -> r[1]==U && r[2]==δ && r[3]==L && r[4]==o && r[5]==1e-2, rows)
            @printf("%9s", i === nothing ? "-" : (isnan(rows[i][6]) ? ">256" : string(Int(rows[i][6]))))
        end
        println()
    end
end

main()
