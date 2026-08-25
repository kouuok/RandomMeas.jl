# ============================================================
# prior の対称性を系統的な設計軸にする — 並進対称性の回復
#
# 主張(一般形):
#   ρ が対称群 S を持ち、prior σ がそれを破っているとする。CRM が prior に
#   要求するのは項ごとの期待値 ⟨P_t⟩_σ だけなので、群平均した混合状態
#       σ̄ = ∫dg U_g σ U_g†      ⟹  ⟨P_t⟩_σ̄ = ∫dg ⟨U_g† P_t U_g⟩_σ
#   を prior として使える。σ̄ は S を持つので、
#     (A) S が ⟨P⟩_ρ = 0 を強制する観測量では ⟨P⟩_σ̄ = 0 となり Δ = 0、
#         したがって G = G_max(損が厳密に消える)
#     (B) 多重項をなす観測量(S^xS^x, S^yS^y, S^zS^z など)では重みが
#         等方的に配分される
#   (A) は §2f の「UHF-sym で損が 0/49 になる」、(B) は §2i の
#   「横縦比が 7.14 → 1.000 になる」として既に確認済み(どちらもスピン SU(2))。
#
# ここで確かめるのは **並進対称性** の版である:
#   周期境界の ρ は並進対称だが、χ 切断MPS prior は切断誤差が位置に依存する
#   ので並進対称でない。§2h(5) で「PBC でも G のサイト依存が 2–5倍残るのは
#   物理ではなく prior の構成に由来する」と結論した部分である。
#   並進平均 σ̄ = (1/L) Σ_r T_r σ T_r† を prior にすればこれが消えるはず。
#   消えなければ §2h(5) の結論が誤りだったことになる。
#
#   CRM 側の実装は「同じ観測量を全サイトで測った prior の値を平均する」だけで、
#   窓RDMは既に全サイト分あるので追加コストはほぼゼロである。
#
# 実行: julia --project=Hubbard_MPS_Env_v2 crm_symmetry_prior.jl
# 環境変数: L_LIST(既定 "16"), U(既定 4.0), CHI_EXP(既定 512), N_REPEAT(既定 300)
# ============================================================
include(joinpath(@__DIR__, "crm_chain_common.jl"))

function chain_mpo_pbc(sites, L, t, U, mu)
    os = OpSum()
    for i in 1:L, q in (qup, qdn)
        j = mod1(i+1, L)
        os += -t, "Cdag", q(i), "C", q(j)
        os += -t, "Cdag", q(j), "C", q(i)
    end
    for i in 1:L
        os += U, "N", qup(i), "N", qdn(i)
        os += -mu, "N", qup(i); os += -mu, "N", qdn(i)
    end
    MPO(os, sites)
end

function gs_pbc(L, t, U, mu; chi_max, nsweeps=30)
    sites = siteinds("Fermion", 2L; conserve_qns=true)
    H = chain_mpo_pbc(sites, L, t, U, mu)
    init, _ = initial_config(L, 0.0)
    ramp = vcat([50,100,200,400], fill(chi_max, max(0, nsweeps-4)))
    noise = vcat([1e-5,1e-6,1e-7,1e-8,1e-9], zeros(max(0, nsweeps-5)))[1:nsweeps]
    E, ψ = dmrg(H, productMPS(sites, init); nsweeps,
                maxdim=min.(chi_max, ramp)[1:nsweeps], cutoff=1e-12, noise, outputlevel=0)
    E, ψ
end

function main()
    t = 1.0
    U  = parse(Float64, get(ENV, "U", "4.0"))
    chi_exp  = parse(Int, get(ENV, "CHI_EXP", "512"))
    n_repeat = parse(Int, get(ENV, "N_REPEAT", "300"))
    nu, nm = 50, 100
    chi_priors = [2, 4, 8]
    L_list = parse.(Int, split(get(ENV, "L_LIST", "16"), ","))

    println("=== validation ==="); flush(stdout); dense_chain_check()

    rows = []
    for L in L_list
        E, ψ = gs_pbc(L, t, U, U/2; chi_max=chi_exp)
        @printf("\n%s\nL=%d PBC  U=%.1f  E0=%.6f  chi=%d\n", "="^92, L, U, E, maxlinkdim(ψ))
        flush(stdout)

        sites = 1:L-1
        obs_by_site = [site_observables(i; bond=true) for i in sites]
        windows = [obs_support(o) for o in obs_by_site]
        rdm_ρ = window_rdms(ψ, windows)
        nobs = length(obs_by_site[1])

        priors = MPS[]
        for cp in chi_priors
            σ = truncate(ψ; maxdim=cp); normalize!(σ); push!(priors, σ)
        end
        push!(priors, ψ)
        rdm_p = [window_rdms(σ, windows) for σ in priors]

        # 各サイトの Mats を作り、prior の項ごとの期待値を全サイトぶん集める
        Mats_all = Vector{Vector{Vector{Matrix{ComplexF64}}}}(undef, length(sites))
        ow_all   = Vector{Vector{Obs}}(undef, length(sites))
        Pρ_all   = Vector{Vector{Vector{Float64}}}(undef, length(sites))
        Pσ_site  = [Vector{Vector{Vector{Float64}}}(undef, length(sites)) for _ in priors]
        for (si, i) in enumerate(sites)
            obs = obs_by_site[si]; q1,q2 = windows[si]; Nw = q2-q1+1
            ow = shift_obs(obs, q1-1); ow_all[si] = ow
            Mats = [[term_window_matrix(tm, Nw) for tm in o.terms] for o in ow]
            Mats_all[si] = Mats
            Pρ_all[si] = [[expect_rdm(rdm_ρ[si], M) for M in Ms] for Ms in Mats]
            for p in 1:length(priors)
                Pσ_site[p][si] = [[expect_rdm(rdm_p[p][si], M) for M in Ms] for Ms in Mats]
            end
        end

        # 並進平均: 同じ観測量・同じ項を全サイトで平均する
        Pσ_avg = [[[mean(Pσ_site[p][si][k][j] for si in 1:length(sites))
                    for j in 1:length(obs_by_site[1][k].terms)] for k in 1:nobs]
                  for p in 1:length(priors)]

        labels = vcat(["chi$c" for c in chi_priors], ["exact"],
                      ["chi$(c)_tavg" for c in chi_priors])
        @printf("\n  %-14s %10s | %s\n", "観測量", "<P>_ρ",
                "Δ のサイト依存 (max/min):  素の prior → 並進平均")
        for k in 1:nobs
            o = obs_by_site[1][k]
            Pr = [sum(tm.coeff*Pρ_all[si][k][ti] for (ti,tm) in enumerate(o.terms))
                  for si in 1:length(sites)]
            cp_i = 2   # chi_p = 4 を代表に
            Δ_plain = [abs(Pr[si] - sum(tm.coeff*Pσ_site[cp_i][si][k][ti]
                                        for (ti,tm) in enumerate(o.terms)))
                       for si in 1:length(sites)]
            Δ_avg = abs(mean(Pr) - sum(tm.coeff*Pσ_avg[cp_i][k][ti]
                                       for (ti,tm) in enumerate(o.terms)))
            r_plain = maximum(Δ_plain)/max(minimum(Δ_plain), 1e-300)
            @printf("  %-14s %10.5f | %.2e–%.2e (%.1f倍) → %.2e (一定)\n",
                    o.name, mean(Pr), minimum(Δ_plain), maximum(Δ_plain), r_plain, Δ_avg)
        end
        flush(stdout)

        # 利得の測定
        for (si, i) in enumerate(sites)
            obs = obs_by_site[si]; q1,q2 = windows[si]; Nw = q2-q1+1
            ow = ow_all[si]
            Pσ_all = Vector{Vector{Vector{Float64}}}(); trOσ_all = Vector{Vector{Float64}}()
            for p in 1:length(priors)
                push!(Pσ_all, Pσ_site[p][si])
                push!(trOσ_all, [sum(tm.coeff*Pσ_site[p][si][k][ti]
                                     for (ti,tm) in enumerate(obs[k].terms)) for k in 1:nobs])
            end
            for p in 1:length(chi_priors)          # 並進平均版
                push!(Pσ_all, Pσ_avg[p])
                push!(trOσ_all, [sum(tm.coeff*Pσ_avg[p][k][ti]
                                     for (ti,tm) in enumerate(obs[k].terms)) for k in 1:nobs])
            end
            es, ec = run_locals(rdm_ρ[si], Nw, ow, Pσ_all, trOσ_all;
                                nu, nm, n_repeat, seed = 660_000 + 100L + i)
            for k in 1:nobs
                v = var(es[:,k]); gmx = v/var(ec[:,k,length(chi_priors)+1])
                Otrue = sum(tm.coeff*Pρ_all[si][k][ti] for (ti,tm) in enumerate(obs[k].terms))
                for (p,lb) in enumerate(labels)
                    push!(rows, (L, i, obs[k].name, lb, Otrue,
                                 Otrue - trOσ_all[p][k], v/var(ec[:,k,p]), gmx))
                end
            end
        end

        @printf("\n  %-14s", "観測量")
        for lb in labels; @printf("%14s", lb); end; println()
        for o in unique(r[3] for r in rows if r[1]==L)
            @printf("  %-14s", o)
            for lb in labels
                h = [r[7] for r in rows if r[1]==L && r[3]==o && r[4]==lb]
                @printf("%14.2f", isempty(h) ? NaN : median(h))
            end
            println()
        end
        @printf("  %-14s", "最悪サイトG")
        for lb in labels
            h = [r[7] for r in rows if r[1]==L && r[4]==lb]
            @printf("%14.2f", isempty(h) ? NaN : minimum(h))
        end
        println(); flush(stdout)
    end

    out = joinpath(@__DIR__, "crm_symmetry_prior_results.tsv")
    open(out,"w") do io
        println(io, "L\tsite\tobservable\tprior\ttrue\tDelta\tG\tG_max")
        for r in rows; println(io, join(r,"\t")); end
    end
    println("\nresults saved: $out")
end
main()
