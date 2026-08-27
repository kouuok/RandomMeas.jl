# ============================================================
# 設定数 n_u を振ると利得は動くか — §0.5b の結論の実測検証
#
# §0.5b では「n_u は利得 G に現れない(比で約分される)ので n_u=50 で足りる」
# と結論した。根拠は
#   (i) 分散分解 Var = (1/n_u)[...] から n_u が比で消えること
#   (ii) 効くのは 3^N ではなく 3^|A| であること(台の外の基底は周辺化される)
#   (iii) 測定された G の CI 幅が sqrt((2 + 3^|A|/n_u)/n_repeat) と合うこと
# だが、**n_u を実際に振って G が動かないことは確かめていなかった**。
# 動くなら (i)(ii) のどこかが誤っている。直接確かめる。
#
# あわせて |A|=4(横スピン相関)を入れる。n_u=50 では基底一致の期待数が
# 50/81 = 0.62 で、多くの実験で一度も一致しない領域である。
#
# 実行: julia --project=Hubbard_MPS_Env_v2 crm_nu_scan.jl
# 環境変数: L(既定 16), U(既定 4.0), NU_LIST(既定 "50,100,200,400,800,1600,3200"),
#           N_REPEAT(既定 200), NM(既定 100)
# ============================================================
include(joinpath(@__DIR__, "crm_chain_common.jl"))
include(joinpath(@__DIR__, "crm_wick_pauli.jl"))

"""|A| が 1〜4 にわたる観測量を1つの窓に揃える。"""
function scan_observables(i::Int)
    u,  d  = qup(i),   qdn(i)
    u2, d2 = qup(i+1), qdn(i+1)
    [Obs("Sz  |A|=1", [Term(-0.25,[(u,3)]), Term(0.25,[(d,3)])], false),
     Obs("ZZ onsite |A|=2", [Term(1.0,[(u,3),(d,3)])], true),
     Obs("ZZ up-up |A|=2", [Term(1.0,[(u,3),(u2,3)])], true),
     Obs("DoubleOcc |A|=2", [Term(0.25,Tuple{Int,Int}[]), Term(-0.25,[(u,3)]),
                             Term(-0.25,[(d,3)]), Term(0.25,[(u,3),(d,3)])], false),
     Obs("SzSz |A|=2", [Term(1/16,[(u,3),(u2,3)]), Term(-1/16,[(u,3),(d2,3)]),
                        Term(-1/16,[(d,3),(u2,3)]), Term(1/16,[(d,3),(d2,3)])], false),
     Obs("hop |A|=3", [Term(0.5,[(u,1),(d,3),(u2,1)]), Term(0.5,[(u,2),(d,3),(u2,2)])], false),
     Obs("ZZZZ |A|=4", [Term(1.0,[(u,3),(d,3),(u2,3),(d2,3)])], true),
     Obs("SxSx |A|=4", [Term(1/16,[(u,1),(d,1),(u2,1),(d2,1)]),
                        Term(1/16,[(u,1),(d,1),(u2,2),(d2,2)]),
                        Term(1/16,[(u,2),(d,2),(u2,1),(d2,1)]),
                        Term(1/16,[(u,2),(d,2),(u2,2),(d2,2)])], false)]
end

absA_of(o) = o.pure ? length(o.terms[1].sup) :
             maximum(length(tm.sup) for tm in o.terms)

function main()
    t = 1.0
    U  = parse(Float64, get(ENV, "U", "4.0"))
    L  = parse(Int, get(ENV, "L", "16"))
    nm = parse(Int, get(ENV, "NM", "100"))
    n_repeat = parse(Int, get(ENV, "N_REPEAT", "200"))
    nu_list = parse.(Int, split(get(ENV, "NU_LIST", "50,100,200,400,800,1600,3200"), ","))
    chi_priors = [2, 4, 8, 32]

    println("=== validation ==="); flush(stdout); dense_chain_check()

    E, ψ, _, _, nel = ground_state(L, t, U, U/2; chi_max=256, nsweeps=20)
    @printf("\nL=%d U=%.1f E0=%.6f chi=%d\n", L, U, E, maxlinkdim(ψ)); flush(stdout)

    i0 = L ÷ 2
    obs = scan_observables(i0)
    q1, q2 = obs_support(obs); Nw = q2-q1+1
    ow = shift_obs(obs, q1-1)
    Mats = [[term_window_matrix(tm, Nw) for tm in o.terms] for o in ow]
    ρw = window_rdms(ψ, [(q1,q2)])[1]
    Otrue = [sum(tm.coeff*expect_rdm(ρw, Mats[k][ti]) for (ti,tm) in enumerate(ow[k].terms))
             for k in 1:length(obs)]

    priors = MPS[]
    for cp in chi_priors
        σ = truncate(ψ; maxdim=cp); normalize!(σ); push!(priors, σ)
    end
    push!(priors, ψ)
    Pσ_all = Vector{Vector{Vector{Float64}}}(); trOσ_all = Vector{Vector{Float64}}()
    for σ in priors
        rσ = window_rdms(σ, [(q1,q2)])[1]
        Pσ = [[expect_rdm(rσ, M) for M in Ms] for Ms in Mats]
        push!(Pσ_all, Pσ)
        push!(trOσ_all, [sum(tm.coeff*Pσ[k][ti] for (ti,tm) in enumerate(obs[k].terms))
                         for k in 1:length(obs)])
    end
    labels = vcat(["chi$c" for c in chi_priors], ["exact"])

    @printf("\n%-18s %5s %10s %12s %14s\n", "観測量", "|A|", "3^|A|", "<P>", "Δ(χ_p=8)")
    for (k,o) in enumerate(obs)
        a = absA_of(o)
        @printf("%-18s %5d %10d %12.5f %14.3e\n", o.name, a, 3^a, Otrue[k],
                Otrue[k]-trOσ_all[3][k])
    end
    flush(stdout)

    rows = []
    for nu in nu_list
        t0 = time()
        es, ec = run_locals(ρw, Nw, ow, Pσ_all, trOσ_all;
                            nu, nm, n_repeat, seed = 424242 + nu)
        @printf("\n-- n_u=%-5d (n_repeat=%d, %.0fs) --\n", nu, n_repeat, time()-t0)
        @printf("  %-18s %8s", "観測量", "一致数")
        for lb in labels; @printf("%12s", lb); end
        @printf("%12s%12s\n", "G理論(χ8)", "CI相対幅")
        for (k,o) in enumerate(obs)
            a = absA_of(o); v = var(es[:,k])
            gs = [v/var(ec[:,k,p]) for p in 1:length(priors)]
            lo, hi = boot_ratio(es[:,k], ec[:,k,3])
            ciw = (hi-lo)/2/max(gs[3], 1e-12)
            Δ8 = Otrue[k] - trOσ_all[3][k]
            gth = o.pure ? ((3.0^a-1)*Otrue[k]^2 + 3.0^a*(1-Otrue[k]^2)/nm) /
                           ((3.0^a-1)*Δ8^2      + 3.0^a*(1-Otrue[k]^2)/nm) : NaN
            @printf("  %-18s %8.2f", o.name, nu/3.0^a)
            for g in gs; @printf("%12.3f", g); end
            @printf("%12s%12.4f\n", isnan(gth) ? "—" : @sprintf("%.3f", gth), ciw)
            for (p,lb) in enumerate(labels)
                push!(rows, (L, nu, nm, n_repeat, o.name, a, lb, Otrue[k],
                             Otrue[k]-trOσ_all[p][k], gs[p], gth, ciw))
            end
        end
        flush(stdout)
    end

    out = joinpath(@__DIR__, "crm_nu_scan_results.tsv")
    open(out,"w") do io
        println(io, "L\tn_u\tn_m\tn_repeat\tobservable\tabsA\tprior\ttrue\tDelta\tG\tG_theory\tCI_relwidth")
        for r in rows; println(io, join(r,"\t")); end
    end
    println("\nresults saved: $out")
end
main()
