# ============================================================
# PBC で主張は成り立つか — 「χ_p≥4 なら損をしない」の境界条件依存性
#
# 動機:
#   §2b の主張は OBC で検証した。PBC はエンタングルメントが約2倍
#   (S=(c/3)lnL vs (c/6)lnL)なので、同じ χ_p の prior は PBC の方が
#   悪いはずである。したがって「χ_p≥4 なら損なサイトはゼロ」が PBC でも
#   成り立つかは自明でない。これを直接確かめる。
#
#   PBC では全サイトが等価なので、サイト分解は不要（中央1サイトで十分）。
#   代わりに OBC/PBC を並べて、同じ χ_p での忠実度・Δ・G を比べる。
#
# 実行: julia --project=Hubbard_MPS_Env_v2 crm_pbc_gain.jl
# 環境変数: L_LIST(既定 "16,32"), CHI_EXP(既定 512), N_REPEAT(既定 500)
# ============================================================
include(joinpath(@__DIR__, "crm_chain_common.jl"))

function chain_mpo_bc(sites, L, t, U, mu; pbc::Bool)
    os = OpSum()
    bonds = pbc ? [(i, mod1(i+1, L)) for i in 1:L] : [(i, i+1) for i in 1:L-1]
    for (i, j) in bonds, q in (qup, qdn)
        os += -t, "Cdag", q(i), "C", q(j)
        os += -t, "Cdag", q(j), "C", q(i)
    end
    for i in 1:L
        os += U, "N", qup(i), "N", qdn(i)
        os += -mu, "N", qup(i); os += -mu, "N", qdn(i)
    end
    return MPO(os, sites)
end

function gs_bc(L, t, U, mu; chi_max, nsweeps, pbc)
    sites = siteinds("Fermion", 2L; conserve_qns=true)
    H = chain_mpo_bc(sites, L, t, U, mu; pbc)
    init, nel = initial_config(L, 0.0)
    ramp = vcat([50,100,200,400], fill(chi_max, max(0, nsweeps-4)))
    noise = vcat([1e-5,1e-6,1e-7,1e-8,1e-9], zeros(max(0, nsweeps-5)))[1:nsweeps]
    E, ψ = dmrg(H, productMPS(sites, init); nsweeps,
                maxdim=min.(chi_max, ramp)[1:nsweeps], cutoff=1e-11, noise, outputlevel=0)
    return E, ψ, nel
end

function center_S(ψ)
    N=length(ψ); ψo=copy(ψ); orthogonalize!(ψo, N÷2)
    _,Sv,_ = svd(ψo[N÷2], (linkinds(ψo)[N÷2-1], siteinds(ψo)[N÷2]))
    p=diag(Array(Sv,inds(Sv)...)).^2; p=p[p.>1e-14]
    return -sum(p.*log.(p))
end

function main()
    t, U, nu, nm = 1.0, 4.0, 50, 100
    L_list = parse.(Int, split(get(ENV, "L_LIST", "16,32"), ","))
    chi_exp = parse(Int, get(ENV, "CHI_EXP", "512"))
    n_repeat = parse(Int, get(ENV, "N_REPEAT", "500"))
    chi_priors = [2, 4, 8, 16, 32]

    println("=== validation (L=4, OBC) ==="); flush(stdout)
    dense_chain_check()

    rows = []
    for L in L_list, pbc in (false, true)
        E, ψ, nel = gs_bc(L, t, U, U/2; chi_max=chi_exp, nsweeps=30, pbc)
        Sc = center_S(ψ)
        @printf("\n=== L=%d %s : E0=%.6f chi=%d S_center=%.4f ===\n",
                L, pbc ? "PBC" : "OBC", E, maxlinkdim(ψ), Sc); flush(stdout)

        i0 = L ÷ 2                       # PBCではどこでも同じ。OBCは中央を取る
        obs = site_observables(i0; bond=true)
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
        fids = [abs2(inner(σ, ψ)) for σ in priors]
        Pσ_all=Vector{Vector{Vector{Float64}}}(); trOσ_all=Vector{Vector{Float64}}()
        for σ in priors
            rσ = window_rdms(σ, [(q1,q2)])[1]
            Pσ = [[expect_rdm(rσ, M) for M in Ms] for Ms in Mats]
            push!(Pσ_all, Pσ)
            push!(trOσ_all, [sum(tm.coeff*Pσ[k][ti] for (ti,tm) in enumerate(obs[k].terms))
                             for k in 1:length(obs)])
        end
        @printf("  prior忠実度: %s\n", join([@sprintf("χ=%d:%.3e", c, f)
                for (c,f) in zip(chi_priors, fids)], "  ")); flush(stdout)

        es, ec = run_locals(ρw, Nw, ow, Pσ_all, trOσ_all; nu, nm, n_repeat, seed=4242+L)
        @printf("  %-14s %9s", "observable", "<P>")
        for c in chi_priors; @printf("%11s", "G(χ=$c)"); end
        @printf("%11s\n", "G_max")
        for (k,o) in enumerate(obs)
            v = var(es[:,k]); gmx = v/var(ec[:,k,end])
            @printf("  %-14s %9.4f", o.name, Otrue[k])
            for p in 1:length(chi_priors)
                g = v/var(ec[:,k,p]); @printf("%11.2f", g)
                push!(rows, (L, pbc ? "PBC" : "OBC", chi_priors[p], o.name,
                             Otrue[k], Otrue[k]-trOσ_all[p][k], fids[p], g, gmx, Sc))
            end
            @printf("%11.2f\n", gmx)
        end
        flush(stdout)
    end
    out = joinpath(@__DIR__, "crm_pbc_gain_results.tsv")
    open(out,"w") do io
        println(io, "L\tbc\tchi_prior\tobservable\ttrue\tDelta\tprior_fid\tG\tG_max\tS_center")
        for r in rows; println(io, join(r,"\t")); end
    end
    println("\nresults saved: $out")
end
main()
