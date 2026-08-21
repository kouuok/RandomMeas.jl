# ============================================================
# §2 の局所性検証に HF(平均場)prior を加える
#
# 動機:
#   §2 は prior として「ρ 自身の χ 切断」だけを使った。これは ρ から作るので
#   構成上 ρ と整合しており、いわば有利な prior である。実務で使いたいのは
#   O(N^3) で作れる HF のような**独立に構成された** prior だが、そちらは
#   質的に違う間違い方をしうる。
#
#   特に1次元では HF は質的に誤る: 半充填のハバード鎖の真の基底状態は
#   長距離磁気秩序を持たない(スピン相関は代数的減衰)のに、UHF は π での
#   nesting により任意の U>0 で**Néel秩序を持つ状態**に落ちる。したがって
#     - <S^z_i> は真値 0 に対し UHF は ≠ 0  → Δ≠0 かつ <P>=0 ⇒ G<1(損)
#     - サイト間スピン相関は凍結秩序で過大評価
#   が予想される。局所性の主張がこの「質的に誤った prior」でも成り立つかを見る。
#
# 比較する prior:
#   UHF      : 共線(z軸)自己無撞着平均場。Slater行列式 ⇒ Wick で閉形式
#   UHF-sym  : スピン回転平均した混合状態(対称性回復)
#   χ切断MPS : §2 と同じ(比較対象)
#
# 実行: julia --project=Hubbard_MPS_Env_v2 crm_1d_hf.jl
# 環境変数: L_LIST(既定 "8,16,32"), CHI_EXP(既定 256), N_REPEAT(既定 500),
#           DOPE_LIST(既定 "0.0"), ALLSITES(既定 0; 1なら全サイト)
# ============================================================
include(joinpath(@__DIR__, "crm_chain_common.jl"))

# ---------------- 1D鎖の UHF ----------------
chain_hop_matrix(L, t) = (T = zeros(L, L);
    for i in 1:L-1; T[i,i+1] = T[i+1,i] = -t; end; T)

function solve_uhf_chain(L, t, U; Nup, Ndn, iters=4000, tol=1e-12, init=:neel)
    T = chain_hop_matrix(L, t)
    fill_avg = (Nup + Ndn) / (2L)
    nup = zeros(L); ndn = zeros(L)
    for i in 1:L
        if init == :neel
            nup[i] = fill_avg + 0.4*(-1)^i;  ndn[i] = fill_avg - 0.4*(-1)^i
        else
            nup[i] = fill_avg + 0.3*(rand()-0.5); ndn[i] = fill_avg + 0.3*(rand()-0.5)
        end
    end
    clamp!(nup, 0.02, 0.98); clamp!(ndn, 0.02, 0.98)
    nup .*= Nup/sum(nup); ndn .*= Ndn/sum(ndn)
    local Fu, Fd
    for _ in 1:iters
        Fu = eigen(Symmetric(T + diagm(U .* ndn)))
        Fd = eigen(Symmetric(T + diagm(U .* nup)))
        nu2 = vec(sum(abs2, Fu.vectors[:, 1:Nup]; dims=2))
        nd2 = vec(sum(abs2, Fd.vectors[:, 1:Ndn]; dims=2))
        if max(maximum(abs.(nu2 .- nup)), maximum(abs.(nd2 .- ndn))) < tol
            nup, ndn = nu2, nd2; break
        end
        nup = 0.5.*nu2 .+ 0.5.*nup; ndn = 0.5.*nd2 .+ 0.5.*ndn
    end
    Φu = Fu.vectors[:, 1:Nup]; Φd = Fd.vectors[:, 1:Ndn]
    Cu = Φu*Φu'; Cd = Φd*Φd'
    E = tr(T*Cu) + tr(T*Cd) + U*sum(diag(Cu).*diag(Cd))
    return (Cu=Cu, Cd=Cd, E=E, m=mean((-1)^i*(Cu[i,i]-Cd[i,i])/2 for i in 1:L))
end

function best_uhf_chain(L, t, U; Nup, Ndn)
    Random.seed!(17); best = nothing
    for ini in (:neel, :random, :random)
        u = solve_uhf_chain(L, t, U; Nup, Ndn, init=ini)
        (best === nothing || u.E < best.E) && (best = u)
    end
    return best
end

# ---------------- Wick(共線 / 一般化) ----------------
site_spin(q) = (div(q+1, 2), isodd(q) ? 1 : 2)

full_C(uhf, L) = (C = zeros(2L, 2L);
    for a in 1:L, b in 1:L
        C[qup(a), qup(b)] = uhf.Cu[a,b]; C[qdn(a), qdn(b)] = uhf.Cd[a,b]
    end; C)

function rotate_C(C0, L, θ)
    c, s = cos(θ/2), sin(θ/2); Um = zeros(2L, 2L)
    for i in 1:L
        a, b = qup(i), qdn(i)
        Um[a,a] = c; Um[a,b] = -s; Um[b,a] = s; Um[b,b] = c
    end
    return Um * C0 * Um'
end

"""一般化相関行列 C に対する Pauli 項の Wick 期待値(Z列と (XZX+YZY)/2 のみ)。"""
function wick_gen(sup, C)
    isempty(sup) && return 1.0
    ps = [a for (_, a) in sup]
    if all(==(3), ps)
        qs = [q for (q, _) in sup]
        length(qs) == 1 && return 1 - 2C[qs[1], qs[1]]
        a, b = qs
        return 1 - 2C[a,a] - 2C[b,b] + 4*(C[a,a]*C[b,b] - C[a,b]*C[b,a])
    elseif length(sup) == 3 && (ps == [1,3,1] || ps == [2,3,2])
        return 2 * C[sup[1][1], sup[3][1]]
    end
    error("wick_gen: 未対応の項 $(sup)")
end

symavg(obs, C0, L; nθ=64) = begin
    P = [[0.0 for _ in o.terms] for o in obs]
    for k in 1:nθ
        θ = acos(-1 + (k-0.5)*2/nθ); C = rotate_C(C0, L, θ)
        for (i,o) in enumerate(obs), (j,tm) in enumerate(o.terms)
            P[i][j] += wick_gen(tm.sup, C)/nθ
        end
    end
    P
end

# ---------------- メイン ----------------
function main()
    t, U, nu, nm = 1.0, 4.0, 50, 100
    L_list = parse.(Int, split(get(ENV, "L_LIST", "8,16,32"), ","))
    dope_list = parse.(Float64, split(get(ENV, "DOPE_LIST", "0.0"), ","))
    chi_exp = parse(Int, get(ENV, "CHI_EXP", "256"))
    n_repeat = parse(Int, get(ENV, "N_REPEAT", "500"))
    allsites = get(ENV, "ALLSITES", "0") == "1"
    chi_priors = [2, 4, 8, 16, 32]

    println("=== validation (L=4) ==="); flush(stdout)
    dense_chain_check()

    rows = []
    for L in L_list, δ in dope_list
        E, ψ, _, _, nel = ground_state(L, t, U, U/2; chi_max=chi_exp, nsweeps=20, δ)
        Nup = (nel+1)÷2; Ndn = nel÷2
        uhf = best_uhf_chain(L, t, U; Nup, Ndn)
        C0 = full_C(uhf, L)
        @printf("\n=== L=%d δ=%.3f : E0(DMRG)=%.6f  E(UHF)=%.6f  UHFのスタッガード磁化 m=%.4f ===\n",
                L, δ, E, uhf.E, uhf.m); flush(stdout)

        sites = allsites ? (1:L-1) : (L÷2:L÷2)
        obs_by_site = [site_observables(i; bond=true) for i in sites]
        windows = [obs_support(o) for o in obs_by_site]
        rdm_ρ = window_rdms(ψ, windows)
        priors = MPS[]
        for cp in chi_priors
            σ = truncate(ψ; maxdim=cp); normalize!(σ); push!(priors, σ)
        end
        push!(priors, ψ)
        fids = [abs2(inner(σ, ψ)) for σ in priors]
        rdm_p = [window_rdms(σ, windows) for σ in priors]
        labels = vcat(["chi$c" for c in chi_priors], ["exact", "UHF", "UHFsym"])

        for (si, i) in enumerate(sites)
            obs = obs_by_site[si]; q1, q2 = windows[si]; Nw = q2-q1+1
            ow = shift_obs(obs, q1-1)
            Mats = [[term_window_matrix(tm, Nw) for tm in o.terms] for o in ow]
            Otrue = [sum(tm.coeff*expect_rdm(rdm_ρ[si], Mats[k][ti])
                         for (ti,tm) in enumerate(ow[k].terms)) for k in 1:length(obs)]
            Pσ_all = Vector{Vector{Vector{Float64}}}(); trOσ_all = Vector{Vector{Float64}}()
            for p in 1:length(priors)
                Pσ = [[expect_rdm(rdm_p[p][si], M) for M in Ms] for Ms in Mats]
                push!(Pσ_all, Pσ)
                push!(trOσ_all, [sum(tm.coeff*Pσ[k][ti] for (ti,tm) in enumerate(obs[k].terms))
                                 for k in 1:length(obs)])
            end
            # UHF(共線)と UHF-sym(回転平均)
            Pu = [[wick_gen(tm.sup, C0) for tm in o.terms] for o in obs]
            Ps = symavg(obs, C0, L)
            for P in (Pu, Ps)
                push!(Pσ_all, P)
                push!(trOσ_all, [sum(tm.coeff*P[k][ti] for (ti,tm) in enumerate(obs[k].terms))
                                 for k in 1:length(obs)])
            end

            es, ec = run_locals(rdm_ρ[si], Nw, ow, Pσ_all, trOσ_all;
                                nu, nm, n_repeat, seed = 8100 + 10L + i)
            for (k,o) in enumerate(obs)
                v = var(es[:,k]); gmx = v/var(ec[:,k,length(chi_priors)+1])
                for (p, lb) in enumerate(labels)
                    push!(rows, (L, δ, i, lb, o.name, Otrue[k], Otrue[k]-trOσ_all[p][k],
                                 p<=length(chi_priors) ? fids[p] : NaN,
                                 v/var(ec[:,k,p]), gmx))
                end
            end
        end

        # 画面表示: 中央サイトの比較表
        ic = allsites ? L÷2 : L÷2
        @printf("  %-14s %9s | %8s %8s %8s | %8s %8s | %8s\n",
                "observable", "<P>", "χ=2", "χ=8", "χ=32", "UHF", "UHFsym", "天井")
        for o in unique(r[5] for r in rows if r[1]==L && r[3]==ic)
            get1(lb) = (h=[r for r in rows if r[1]==L && r[3]==ic && r[4]==lb && r[5]==o];
                        isempty(h) ? NaN : h[1][9])
            gm = [r[10] for r in rows if r[1]==L && r[3]==ic && r[5]==o][1]
            tv = [r[6] for r in rows if r[1]==L && r[3]==ic && r[5]==o][1]
            @printf("  %-14s %9.4f | %8.2f %8.2f %8.2f | %8.2f %8.2f | %8.2f\n",
                    o, tv, get1("chi2"), get1("chi8"), get1("chi32"),
                    get1("UHF"), get1("UHFsym"), gm)
        end
        flush(stdout)
    end

    out = joinpath(@__DIR__, "crm_1d_hf_results.tsv")
    open(out, "w") do io
        println(io, "L\tdoping\tsite\tprior\tobservable\ttrue\tDelta\tprior_fid\tG\tG_max")
        for r in rows; println(io, join(r, "\t")); end
    end
    println("\nresults saved: $out")
end
main()
