# ============================================================
# 横スピン相関 <S^x_i S^x_j> — UHF の対称性の破れを直接見る
#
# 動機:
#   これまで測ってきたのは密度・縦スピン(Z型)と運動項だけだった。
#   半充填ハバード鎖の基底状態は SU(2) 一重項なので
#       <S^x_i S^x_j> = <S^y_i S^y_j> = <S^z_i S^z_j>
#   が厳密に成り立つ。ところが UHF は z 軸に共線な Néel 解なので、
#   縦は大きく横はほぼゼロを返す。**同じ prior が同じ物理量に対して
#   成分によって桁違いに違う誤差を持つ**わけで、Z 型しか測っていない
#   限りこの破れは一度も観測にかからない。
#   対称性回復(UHF-sym)は回転平均なので等方性を回復するはずで、
#   §2f/§4 の主張に対する定量的で反証可能な検証になる。
#
# JW での表現(a = qup(i), b = qdn(i) は隣接なので弦が付かない):
#       S^x_i = (X_a X_b + Y_a Y_b)/4      → |A| = 2
#       S^x_i S^x_j                        → |A| = 4 (4項)
#   |A|=4 なので基底一致率は 3^-4 = 1/81。§0.5b の指針に従い n_u = 240
#   (一致設定の期待数 ≈ 3)を使う。n_u=50 では期待数 0.62 で測定にならない。
#
# prior 側は crm_wick_pauli.jl の一般 Wick(Majorana + Pfaffian)で評価する。
# 既存の wick_gen は全Z列とホッピング1種しか扱えず、この実験ができなかった。
#
# 実行: julia --project=Hubbard_MPS_Env_v2 crm_transverse.jl
# 環境変数: L_LIST(既定 "8,16"), U(既定 4.0), NU(既定 240), N_REPEAT(既定 400)
# ============================================================
include(joinpath(@__DIR__, "crm_chain_common.jl"))
include(joinpath(@__DIR__, "crm_wick_pauli.jl"))

"""サイト i, i+1 の縦・横スピン相関と対照のオンサイト量。"""
function spin_observables(i::Int)
    u,  d  = qup(i),   qdn(i)
    u2, d2 = qup(i+1), qdn(i+1)
    obs = Obs[]
    push!(obs, Obs("ZZ onsite", [Term(1.0, [(u,3),(d,3)])], true))
    # S^z_i S^z_j = (Z_d - Z_u)(Z_d2 - Z_u2)/16
    push!(obs, Obs("SzSz r=1",
        [Term( 1/16, [(u,3),(u2,3)]), Term(-1/16, [(u,3),(d2,3)]),
         Term(-1/16, [(d,3),(u2,3)]), Term( 1/16, [(d,3),(d2,3)])], false))
    # S^x_i S^x_j = (X_uX_d + Y_uY_d)(X_u2X_d2 + Y_u2Y_d2)/16
    push!(obs, Obs("SxSx r=1",
        [Term(1/16, [(u,1),(d,1),(u2,1),(d2,1)]), Term(1/16, [(u,1),(d,1),(u2,2),(d2,2)]),
         Term(1/16, [(u,2),(d,2),(u2,1),(d2,1)]), Term(1/16, [(u,2),(d,2),(u2,2),(d2,2)])], false))
    return obs
end

function chain_T(L, t)
    T = zeros(L, L)
    for i in 1:L-1; T[i,i+1] = T[i+1,i] = -t; end
    T
end

function solve_uhf_chain(L, t, U; Nup, Ndn, iters=4000, tol=1e-12)
    T = chain_T(L, t); fill_avg = (Nup+Ndn)/(2L)
    nup = [fill_avg + 0.4*(-1)^i for i in 1:L]
    ndn = [fill_avg - 0.4*(-1)^i for i in 1:L]
    clamp!(nup, 0.02, 0.98); clamp!(ndn, 0.02, 0.98)
    nup .*= Nup/sum(nup); ndn .*= Ndn/sum(ndn)
    local Fu, Fd
    for _ in 1:iters
        Fu = eigen(Symmetric(T + diagm(U .* ndn))); Fd = eigen(Symmetric(T + diagm(U .* nup)))
        nu2 = vec(sum(abs2, Fu.vectors[:, 1:Nup]; dims=2))
        nd2 = vec(sum(abs2, Fd.vectors[:, 1:Ndn]; dims=2))
        if max(maximum(abs.(nu2 .- nup)), maximum(abs.(nd2 .- ndn))) < tol
            nup, ndn = nu2, nd2; break
        end
        nup = 0.5.*nu2 .+ 0.5.*nup; ndn = 0.5.*nd2 .+ 0.5.*ndn
    end
    Φu = Fu.vectors[:, 1:Nup]; Φd = Fd.vectors[:, 1:Ndn]
    (Cu = Φu*Φu', Cd = Φd*Φd',
     m = mean((-1)^i*(Φu[i,:]'Φu[i,:] - Φd[i,:]'Φd[i,:])/2 for i in 1:L))
end

function full_C(uhf, L)
    C = zeros(2L, 2L)
    for a in 1:L, b in 1:L
        C[qup(a), qup(b)] = uhf.Cu[a,b]; C[qdn(a), qdn(b)] = uhf.Cd[a,b]
    end
    C
end

"""スピン量子化軸を n̂(θ,φ) に向ける回転。φ を含むので一般には複素。

従来実装は φ=0 に固定した「x–z 大円上の平均」だった。⟨n_z²⟩ は大円でも
球面でも 1/3 なので Z 型観測量には正しい答えを与えるが、⟨n_x²⟩ は
大円 2/3 / 球面 1/3 と2倍ずれる。横スピン相関を測って初めて表に出る。"""
function rotate_C(C0, L, θ, φ=0.0)
    c, s = cos(θ/2), sin(θ/2)
    Um = zeros(ComplexF64, 2L, 2L)
    for i in 1:L
        a, b = qup(i), qdn(i)
        Um[a,a] = c;                 Um[a,b] = -exp(-im*φ)*s
        Um[b,a] = exp(im*φ)*s;       Um[b,b] = c
    end
    Um*C0*Um'
end

"""スピン回転平均した混合状態 prior の、項ごとの期待値。

球面上の一様平均 ∫ dΩ/4π を極角 nθ 点 × 方位角 nφ 点で取る。
nφ=1(φ=0 固定)にすると従来実装(大円平均)を再現する。"""
function symavg_P(obs, C0, L; nθ=48, nφ=16)
    P = [[0.0 for _ in o.terms] for o in obs]
    w = 1.0/(nθ*nφ)
    for k in 1:nθ, l in 1:nφ
        θ = acos(-1 + (k-0.5)*2/nθ)          # cosθ が [-1,1] 上一様
        φ = 2π*(l-0.5)/nφ
        M = majorana_M(rotate_C(C0, L, θ, φ))
        for (i,o) in enumerate(obs), (j,tm) in enumerate(o.terms)
            P[i][j] += gauss_pauli_expect(tm.sup, M)*w
        end
    end
    return P
end

"""S^x の Pauli 表現が本当に (c†_up c_dn + h.c.)/2 かを密行列で確かめる。"""
function validate_sx()
    L = 3; N = 2L
    cdag = [0.0 0.0; 1.0 0.0]; I2 = Matrix(1.0I,2,2); Fm = [1.0 0.0; 0.0 -1.0]
    Cd = [reduce(kron, [j<k ? Fm : (j==k ? cdag : I2) for j in 1:N]) for k in 1:N]
    C  = [Matrix(m') for m in Cd]
    I2c = Matrix{ComplexF64}(I,2,2)
    dense(t::Term) = t.coeff * reduce(kron,
        [(a = findfirst(x->x[1]==q, t.sup); a===nothing ? I2c : SIGMA[t.sup[a][2]]) for q in 1:N])
    maxerr = 0.0
    for i in 1:L-1
        # S^x_i = (c†_{i↑}c_{i↓} + c†_{i↓}c_{i↑})/2
        Sx(k) = 0.5*(Cd[qup(k)]*C[qdn(k)] + Cd[qdn(k)]*C[qup(k)])
        Sz(k) = 0.5*(Cd[qup(k)]*C[qup(k)] - Cd[qdn(k)]*C[qdn(k)])
        for (name, ref) in (("SxSx", Sx(i)*Sx(i+1)), ("SzSz", Sz(i)*Sz(i+1)))
            o = spin_observables(i)[name == "SxSx" ? 3 : 2]
            got = sum(dense(tm) for tm in o.terms)
            maxerr = max(maxerr, maximum(abs.(got .- ref)))
        end
    end
    @printf("  [check] S^xS^x / S^zS^z の Pauli 表現 vs フェルミオン演算子: max差 = %.2e\n", maxerr)
    @assert maxerr < 1e-12
end

function main()
    t = 1.0
    U  = parse(Float64, get(ENV, "U", "4.0"))
    nu = parse(Int, get(ENV, "NU", "240"))
    nm = parse(Int, get(ENV, "NM", "100"))
    n_repeat = parse(Int, get(ENV, "N_REPEAT", "400"))
    L_list = parse.(Int, split(get(ENV, "L_LIST", "8,16"), ","))
    chi_priors = [2, 4, 8, 32]

    println("=== 検証 ==="); flush(stdout)
    validate_sx()
    err, nt, _ = crm_wick_pauli_check()
    @printf("  [check] 一般Wick vs 密行列 Slater: max差 = %.2e (%d項)\n", err, nt)
    @assert err < 1e-10

    rows = []
    for L in L_list
        E, ψ, _, _, nel = ground_state(L, t, U, U/2; chi_max=256, nsweeps=20)
        Nup = (nel+1) ÷ 2; Ndn = nel ÷ 2
        uhf = solve_uhf_chain(L, t, U; Nup, Ndn); C0 = full_C(uhf, L)
        @printf("\n%s\nL=%d U=%.1f  E0=%.6f chi=%d | UHF スタッガード磁化 m=%.4f\n",
                "="^92, L, U, E, maxlinkdim(ψ), uhf.m); flush(stdout)

        sites = 1:L-1
        obs_by_site = [spin_observables(i) for i in sites]
        windows = [obs_support(o) for o in obs_by_site]
        rdm_ρ = window_rdms(ψ, windows)
        priors = MPS[]
        for cp in chi_priors
            σ = truncate(ψ; maxdim=cp); normalize!(σ); push!(priors, σ)
        end
        push!(priors, ψ)
        rdm_p = [window_rdms(σ, windows) for σ in priors]
        labels = vcat(["chi$c" for c in chi_priors], ["exact", "UHF", "UHFsym"])
        Muhf = majorana_M(C0)

        @printf("\n  %-12s %10s %10s %10s | %s\n", "サイト", "<SzSz>", "<SxSx>", "差",
                "prior の <SzSz> / <SxSx>  (UHF, UHFsym)")
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
            Pu = [[gauss_pauli_expect(tm.sup, Muhf) for tm in o.terms] for o in obs]
            Ps = symavg_P(obs, C0, L)
            for P in (Pu, Ps)
                push!(Pσ_all, P)
                push!(trOσ_all, [sum(tm.coeff*P[k][ti] for (ti,tm) in enumerate(obs[k].terms))
                                 for k in 1:length(obs)])
            end
            nlab = length(labels)
            es, ec = run_locals(rdm_ρ[si], Nw, ow, Pσ_all, trOσ_all;
                                nu, nm, n_repeat, seed = 880_000 + 100L + i)
            for (k,o) in enumerate(obs)
                v = var(es[:,k]); gmx = v/var(ec[:,k,length(chi_priors)+1])
                for (p,lb) in enumerate(labels)
                    Δ = Otrue[k] - trOσ_all[p][k]
                    push!(rows, (L, i, o.name, lb, Otrue[k], trOσ_all[p][k], Δ,
                                 abs(Otrue[k])>1e-6 ? abs(Δ/Otrue[k]) : NaN,
                                 v/var(ec[:,k,p]), gmx))
                end
            end
            iu = findfirst(==("UHF"), labels); isy = findfirst(==("UHFsym"), labels)
            @printf("  %-12d %10.5f %10.5f %10.2e | UHF %8.5f/%8.5f  sym %8.5f/%8.5f\n",
                    i, Otrue[2], Otrue[3], abs(Otrue[2]-Otrue[3]),
                    trOσ_all[iu][2], trOσ_all[iu][3], trOσ_all[isy][2], trOσ_all[isy][3])
            flush(stdout)
        end
    end

    out = joinpath(@__DIR__, "crm_transverse_results.tsv")
    open(out,"w") do io
        println(io, "L\tsite\tobservable\tprior\ttrue\tprior_val\tDelta\trelerr\tG\tG_max")
        for r in rows; println(io, join(r,"\t")); end
    end
    println("\nresults saved: $out")
end
get(ENV, "NO_MAIN", "0") == "1" || main()
