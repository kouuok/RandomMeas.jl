# ============================================================
# 2Dシリンダーのサイト分解: 安価な prior はどこで破綻するか
#
# 動機:
#   1Dのサイト分解(crm_site_resolved.jl)では、切断MPS prior の局所誤差
#   Δ(i) はサイトにも系サイズにもほとんど依存しなかった。一方 §9 の結果は
#   「安い prior(UHF)が破綻するのはドープ系の電荷テクスチャ」を示唆する。
#   ならば2Dドープ系では Δ(i) が空間的に大きく変化するはずである。
#   これを直接地図にする。
#
# 問い:
#   Q1. UHF prior の局所誤差はストライプの壁とドメイン内部でどう違うか?
#   Q2. その空間構造は電荷密度・スタッガード磁化のどちらと相関するか?
#   Q3. 切断MPS prior は同じ空間構造を持つか(それとも一様か)?
#
# 設計:
#   - 観測量の台が連続窓に載る量に限る(スネーク順で同一サイトの上下スピンは
#     隣接、y方向ボンドも隣接)。窓RDMサンプリングでボンド次元に依存しない。
#   - Δ の計算に Monte Carlo は不要(窓RDMから厳密)。まず全サイトの Δ 地図を
#     作り、その後 CRM 利得 G を同じ窓からサンプルする。
#
# 実行:
#   JULIA_LOAD_PATH="@:@v#.#:@stdlib" julia --project=Hubbard_MPS_Env_v2 crm_2d_site.jl
# 環境変数: UVAL(既定 8.0), NHOLE(既定 0), CHI_EXP(既定 512), LX(既定 8),
#           N_REPEAT(既定 300), DO_SAMPLE(既定 1)
# ============================================================

include(joinpath(@__DIR__, "crm_chain_common.jl"))   # 窓RDM・観測量・サンプラー

const W  = 4
const LX = parse(Int, get(ENV, "LX", "8"))
sidx(x, y) = (x - 1) * W + y

function cyl_edges(Lx, W)
    edges = Tuple{Int,Int}[]
    for x in 1:Lx, y in 1:W
        s = sidx(x, y)
        x < Lx && push!(edges, (s, sidx(x+1, y)))
        push!(edges, (s, sidx(x, y == W ? 1 : y + 1)))
    end
    return edges
end

function hubbard_cyl_mpo(sites, Lx, W, t, U, mu)
    os = OpSum()
    for (s, s2) in cyl_edges(Lx, W)
        for q in (qup, qdn)
            os += -t, "Cdag", q(s), "C", q(s2)
            os += -t, "Cdag", q(s2), "C", q(s)
        end
    end
    for s in 1:Lx*W
        os += U, "N", qup(s), "N", qdn(s)
        os += -mu, "N", qup(s)
        os += -mu, "N", qdn(s)
    end
    return MPO(os, sites)
end

function ground_state_cyl(Lx, W, t, U, mu; chi_max=512, nsweeps=14, nhole=0)
    n = Lx * W; N = 2n
    sites = siteinds("Fermion", N; conserve_qns=true)
    H = hubbard_cyl_mpo(sites, Lx, W, t, U, mu)
    occ = falses(N)
    for x in 1:Lx, y in 1:W
        s = sidx(x, y)
        isodd(x + y) ? (occ[qup(s)] = true) : (occ[qdn(s)] = true)
    end
    # ホールは中央の列から等間隔に抜く（ストライプが立ちやすい初期条件）
    if nhole > 0
        cand = [sidx(x, y) for x in 1:Lx for y in 1:W]
        pick = [cand[round(Int, (k - 0.5) * length(cand) / nhole)] for k in 1:nhole]
        for s in unique(pick)
            occ[qup(s)] = false; occ[qdn(s)] = false
        end
    end
    ψ0 = productMPS(sites, [occ[q] ? "Occ" : "Emp" for q in 1:N])
    ramp = vcat([50, 100, 200, 400], fill(chi_max, max(0, nsweeps - 4)))
    maxdims = min.(chi_max, ramp)[1:nsweeps]
    noise = vcat([1e-5, 1e-6, 1e-7, 1e-8, 1e-9], zeros(max(0, nsweeps - 5)))[1:nsweeps]
    E, ψ = dmrg(H, ψ0; nsweeps, maxdim=maxdims, cutoff=1e-10, noise, outputlevel=0)
    return E, ψ, count(occ)
end

# ------------------------------------------------------------
# UHF平均場 prior (crm_2d_doped.jl と同じ実装)
# ------------------------------------------------------------
function hopping_matrix(Lx, W, t)
    n = Lx * W; T = zeros(n, n)
    for (s, s2) in cyl_edges(Lx, W); T[s, s2] = T[s2, s] = -t; end
    return T
end

function solve_uhf(Lx, W, t, U; Nup, Ndn, iters=3000, tol=1e-11, init=:neel)
    n = Lx * W; T = hopping_matrix(Lx, W, t)
    nup = zeros(n); ndn = zeros(n); fill_avg = (Nup + Ndn) / (2n)
    for x in 1:Lx, y in 1:W
        s = sidx(x, y)
        if init == :neel
            nup[s] = fill_avg + 0.4 * (-1)^(x + y); ndn[s] = fill_avg - 0.4 * (-1)^(x + y)
        elseif init == :stripe
            wall = Lx / 2 + 0.5; dom = x <= Lx ÷ 2 ? 1.0 : -1.0
            hole = 0.35 * exp(-abs(x - wall))
            nup[s] = fill_avg - hole + dom * 0.4 * (-1)^(x + y)
            ndn[s] = fill_avg - hole - dom * 0.4 * (-1)^(x + y)
        else
            nup[s] = fill_avg + 0.3 * (rand() - 0.5); ndn[s] = fill_avg + 0.3 * (rand() - 0.5)
        end
    end
    clamp!(nup, 0.02, 0.98); clamp!(ndn, 0.02, 0.98)
    nup .*= Nup / sum(nup); ndn .*= Ndn / sum(ndn)
    local Fu, Fd
    for _ in 1:iters
        Fu = eigen(Symmetric(T + diagm(U .* ndn))); Fd = eigen(Symmetric(T + diagm(U .* nup)))
        new_up = vec(sum(abs2, Fu.vectors[:, 1:Nup]; dims=2))
        new_dn = vec(sum(abs2, Fd.vectors[:, 1:Ndn]; dims=2))
        if max(maximum(abs.(new_up .- nup)), maximum(abs.(new_dn .- ndn))) < tol
            nup, ndn = new_up, new_dn; break
        end
        nup = 0.5 .* new_up .+ 0.5 .* nup; ndn = 0.5 .* new_dn .+ 0.5 .* ndn
    end
    Φu = Fu.vectors[:, 1:Nup]; Φd = Fd.vectors[:, 1:Ndn]
    Cu = Φu * Φu'; Cd = Φd * Φd'
    E = tr(T * Cu) + tr(T * Cd) + U * sum(diag(Cu) .* diag(Cd))
    return (nup=diag(Cu), ndn=diag(Cd), Cu=Cu, Cd=Cd, E=E)
end

function solve_uhf_best(Lx, W, t, U; Nup, Ndn)
    Random.seed!(31); best = nothing
    for ini in (:neel, :stripe, :random)
        uhf = solve_uhf(Lx, W, t, U; Nup, Ndn, init=ini)
        @printf("  UHF init=%-7s E=%.6f\n", ini, uhf.E)
        (best === nothing || uhf.E < best[2].E) && (best = (ini, uhf))
    end
    @printf("  -> selected init=%s (E=%.6f)\n", best[1], best[2].E)
    return best[2]
end

site_spin(q) = (div(q + 1, 2), isodd(q) ? 1 : 2)

"""Slater行列式 UHF に対する Pauli 項の期待値（Wickの定理）。
   Z_q = 1-2n_q、および同一スピン間の 2 演算子項のみ扱う。"""
function wick_term(sup::Vector{Tuple{Int,Int}}, uhf)
    isempty(sup) && return 1.0
    all(a == 3 for (_, a) in sup) || return 0.0    # 本スクリプトはZ列のみ
    Cs = (uhf.Cu, uhf.Cd)
    if length(sup) == 1
        s, σ = site_spin(sup[1][1]); return 1 - 2 * Cs[σ][s, s]
    elseif length(sup) == 2
        (qa, _), (qb, _) = sup[1], sup[2]
        sa, σa = site_spin(qa); sb, σb = site_spin(qb)
        Ca = Cs[σa][sa, sa]; Cb = Cs[σb][sb, sb]
        cross = σa == σb ? Cs[σa][sa, sb] * Cs[σa][sb, sa] : 0.0
        return 1 - 2Ca - 2Cb + 4 * (Ca * Cb - cross)
    end
    error("wick_term: 3体以上は未対応")
end

# ------------------------------------------------------------
# 2Dサイトの観測量（窓が連続になる量に限る）
# ------------------------------------------------------------
"""サイト s のオンサイト量と、スネーク順で隣接するサイト s+1 との相関。"""
function site_observables_2d(s::Int; bond::Bool)
    u, d = qup(s), qdn(s)
    obs = Obs[]
    push!(obs, Obs("n", [Term(1.0, Tuple{Int,Int}[]), Term(-0.5, [(u,3)]), Term(-0.5, [(d,3)])], false))
    push!(obs, Obs("Sz", [Term(-0.25, [(u,3)]), Term(0.25, [(d,3)])], false))
    push!(obs, Obs("DoubleOcc", [Term(0.25, Tuple{Int,Int}[]), Term(-0.25, [(u,3)]),
                                 Term(-0.25, [(d,3)]), Term(0.25, [(u,3),(d,3)])], false))
    push!(obs, Obs("ZZ onsite", [Term(1.0, [(u,3),(d,3)])], true))
    if bond
        u2, d2 = qup(s+1), qdn(s+1)
        push!(obs, Obs("ZZ up-up nb", [Term(1.0, [(u,3),(u2,3)])], true))
        push!(obs, Obs("SzSz nb",
            [Term( 1/16, [(u,3),(u2,3)]), Term(-1/16, [(u,3),(d2,3)]),
             Term(-1/16, [(d,3),(u2,3)]), Term( 1/16, [(d,3),(d2,3)])], false))
    end
    return obs
end

# ------------------------------------------------------------
# メイン
# ------------------------------------------------------------
function main()
    t = 1.0
    U = parse(Float64, get(ENV, "UVAL", "8.0"))
    nhole = parse(Int, get(ENV, "NHOLE", "0"))
    chi_exp = parse(Int, get(ENV, "CHI_EXP", "512"))
    n_repeat = parse(Int, get(ENV, "N_REPEAT", "300"))
    do_sample = get(ENV, "DO_SAMPLE", "1") == "1"
    nu, nm = 50, 100
    chi_priors = [8, 32]
    n = LX * W; mu = U/2

    @printf("=== 2D site-resolved: W=%d Lx=%d (N=%d qubits), U=%.1f, nhole=%d ===\n",
            W, LX, 2n, U, nhole); flush(stdout)
    t0 = time()
    E, ψ, nel = ground_state_cyl(LX, W, t, U, mu; chi_max=chi_exp, nhole)
    @printf("DMRG: E0=%.6f  N_el=%d  maxlinkdim=%d  (%.0fs)\n",
            E, nel, maxlinkdim(ψ), time()-t0); flush(stdout)

    Nup = (nel + 1) ÷ 2; Ndn = nel ÷ 2
    uhf = solve_uhf_best(LX, W, t, U; Nup, Ndn); flush(stdout)

    obs_by_site = [site_observables_2d(s; bond = s < n) for s in 1:n]
    windows = [obs_support(o) for o in obs_by_site]
    rdm_ρ = window_rdms(ψ, windows)

    priors = MPS[]
    for cp in chi_priors
        σ = truncate(ψ; maxdim=cp); normalize!(σ); push!(priors, σ)
    end
    push!(priors, ψ)
    prior_fid = [abs2(inner(σ, ψ)) for σ in priors]
    rdm_p = [window_rdms(σ, windows) for σ in priors]
    plabels = vcat(["chi$(c)" for c in chi_priors], ["exact"])
    @printf("prior fidelities: %s\n",
            join([@sprintf("%s=%.3e", l, f) for (l,f) in zip(plabels, prior_fid)], ", ")); flush(stdout)

    rows = []
    for s in 1:n
        obs = obs_by_site[s]; q1, q2 = windows[s]; Nw = q2 - q1 + 1
        ow = shift_obs(obs, q1 - 1)
        Mats = [[term_window_matrix(tm, Nw) for tm in o.terms] for o in ow]
        Pρ = [[expect_rdm(rdm_ρ[s], M) for M in Ms] for Ms in Mats]
        Otrue = [sum(tm.coeff * Pρ[k][ti] for (ti, tm) in enumerate(obs[k].terms))
                 for k in 1:length(obs)]
        # prior 群: 切断MPS + 完全 + UHF
        Pσ_all = Vector{Vector{Vector{Float64}}}(); trOσ_all = Vector{Vector{Float64}}()
        for p in 1:length(priors)
            Pσ = [[expect_rdm(rdm_p[p][s], M) for M in Ms] for Ms in Mats]
            push!(Pσ_all, Pσ)
            push!(trOσ_all, [sum(tm.coeff * Pσ[k][ti] for (ti, tm) in enumerate(obs[k].terms))
                             for k in 1:length(obs)])
        end
        Pu = [[wick_term(tm.sup, uhf) for tm in o.terms] for o in obs]
        push!(Pσ_all, Pu)
        push!(trOσ_all, [sum(tm.coeff * Pu[k][ti] for (ti, tm) in enumerate(obs[k].terms))
                         for k in 1:length(obs)])
        labels = vcat(plabels, ["UHF"])

        Gs = fill(NaN, length(obs), length(labels))
        if do_sample
            es, ec = run_locals(rdm_ρ[s], Nw, ow, Pσ_all, trOσ_all;
                                nu, nm, n_repeat, seed = 500_000 + 100s)
            for k in 1:length(obs)
                v = var(es[:, k])
                for p in 1:length(labels); Gs[k, p] = v / var(ec[:, k, p]); end
            end
        end
        x = (s - 1) ÷ W + 1; y = (s - 1) % W + 1
        for (k, o) in enumerate(obs), (p, lb) in enumerate(labels)
            Δ = Otrue[k] - trOσ_all[p][k]
            push!(rows, (U, nhole, s, x, y, o.name, lb, Otrue[k], Δ,
                         abs(Otrue[k]) > 1e-3 ? abs(Δ/Otrue[k]) : NaN, Gs[k, p],
                         Otrue[1], Otrue[2]))   # <n>(s), <Sz>(s) を毎行に付ける
        end
        s % max(1, n ÷ 8) == 0 && (@printf("  site %d/%d (%.0fs)\n", s, n, time()-t0); flush(stdout))
    end

    out = joinpath(@__DIR__, @sprintf("crm_2d_site_U%s_h%d.tsv",
                                      replace(@sprintf("%.1f", U), "."=>"p"), nhole))
    open(out, "w") do io
        println(io, "U\tnhole\tsite\tx\ty\tobservable\tprior\ttrue\tDelta\trelerr\tG\tdens\tSz")
        for r in rows; println(io, join(r, "\t")); end
    end
    println("results saved: $out"); flush(stdout)

    # ---- 画面サマリ: prior ごとの相対誤差の空間構造 ----
    for o in ["ZZ onsite", "DoubleOcc", "SzSz nb"]
        println("\n[$o] 列 x ごとの median 相対誤差")
        @printf("%-8s", "prior"); for x in 1:LX; @printf("%9s", "x=$x"); end
        @printf("%11s\n", "max/min")
        for lb in vcat(plabels, ["UHF"])
            lb == "exact" && continue
            @printf("%-8s", lb); v = Float64[]
            for x in 1:LX
                sel = [r[10] for r in rows if r[6]==o && r[7]==lb && r[4]==x && !isnan(r[10])]
                m = isempty(sel) ? NaN : median(sel); push!(v, m)
                @printf("%9.2e", m)
            end
            vv = filter(!isnan, v)
            @printf("%11.1f\n", isempty(vv) ? NaN : maximum(vv)/minimum(vv))
        end
    end
    println("\n列ごとの電子密度 <n> と スタッガード磁化")
    @printf("%-10s", "x"); for x in 1:LX; @printf("%9d", x); end; println()
    @printf("%-10s", "<n>")
    for x in 1:LX
        v = [r[12] for r in rows if r[4]==x]; @printf("%9.3f", isempty(v) ? NaN : mean(v))
    end
    println()
    @printf("%-10s", "m_stag")
    for x in 1:LX
        v = [(-1)^(x+r[5]) * r[13] for r in rows if r[4]==x]
        @printf("%9.3f", isempty(v) ? NaN : mean(v))
    end
    println()
    return rows
end

rows = main()
