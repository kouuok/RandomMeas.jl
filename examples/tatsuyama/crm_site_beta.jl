# ============================================================
# 実験A: 「損の領域」は最適係数CRMで修復できるか
#
# 動機:
#   crm_site_resolved.jl で、χ_p=2 の粗い prior を使うと δ=1/8 のドープ鎖の
#   2〜3割のサイトで CRM が有意に損をする（最悪 G=0.039、標準シャドウの25倍悪い）
#   ことが分かった。制御変量法の理論では、係数 β を最適化すれば漸近的に
#   G≥1 が保証される。しかし β は同じ n_u 設定から推定する（プラグイン）ので、
#   有限標本ノイズが修復分を食い潰さない保証はない。
#
# 問い:
#   Q1. プラグイン β̂ は本当に損を消すか? どこまで戻せるか?
#   Q2. そのために何設定 (n_u) 必要か?
#   Q3. 損をしていない領域で β̂ を使うと、逆にどれだけ損するか (保険料)?
#
# 設計:
#   - 各設定 u について m_ρ(u)（実験側ショット平均）と Y(u)=E[X_σ|u]（厳密）を保持
#   - 3つの推定量を同一データで比較:
#       標準   : m̄_ρ
#       CRM    : m̄_ρ - (Ȳ - Tr[Oσ])                     （β=1）
#       最適β  : m̄_ρ - β̂(Ȳ - Tr[Oσ]),  β̂ = Cov(m_ρ,Y)/Var(Y)
#   - n_u は入れ子（先頭 50 / 100 / 200 設定）で切り出すので 1 回のサンプリングで済む
#
# 実行:
#   JULIA_LOAD_PATH="@:@v#.#:@stdlib" julia --project=Hubbard_MPS_Env_v2 crm_site_beta.jl
# 環境変数: L_LIST(既定 "16,64"), DOPE_LIST(既定 "0.0,0.125"),
#           N_REPEAT(既定 500), NU_CUTS(既定 "50,100,200"), CHI_EXP(既定 256)
# ============================================================

include(joinpath(@__DIR__, "crm_chain_common.jl"))


"""prior側の条件付き平均 Y(u)=E[X_σ|u] の分散を、標本ではなく厳密に計算する。

Y(u) = Σ_t c_t · 1[基底が項tの台で一致] · 3^{|A_t|} · <P_t>_σ であり、
基底は一様独立なので Pr[項t が一致] = 3^{-|A_t|}、
Pr[項t,t' が同時に一致] = 3^{-|A_t ∪ A_t'|}（重なり上でPauli指定が矛盾しなければ、
矛盾すれば 0）である。したがって Var(Y) は prior の期待値だけから閉形式で出る。

これが要る理由: β̂ = Cov/Var(Y) を標本 Var(Y) で割ると、prior がその観測量を
ほぼ 0 と予言する場合（例: χ_p=2 の prior は二重占有を厳密に 0 と予言する）に
Var(Y)→0 となり、β̂ が 0/0 で発散する。厳密な Var(Y) を使えばこの発散は起きない。
元の CRM 推定器を「prior 側はサンプルせず厳密に計算する」ことで直したのと同じ原理。"""
function exact_var_Y(o::Obs, Pσ::Vector{Float64})
    nt = length(o.terms)
    a = [o.terms[i].coeff * 3.0^length(o.terms[i].sup) * Pσ[i] for i in 1:nt]
    m1 = sum(a[i] * 3.0^(-length(o.terms[i].sup)) for i in 1:nt)
    m2 = 0.0
    for i in 1:nt, j in 1:nt
        si = o.terms[i].sup; sj = o.terms[j].sup
        compatible = true
        for (q, pa) in si
            k = findfirst(x -> x[1] == q, sj)
            if k !== nothing && sj[k][2] != pa; compatible = false; break; end
        end
        compatible || continue
        qs = Set(q for (q, _) in si)
        for (q, _) in sj; push!(qs, q); end
        m2 += a[i] * a[j] * 3.0^(-length(qs))
    end
    return max(m2 - m1^2, 0.0)
end

# ------------------------------------------------------------
# サンプリング本体: 設定ごとの (m_ρ, Y) を保持して3推定量を作る
# ------------------------------------------------------------
function run_locals_beta(ρw, Nw, obs_w, Pσ_all, trOσ_all;
                         nu_max, nu_cuts, nm, n_repeat, seed)
    # Y の厳密分散は prior だけで決まるので事前に一度計算しておく
    vY_exact = [exact_var_Y(obs_w[k], Pσ_all[p][k])
                for k in 1:length(obs_w), p in 1:length(Pσ_all)]
    Random.seed!(seed)
    nobs = length(obs_w); np = length(Pσ_all); nc = length(nu_cuts)
    bits = zeros(Int, Nw)
    cum_cache = Dict{Vector{Int},Vector{Float64}}()

    mρ = zeros(nu_max, nobs)              # 設定ごとの実験側ショット平均
    Yv = zeros(nu_max, nobs, np)          # 設定ごとの prior 側厳密条件付き平均
    est_std  = zeros(n_repeat, nobs, nc)
    est_crm  = zeros(n_repeat, nobs, np, nc)
    est_beta = zeros(n_repeat, nobs, np, nc)
    beta_hat = zeros(n_repeat, nobs, np, nc)
    est_shr   = zeros(n_repeat, nobs, np, nc)   # 縮小推定
    beta_shr  = zeros(n_repeat, nobs, np, nc)
    est_split = zeros(n_repeat, nobs, np, nc)   # 標本分割（交差適用）
    est_ev    = zeros(n_repeat, nobs, np, nc)   # 厳密Var(Y)版
    beta_ev   = zeros(n_repeat, nobs, np, nc)
    est_sev   = zeros(n_repeat, nobs, np, nc)   # 標本分割×厳密Var(Y)
    xs = zeros(nobs)

    for rep in 1:n_repeat
        for u in 1:nu_max
            basis = rand(1:3, Nw)
            cum = get!(cum_cache, copy(basis)) do
                cumsum(basis_probs(ρw, basis))
            end
            fill!(xs, 0.0)
            for _ in 1:nm
                draw_bits!(bits, cum)
                for k in 1:nobs
                    xs[k] += estimate_obs(obs_w[k], basis, bits)
                end
            end
            for k in 1:nobs
                mρ[u, k] = xs[k] / nm
                for p in 1:np
                    Yv[u, k, p] = exact_prior_mean(obs_w[k], basis, Pσ_all[p][k])
                end
            end
        end
        for (ci, nu) in enumerate(nu_cuts)
            for k in 1:nobs
                mv = @view mρ[1:nu, k]
                mbar = mean(mv)
                est_std[rep, k, ci] = mbar
                vm = var(mv)
                for p in 1:np
                    yv = @view Yv[1:nu, k, p]
                    ybar = mean(yv)
                    est_crm[rep, k, p, ci] = mbar - ybar + trOσ_all[p][k]
                    vY = var(yv)
                    # Var(Y)=0 は「prior の条件付き平均が設定に依らない」ケース。
                    # 制御変量として使えないので β̂=0（=標準シャドウ）に落とす。
                    cv = (vY > 1e-14 && vm > 1e-14) ? cov(mv, yv) : 0.0
                    β = vY > 1e-14 ? cv / vY : 0.0
                    beta_hat[rep, k, p, ci] = β
                    est_beta[rep, k, p, ci] = mbar - β * (ybar - trOσ_all[p][k])

                    # --- 防御版 ---
                    # 素朴な β̂ は Var(Y)→0 で発散する。標本相関 r² は
                    # 「制御変量として使う価値」そのものなので、これで縮小する:
                    #   β_shrunk = β̂ · r²/(r² + 1/n_u)
                    # r² ≫ 1/n_u（有意な相関）なら β̂ をそのまま、
                    # r² ≲ 1/n_u（相関がノイズと区別できない）なら 0 に落ちる。
                    r2 = (vY > 1e-14 && vm > 1e-14) ? cv^2 / (vY * vm) : 0.0
                    βs = β * r2 / (r2 + 1/nu)
                    beta_shr[rep, k, p, ci] = βs
                    est_shr[rep, k, p, ci] = mbar - βs * (ybar - trOσ_all[p][k])

                    # 標本分割: 前半で β を推定し後半に適用（自己相関バイアスの除去）。
                    # 推定と適用が独立になる代わり、各半分の標本数は半減する。
                    h = nu ÷ 2
                    m1 = @view mρ[1:h, k];      y1 = @view Yv[1:h, k, p]
                    m2 = @view mρ[h+1:nu, k];   y2 = @view Yv[h+1:nu, k, p]
                    v1 = var(y1); v2 = var(y2)
                    β1 = v1 > 1e-14 ? cov(m1, y1) / v1 : 0.0
                    β2 = v2 > 1e-14 ? cov(m2, y2) / v2 : 0.0
                    # 前半の β を後半に、後半の β を前半に当てて平均（交差適用）
                    est_split[rep, k, p, ci] = 0.5 * (
                        (mean(m2) - β1 * (mean(y2) - trOσ_all[p][k])) +
                        (mean(m1) - β2 * (mean(y1) - trOσ_all[p][k])))

                    # 厳密Var(Y)版: 分母をサンプルから推定しない
                    ve = vY_exact[k, p]
                    βe = ve > 1e-14 ? cv / ve : 0.0
                    beta_ev[rep, k, p, ci] = βe
                    est_ev[rep, k, p, ci] = mbar - βe * (ybar - trOσ_all[p][k])

                    # 標本分割 × 厳密Var(Y): 分散の暴走（厳密Varで解消）と
                    # プラグインのバイアス（分割で解消）を同時に潰す
                    β1e = ve > 1e-14 ? cov(m1, y1) / ve : 0.0
                    β2e = ve > 1e-14 ? cov(m2, y2) / ve : 0.0
                    est_sev[rep, k, p, ci] = 0.5 * (
                        (mean(m2) - β1e * (mean(y2) - trOσ_all[p][k])) +
                        (mean(m1) - β2e * (mean(y1) - trOσ_all[p][k])))
                end
            end
        end
    end
    return est_std, est_crm, est_beta, beta_hat, est_shr, beta_shr, est_split, est_ev, beta_ev, est_sev
end

# ------------------------------------------------------------
# メイン
# ------------------------------------------------------------
function main()
    t, U = 1.0, 4.0; mu = U/2
    L_list = parse.(Int, split(get(ENV, "L_LIST", "16,64"), ","))
    dope_list = parse.(Float64, split(get(ENV, "DOPE_LIST", "0.0,0.125"), ","))
    nu_cuts = parse.(Int, split(get(ENV, "NU_CUTS", "50,100,200"), ","))
    chi_priors = [2, 4, 8, 16, 32]
    nm = 100
    n_repeat = parse(Int, get(ENV, "N_REPEAT", "500"))
    chi_exp = parse(Int, get(ENV, "CHI_EXP", "256"))
    nu_max = maximum(nu_cuts)

    println("=== validation (L=4) ==="); flush(stdout)
    dense_chain_check()

    rows = []
    for L in L_list, δ in dope_list
        println("\n", "="^72)
        @printf("L=%d, U=%.1f, δ=%.3f, chi_exp=%d, n_repeat=%d, nu_cuts=%s\n",
                L, U, δ, chi_exp, n_repeat, string(nu_cuts)); flush(stdout)
        tstart = time()
        E, ψ, _, _, nel = ground_state(L, t, U, mu; chi_max=chi_exp, δ)
        @printf("  DMRG: E0=%.6f  N_el=%d  maxlinkdim=%d  (%.0fs)\n",
                E, nel, maxlinkdim(ψ), time()-tstart); flush(stdout)

        obs_by_site = [site_observables(i; bond = i < L) for i in 1:L]
        windows = [obs_support(o) for o in obs_by_site]
        priors = MPS[]
        for chi_p in chi_priors
            σ = truncate(ψ; maxdim=chi_p); normalize!(σ); push!(priors, σ)
        end
        push!(priors, ψ)                                   # 完全prior = 天井
        prior_fid = [abs2(inner(σ, ψ)) for σ in priors]
        rdm_ρ = window_rdms(ψ, windows)
        rdm_p = [window_rdms(σ, windows) for σ in priors]
        @printf("  prior fidelities: %s\n",
                join([@sprintf("%.2e", f) for f in prior_fid], ", ")); flush(stdout)

        tstart = time()
        for i in 1:L
            obs = obs_by_site[i]; q1, q2 = windows[i]; Nw = q2 - q1 + 1
            obs_w = shift_obs(obs, q1 - 1)
            Mats = [[term_window_matrix(tm, Nw) for tm in o.terms] for o in obs_w]
            Pρ = [[expect_rdm(rdm_ρ[i], M) for M in Ms] for Ms in Mats]
            Otrue = [sum(tm.coeff * Pρ[k][ti] for (ti, tm) in enumerate(obs[k].terms))
                     for k in 1:length(obs)]
            Pσ_all = Vector{Vector{Vector{Float64}}}(); trOσ_all = Vector{Vector{Float64}}()
            for p in 1:length(priors)
                Pσ = [[expect_rdm(rdm_p[p][i], M) for M in Ms] for Ms in Mats]
                push!(Pσ_all, Pσ)
                push!(trOσ_all, [sum(tm.coeff * Pσ[k][ti] for (ti, tm) in enumerate(obs[k].terms))
                                 for k in 1:length(obs)])
            end

            es, ec, eb, bh, esh, bsh, esp, eev, bev, esv = run_locals_beta(rdm_ρ[i], Nw, obs_w, Pσ_all, trOσ_all;
                                             nu_max, nu_cuts, nm, n_repeat,
                                             seed = 90_000 + 1000L + i)

            for (ci, nu) in enumerate(nu_cuts), (k, o) in enumerate(obs)
                a = es[:, k, ci]; vstd = var(a)
                # G は分散比なのでバイアスを見ない。プラグイン β̂ は同一データで
                # 係数を推定するため原理的にバイアスを持つので、MSE でも比較する。
                mse(v) = var(v) + (mean(v) - Otrue[k])^2
                mse_std = mse(a)
                Gmax = vstd / var(ec[:, k, end, ci])
                for p in 1:length(chi_priors)
                    bc = ec[:, k, p, ci]; bb = eb[:, k, p, ci]
                    bs = esh[:, k, p, ci]; bp = esp[:, k, p, ci]; be = eev[:, k, p, ci]
                    bv = esv[:, k, p, ci]
                    Gc = vstd / var(bc); Gb = vstd / var(bb)
                    Gs = vstd / var(bs); Gp = vstd / var(bp); Ge = vstd / var(be)
                    Mc = mse_std/mse(bc); Mb = mse_std/mse(bb)
                    Ms = mse_std/mse(bs); Mp = mse_std/mse(bp); Me = mse_std/mse(be)
                    Gv = vstd / var(bv); Mv = mse_std/mse(bv)
                    vlo, vhi = boot_ratio(a, bv; rng = MersenneTwister(41 + i))
                    # バイアスを推定量の標準偏差で規格化（1を超えると無視できない）
                    bias_n(v) = (mean(v) - Otrue[k]) / sqrt(var(v))
                    clo, chi_ = boot_ratio(a, bc; rng = MersenneTwister(41 + i))
                    blo, bhi  = boot_ratio(a, bb; rng = MersenneTwister(41 + i))
                    slo, shi  = boot_ratio(a, bs; rng = MersenneTwister(41 + i))
                    plo, phi  = boot_ratio(a, bp; rng = MersenneTwister(41 + i))
                    elo, ehi  = boot_ratio(a, be; rng = MersenneTwister(41 + i))
                    Δ = Otrue[k] - trOσ_all[p][k]
                    push!(rows, (L, δ, i, chi_priors[p], o.name, nu, Otrue[k], Δ,
                                 prior_fid[p], Gc, clo, chi_, Gb, blo, bhi,
                                 Gs, slo, shi, Gp, plo, phi, Ge, elo, ehi,
                                 Gv, vlo, vhi, Mc, Mb, Ms, Mp, Me, Mv,
                                 bias_n(bc), bias_n(bb), bias_n(bs), bias_n(bp),
                                 bias_n(be), bias_n(bv),
                                 median(bh[:, k, p, ci]), std(bh[:, k, p, ci]),
                                 median(bsh[:, k, p, ci]), median(bev[:, k, p, ci]),
                                 std(bev[:, k, p, ci]), Gmax))
                end
            end
            if i % max(1, L ÷ 8) == 0 || i == L
                @printf("    site %2d/%d (%.0fs)\n", i, L, time()-tstart); flush(stdout)
            end
        end
        @printf("  run: %.0fs\n", time()-tstart); flush(stdout)

        # 画面サマリ: 損をしていた領域 (χ_p=2) の修復状況
        for nu in nu_cuts
            sel = [r for r in rows if r[1]==L && r[2]==δ && r[4]==2 && r[6]==nu]
            isempty(sel) && continue
            n = length(sel)
            @printf("  [χ_p=2, n_u=%3d]  %-11s %-11s %-11s %-11s %-11s %-11s\n", nu,
                    "CRM(b=1)", "plain-b", "shrunk", "split", "exactVar", "split+eV")
            gcol = (10, 13, 16, 19, 22, 25); hcol = (12, 15, 18, 21, 24, 27)
            mcol = (28, 29, 30, 31, 32, 33); bcol = (34, 35, 36, 37, 38, 39)
            @printf("      min G      "); for j in gcol; @printf("%-11.3g ", minimum(r[j] for r in sel)); end; println()
            @printf("      lost       "); for j in hcol; @printf("%-11s ", "$(count(r -> r[j] < 0.9, sel))/$n"); end; println()
            @printf("      median G   "); for j in gcol; @printf("%-11.2f ", median(r[j] for r in sel)); end; println()
            @printf("      median MSE "); for j in mcol; @printf("%-11.2f ", median(r[j] for r in sel)); end; println()
            @printf("      max|bias/sd|"); for j in bcol; @printf("%-11.2f ", maximum(abs(r[j]) for r in sel)); end; println()
        end
        flush(stdout)
    end

    out = joinpath(@__DIR__, "crm_site_beta_results.tsv")
    open(out, "w") do io
        println(io, "L\tdoping\tsite\tchi_prior\tobservable\tnu\ttrue\tDelta\tprior_fid\t" *
                    "G_crm\tG_crm_lo\tG_crm_hi\tG_beta\tG_beta_lo\tG_beta_hi\t" *
                    "G_shrunk\tG_shrunk_lo\tG_shrunk_hi\tG_split\tG_split_lo\tG_split_hi\t" *
                    "G_exactvar\tG_exactvar_lo\tG_exactvar_hi\t" *
                    "G_splitev\tG_splitev_lo\tG_splitev_hi\t" *
                    "M_crm\tM_beta\tM_shrunk\tM_split\tM_exactvar\tM_splitev\t" *
                    "bias_crm\tbias_beta\tbias_shrunk\tbias_split\tbias_exactvar\tbias_splitev\t" *
                    "beta_median\tbeta_std\tbeta_shrunk_median\tbeta_ev_median\tbeta_ev_std\tG_max")
        for r in rows; println(io, join(r, "\t")); end
    end
    println("\nresults saved: $out"); flush(stdout)
    return rows
end

rows = main()
