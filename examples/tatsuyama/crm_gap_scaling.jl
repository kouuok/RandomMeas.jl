# ============================================================
# 実験B: 局所prior誤差 Δ の L 依存性は「ギャップの有無」で決まるか
#
# 動機:
#   crm_site_resolved.jl の結果を精査すると、half-filling (U=4) では
#   median|Δ| が L に依存しないのに、δ=1/8 にドープすると全観測量で
#   L とともに増加した (ZZ onsite: 5.8e-3 → 3.4e-2, L=8→128)。
#   これは「ドープにより電荷セクターが gapless (朝永-Luttinger液体) になり、
#   エンタングルメントが ~(c/3)log L で増えるので、固定χの prior が
#   捉えられる割合が落ちる」と解釈できる。局所性の主張に適用条件を
#   与える重要な限界なので、直接検証する。
#
# 予言:
#   ギャップのある系  : median|Δ|(L) → 定数（Lに非依存）
#   gapless な系      : median|Δ|(L) ~ log L で増加
#
# 1Dハバード鎖の電荷ギャップは half-filling で任意の U>0 に対し開くが、
# 小さい U では指数的に小さい (Δ_c ~ e^{-2πt/U}, 相関長 ξ_c ~ e^{+2πt/U})。
# したがって U を振ると:
#   U=0      : 厳密に gapless          → 増加
#   U=1      : ξ_c が巨大 → 実効的に gapless → L≲ξ_c では増加
#   U=4, 8   : ξ_c が O(1)             → すぐ飽和
#   U=4,8 δ=1/8 : 電荷セクター gapless  → 増加
# という階段が見えるはずで、これが見えれば因果が確定する。
#
# 設計上の要点:
#   Δ = <P>_ρ - <P>_σ は窓の縮約密度行列から厳密に計算できるので、
#   この実験に Monte Carlo は一切不要（DMRG と RDM 縮約のみ）。
#   あわせて鎖中央のボンドエントロピー S_c(L) も測る。開放端では
#   gapless なら S_c = (c/6)log L + const、gapped なら飽和するので、
#   Δ の振る舞いと中心電荷を同時に読める独立なチェックになる。
#
# 実行:
#   JULIA_LOAD_PATH="@:@v#.#:@stdlib" julia --project=Hubbard_MPS_Env_v2 crm_gap_scaling.jl
# 環境変数: L_LIST(既定 "8,16,32,64,128"), CHI_EXP(既定 256), NSWEEPS(既定 16)
# ============================================================

include(joinpath(@__DIR__, "crm_chain_common.jl"))

"""鎖中央付近のボンドエントロピーの平均。

単一ボンドの値は §2c で見たのと同じ偶奇交替を拾ってしまうので、
中央の 4 ボンド（=2 サイト分、交替の1周期）で平均して平滑化する。"""
function center_entropy(ψ::MPS)
    N = length(ψ)
    ψo = copy(ψ); orthogonalize!(ψo, 1)
    vals = Float64[]
    for b in (N ÷ 2 - 1):(N ÷ 2 + 2)
        (2 <= b <= N - 1) || continue
        orthogonalize!(ψo, b)
        _, Sv, _ = svd(ψo[b], (linkinds(ψo)[b-1], siteinds(ψo)[b]))
        p = diag(Array(Sv, inds(Sv)...)).^2
        p = p[p .> 1e-14]
        push!(vals, -sum(p .* log.(p)))
    end
    return mean(vals)
end

function main()
    t = 1.0
    L_list = parse.(Int, split(get(ENV, "L_LIST", "8,16,32,64,128"), ","))
    chi_exp = parse(Int, get(ENV, "CHI_EXP", "256"))
    nsweeps = parse(Int, get(ENV, "NSWEEPS", "16"))
    chi_priors = [4, 8, 16, 32]

    # (U, δ, ラベル, 予想されるギャップ構造)
    configs = [(0.0, 0.0,   "U=0  half",      "charge gapless"),
               (1.0, 0.0,   "U=1  half",      "charge gap ~ e^{-2pi t/U} (huge xi)"),
               (2.0, 0.0,   "U=2  half",      "charge gap small"),
               (4.0, 0.0,   "U=4  half",      "charge gapped"),
               (8.0, 0.0,   "U=8  half",      "charge strongly gapped"),
               (4.0, 0.125, "U=4  delta=1/8", "charge gapless (doped)"),
               (8.0, 0.125, "U=8  delta=1/8", "charge gapless (doped)")]

    println("=== validation (L=4) ==="); flush(stdout)
    dense_chain_check()

    rows = []; srows = []
    for (U, δ, lab, expectation) in configs
        mu = U/2
        println("\n", "="^76)
        @printf("%s   [%s]\n", lab, expectation); flush(stdout)
        for L in L_list
            tstart = time()
            E, ψ, _, _, nel = ground_state(L, t, U, mu; chi_max=chi_exp, nsweeps, δ)
            Sc = center_entropy(ψ)
            @printf("  L=%3d: E0=%12.6f  N_el=%3d  chi=%3d  S_center=%.4f  (%.0fs)\n",
                    L, E, nel, maxlinkdim(ψ), Sc, time()-tstart); flush(stdout)
            push!(srows, (U, δ, L, E, nel, maxlinkdim(ψ), Sc))

            obs_by_site = [site_observables(i; bond = i < L) for i in 1:L]
            windows = [obs_support(o) for o in obs_by_site]
            rdm_ρ = window_rdms(ψ, windows)
            for chi_p in chi_priors
                σ = truncate(ψ; maxdim=chi_p); normalize!(σ)
                fid = abs2(inner(σ, ψ))
                rdm_σ = window_rdms(σ, windows)
                acc = Dict{String,Vector{Float64}}()
                accP = Dict{String,Vector{Float64}}()
                accE = Dict{String,Vector{Float64}}()
                for i in 1:L
                    obs = obs_by_site[i]; q1, q2 = windows[i]; Nw = q2 - q1 + 1
                    obs_w = shift_obs(obs, q1 - 1)
                    for (o, ow) in zip(obs, obs_w)
                        Ms = [term_window_matrix(tm, Nw) for tm in ow.terms]
                        vρ = sum(tm.coeff * expect_rdm(rdm_ρ[i], Ms[ti])
                                 for (ti, tm) in enumerate(ow.terms))
                        vσ = sum(tm.coeff * expect_rdm(rdm_σ[i], Ms[ti])
                                 for (ti, tm) in enumerate(ow.terms))
                        push!(get!(acc, o.name, Float64[]), abs(vρ - vσ))
                        push!(get!(accP, o.name, Float64[]), abs(vρ))
                        # 相対誤差 ε=Δ/<P>。<P> が対称性で 0 になる観測量
                        # (half-filling の n, Sz) では ε が発散するので除外する。
                        abs(vρ) > 1e-3 && push!(get!(accE, o.name, Float64[]),
                                                abs((vρ - vσ) / vρ))
                    end
                end
                for (name, v) in acc
                    e = get(accE, name, Float64[])
                    push!(rows, (U, δ, L, chi_p, name, fid, Sc, median(v), mean(v), maximum(v),
                                 median(accP[name]),
                                 isempty(e) ? NaN : median(e),
                                 isempty(e) ? NaN : maximum(e), length(e)))
                end
            end
        end
    end

    out = joinpath(@__DIR__, "crm_gap_scaling_results.tsv")
    open(out, "w") do io
        println(io, "U\tdoping\tL\tchi_prior\tobservable\tprior_fid\tS_center\t" *
                    "median_absDelta\tmean_absDelta\tmax_absDelta\t" *
                    "median_absP\tmedian_relerr\tmax_relerr\tn_relerr")
        for r in rows; println(io, join(r, "\t")); end
    end
    outs = joinpath(@__DIR__, "crm_gap_states.tsv")
    open(outs, "w") do io
        println(io, "U\tdoping\tL\tE0\tN_el\tmaxlinkdim\tS_center")
        for r in srows; println(io, join(r, "\t")); end
    end
    println("\nresults saved: $out, $outs"); flush(stdout)

    # ---- 画面サマリ ----
    println("\n", "="^76)
    println("median 相対誤差 ε=|Δ/<P>| の L 依存 (χ_p=8, ZZ onsite) — こちらが利得を支配する量")
    @printf("%-16s", "config"); for L in L_list; @printf("%11s", "L=$L"); end
    @printf("%16s\n", "L=max/L=min")
    for (U, δ, lab, _) in configs
        @printf("%-16s", lab)
        v = Float64[]
        for L in L_list
            idx = [r for r in rows if r[1]==U && r[2]==δ && r[3]==L && r[4]==8 && r[5]=="ZZ onsite"]
            x = isempty(idx) ? NaN : idx[1][12]; push!(v, x)
            @printf("%11.3e", x)
        end
        @printf("%16.3g\n", v[end]/v[1])
    end

    println("\n中央ボンドエントロピー S_c(L) — gapless なら (c/6)log L で増加、gapped なら飽和")
    @printf("%-16s", "config"); for L in L_list; @printf("%9s", "L=$L"); end
    @printf("%22s\n", "fit c (S=c/6*lnL+a)")
    for (U, δ, lab, _) in configs
        @printf("%-16s", lab)
        xs = Float64[]; ys = Float64[]
        for L in L_list
            idx = [r for r in srows if r[1]==U && r[2]==δ && r[3]==L]
            s = isempty(idx) ? NaN : idx[1][7]
            @printf("%9.4f", s)
            if !isnan(s); push!(xs, log(L)); push!(ys, s); end
        end
        # 最小二乗で傾き → c = 6*slope
        if length(xs) >= 2
            sl = (mean(xs .* ys) - mean(xs)*mean(ys)) / (mean(xs.^2) - mean(xs)^2)
            @printf("%22.3f\n", 6sl)
        else
            println()
        end
    end
    return rows, srows
end

rows, srows = main()
