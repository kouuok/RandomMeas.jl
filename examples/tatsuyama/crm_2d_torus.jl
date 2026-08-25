# ============================================================
# 2D でも周期境界にするとどうなるか — シリンダー vs トーラス
#
# 現行の2D計算(§4, §4b, §8, §9)は「シリンダー」で、円周(y)方向だけが周期、
# 軸(x)方向は開放だった。ここで x 方向も周期にして完全な2次元トーラスにする。
#
# なぜ意味があるか:
#   MPS は列順に1次元に並べた表現なので、切り口はシリンダーだと W 本の結合を
#   横切る(S ~ cW)。トーラスでは巻き付き結合も同時に横切るので S ~ 2cW となり、
#   必要ボンド次元が chi ~ chi_cyl^2 に跳ね上がる。一方 UHF は跳び移り行列を
#   対角化して Wick を使うだけなので、幾何にまったく依存しない。
#   → §8 で見た「幅を増やすと MPS だけが苦しくなる」の、より強い版になるはず。
#
# 1D の §2h と同じ姿勢で、まず DMRG が収束しているかを chi 走査で確かめてから
# 利得を測る。トーラスは不利なので、収束していなければそう明記する。
#
# 実行: julia --project=Hubbard_MPS_Env_v2 crm_2d_torus.jl
# 環境変数: LX(既定 4), W(既定 4), UVAL(既定 8.0), CHI_EXP(既定 1024),
#           CHI_SCAN(既定 "128,256,512,1024"), N_REPEAT(既定 300), DO_SAMPLE(既定 1)
# ============================================================
include(joinpath(@__DIR__, "crm_chain_common.jl"))

const W  = parse(Int, get(ENV, "W", "4"))
const LX = parse(Int, get(ENV, "LX", "4"))
sidx(x, y) = (x - 1) * W + y

"""格子の結合。pbc_x=false ならシリンダー(y のみ周期)、true ならトーラス。"""
function lat_edges(Lx, W; pbc_x::Bool)
    edges = Tuple{Int,Int}[]
    for x in 1:Lx, y in 1:W
        s = sidx(x, y)
        if x < Lx
            push!(edges, (s, sidx(x+1, y)))
        elseif pbc_x && Lx > 2
            push!(edges, (s, sidx(1, y)))          # x 方向の巻き付き
        end
        W > 2 && push!(edges, (s, sidx(x, y == W ? 1 : y + 1)))
    end
    return edges
end

function hubbard_mpo(sites, Lx, W, t, U, mu; pbc_x)
    os = OpSum()
    for (s, s2) in lat_edges(Lx, W; pbc_x), q in (qup, qdn)
        os += -t, "Cdag", q(s), "C", q(s2)
        os += -t, "Cdag", q(s2), "C", q(s)
    end
    for s in 1:Lx*W
        os += U, "N", qup(s), "N", qdn(s)
        os += -mu, "N", qup(s); os += -mu, "N", qdn(s)
    end
    return MPO(os, sites)
end

function ground_state_2d(Lx, W, t, U, mu; chi_max, nsweeps=30, pbc_x)
    n = Lx * W; N = 2n
    sites = siteinds("Fermion", N; conserve_qns=true)
    H = hubbard_mpo(sites, Lx, W, t, U, mu; pbc_x)
    occ = falses(N)
    for x in 1:Lx, y in 1:W
        s = sidx(x, y)
        isodd(x + y) ? (occ[qup(s)] = true) : (occ[qdn(s)] = true)
    end
    ψ0 = productMPS(sites, [occ[q] ? "Occ" : "Emp" for q in 1:N])
    ramp = vcat([64, 128, 256, 512], fill(chi_max, max(0, nsweeps - 4)))
    noise = vcat([1e-5, 1e-6, 1e-7, 1e-8, 1e-9], zeros(max(0, nsweeps - 5)))[1:nsweeps]
    E, ψ = dmrg(H, ψ0; nsweeps, maxdim=min.(chi_max, ramp)[1:nsweeps],
                cutoff=1e-11, noise, outputlevel=0)
    return E, ψ, count(occ)
end

center_S(ψ) = begin
    N = length(ψ); ψo = copy(ψ); orthogonalize!(ψo, N÷2)
    _, Sv, _ = svd(ψo[N÷2], (linkinds(ψo)[N÷2-1], siteinds(ψo)[N÷2]))
    p = diag(Array(Sv, inds(Sv)...)).^2; p = p[p .> 1e-14]
    -sum(p .* log.(p))
end

# ---------------- UHF(幾何に依存しない) ----------------
function hopping_matrix(Lx, W, t; pbc_x)
    n = Lx * W; T = zeros(n, n)
    for (s, s2) in lat_edges(Lx, W; pbc_x); T[s, s2] = T[s2, s] = -t; end
    return T
end

function solve_uhf(Lx, W, t, U; Nup, Ndn, pbc_x, iters=4000, tol=1e-11, init=:neel)
    n = Lx*W; T = hopping_matrix(Lx, W, t; pbc_x); fill_avg = (Nup+Ndn)/(2n)
    nup = zeros(n); ndn = zeros(n)
    for x in 1:Lx, y in 1:W
        s = sidx(x,y); sgn = isodd(x+y) ? 1 : -1
        if init == :neel
            nup[s] = fill_avg + 0.4sgn; ndn[s] = fill_avg - 0.4sgn
        else
            nup[s] = fill_avg + 0.3*(rand()-0.5); ndn[s] = fill_avg + 0.3*(rand()-0.5)
        end
    end
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
    Cu = Φu*Φu'; Cd = Φd*Φd'
    (Cu=Cu, Cd=Cd, E=tr(T*Cu)+tr(T*Cd)+U*sum(diag(Cu).*diag(Cd)),
     m=mean(abs(Cu[s,s]-Cd[s,s])/2 for s in 1:n))
end

function solve_uhf_best(Lx, W, t, U; Nup, Ndn, pbc_x)
    Random.seed!(31); best = nothing
    for ini in (:neel, :rand, :rand, :rand)
        u = solve_uhf(Lx, W, t, U; Nup, Ndn, pbc_x, init=ini)
        (best === nothing || u.E < best.E) && (best = u)
    end
    best
end

full_C(uhf, n) = (C = zeros(2n, 2n);
    for a in 1:n, b in 1:n
        C[qup(a), qup(b)] = uhf.Cu[a,b]; C[qdn(a), qdn(b)] = uhf.Cd[a,b]
    end; C)

function rotate_C(C0, n, θ)
    c, s = cos(θ/2), sin(θ/2); Um = zeros(2n, 2n)
    for i in 1:n
        a, b = qup(i), qdn(i)
        Um[a,a]=c; Um[a,b]=-s; Um[b,a]=s; Um[b,b]=c
    end
    Um*C0*Um'
end

function wick_gen(sup, C)
    isempty(sup) && return 1.0
    ps = [a for (_,a) in sup]
    all(==(3), ps) || error("wick_gen: Z列のみ対応")
    qs = [q for (q,_) in sup]
    length(qs)==1 && return 1 - 2C[qs[1],qs[1]]
    a,b = qs
    1 - 2C[a,a] - 2C[b,b] + 4*(C[a,a]*C[b,b] - C[a,b]*C[b,a])
end

symavg_P(obs, C0, n; nθ=64) = begin
    P = [[0.0 for _ in o.terms] for o in obs]
    for k in 1:nθ
        θ = acos(-1 + (k-0.5)*2/nθ); C = rotate_C(C0, n, θ)
        for (i,o) in enumerate(obs), (j,tm) in enumerate(o.terms)
            P[i][j] += wick_gen(tm.sup, C)/nθ
        end
    end
    P
end

"""サイト s の観測量。bond=true のときの隣は s+1(= 同じ列の y+1)なので、
   列の下端(y=W)では物理的な最近接でなくなる。トーラスでは全サイトが等価な
   はずなので、その判定を汚さないよう y<W のサイトだけボンド量を測る。"""
function site_obs_2d(s::Int; bond::Bool)
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
    obs
end

# ---------------- メイン ----------------
function main()
    t = 1.0
    U = parse(Float64, get(ENV, "UVAL", "8.0"))
    chi_exp   = parse(Int, get(ENV, "CHI_EXP", "1024"))
    chi_scan  = parse.(Int, split(get(ENV, "CHI_SCAN", "128,256,512,1024"), ","))
    n_repeat  = parse(Int, get(ENV, "N_REPEAT", "300"))
    do_sample = get(ENV, "DO_SAMPLE", "1") == "1"
    nu, nm = 50, 100
    chi_priors = [4, 8, 32]
    n = LX * W; mu = U/2

    @printf("=== 2D シリンダー vs トーラス: W=%d Lx=%d (%d量子ビット) U=%.1f ===\n",
            W, LX, 2n, U); flush(stdout)

    # ---- (1) DMRG の収束を chi 走査で確認 ----
    println("\n[1] DMRG 収束チェック")
    @printf("%-10s %8s %18s %10s %10s %8s\n", "幾何", "chi", "E0", "S_center", "到達chi", "秒")
    states = Dict{Bool,Any}()
    for pbc_x in (false, true)
        geo = pbc_x ? "トーラス" : "シリンダー"
        Es = Float64[]; Ss = Float64[]; Ts = Float64[]; Ms = Int[]
        for cx in chi_scan
            t0 = time()
            E, ψ, nel = ground_state_2d(LX, W, t, U, mu; chi_max=cx, pbc_x)
            push!(Es, E); push!(Ss, center_S(ψ)); push!(Ts, time()-t0); push!(Ms, maxlinkdim(ψ))
            cx == chi_scan[end] && (states[pbc_x] = (E=E, ψ=ψ, nel=nel))
            @printf("%-10s %8d %18.8f %10.4f %10d %8.0f\n", geo, cx, E, Ss[end], Ms[end], Ts[end])
            flush(stdout)
        end
        # 最大 chi の値を基準にした収束(chi を上げてもエネルギーが動かなくなったか)
        Eref = Es[end]
        @printf("   → %s の収束: %s\n", geo,
                join([@sprintf("chi=%d:%.2e", c, abs(e-Eref)) for (c,e) in zip(chi_scan, Es)], "  "))
        @printf("      chi が上限に張り付いたか: %s\n",
                Ms[end] >= chi_scan[end] ? "**はい(未収束の可能性)**" : "いいえ(切断が厳密)")
        flush(stdout)
    end

    # ---- (2) 利得の測定 ----
    rows = []
    for pbc_x in (false, true)
        geo = pbc_x ? "torus" : "cylinder"
        E, ψ, nel = states[pbc_x].E, states[pbc_x].ψ, states[pbc_x].nel
        Sc = center_S(ψ)
        @printf("\n[2] %s: E0=%.6f  N_el=%d  chi=%d  S_center=%.4f\n",
                geo, E, nel, maxlinkdim(ψ), Sc); flush(stdout)

        Nup = (nel+1) ÷ 2; Ndn = nel ÷ 2
        uhf = solve_uhf_best(LX, W, t, U; Nup, Ndn, pbc_x)
        C0  = full_C(uhf, n)
        @printf("   UHF: E=%.6f  スタッガード磁化 m=%.4f\n", uhf.E, uhf.m); flush(stdout)

        sites_meas = [s for s in 1:n if (s-1) % W < W-1]     # y<W のサイトのみ
        obs_by_site = [site_obs_2d(s; bond=true) for s in sites_meas]
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
        @printf("   prior忠実度: %s\n",
                join([@sprintf("chi=%d:%.3e", c, f) for (c,f) in zip(chi_priors, fids)], "  "))
        flush(stdout)

        for (si, s) in enumerate(sites_meas)
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
            Pu = [[wick_gen(tm.sup, C0) for tm in o.terms] for o in obs]
            Ps = symavg_P(obs, C0, n)
            for P in (Pu, Ps)
                push!(Pσ_all, P)
                push!(trOσ_all, [sum(tm.coeff*P[k][ti] for (ti,tm) in enumerate(obs[k].terms))
                                 for k in 1:length(obs)])
            end
            Gs = fill(NaN, length(obs), length(labels)); Gmx = fill(NaN, length(obs))
            if do_sample
                es, ec = run_locals(rdm_ρ[si], Nw, ow, Pσ_all, trOσ_all;
                                    nu, nm, n_repeat, seed = 770_000 + 100s + (pbc_x ? 1 : 0))
                for k in 1:length(obs)
                    v = var(es[:,k]); Gmx[k] = v/var(ec[:,k,length(chi_priors)+1])
                    for p in 1:length(labels); Gs[k,p] = v/var(ec[:,k,p]); end
                end
            end
            x = (s-1) ÷ W + 1; y = (s-1) % W + 1
            for (k,o) in enumerate(obs), (p,lb) in enumerate(labels)
                Δ = Otrue[k] - trOσ_all[p][k]
                push!(rows, (geo, U, s, x, y, o.name, lb, Otrue[k], Δ,
                             abs(Otrue[k])>1e-3 ? abs(Δ/Otrue[k]) : NaN,
                             Gs[k,p], Gmx[k], fids[min(p,length(fids))], Sc))
            end
        end

        # サイト依存性(トーラスなら並進対称で消えるはず)
        @printf("\n   %s の <P> のサイト依存性(全 %d サイト)\n", geo, length(sites_meas))
        @printf("   %-14s %12s %12s %12s\n", "観測量", "最小", "最大", "変動係数")
        for o in unique(r[6] for r in rows if r[1]==geo)
            v = [r[8] for r in rows if r[1]==geo && r[6]==o && r[7]=="chi4"]
            isempty(v) && continue
            cv = abs(mean(v)) > 1e-9 ? std(v)/abs(mean(v)) : NaN
            @printf("   %-14s %12.6f %12.6f %12.2e\n", o, minimum(v), maximum(v), cv)
        end
        flush(stdout)
    end

    out = joinpath(@__DIR__, @sprintf("crm_2d_torus_U%s.tsv", replace(@sprintf("%.1f",U), "."=>"p")))
    open(out, "w") do io
        println(io, "geometry\tU\tsite\tx\ty\tobservable\tprior\ttrue\tDelta\trelerr\tG\tG_max\tprior_fid\tS_center")
        for r in rows; println(io, join(r,"\t")); end
    end
    println("\nresults saved: $out")
end
main()
