# ============================================================
# §2g の ED 参照状態チェックを周期境界条件(PBC)で行う
#
# 動機:
#   §2g は OBC で「ρ を厳密対角化の基底状態に置き換えても利得は 1 点も
#   変わらない」ことを示した。しかし OBC では DMRG(χ=256)が事実上厳密
#   (|ΔE|~1e-9)だったので、「ρ が切断MPSであること」の影響を見るには
#   そもそも切断が効いていなかった、という弱点がある。
#
#   PBC はエンタングルメントが約2倍 (S=(c/3)lnL vs (c/6)lnL) で、
#   さらに JW 順序では巻き付きホッピングが全鎖にまたがる文字列になるため
#   MPS にとって著しく不利である。つまり「ρ が切断MPSであること」が
#   結論に効く余地が OBC より大きい。ここが本節の実験的な狙い。
#
#   ED 側は境界条件に一切影響を受けない（ボンド次元という概念がない）。
#   巻き付き項は JW 文字列が長くなるだけで、符号を正しく入れれば済む。
#
# 巻き付きホッピングの JW 符号 (qup(i)=2i-1, qdn(i)=2i):
#   up  1↔L: モード qup(1)=1 と qup(L)=2L-1 の間 = qup(2..L-1), qdn(1..L-1)
#   dn  1↔L: モード qdn(1)=2 と qdn(L)=2L   の間 = qup(2..L),   qdn(2..L-1)
#   (符号は「厳密に間にあるモード」の占有パリティ。両端は変化するが
#    間の集合には入らないので、始状態で評価してよい。)
#
# 妥当性検証:
#   小さい L で ED(PBC) の基底エネルギーを DMRG(PBC, ITensor の OpSum が
#   JW を内部処理する完全に独立な経路)と突き合わせる。巻き付き符号を
#   間違えていればここで必ず落ちる。
#
# 実行: julia --project=Hubbard_MPS_Env_v2 crm_ed_pbc.jl
# 環境変数: L_LIST(既定 "8,12,14"), BC_LIST(既定 "obc,pbc"),
#           CHI_EXP(既定 1024), N_REPEAT(既定 500)
# ============================================================
include(joinpath(@__DIR__, "crm_chain_common.jl"))

function combos(L, k)
    out = Int[]
    for x in 0:(1<<L - 1); count_ones(x) == k && push!(out, x); end
    return out
end

struct Sector
    L::Int; ups::Vector{Int}; dns::Vector{Int}
    upidx::Dict{Int,Int}; dnidx::Dict{Int,Int}
end
function Sector(L, Nup, Ndn)
    ups = combos(L, Nup); dns = combos(L, Ndn)
    Sector(L, ups, dns, Dict(a=>i for (i,a) in enumerate(ups)),
           Dict(b=>i for (i,b) in enumerate(dns)))
end
Base.length(s::Sector) = length(s.ups)*length(s.dns)
bit(x, i) = (x >> (i-1)) & 1
maskrange(lo, hi) = lo > hi ? 0 : ((1 << (hi-lo+1)) - 1) << (lo-1)

"""ホッピング先を前計算しておく(内側ループから Dict 参照を除く)。

近接項は crm_ed_reference.jl と同一。pbc=true のときだけ巻き付き項の
テーブルを作る。巻き付きの符号は up マスクと dn マスクの両方に依存する
ので、それぞれのパリティを別々に持ち、内側で XOR して合成する。"""
struct HopTables
    up::Vector{Vector{Tuple{Int,Int}}}   # ia -> [(bond i, 行き先 ia2)]
    dn::Vector{Vector{Tuple{Int,Int}}}
    pbc::Bool
    wup::Vector{Tuple{Int,Int,Int}}      # (ia, ia2, パリティ from a)
    wdn::Vector{Tuple{Int,Int,Int}}      # (ib, ib2, パリティ from b)
    pb_all::Vector{Int}                  # up 巻き付きで使う b のパリティ
    pa_2L::Vector{Int}                   # dn 巻き付きで使う a のパリティ
end
function HopTables(S::Sector; pbc::Bool)
    L = S.L
    mk(cfgs, idx) = [begin
        h = Tuple{Int,Int}[]
        for i in 1:L-1
            bit(x,i) == bit(x,i+1) && continue
            push!(h, (i, idx[x ⊻ ((1<<(i-1)) | (1<<i))]))
        end
        h
    end for x in cfgs]
    up = mk(S.ups, S.upidx); dn = mk(S.dns, S.dnidx)
    if !pbc
        return HopTables(up, dn, false, Tuple{Int,Int,Int}[], Tuple{Int,Int,Int}[], Int[], Int[])
    end
    flip = (1 << 0) | (1 << (L-1))
    m_up_between = maskrange(2, L-1)   # up 巻き付きの「間」にある up サイト
    m_dn_all     = maskrange(1, L-1)   # 同上の dn サイト
    m_up_2L      = maskrange(2, L)     # dn 巻き付きの「間」にある up サイト
    m_dn_between = maskrange(2, L-1)   # 同上の dn サイト
    wup = Tuple{Int,Int,Int}[]
    for (ia, a) in enumerate(S.ups)
        bit(a,1) == bit(a,L) && continue
        push!(wup, (ia, S.upidx[a ⊻ flip], count_ones(a & m_up_between) & 1))
    end
    wdn = Tuple{Int,Int,Int}[]
    for (ib, b) in enumerate(S.dns)
        bit(b,1) == bit(b,L) && continue
        push!(wdn, (ib, S.dnidx[b ⊻ flip], count_ones(b & m_dn_between) & 1))
    end
    pb_all = [count_ones(b & m_dn_all) & 1 for b in S.dns]
    pa_2L  = [count_ones(a & m_up_2L)  & 1 for a in S.ups]
    return HopTables(up, dn, true, wup, wdn, pb_all, pa_2L)
end

function apply_H!(w, v, S::Sector, HT::HopTables, t, U, mu)
    nu = length(S.ups); nd = length(S.dns)
    fill!(w, 0)
    @inbounds for ia in 1:nu
        a = S.ups[ia]; base = (ia-1)*nd
        na = count_ones(a)
        for ib in 1:nd                                   # 対角
            b = S.dns[ib]
            w[base+ib] += (U*count_ones(a & b) - mu*(na + count_ones(b))) * v[base+ib]
        end
        for (i, ia2) in HT.up[ia]                        # up 近接
            tb = (ia2-1)*nd
            for ib in 1:nd
                sgn = iseven(bit(S.dns[ib], i)) ? 1.0 : -1.0
                w[tb+ib] += -t*sgn*v[base+ib]
            end
        end
        for ib in 1:nd                                   # dn 近接
            amp = v[base+ib]
            amp == 0 && continue
            for (i, ib2) in HT.dn[ib]
                sgn = iseven(bit(a, i+1)) ? 1.0 : -1.0
                w[base+ib2] += -t*sgn*amp
            end
        end
    end
    HT.pbc || return w
    @inbounds for (ia, ia2, pa) in HT.wup                # up 巻き付き 1↔L
        base = (ia-1)*nd; tb = (ia2-1)*nd
        for ib in 1:nd
            sgn = iseven(pa ⊻ HT.pb_all[ib]) ? 1.0 : -1.0
            w[tb+ib] += -t*sgn*v[base+ib]
        end
    end
    @inbounds for ia in 1:nu                             # dn 巻き付き 1↔L
        base = (ia-1)*nd; pa = HT.pa_2L[ia]
        for (ib, ib2, pb) in HT.wdn
            amp = v[base+ib]
            amp == 0 && continue
            sgn = iseven(pa ⊻ pb) ? 1.0 : -1.0
            w[base+ib2] += -t*sgn*amp
        end
    end
    return w
end

"""Lanczos(完全再直交化・thick restart)。第2固有値も返す(縮退の検知用)。"""
function lanczos_gs(S::Sector, t, U, mu; pbc, m=50, maxrestart=80, tol=1e-11, seed=1234)
    D = length(S); HT = HopTables(S; pbc); Random.seed!(seed)
    v = randn(D); v ./= norm(v); Eprev = Inf; w = similar(v); E1 = NaN
    for restart in 1:maxrestart
        Q = [copy(v)]; α = Float64[]; β = Float64[]
        for j in 1:m
            apply_H!(w, Q[j], S, HT, t, U, mu)
            a = dot(Q[j], w); push!(α, a)
            w .-= a .* Q[j]
            j > 1 && (w .-= β[j-1] .* Q[j-1])
            for q in Q; w .-= dot(q, w) .* q; end
            b = norm(w)
            if b < 1e-12 || j == m; push!(β, b); break; end
            push!(β, b); push!(Q, w ./ b)
        end
        k = length(α)
        F = eigen(SymTridiagonal(α, β[1:k-1]))
        E0 = F.values[1]; c = F.vectors[:, 1]
        k > 1 && (E1 = F.values[2])
        v = zeros(D); for j in 1:k; v .+= c[j] .* Q[j]; end
        v ./= norm(v)
        abs(E0 - Eprev) < tol && return E0, E1, v, restart
        Eprev = E0
    end
    return Eprev, E1, v, maxrestart
end

function expand_full(v, S::Sector)
    L = S.L; nd = length(S.dns)
    full = zeros(ComplexF64, 1 << (2L))
    @inbounds for (ia, a) in enumerate(S.ups), (ib, b) in enumerate(S.dns)
        q = 0
        for i in 1:L
            q |= bit(a,i) << (2i-2)
            q |= bit(b,i) << (2i-1)
        end
        full[q+1] = v[(ia-1)*nd + ib]
    end
    return full
end

"""MPS を全ベクトルに展開する(expand_full と同じ添字規約: 量子ビット n が bit n-1)。
中間配列は 2^k × χ_k ≤ 2^N に抑えられるので、L≤12 なら数百MBで済む。"""
function mps_full(ψ::MPS)
    tens = extract_tensors(ψ)
    Ψ = ComplexF64[1.0;;]
    for (A0, A1) in tens
        Ψ = vcat(Ψ*ComplexF64.(A0), Ψ*ComplexF64.(A1))
    end
    return vec(Ψ)
end

"""全ベクトルから窓 [q1,q2] の縮約密度行列。窓内は b_{q1} が最上位になるよう並べ替える。"""
function window_rdm_full(full, N, q1, q2)
    Nw = q2 - q1 + 1
    dl = 1 << (q1-1); dw = 1 << Nw; dr = 1 << (N - q2)
    A = reshape(full, dl, dw, dr)
    ρ = zeros(ComplexF64, dw, dw)
    @inbounds for r in 1:dr, s2 in 1:dw, s1 in 1:dw
        acc = zero(ComplexF64)
        for l in 1:dl; acc += A[l, s1, r] * conj(A[l, s2, r]); end
        ρ[s1, s2] += acc
    end
    perm = [begin
                x = s - 1; y = 0
                for k in 0:Nw-1; y |= ((x >> k) & 1) << (Nw-1-k); end
                y + 1
            end for s in 1:dw]
    R = ρ[perm, perm]
    return Hermitian((R + R')/2)
end

"""窓RDM が ED と DMRG でどれだけ違うかを、実際にサンプリングに使う量で測る。

「利得が厳密に一致した」理由を推測でなく実測にするためのもの。サンプラーは
basis_probs(ρ_w, basis) の累積分布からビット列を引くので、効くのは ρ の
行列要素そのものではなく **各基底での測定確率** である。3^4=81 通りの基底
すべてについて、確率の最大差と全変動距離 (TV) を返す。

1 回の抽選でビットが変わる確率は高々 TV なので、総抽選数 × max TV が
「ビットが 1 つでも反転する回数」の期待値の上限になる。"""
function prob_gap(ρa, ρb, Nw)
    maxdp = 0.0; maxtv = 0.0
    for idx in 0:(3^Nw - 1)
        basis = Vector{Int}(undef, Nw); x = idx
        for j in 1:Nw; basis[j] = x % 3 + 1; x ÷= 3; end
        pa = basis_probs(ρa, basis); pb = basis_probs(ρb, basis)
        d = abs.(pa .- pb)
        maxdp = max(maxdp, maximum(d)); maxtv = max(maxtv, 0.5*sum(d))
    end
    return maxdp, maxtv
end

# ---------------- DMRG (境界条件つき) ----------------
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
    init, _ = initial_config(L, 0.0)
    ramp = vcat([50,100,200,400], fill(chi_max, max(0, nsweeps-4)))
    noise = vcat([1e-5,1e-6,1e-7,1e-8,1e-9], zeros(max(0, nsweeps-5)))[1:nsweeps]
    E, ψ = dmrg(H, productMPS(sites, init); nsweeps,
                maxdim=min.(chi_max, ramp)[1:nsweeps], cutoff=1e-12, noise, outputlevel=0)
    return E, ψ
end

function center_S(ψ)
    N=length(ψ); ψo=copy(ψ); orthogonalize!(ψo, N÷2)
    _,Sv,_ = svd(ψo[N÷2], (linkinds(ψo)[N÷2-1], siteinds(ψo)[N÷2]))
    p=diag(Array(Sv,inds(Sv)...)).^2; p=p[p.>1e-14]
    return -sum(p.*log.(p))
end

# ---------------- 1D鎖 UHF (境界条件つき) ----------------
function chain_hop_matrix(L, t; pbc::Bool)
    T = zeros(L, L)
    for i in 1:L-1; T[i,i+1] = T[i+1,i] = -t; end
    pbc && (T[1,L] = T[L,1] = -t)
    return T
end

function solve_uhf_chain(L, t, U; Nup, Ndn, pbc, iters=4000, tol=1e-12)
    T = chain_hop_matrix(L, t; pbc); fill_avg = (Nup+Ndn)/(2L)
    nup = zeros(L); ndn = zeros(L)
    for i in 1:L
        nup[i] = fill_avg + 0.4*(-1)^i; ndn[i] = fill_avg - 0.4*(-1)^i
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
    return (Cu=Cu, Cd=Cd, E=tr(T*Cu)+tr(T*Cd)+U*sum(diag(Cu).*diag(Cd)))
end

full_C(uhf, L) = (C = zeros(2L, 2L);
    for a in 1:L, b in 1:L
        C[qup(a), qup(b)] = uhf.Cu[a,b]; C[qdn(a), qdn(b)] = uhf.Cd[a,b]
    end; C)

function rotate_C(C0, L, θ)
    c, s = cos(θ/2), sin(θ/2); Um = zeros(2L, 2L)
    for i in 1:L
        a, b = qup(i), qdn(i)
        Um[a,a]=c; Um[a,b]=-s; Um[b,a]=s; Um[b,b]=c
    end
    return Um*C0*Um'
end

function wick_gen(sup, C)
    isempty(sup) && return 1.0
    ps = [a for (_,a) in sup]
    if all(==(3), ps)
        qs = [q for (q,_) in sup]
        length(qs)==1 && return 1 - 2C[qs[1],qs[1]]
        a,b = qs
        return 1 - 2C[a,a] - 2C[b,b] + 4*(C[a,a]*C[b,b] - C[a,b]*C[b,a])
    elseif length(sup)==3 && (ps==[1,3,1] || ps==[2,3,2])
        return 2*C[sup[1][1], sup[3][1]]
    end
    error("wick_gen: 未対応 $(sup)")
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

# ---------------- 巻き付き符号の妥当性検証 ----------------
"""小さい L で ED(PBC) と DMRG(PBC) の基底エネルギーを突き合わせる。
DMRG 側は ITensor の OpSum が JW 変換を内部で行う独立実装なので、
巻き付き符号を間違えていればここで必ず食い違う。"""
function wrap_sign_check(t, U)
    for L in (4, 6)
        S = Sector(L, L÷2, L - L÷2)
        for pbc in (false, true)
            E_ed, _, _, _ = lanczos_gs(S, t, U, U/2; pbc)
            E_dm, _ = gs_bc(L, t, U, U/2; chi_max=256, nsweeps=20, pbc)
            @printf("  [check] L=%d %s: E_ED=%.10f  E_DMRG=%.10f  diff=%.2e\n",
                    L, pbc ? "PBC" : "OBC", E_ed, E_dm, abs(E_ed-E_dm))
            @assert abs(E_ed - E_dm) < 1e-7
        end
    end
    println("  [check] 巻き付き JW 符号 OK"); flush(stdout)
end

# ---------------- メイン ----------------
function main()
    t, U, nu, nm = 1.0, 4.0, 50, 100
    L_list  = parse.(Int, split(get(ENV, "L_LIST",  "8,12,14"), ","))
    bc_list = [strip(s) == "pbc" for s in split(get(ENV, "BC_LIST", "obc,pbc"), ",")]
    n_repeat = parse(Int, get(ENV, "N_REPEAT", "500"))
    seeds    = parse.(Int, split(get(ENV, "SEED_LIST", "0,1,2"), ","))
    chi_exp  = parse(Int, get(ENV, "CHI_EXP", "1024"))
    fid_lmax = parse(Int, get(ENV, "FID_LMAX", "12"))
    chi_priors = [2, 4, 8, 16, 32]

    println("=== validation (L=4, OBC 既存系) ==="); flush(stdout)
    dense_chain_check()
    println("\n=== validation (巻き付き符号: ED(PBC) vs DMRG(PBC)) ==="); flush(stdout)
    wrap_sign_check(t, U)

    rows = []; meta = []
    for L in L_list, pbc in bc_list
        bclab = pbc ? "PBC" : "OBC"
        N = 2L; Nup = L ÷ 2; Ndn = L - Nup
        S = Sector(L, Nup, Ndn)
        @printf("\n%s\nL=%d %s (N=%d qubits)  セクター次元=%d  全ベクトル=%.2f GB\n",
                "="^76, L, bclab, N, length(S), (1<<N)*16/1e9); flush(stdout)

        t0 = time()
        E_ed, E1, vsec, nres = lanczos_gs(S, t, U, U/2; pbc)
        @printf("  ED(Lanczos): E0=%.10f  E1-E0=%.3e  restart=%d  (%.0fs)\n",
                E_ed, E1-E_ed, nres, time()-t0); flush(stdout)
        full = expand_full(vsec, S); vsec = nothing; GC.gc()

        t1 = time()
        E_dm, ψ = gs_bc(L, t, U, U/2; chi_max=chi_exp, nsweeps=40, pbc)
        Sc = center_S(ψ)
        @printf("  DMRG(χ≤%d): E0=%.10f  maxlinkdim=%d  |ΔE|=%.2e  S_center=%.4f  (%.0fs)\n",
                chi_exp, E_dm, maxlinkdim(ψ), abs(E_dm-E_ed), Sc, time()-t1); flush(stdout)

        fid = NaN
        if L <= fid_lmax
            fid = abs2(dot(full, mps_full(ψ))); GC.gc()
            @printf("  |⟨ED|DMRG⟩|² = %.12f  (1-F = %.2e)\n", fid, 1-fid); flush(stdout)
        end

        uhf = solve_uhf_chain(L, t, U; Nup, Ndn, pbc); C0 = full_C(uhf, L)

        sites = 1:L-1
        obs_by_site = [site_observables(i; bond=true) for i in sites]
        windows = [obs_support(o) for o in obs_by_site]
        rdm_dm = window_rdms(ψ, windows)
        priors = MPS[]
        for cp in chi_priors
            σ = truncate(ψ; maxdim=cp); normalize!(σ); push!(priors, σ)
        end
        push!(priors, ψ)
        fids = [abs2(inner(σ, ψ)) for σ in priors]
        @printf("  prior忠実度(vs DMRG): %s\n", join([@sprintf("χ=%d:%.3e", c, f)
                for (c,f) in zip(chi_priors, fids)], "  ")); flush(stdout)
        rdm_p = [window_rdms(σ, windows) for σ in priors]
        labels = vcat(["chi$c" for c in chi_priors], ["exact", "UHF", "UHFsym"])

        maxdiff = 0.0; maxdp = 0.0; maxtv = 0.0
        for (si, i) in enumerate(sites)
            obs = obs_by_site[si]; q1, q2 = windows[si]; Nw = q2-q1+1
            ow = shift_obs(obs, q1-1)
            Mats = [[term_window_matrix(tm, Nw) for tm in o.terms] for o in ow]
            ρ_ed = window_rdm_full(full, N, q1, q2)
            for (name, ρref) in (("ED", ρ_ed), ("DMRG", rdm_dm[si]))
                Otrue = [sum(tm.coeff*expect_rdm(ρref, Mats[k][ti])
                             for (ti,tm) in enumerate(ow[k].terms)) for k in 1:length(obs)]
                if name == "ED"
                    Odm = [sum(tm.coeff*expect_rdm(rdm_dm[si], Mats[k][ti])
                               for (ti,tm) in enumerate(ow[k].terms)) for k in 1:length(obs)]
                    maxdiff = max(maxdiff, maximum(abs.(Otrue .- Odm)))
                    dp, tv = prob_gap(ρ_ed, rdm_dm[si], Nw)
                    maxdp = max(maxdp, dp); maxtv = max(maxtv, tv)
                end
                Pσ_all = Vector{Vector{Vector{Float64}}}(); trOσ_all = Vector{Vector{Float64}}()
                for p in 1:length(priors)
                    Pσ = [[expect_rdm(rdm_p[p][si], M) for M in Ms] for Ms in Mats]
                    push!(Pσ_all, Pσ)
                    push!(trOσ_all, [sum(tm.coeff*Pσ[k][ti] for (ti,tm) in enumerate(obs[k].terms))
                                     for k in 1:length(obs)])
                end
                Pu = [[wick_gen(tm.sup, C0) for tm in o.terms] for o in obs]
                Ps = symavg(obs, C0, L)
                for P in (Pu, Ps)
                    push!(Pσ_all, P)
                    push!(trOσ_all, [sum(tm.coeff*P[k][ti] for (ti,tm) in enumerate(obs[k].terms))
                                     for k in 1:length(obs)])
                end
                for sd in seeds
                    es, ec = run_locals(ρref, Nw, ow, Pσ_all, trOσ_all;
                                        nu, nm, n_repeat, seed = 9100 + 100L + i + 7919*sd)
                    for (k,o) in enumerate(obs)
                        v = var(es[:,k]); gmx = v/var(ec[:,k,length(chi_priors)+1])
                        for (p, lb) in enumerate(labels)
                            lo, hi = boot_ratio(es[:,k], ec[:,k,p])
                            push!(rows, (L, bclab, sd, name, i, lb, o.name, Otrue[k],
                                         Otrue[k]-trOσ_all[p][k], v/var(ec[:,k,p]), lo, hi, gmx))
                        end
                    end
                end
            end
        end
        full = nothing; GC.gc()
        ndraw = length(sites) * n_repeat * nu * nm
        @printf("  局所観測量の ED vs DMRG 最大差: %.3e\n", maxdiff)
        @printf("  測定確率の ED vs DMRG 差: max|Δp|=%.3e  max TV=%.3e\n", maxdp, maxtv)
        @printf("  → 総抽選数 %.2e に対し、ビットが反転する期待回数の上限 = %.3e\n",
                Float64(ndraw), ndraw*maxtv); flush(stdout)
        push!(meta, (L, bclab, E_ed, E_dm, abs(E_dm-E_ed), maxlinkdim(ψ), Sc, fid,
                     maxdiff, maxdp, maxtv, ndraw, ndraw*maxtv))

        @printf("\n  %-14s %-6s %8s %8s %8s %8s %8s | %8s\n",
                "observable","ρ","χ=2","χ=4","χ=8","χ=32","UHFsym","天井")
        for o in unique(r[7] for r in rows if r[1]==L && r[2]==bclab)
            for nmr in ("ED","DMRG")
                f(lb)=(h=[r[10] for r in rows if r[1]==L&&r[2]==bclab&&r[4]==nmr&&r[6]==lb&&r[7]==o];
                       isempty(h) ? NaN : median(h))
                gm=median([r[13] for r in rows if r[1]==L&&r[2]==bclab&&r[4]==nmr&&r[7]==o])
                @printf("  %-14s %-6s %8.2f %8.2f %8.2f %8.2f %8.2f | %8.2f\n",
                        nmr=="ED" ? o : "", nmr, f("chi2"),f("chi4"),f("chi8"),
                        f("chi32"),f("UHFsym"),gm)
            end
        end
        flush(stdout)
    end

    # ---- ED と DMRG の一致度サマリ ----
    println("\n", "="^76)
    @printf("%-4s %-5s %10s %10s %10s %12s\n", "L","bc","点数","max|ΔG|","max|ΔΔ|","G が完全一致")
    for (L, bclab) in unique((r[1], r[2]) for r in rows)
        gd = Float64[]; dd = Float64[]
        for r in rows
            (r[1]==L && r[2]==bclab && r[4]=="ED") || continue
            m = findfirst(s -> s[1]==L&&s[2]==bclab&&s[3]==r[3]&&s[4]=="DMRG"&&
                               s[5]==r[5]&&s[6]==r[6]&&s[7]==r[7], rows)
            m === nothing && continue
            push!(gd, abs(r[10]-rows[m][10])); push!(dd, abs(r[9]-rows[m][9]))
        end
        @printf("%-4d %-5s %10d %10.3e %10.3e %12s\n", L, bclab, length(gd),
                maximum(gd), maximum(dd), all(iszero, gd) ? "はい" : "いいえ")
    end

    # 損をする点の数 (プロジェクト共通の基準: ブートストラップCI上端 < 0.9)
    println("\n損をする点の数 (CI上端<0.9、ρ=ED、seed=$(seeds[1]))")
    plist = vcat(["chi$c" for c in chi_priors], ["UHF","UHFsym"])
    @printf("%-4s %-5s %8s", "L","bc","総点数")
    for p in plist; @printf("%10s", p); end; println()
    for (L, bclab) in unique((r[1], r[2]) for r in rows)
        tot = count(r -> r[1]==L&&r[2]==bclab&&r[3]==seeds[1]&&r[4]=="ED"&&r[6]=="chi2", rows)
        @printf("%-4d %-5s %8d", L, bclab, tot)
        for p in plist
            nb = count(r -> r[1]==L&&r[2]==bclab&&r[3]==seeds[1]&&r[4]=="ED"&&r[6]==p&&r[12]<0.9, rows)
            @printf("%10s", "$nb/$tot")
        end
        println()
    end

    out = joinpath(@__DIR__, "crm_ed_pbc_results.tsv")
    open(out,"w") do io
        println(io, "L\tbc\tseed\trho\tsite\tprior\tobservable\ttrue\tDelta\tG\tG_lo\tG_hi\tG_max")
        for r in rows; println(io, join(r,"\t")); end
    end
    out2 = joinpath(@__DIR__, "crm_ed_pbc_meta.tsv")
    open(out2,"w") do io
        println(io, "L\tbc\tE_ed\tE_dmrg\tdE\tmaxlinkdim\tS_center\tfid\tmax_obs_diff\tmax_dp\tmax_tv\tndraw\texp_flips")
        for r in meta; println(io, join(r,"\t")); end
    end
    println("\nresults saved: $out , $out2")
end
main()
