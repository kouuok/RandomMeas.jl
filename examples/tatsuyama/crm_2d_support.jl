# ============================================================
# 2D シリンダー vs トーラス — 参照状態を厳密対角化にして決着させる
#
# 経緯:
#   crm_2d_torus.jl で 4×4 を回したが結論が出なかった。chi=256 でも
#   両幾何とも未収束(chi を倍にして誤差が1/3にしかならない)で、
#   S_center が切断で頭打ちになり比較が成立しなかった。さらに悪いことに、
#   トーラスでは ρ 自身が質の悪い MPS なのに prior はその切断なので、
#   **MPS prior に構造的に有利な比較**になっている(§2g で1Dについて
#   潰した疑いが2Dで復活する)。
#
#   4×3 なら 12 サイト・24 量子ビット、半充填セクター次元 C(12,6)^2 = 853,776 で
#   厳密対角化が確実に届く(§2g の 1D L=12 と同規模)。収束の疑いを完全に
#   消したうえで幾何を比べる。
#
# JW 符号(任意のボンド s1 < s2 に一般化):
#   qup(s)=2s-1, qdn(s)=2s。符号は「厳密に間にあるモード」の占有パリティ。
#     up  hop: 間 = qdn(s1..s2-1) と qup(s1+1..s2-1)
#     dn  hop: 間 = qup(s1+1..s2) と qdn(s1+1..s2-1)
#   1D 最近接・巻き付きの両方がこの特別な場合になっている(crm_ed_pbc.jl と一致)。
#
# 実行: julia --project=Hubbard_MPS_Env_v2 crm_2d_torus_ed.jl
# 環境変数: LX(既定 3), W(既定 4), UVAL(既定 8.0), CHI_EXP(既定 2048),
#           N_REPEAT(既定 300), NU(既定 50)
# ============================================================
include(joinpath(@__DIR__, "crm_chain_common.jl"))
include(joinpath(@__DIR__, "crm_wick_pauli.jl"))

const W  = parse(Int, get(ENV, "W", "4"))
const LX = parse(Int, get(ENV, "LX", "3"))
sidx(x, y) = (x - 1) * W + y

function lat_edges(Lx, W; pbc_x::Bool)
    edges = Tuple{Int,Int}[]
    for x in 1:Lx, y in 1:W
        s = sidx(x, y)
        if x < Lx
            push!(edges, minmax(s, sidx(x+1, y)))
        elseif pbc_x && Lx > 2
            push!(edges, minmax(s, sidx(1, y)))
        end
        if W > 2
            push!(edges, minmax(s, sidx(x, y == W ? 1 : y + 1)))
        elseif W == 2 && y == 1
            push!(edges, minmax(s, sidx(x, 2)))      # ラダーの桟は1本だけ
        end
    end
    return unique(edges)
end

# ---------------- 厳密対角化 ----------------
combos(L, k) = [x for x in 0:(1<<L - 1) if count_ones(x) == k]
bit(x, i) = (x >> (i-1)) & 1
maskrange(lo, hi) = lo > hi ? 0 : ((1 << (hi-lo+1)) - 1) << (lo-1)

struct Sector
    n::Int; ups::Vector{Int}; dns::Vector{Int}
    upidx::Dict{Int,Int}; dnidx::Dict{Int,Int}
end
function Sector(n, Nup, Ndn)
    ups = combos(n, Nup); dns = combos(n, Ndn)
    Sector(n, ups, dns, Dict(a=>i for (i,a) in enumerate(ups)),
           Dict(b=>i for (i,b) in enumerate(dns)))
end
Base.length(s::Sector) = length(s.ups)*length(s.dns)

"""1つのボンドについて、上下スピンそれぞれの飛び先とパリティを前計算する。"""
struct BondTab
    up::Vector{Tuple{Int,Int,Int}}   # (ia, ia2, a 由来のパリティ)
    dn::Vector{Tuple{Int,Int,Int}}   # (ib, ib2, b 由来のパリティ)
    pb_up::Vector{Int}               # up hop で使う b 由来のパリティ
    pa_dn::Vector{Int}               # dn hop で使う a 由来のパリティ
end
function BondTab(S::Sector, s1::Int, s2::Int)
    flip = (1 << (s1-1)) | (1 << (s2-1))
    m_up_u = maskrange(s1+1, s2-1); m_dn_u = maskrange(s1, s2-1)   # up hop の「間」
    m_up_d = maskrange(s1+1, s2);   m_dn_d = maskrange(s1+1, s2-1) # dn hop の「間」
    up = Tuple{Int,Int,Int}[]
    for (ia, a) in enumerate(S.ups)
        bit(a,s1) == bit(a,s2) && continue
        push!(up, (ia, S.upidx[a ⊻ flip], count_ones(a & m_up_u) & 1))
    end
    dn = Tuple{Int,Int,Int}[]
    for (ib, b) in enumerate(S.dns)
        bit(b,s1) == bit(b,s2) && continue
        push!(dn, (ib, S.dnidx[b ⊻ flip], count_ones(b & m_dn_d) & 1))
    end
    BondTab(up, dn,
            [count_ones(b & m_dn_u) & 1 for b in S.dns],
            [count_ones(a & m_up_d) & 1 for a in S.ups])
end

function apply_H!(w, v, S::Sector, bonds, tabs, t, U, mu)
    nu = length(S.ups); nd = length(S.dns)
    fill!(w, 0)
    @inbounds for ia in 1:nu                       # 対角項
        a = S.ups[ia]; base = (ia-1)*nd; na = count_ones(a)
        for ib in 1:nd
            b = S.dns[ib]
            w[base+ib] += (U*count_ones(a & b) - mu*(na + count_ones(b))) * v[base+ib]
        end
    end
    @inbounds for tb in tabs
        for (ia, ia2, pa) in tb.up                 # up ホッピング
            base = (ia-1)*nd; tgt = (ia2-1)*nd
            for ib in 1:nd
                sgn = iseven(pa ⊻ tb.pb_up[ib]) ? 1.0 : -1.0
                w[tgt+ib] += -t*sgn*v[base+ib]
            end
        end
        for ia in 1:nu                             # dn ホッピング
            base = (ia-1)*nd; pa = tb.pa_dn[ia]
            for (ib, ib2, pb) in tb.dn
                amp = v[base+ib]
                amp == 0 && continue
                sgn = iseven(pa ⊻ pb) ? 1.0 : -1.0
                w[base+ib2] += -t*sgn*amp
            end
        end
    end
    return w
end

function lanczos_gs(S::Sector, bonds, t, U, mu;
                    m=parse(Int, get(ENV, "LANCZOS_M", "60")),
                    maxrestart=parse(Int, get(ENV, "LANCZOS_RESTART", "80")),
                    tol=1e-11, seed=1234)
    D = length(S); tabs = [BondTab(S, s1, s2) for (s1,s2) in bonds]
    Random.seed!(seed)
    v = randn(D); v ./= norm(v); Eprev = Inf; w = similar(v); E1 = NaN
    for restart in 1:maxrestart
        Q = [copy(v)]; α = Float64[]; β = Float64[]
        for j in 1:m
            apply_H!(w, Q[j], S, bonds, tabs, t, U, mu)
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
        E0 = F.values[1]; c = F.vectors[:,1]
        k > 1 && (E1 = F.values[2])
        v = zeros(D); for j in 1:k; v .+= c[j].*Q[j]; end
        v ./= norm(v)
        abs(E0 - Eprev) < tol && return E0, E1, v, restart
        Eprev = E0
    end
    return Eprev, E1, v, maxrestart
end

function expand_full(v, S::Sector)
    n = S.n; nd = length(S.dns)
    full = zeros(ComplexF64, 1 << (2n))
    @inbounds for (ia, a) in enumerate(S.ups), (ib, b) in enumerate(S.dns)
        q = 0
        for i in 1:n
            q |= bit(a,i) << (2i-2); q |= bit(b,i) << (2i-1)
        end
        full[q+1] = v[(ia-1)*nd + ib]
    end
    full
end

function window_rdm_full(full, N, q1, q2)
    Nw = q2 - q1 + 1
    dl = 1 << (q1-1); dw = 1 << Nw; dr = 1 << (N - q2)
    A = reshape(full, dl, dw, dr)
    ρ = zeros(ComplexF64, dw, dw)
    @inbounds for r in 1:dr, s2 in 1:dw, s1 in 1:dw
        acc = zero(ComplexF64)
        for l in 1:dl; acc += A[l,s1,r]*conj(A[l,s2,r]); end
        ρ[s1,s2] += acc
    end
    perm = [begin x=s-1; y=0; for k in 0:Nw-1; y |= ((x>>k)&1)<<(Nw-1-k); end; y+1 end for s in 1:dw]
    R = ρ[perm,perm]
    Hermitian((R + R')/2)
end

"""セクター表現のまま窓の縮約密度行列を作る(全ベクトル展開を避ける)。

4x4 (32量子ビット) では全ベクトルが 2^32 * 16 バイト = 68 GB で載らない。
窓はサイト i, i+1 の up/dn 占有 4 量子ビットなので、状態を
ψ[ia, ib](up配位 × dn配位の行列)として持てば
    ρ_w[(au,bu),(au',bu')] = Σ_{aenv} Σ_{benv}
        conj(ψ[ia(aenv,au), ib(benv,bu)]) ψ[ia(aenv,au'), ib(benv,bu')]
と書け、部分行列2枚のフロベニウス内積になる。粒子数保存から
popcount(au)=popcount(au') かつ popcount(bu)=popcount(bu') のブロックだけ残る。

窓内のビット順は window_rdm_full と同じ規約に合わせる:
    s = a_i*8 + b_i*4 + a_{i+1}*2 + b_{i+1}
"""
function window_rdm_sector(v::AbstractVector, S::Sector, i::Int)
    n = S.n; nd = length(S.dns)
    @inline ψ(ia, ib) = @inbounds v[(ia-1)*nd + ib]
    mask = (1 << (i-1)) | (1 << i)                    # サイト i, i+1
    # 各配位を (窓の2ビット, 環境ランク) に分解する
    function split(cfgs)
        wv = Vector{Int}(undef, length(cfgs))
        env = Vector{Int}(undef, length(cfgs))
        seen = Dict{Int,Int}()
        for (k, x) in enumerate(cfgs)
            wv[k] = (bit(x,i) << 1) | bit(x, i+1)
            e = x & ~mask
            env[k] = get!(seen, e, length(seen)+1)
        end
        wv, env, length(seen)
    end
    wa, ea, nea = split(S.ups)
    wb, eb, neb = split(S.dns)
    # (窓ビット値) -> 環境ランクごとの行/列番号
    rows = [fill(0, nea) for _ in 0:3]; cols = [fill(0, neb) for _ in 0:3]
    for k in 1:length(wa); rows[wa[k]+1][ea[k]] = k; end
    for k in 1:length(wb); cols[wb[k]+1][eb[k]] = k; end

    sidxw(au, bu) = ((au >> 1) << 3) | ((bu >> 1) << 2) | ((au & 1) << 1) | (bu & 1)
    ρ = zeros(ComplexF64, 16, 16)
    for au in 0:3, au2 in 0:3
        count_ones(au) == count_ones(au2) || continue
        ra, ra2 = rows[au+1], rows[au2+1]
        for bu in 0:3, bu2 in 0:3
            count_ones(bu) == count_ones(bu2) || continue
            cb, cb2 = cols[bu+1], cols[bu2+1]
            acc = zero(ComplexF64)
            @inbounds for f in 1:neb
                c1 = cb[f]; c2 = cb2[f]
                (c1 == 0 || c2 == 0) && continue
                for e in 1:nea
                    r1 = ra[e]; r2 = ra2[e]
                    (r1 == 0 || r2 == 0) && continue
                    acc += conj(ψ(r1, c1)) * ψ(r2, c2)
                end
            end
            ρ[sidxw(au,bu)+1, sidxw(au2,bu2)+1] = acc
        end
    end
    Hermitian((ρ + ρ')/2)
end

"""窓RDM を、窓内の部分集合 keep(1始まり、q1 が位置1)に部分トレースする。

観測量の台に含まれない量子ビットは Δ に一切寄与しないので、
「部分系の忠実度が利得を決めるか」を検証するには窓を台そのものに絞る必要がある。
"""
function ptrace_window(ρw::AbstractMatrix, Nw::Int, keep::Vector{Int})
    ks = sort(keep); nk = length(ks); dk = 1 << nk
    rest = [q for q in 1:Nw if !(q in ks)]
    out = zeros(ComplexF64, dk, dk)
    bitof(s, pos) = (s >> (Nw - pos)) & 1
    kidx(s) = (v = 0; for (j,q) in enumerate(ks); v |= bitof(s,q) << (nk - j); end; v)
    ridx(s) = (v = 0; for (j,q) in enumerate(rest); v |= bitof(s,q) << (length(rest) - j); end; v)
    for s1 in 0:(1<<Nw)-1, s2 in 0:(1<<Nw)-1
        ridx(s1) == ridx(s2) || continue
        out[kidx(s1)+1, kidx(s2)+1] += ρw[s1+1, s2+1]
    end
    Hermitian((out + out')/2)
end

"""Uhlmann 忠実度とトレース距離。"""
function fid_trdist(ρ, σ)
    A = Matrix(Hermitian((ρ+ρ')/2)); B = Matrix(Hermitian((σ+σ')/2))
    Fa = eigen(Hermitian(A)); sq = Fa.vectors*Diagonal(sqrt.(max.(Fa.values,0)))*Fa.vectors'
    F = (sum(sqrt.(max.(eigen(Hermitian(sq*B*sq)).values, 0))))^2
    D = 0.5*sum(abs.(eigen(Hermitian(A-B)).values))
    real(F), real(D)
end

# ---------------- DMRG / UHF ----------------
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
    MPO(os, sites)
end

function ground_state_2d(Lx, W, t, U, mu; chi_max, nsweeps=40, pbc_x)
    n = Lx*W; N = 2n
    sites = siteinds("Fermion", N; conserve_qns=true)
    H = hubbard_mpo(sites, Lx, W, t, U, mu; pbc_x)
    occ = falses(N)
    for x in 1:Lx, y in 1:W
        s = sidx(x,y); isodd(x+y) ? (occ[qup(s)] = true) : (occ[qdn(s)] = true)
    end
    ψ0 = productMPS(sites, [occ[q] ? "Occ" : "Emp" for q in 1:N])
    ramp = vcat([64,128,256,512], fill(chi_max, max(0, nsweeps-4)))
    noise = vcat([1e-5,1e-6,1e-7,1e-8,1e-9], zeros(max(0, nsweeps-5)))[1:nsweeps]
    E, ψ = dmrg(H, ψ0; nsweeps, maxdim=min.(chi_max, ramp)[1:nsweeps],
                cutoff=1e-12, noise, outputlevel=0)
    E, ψ, count(occ)
end

center_S(ψ) = begin
    N=length(ψ); ψo=copy(ψ); orthogonalize!(ψo, N÷2)
    _,Sv,_ = svd(ψo[N÷2], (linkinds(ψo)[N÷2-1], siteinds(ψo)[N÷2]))
    p = diag(Array(Sv,inds(Sv)...)).^2; p = p[p .> 1e-14]
    -sum(p.*log.(p))
end

function solve_uhf(Lx, W, t, U; Nup, Ndn, pbc_x, iters=5000, tol=1e-12, init=:neel)
    n = Lx*W; T = zeros(n,n)
    for (s,s2) in lat_edges(Lx,W; pbc_x); T[s,s2] = T[s2,s] = -t; end
    fill_avg = (Nup+Ndn)/(2n)
    nup = zeros(n); ndn = zeros(n)
    for x in 1:Lx, y in 1:W
        s = sidx(x,y); g = isodd(x+y) ? 1 : -1
        if init == :neel
            nup[s] = fill_avg + 0.4g; ndn[s] = fill_avg - 0.4g
        else
            nup[s] = fill_avg + 0.3*(rand()-0.5); ndn[s] = fill_avg + 0.3*(rand()-0.5)
        end
    end
    clamp!(nup,0.02,0.98); clamp!(ndn,0.02,0.98)
    nup .*= Nup/sum(nup); ndn .*= Ndn/sum(ndn)
    local Fu, Fd
    for _ in 1:iters
        Fu = eigen(Symmetric(T + diagm(U .* ndn))); Fd = eigen(Symmetric(T + diagm(U .* nup)))
        nu2 = vec(sum(abs2, Fu.vectors[:,1:Nup]; dims=2))
        nd2 = vec(sum(abs2, Fd.vectors[:,1:Ndn]; dims=2))
        if max(maximum(abs.(nu2 .- nup)), maximum(abs.(nd2 .- ndn))) < tol
            nup, ndn = nu2, nd2; break
        end
        nup = 0.5.*nu2 .+ 0.5.*nup; ndn = 0.5.*nd2 .+ 0.5.*ndn
    end
    Φu = Fu.vectors[:,1:Nup]; Φd = Fd.vectors[:,1:Ndn]
    Cu = Φu*Φu'; Cd = Φd*Φd'
    (Cu=Cu, Cd=Cd, E=tr(T*Cu)+tr(T*Cd)+U*sum(diag(Cu).*diag(Cd)),
     m=mean(abs(Cu[s,s]-Cd[s,s])/2 for s in 1:n))
end

function solve_uhf_best(Lx,W,t,U; Nup,Ndn,pbc_x)
    Random.seed!(31); best=nothing
    for ini in (:neel,:rand,:rand,:rand)
        u = solve_uhf(Lx,W,t,U; Nup,Ndn,pbc_x, init=ini)
        (best===nothing || u.E < best.E) && (best = u)
    end
    best
end

full_C(uhf, n) = begin
    C = zeros(ComplexF64, 2n, 2n)
    for a in 1:n, b in 1:n
        C[qup(a),qup(b)] = uhf.Cu[a,b]; C[qdn(a),qdn(b)] = uhf.Cd[a,b]
    end
    C
end

function rotate_C(C0, n, θ, φ)
    c, s = cos(θ/2), sin(θ/2); Um = zeros(ComplexF64, 2n, 2n)
    for i in 1:n
        a,b = qup(i), qdn(i)
        Um[a,a]=c; Um[a,b]=-exp(-im*φ)*s; Um[b,a]=exp(im*φ)*s; Um[b,b]=c
    end
    Um*C0*Um'
end

function symavg_P(obs, C0, n; nθ=32, nφ=12)
    P = [[0.0 for _ in o.terms] for o in obs]; w = 1.0/(nθ*nφ)
    for k in 1:nθ, l in 1:nφ
        θ = acos(-1 + (k-0.5)*2/nθ); φ = 2π*(l-0.5)/nφ
        M = majorana_M(rotate_C(C0, n, θ, φ))
        for (i,o) in enumerate(obs), (j,tm) in enumerate(o.terms)
            P[i][j] += gauss_pauli_expect(tm.sup, M)*w
        end
    end
    P
end

function site_obs_2d(s::Int)
    u,d = qup(s), qdn(s); u2,d2 = qup(s+1), qdn(s+1)
    [Obs("n", [Term(1.0,Tuple{Int,Int}[]), Term(-0.5,[(u,3)]), Term(-0.5,[(d,3)])], false),
     Obs("Sz", [Term(-0.25,[(u,3)]), Term(0.25,[(d,3)])], false),
     Obs("DoubleOcc", [Term(0.25,Tuple{Int,Int}[]), Term(-0.25,[(u,3)]),
                       Term(-0.25,[(d,3)]), Term(0.25,[(u,3),(d,3)])], false),
     Obs("ZZ onsite", [Term(1.0,[(u,3),(d,3)])], true),
     Obs("ZZ up-up nb", [Term(1.0,[(u,3),(u2,3)])], true),
     Obs("SzSz nb", [Term(1/16,[(u,3),(u2,3)]), Term(-1/16,[(u,3),(d2,3)]),
                     Term(-1/16,[(d,3),(u2,3)]), Term(1/16,[(d,3),(d2,3)])], false),
     Obs("SxSx nb", [Term(1/16,[(u,1),(d,1),(u2,1),(d2,1)]),
                     Term(1/16,[(u,1),(d,1),(u2,2),(d2,2)]),
                     Term(1/16,[(u,2),(d,2),(u2,1),(d2,1)]),
                     Term(1/16,[(u,2),(d,2),(u2,2),(d2,2)])], false)]
end

# ---------------- メイン ----------------
function main()
    t = 1.0
    U = parse(Float64, get(ENV, "UVAL", "8.0"))
    chi_exp  = parse(Int, get(ENV, "CHI_EXP", "2048"))
    n_repeat = parse(Int, get(ENV, "N_REPEAT", "300"))
    nu = parse(Int, get(ENV, "NU", "240")); nm = 100
    chi_priors = [4, 8, 32, 128]
    n = LX*W; N = 2n; mu = U/2
    Nup = n ÷ 2; Ndn = n - Nup

    @printf("=== 2D シリンダー vs トーラス (ED参照): W=%d Lx=%d, %d量子ビット, U=%.1f ===\n",
            W, LX, N, U)
    S = Sector(n, Nup, Ndn)
    @printf("セクター次元 = %d (%.2f GB の全ベクトルに展開)\n", length(S), (1<<N)*16/1e9)
    flush(stdout)

    rows = []; meta = []
    for pbc_x in (false, true)
        geo = pbc_x ? "torus" : "cylinder"
        bonds = lat_edges(LX, W; pbc_x)
        # 二部格子か(奇環がないか)を2彩色で判定する。奇環があると反強磁性が
        # フラストレートし、<P> 自体が潰れて「prior の良し悪し」の比較にならない。
        col = fill(-1, n); col[1] = 0; stack = [1]; bip = true
        adj = [Int[] for _ in 1:n]
        for (u,v) in bonds; push!(adj[u], v); push!(adj[v], u); end
        while !isempty(stack)
            u = pop!(stack)
            for v in adj[u]
                if col[v] == -1; col[v] = 1 - col[u]; push!(stack, v)
                elseif col[v] == col[u]; bip = false; end
            end
        end
        @printf("\n%s\n[%s] ボンド数 %d  二部格子: %s\n", "="^92, geo, length(bonds),
                bip ? "はい" : "**いいえ(奇環あり → 反強磁性がフラストレート)**")
        flush(stdout)

        t0 = time()
        E_ed, E1, vsec, nres = lanczos_gs(S, bonds, t, U, mu)
        @printf("  ED: E0=%.10f  E1-E0=%.3e  restart=%d (%.0fs)\n",
                E_ed, E1-E_ed, nres, time()-t0); flush(stdout)
        @assert E1 - E_ed > 1e-6 "基底状態が縮退している"
        GC.gc()   # Lanczos の Krylov 基底を解放してから測定に入る

        t1 = time()
        E_dm, ψ, nel = ground_state_2d(LX, W, t, U, mu; chi_max=chi_exp, pbc_x)
        Sc = center_S(ψ)
        @printf("  DMRG(χ≤%d): E0=%.10f  到達χ=%d  |ΔE|=%.2e  S_center=%.4f (%.0fs)\n",
                chi_exp, E_dm, maxlinkdim(ψ), abs(E_dm-E_ed), Sc, time()-t1); flush(stdout)

        uhf = solve_uhf_best(LX, W, t, U; Nup, Ndn, pbc_x); C0 = full_C(uhf, n)
        Muhf = majorana_M(C0)
        Muhf_list = [Muhf]                      # 共線UHF: 回転なし1点
        Msym_list = [majorana_M(rotate_C(C0, n, acos(-1+(k-0.5)*2/16), 2π*(l-0.5)/8))
                     for k in 1:16 for l in 1:8]   # 対称性回復: 球面16x8点
        @printf("  UHF: E=%.6f  スタッガード磁化 m=%.4f\n", uhf.E, uhf.m); flush(stdout)

        sites_meas = [s for s in 1:n if (s-1) % W < W-1]
        obs_by_site = [site_obs_2d(s) for s in sites_meas]
        windows = [obs_support(o) for o in obs_by_site]
        rdm_dm = window_rdms(ψ, windows)
        priors = MPS[]
        for cp in chi_priors
            σ = truncate(ψ; maxdim=cp); normalize!(σ); push!(priors, σ)
        end
        push!(priors, ψ)
        fids = [abs2(inner(σ, ψ)) for σ in priors]
        rdm_p = [window_rdms(σ, windows) for σ in priors]
        labels = vcat(["chi$c" for c in chi_priors], ["exact", "UHF", "UHFsym"])
        @printf("  prior忠実度(vs DMRG): %s\n",
                join([@sprintf("χ=%d:%.3e", c, f) for (c,f) in zip(chi_priors, fids)], "  "))
        flush(stdout)

        maxdiff = 0.0
        for (si, s) in enumerate(sites_meas)
            obs = obs_by_site[si]; q1,q2 = windows[si]; Nw = q2-q1+1
            ow = shift_obs(obs, q1-1)
            Mats = [[term_window_matrix(tm, Nw) for tm in o.terms] for o in ow]
            ρ_ed = window_rdm_sector(vsec, S, s)
            Pσ_all = Vector{Vector{Vector{Float64}}}(); trOσ_all = Vector{Vector{Float64}}()
            for p in 1:length(priors)
                Pσ = [[expect_rdm(rdm_p[p][si], M) for M in Ms] for Ms in Mats]
                push!(Pσ_all, Pσ)
                push!(trOσ_all, [sum(tm.coeff*Pσ[k][ti] for (ti,tm) in enumerate(obs[k].terms))
                                 for k in 1:length(obs)])
            end
            Pu = [[gauss_pauli_expect(tm.sup, Muhf) for tm in o.terms] for o in obs]
            Ps = symavg_P(obs, C0, n)
            for P in (Pu, Ps)
                push!(Pσ_all, P)
                push!(trOσ_all, [sum(tm.coeff*P[k][ti] for (ti,tm) in enumerate(obs[k].terms))
                                 for k in 1:length(obs)])
            end
            supp = [sort(unique(q - (q1-1) for tm in o.terms for (q,_) in tm.sup))
                    for o in obs]
            for (rname, ρref) in (("ED", ρ_ed), ("DMRG", rdm_dm[si]))
                Otrue = [sum(tm.coeff*expect_rdm(ρref, Mats[k][ti])
                             for (ti,tm) in enumerate(ow[k].terms)) for k in 1:length(obs)]
                if rname == "ED"
                    Odm = [sum(tm.coeff*expect_rdm(rdm_dm[si], Mats[k][ti])
                               for (ti,tm) in enumerate(ow[k].terms)) for k in 1:length(obs)]
                    maxdiff = max(maxdiff, maximum(abs.(Otrue .- Odm)))
                end
                es, ec = run_locals(ρref, Nw, ow, Pσ_all, trOσ_all;
                                    nu, nm, n_repeat, seed = 990_000 + 100s + (pbc_x ? 1 : 0))
                # prior ごとの窓RDM。平均場は窓内の全 Pauli 列の期待値から
                # σ_w = 2^{-Nw} Σ_P <P>_σ P と再構成する(§2i の一般Wickが要る)。
                σws = Any[rdm_p[p][si] for p in 1:length(priors)]
                for Ms in (Muhf_list, Msym_list)
                    d = 1 << Nw; Mw = zeros(ComplexF64, d, d)
                    for idx in 0:(4^Nw - 1)
                        labs = Vector{Int}(undef, Nw); x = idx
                        for j in Nw:-1:1; labs[j] = x % 4; x ÷= 4; end
                        sup = [(q + q1 - 1, labs[q]) for q in 1:Nw if labs[q] != 0]
                        v = mean(gauss_pauli_expect(sup, M) for M in Ms)
                        v == 0 && continue
                        Mw .+= v .* reduce(kron, [labs[q]==0 ? Matrix{ComplexF64}(I,2,2) :
                                                  SIGMA[labs[q]] for q in 1:Nw])
                    end
                    push!(σws, Hermitian(Mw ./ d))
                end
                for (k,o) in enumerate(obs)
                    v = var(es[:,k]); gmx = v/var(ec[:,k,length(chi_priors)+1])
                    ρA = isempty(supp[k]) ? nothing : ptrace_window(ρref, Nw, supp[k])
                    for (p,lb) in enumerate(labels)
                        Δ = Otrue[k] - trOσ_all[p][k]
                        FA, DA = (ρA === nothing) ? (1.0, 0.0) :
                            fid_trdist(Matrix(ρA), Matrix(ptrace_window(σws[p], Nw, supp[k])))
                        push!(rows, (geo, rname, U, s, o.name, lb, Otrue[k], Δ,
                                     abs(Otrue[k])>1e-6 ? abs(Δ/Otrue[k]) : NaN,
                                     v/var(ec[:,k,p]), gmx, FA, 1-FA, DA, length(supp[k])))
                    end
                end
            end
        end
        vsec = nothing; GC.gc()
        @printf("  局所観測量の ED vs DMRG 最大差: %.3e\n", maxdiff); flush(stdout)
        push!(meta, (geo, U, E_ed, E_dm, abs(E_dm-E_ed), maxlinkdim(ψ), Sc, uhf.E, uhf.m, maxdiff))

        # ρ=ED での利得(全サイト中央値)
        @printf("\n  %-14s" , "観測量")
        for lb in labels; @printf("%10s", lb); end
        @printf("%10s\n", "天井")
        for o in unique(r[5] for r in rows if r[1]==geo)
            @printf("  %-14s", o)
            for lb in labels
                h = [r[10] for r in rows if r[1]==geo && r[2]=="ED" && r[5]==o && r[6]==lb]
                @printf("%10.2f", isempty(h) ? NaN : median(h))
            end
            gm = median([r[11] for r in rows if r[1]==geo && r[2]=="ED" && r[5]==o])
            @printf("%10.2f\n", gm)
        end
        flush(stdout)
    end

    println("\n", "="^92)
    @printf("%-10s %18s %18s %12s %10s %10s %12s\n",
            "幾何","E0(ED)","E0(DMRG)","|ΔE|","到達χ","S_center","局所量max差")
    for r in meta
        @printf("%-10s %18.10f %18.10f %12.2e %10d %10.4f %12.2e\n",
                r[1], r[3], r[4], r[5], r[6], r[7], r[10])
    end

    # 出力名に格子形状を入れる(入れ忘れて 4x3 の結果を 2x6 が上書きした)
    out = joinpath(@__DIR__, @sprintf("crm_2d_supp_W%dL%d_U%s.tsv", W, LX,
                                      replace(@sprintf("%.1f",U),"."=>"p")))
    open(out,"w") do io
        println(io, "geometry\trho\tU\tsite\tobservable\tprior\ttrue\tDelta\trelerr\tG\tG_max\tF_supp\tinfid_supp\tD_supp\tnsupp")
        for r in rows; println(io, join(r,"\t")); end
    end
    println("\nresults saved: $out")
end
get(ENV, "NO_MAIN", "0") == "1" || main()
