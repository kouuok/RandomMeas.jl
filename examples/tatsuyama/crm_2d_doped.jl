# ============================================================
# CRMシャドウの2D拡張 その4: ドープ系 (1/8 hole doping, ストライプUHF prior)
#   W=4x8, U=8, 電子28個 (Nup=Ndn=14, μ=0; 粒子数はQNで固定)。
#   ドープ2DハバードはDMRGが最も苦しむ領域で、UHFはストライプ解
#   (電荷の川 + 反位相AFドメイン) を出す。安いpriorがどこまで働くかを検証。
#   UHFは複数初期値 (Néel / ストライプ) から最低エネルギー解を採用。
#
# crm_2d_cylinder.jl の結果: UHFはオンサイト量ではχ=32-64のMPS priorに
# 匹敵するが、サイト間スピン相関では凍結Néel秩序が |<SzSz>| を過大評価して
# 失敗する (G<=3)。
#
# 本スクリプト: スピン回転で平均した混合状態
#   σ̄ = ∫dΩ R(Ω) |UHF><UHF| R†(Ω)
# を prior に追加する。
#  - CRMのpriorに必要なのは <P>_σ のみなので★混合状態でも使える
#  - 回転した行列式は一般化(非共線)Slater行列式 → フル 2n×2n 相関行列
#    C(θ) = U(θ) C U(θ)† に対する Wick の定理で厳密計算 (コスト無視可能)
#  - 共線Néel秩序はz軸回り回転で不変 → 平均は極角θの1次元求積のみ
#  - 期待: オンサイト量は回転不変(UHFの強み保持)、サイト間スピン相関は
#    平均で凍結秩序の過大評価が補正される → 全観測量で安価priorが機能するか
#
# 追加検証: 一般化Wickがθ=0で共線Wickと一致 / U=0(一重項RHF)で回転平均が不変
#
# 設計:
#  - 4 x Lx シリンダー (y周期, x開放), スネークJW順序 site s=(x-1)W+y,
#    qubit qup(s)=2s-1, qdn(s)=2s
#  - 参照状態 ρ: DMRG chi=192 (固定参照。完全収束基底状態である必要はない)
#  - prior: (a) ρの切断MPS chi=8..64, (b) UHF (Néel初期値の自己無撞着平均場)
#  - ★UHF prior の厳密古典側は Wick の定理で閉形式:
#      <Z_q> = 1-2n_q
#      <Z_q Z_q'> = (1-2n)(1-2n') - 4C_{ss'}²  (同スピン; 異スピンは積)
#      <(XZX+YZY)/2> = 2 C_{ss'}               (JW隣接ボンド)
#    → MPS化不要。局所Pauliの prior 側は基底一致 × 3^|A| <P>_σ のみ
#  - 検証: U=0 で DMRGエネルギー vs 厳密自由フェルミオン、
#    Wick期待値 vs DMRG期待値 (ハミルトニアン/JW/C行列/Wick式を一括検証)
#  - 観測量は中央サイト周りの連続窓 (~10 qubits) に配置し窓サンプリング
#
# 実行方法:
#  JULIA_LOAD_PATH="@:@v#.#:@stdlib" julia --project=Hubbard_MPS_Env_v2 crm_2d_cylinder.jl
# ============================================================

using ITensors, ITensorMPS
using LinearAlgebra
using Statistics
using Printf
using Random

# ------------------------------------------------------------
# 0. 格子
# ------------------------------------------------------------
const W  = 4
const LX = 8
const CHI_EXP = parse(Int, get(ENV, "CHI_EXP", "512"))
const NHOLE = parse(Int, get(ENV, "NHOLE", "4"))     # ホール数 (1/8 doping = 4)
sidx(x, y) = (x - 1) * W + y
qup(s) = 2s - 1
qdn(s) = 2s

# シリンダーのボンドリスト (site index ペア)
function cyl_edges(Lx, W)
    edges = Tuple{Int,Int}[]
    for x in 1:Lx, y in 1:W
        s = sidx(x, y)
        x < Lx && push!(edges, (s, sidx(x+1, y)))
        push!(edges, (s, sidx(x, y == W ? 1 : y + 1)))
    end
    return edges
end

# ------------------------------------------------------------
# 1. DMRG参照状態
# ------------------------------------------------------------
function hubbard_cyl_mpo(sites, Lx, W, t, U, mu)
    os = OpSum()
    for (s, s2) in cyl_edges(Lx, W)
        for (qa, qb) in ((qup(s), qup(s2)), (qdn(s), qdn(s2)))
            os += -t, "Cdag", qa, "C", qb
            os += -t, "Cdag", qb, "C", qa
        end
    end
    for s in 1:Lx*W
        os += U, "N", qup(s), "N", qdn(s)
        os += -mu, "N", qup(s)
        os += -mu, "N", qdn(s)
    end
    return MPO(os, sites)
end

function ground_state_cyl(Lx, W, t, U, mu; chi_max=192, nsweeps=10, nhole=NHOLE)
    N = 2 * Lx * W
    sites = siteinds("Fermion", N; conserve_qns=true)
    H = hubbard_cyl_mpo(sites, Lx, W, t, U, mu)
    # Néel初期状態から NHOLE 個のホールをx方向に等間隔で抜く
    init = fill("Emp", N)
    occ_sites = Tuple{Int,Int}[]   # (qubit, x)
    for x in 1:Lx, y in 1:W
        s = sidx(x, y)
        q = iseven(x + y) ? qup(s) : qdn(s)
        init[q] = "Occ"
        push!(occ_sites, (q, x))
    end
    if nhole > 0
        # up/dn交互に、x位置を分散させて抜く
        ups = [q for (q, x) in occ_sites if isodd(q)]
        dns = [q for (q, x) in occ_sites if iseven(q)]
        for k in 1:(nhole ÷ 2)
            init[ups[k * length(ups) ÷ (nhole ÷ 2 + 1)]] = "Emp"
            init[dns[k * length(dns) ÷ (nhole ÷ 2 + 1)]] = "Emp"
        end
    end
    ψ0 = productMPS(sites, init)
    ramp = [20, 50, 100, 200, 400, chi_max]
    maxdims = min.(chi_max, [ramp; fill(chi_max, max(nsweeps - length(ramp), 0))])[1:nsweeps]
    nz = [1e-5, 1e-5, 1e-6, 1e-6, 1e-7, 1e-8]
    noise = [nz; fill(0.0, max(nsweeps - length(nz), 0))][1:nsweeps]
    E, ψ = dmrg(H, ψ0; nsweeps, maxdim=maxdims, cutoff=1e-9, noise, outputlevel=1)
    return E, ψ
end

# ------------------------------------------------------------
# 2. UHF平均場 + Wick期待値
# ------------------------------------------------------------
function hopping_matrix(Lx, W, t)
    n = Lx * W
    T = zeros(n, n)
    for (s, s2) in cyl_edges(Lx, W)
        T[s, s2] = T[s2, s] = -t
    end
    return T
end

function solve_uhf(Lx, W, t, U; Nup, Ndn, iters=2000, tol=1e-10, init=:neel)
    n = Lx * W
    T = hopping_matrix(Lx, W, t)
    nup = zeros(n); ndn = zeros(n)
    fill_avg = (Nup + Ndn) / (2n)
    for x in 1:Lx, y in 1:W
        s = sidx(x, y)
        if init == :neel
            nup[s] = fill_avg + 0.4 * (-1)^(x + y)
            ndn[s] = fill_avg - 0.4 * (-1)^(x + y)
        elseif init == :stripe
            # 反位相AFドメイン + ドメイン壁 (x=Lx/2+0.5) にホールの川
            wall = Lx / 2 + 0.5
            dom = x <= Lx ÷ 2 ? 1.0 : -1.0
            hole = 0.35 * exp(-abs(x - wall))
            nup[s] = fill_avg - hole + dom * 0.4 * (-1)^(x + y)
            ndn[s] = fill_avg - hole - dom * 0.4 * (-1)^(x + y)
        else   # :random
            nup[s] = fill_avg + 0.3 * (rand() - 0.5)
            ndn[s] = fill_avg + 0.3 * (rand() - 0.5)
        end
    end
    clamp!(nup, 0.02, 0.98); clamp!(ndn, 0.02, 0.98)
    nup .*= Nup / sum(nup); ndn .*= Ndn / sum(ndn)
    local Fu, Fd
    for it in 1:iters
        Fu = eigen(Symmetric(T + diagm(U .* ndn)))
        Fd = eigen(Symmetric(T + diagm(U .* nup)))
        new_up = vec(sum(abs2, Fu.vectors[:, 1:Nup]; dims=2))
        new_dn = vec(sum(abs2, Fd.vectors[:, 1:Ndn]; dims=2))
        if max(maximum(abs.(new_up .- nup)), maximum(abs.(new_dn .- ndn))) < tol
            nup, ndn = new_up, new_dn
            break
        end
        nup = 0.5 .* new_up .+ 0.5 .* nup
        ndn = 0.5 .* new_dn .+ 0.5 .* ndn
    end
    Φu = Fu.vectors[:, 1:Nup]; Φd = Fd.vectors[:, 1:Ndn]
    Cu = Φu * Φu'; Cd = Φd * Φd'
    E = tr(T * Cu) + tr(T * Cd) + U * sum(diag(Cu) .* diag(Cd))
    return (nup=diag(Cu), ndn=diag(Cd), Cu=Cu, Cd=Cd, E=E)
end

# 複数初期値からエネルギー最小のUHF解を採用
function solve_uhf_best(Lx, W, t, U; Nup, Ndn)
    Random.seed!(31)
    best = nothing
    for ini in (:neel, :stripe, :random)
        uhf = solve_uhf(Lx, W, t, U; Nup, Ndn, init=ini)
        @printf("  UHF init=%-7s E=%.6f\n", ini, uhf.E)
        (best === nothing || uhf.E < best[2].E) && (best = (ini, uhf))
    end
    @printf("  -> selected init=%s (E=%.6f)\n", best[1], best[2].E)
    # プロファイル出力 (列ごとの平均密度とスタッガード磁化)
    uhf = best[2]
    for x in 1:Lx
        dens = mean(uhf.nup[sidx(x,y)] + uhf.ndn[sidx(x,y)] for y in 1:W)
        mag  = mean((-1)^(x+y) * (uhf.nup[sidx(x,y)] - uhf.ndn[sidx(x,y)]) / 2 for y in 1:W)
        @printf("  col x=%d: <n>=%.3f  m_stag=%+.3f\n", x, dens, mag)
    end
    return uhf
end

site_spin(q) = (div(q + 1, 2), isodd(q) ? 1 : 2)   # spin 1=up, 2=dn

# 1つのPauli項のUHF期待値 (Wickの定理; 本スクリプトの観測量クラスのみ対応)
function wick_term(sup::Vector{Tuple{Int,Int}}, uhf)
    isempty(sup) && return 1.0
    ps = [a for (_, a) in sup]
    if all(==(3), ps)
        ss = [site_spin(q) for (q, _) in sup]
        nval(s, σ) = σ == 1 ? uhf.nup[s] : uhf.ndn[s]
        if length(ss) == 1
            (s, σ) = ss[1]
            return 1 - 2nval(s, σ)
        elseif length(ss) == 2
            (s1, σ1), (s2, σ2) = ss
            z = (1 - 2nval(s1, σ1)) * (1 - 2nval(s2, σ2))
            if σ1 == σ2
                Cm = σ1 == 1 ? uhf.Cu : uhf.Cd
                return z - 4 * Cm[s1, s2]^2
            else
                return z
            end
        end
        error("wick_term: unsupported Z-string length $(length(ss))")
    elseif length(sup) == 3 && (ps == [1, 3, 1] || ps == [2, 3, 2])
        (s1, σ1) = site_spin(sup[1][1]); (s2, σ2) = site_spin(sup[3][1])
        @assert σ1 == σ2
        Cm = σ1 == 1 ? uhf.Cu : uhf.Cd
        return 2 * Cm[s1, s2]      # 数保存状態ではXZXもYZYも同値
    end
    error("wick_term: unsupported term $(sup)")
end

# ---- 対称性回復UHF: 一般化(非共線)相関行列に対するWick ----

# UHFのブロック対角相関行列をフル 2n×2n (qubit=JWモード順) に展開
function full_C(uhf, n)
    C = zeros(2n, 2n)
    for s1 in 1:n, s2 in 1:n
        C[qup(s1), qup(s2)] = uhf.Cu[s1, s2]
        C[qdn(s1), qdn(s2)] = uhf.Cd[s1, s2]
    end
    return C
end

# 全サイト一斉スピン回転 (y軸回り角θ) をかけた相関行列
function rotate_C(C0::Matrix{Float64}, n::Int, θ::Float64)
    c, s = cos(θ/2), sin(θ/2)
    U = zeros(2n, 2n)
    for site in 1:n
        a, b = qup(site), qdn(site)
        U[a, a] = c;  U[a, b] = -s
        U[b, a] = s;  U[b, b] = c
    end
    return U * C0 * U'
end

# 一般化相関行列 C (実対称) に対するPauli項のWick期待値
function wick_term_gen(sup::Vector{Tuple{Int,Int}}, C::Matrix{Float64})
    isempty(sup) && return 1.0
    ps = [a for (_, a) in sup]
    if all(==(3), ps)
        qs = [q for (q, _) in sup]
        if length(qs) == 1
            return 1 - 2C[qs[1], qs[1]]
        elseif length(qs) == 2
            a, b = qs
            return 1 - 2C[a,a] - 2C[b,b] + 4*(C[a,a]*C[b,b] - C[a,b]*C[b,a])
        end
        error("wick_term_gen: unsupported Z-string length")
    elseif length(sup) == 3 && (ps == [1, 3, 1] || ps == [2, 3, 2])
        a = sup[1][1]; b = sup[3][1]
        return 2 * C[a, b]
    end
    error("wick_term_gen: unsupported term $(sup)")
end

# スピン回転平均 <P>_σ̄ : 極角cosθの中点則求積 (被積分関数はcosθの低次多項式)
function symavg_P(obs, uhf, n; nθ=64)
    C0 = full_C(uhf, n)
    P = [[0.0 for _ in o.terms] for o in obs]
    for k in 1:nθ
        θ = acos(-1 + (k - 0.5) * 2 / nθ)
        C = rotate_C(C0, n, θ)
        for (i, o) in enumerate(obs), (j, tm) in enumerate(o.terms)
            P[i][j] += wick_term_gen(tm.sup, C) / nθ
        end
    end
    return P
end

# ------------------------------------------------------------
# 3. MPS機構 (crm_mps_scaling.jl と同一の検証済みコード)
# ------------------------------------------------------------
function site_matrices(ψo::MPS, j::Int)
    N = length(ψo); ss = siteinds(ψo); ls = linkinds(ψo)
    if j == 1
        A = Array(ψo[j], ss[j], ls[j])
        return (reshape(A[1, :], 1, :), reshape(A[2, :], 1, :))
    elseif j == N
        A = Array(ψo[j], ls[j-1], ss[j])
        return (reshape(A[:, 1], :, 1), reshape(A[:, 2], :, 1))
    else
        A = Array(ψo[j], ls[j-1], ss[j], ls[j])
        return (Matrix(A[:, 1, :]), Matrix(A[:, 2, :]))
    end
end

extract_tensors(ψ::MPS) = (ψo = orthogonalize(ψ, 1); [site_matrices(ψo, j) for j in 1:length(ψo)])

function extract_window(ψ::MPS, q1::Int, q2::Int)
    ψo = orthogonalize(ψ, q1)
    tens = [site_matrices(ψo, j) for j in q1:q2]
    C1, C2 = tens[1]
    pl = vec(sum(abs2, C1; dims=2) .+ sum(abs2, C2; dims=2))
    pl ./= sum(pl)
    return tens, cumsum(pl)
end

function product_op_expect(tens, ops::Dict{Int,Matrix{ComplexF64}})
    Lenv = ones(ComplexF64, 1, 1)
    for j in 1:length(tens)
        O = get(ops, j, nothing)
        A1, A2 = tens[j]
        if O === nothing
            Lenv = A1' * Lenv * A1 .+ A2' * Lenv * A2
        else
            Lenv = O[1,1] .* (A1' * Lenv * A1) .+ O[1,2] .* (A1' * Lenv * A2) .+
                   O[2,1] .* (A2' * Lenv * A1) .+ O[2,2] .* (A2' * Lenv * A2)
        end
    end
    return real(Lenv[1,1])
end

const SIGMA  = (ComplexF64[0 1; 1 0], ComplexF64[0 -im; im 0], ComplexF64[1 0; 0 -1])
const UBASIS = (ComplexF64[1 1; 1 -1]/sqrt(2), ComplexF64[1 -im; 1 im]/sqrt(2),
                ComplexF64[1 0; 0 1])
const DVAL = ntuple(β -> ntuple(a -> real.(diag(UBASIS[β]*SIGMA[a]*UBASIS[β]')), 3), 3)

make_rot_buffers(tens) = [(zeros(ComplexF64, size(t[1])), zeros(ComplexF64, size(t[2]))) for t in tens]

function rotate!(rot, tens, basis)
    for j in 1:length(tens)
        u = UBASIS[basis[j]]
        A1, A2 = tens[j]
        rot[j][1] .= u[1,1] .* A1 .+ u[1,2] .* A2
        rot[j][2] .= u[2,1] .* A1 .+ u[2,2] .* A2
    end
end

function sample_bits!(bits::Vector{Int}, rot, vbuf, wbuf, cum_pl)
    Nw = length(rot)
    chil = length(cum_pl)
    l = chil == 1 ? 1 : searchsortedfirst(cum_pl, rand())
    v = vbuf
    @inbounds for k in 1:chil; v[k] = 0; end
    v[l] = 1.0 + 0im
    vlen = chil
    for j in 1:Nw
        wlen = size(rot[j][1], 2)
        w0 = view(wbuf, 1:wlen); w1 = view(wbuf, wlen+1:2wlen)
        mul!(w0, transpose(rot[j][1]), view(v, 1:vlen))
        mul!(w1, transpose(rot[j][2]), view(v, 1:vlen))
        p0 = real(dot(w0, w0)); p1 = real(dot(w1, w1))
        b = rand() * (p0 + p1) < p0 ? 0 : 1
        bits[j] = b
        src = b == 0 ? w0 : w1
        nrm = sqrt(b == 0 ? p0 : p1)
        @inbounds for k in 1:wlen; v[k] = src[k] / nrm; end
        vlen = wlen
    end
    return bits
end

# ------------------------------------------------------------
# 4. 観測量 (中央サイト s0=(Lx/2, 2) 周りの窓に配置)
# ------------------------------------------------------------
struct Term
    coeff::Float64
    sup::Vector{Tuple{Int,Int}}
end
struct Obs
    name::String
    terms::Vector{Term}
    pure::Bool
end

term_matrixdict(t::Term) = Dict{Int,Matrix{ComplexF64}}(q => SIGMA[a] for (q,a) in t.sup)

function build_observables_2d(Lx, W)
    s0 = sidx(Lx ÷ 2, 2)
    sy1 = s0 + 1          # y方向隣 (同列)
    sy2 = s0 + 2          # y方向距離2 (リング対面)
    sx1 = s0 + W          # x方向隣 (次列)
    obs = Obs[]
    push!(obs, Obs("ZZ onsite",   [Term(1.0, [(qup(s0),3),(qdn(s0),3)])], true))
    push!(obs, Obs("ZZ up-up y1", [Term(1.0, [(qup(s0),3),(qup(sy1),3)])], true))
    push!(obs, Obs("ZZ up-up y2", [Term(1.0, [(qup(s0),3),(qup(sy2),3)])], true))
    push!(obs, Obs("ZZ up-up x1", [Term(1.0, [(qup(s0),3),(qup(sx1),3)])], true))
    push!(obs, Obs("DoubleOcc",
        [Term(0.25, Tuple{Int,Int}[]), Term(-0.25, [(qup(s0),3)]),
         Term(-0.25, [(qdn(s0),3)]),  Term(0.25, [(qup(s0),3),(qdn(s0),3)])], false))
    push!(obs, Obs("SzSz y-bond",
        [Term( 1/16, [(qup(s0),3),(qup(sy1),3)]), Term(-1/16, [(qup(s0),3),(qdn(sy1),3)]),
         Term(-1/16, [(qdn(s0),3),(qup(sy1),3)]), Term( 1/16, [(qdn(s0),3),(qdn(sy1),3)])], false))
    push!(obs, Obs("SzSz x-bond",
        [Term( 1/16, [(qup(s0),3),(qup(sx1),3)]), Term(-1/16, [(qup(s0),3),(qdn(sx1),3)]),
         Term(-1/16, [(qdn(s0),3),(qup(sx1),3)]), Term( 1/16, [(qdn(s0),3),(qdn(sx1),3)])], false))
    push!(obs, Obs("hop y-bond",
        [Term(0.5, [(qup(s0),1),(qdn(s0),3),(qup(sy1),1)]),
         Term(0.5, [(qup(s0),2),(qdn(s0),3),(qup(sy1),2)])], false))
    # サイト分解密度 (ドープ系の電荷変調 = ストライプ秩序の直接測定)
    for st in (s0, sx1)
        push!(obs, Obs("n(site $st)",
            [Term(1.0, Tuple{Int,Int}[]), Term(-0.5, [(qup(st),3)]),
             Term(-0.5, [(qdn(st),3)])], false))
    end
    return obs
end

obs_support(obs) = (minimum(q for o in obs for t in o.terms for (q,_) in t.sup),
                    maximum(q for o in obs for t in o.terms for (q,_) in t.sup))

shift_obs(obs, off) = [Obs(o.name,
    [Term(t.coeff, [(q - off, a) for (q, a) in t.sup]) for t in o.terms], o.pure) for o in obs]

function estimate_obs(o::Obs, basis::Vector{Int}, bits::Vector{Int})
    x = 0.0
    for t in o.terms
        v = t.coeff
        for (q, a) in t.sup
            v *= 3.0 * DVAL[basis[q]][a][bits[q]+1]
            v == 0.0 && break
        end
        x += v
    end
    return x
end

function exact_prior_mean(o::Obs, basis::Vector{Int}, Pσ::Vector{Float64})
    m = 0.0
    for (ti, t) in enumerate(o.terms)
        ok = all(basis[q] == a for (q, a) in t.sup)
        ok || continue
        m += t.coeff * 3.0^length(t.sup) * Pσ[ti]
    end
    return m
end

theory_var(absA, P, Δ, nu, nm) = (((3.0^absA - 1)*P^2 + 3.0^absA*(1-P^2)/nm)/nu,
                                  ((3.0^absA - 1)*Δ^2 + 3.0^absA*(1-P^2)/nm)/nu)

function run_locals(wtens, cum_pl, obs_w, Pσ_all, trOσ_all; nu, nm, n_repeat, seed)
    Random.seed!(seed)
    Nw = length(wtens); nobs = length(obs_w); np = length(Pσ_all)
    chimax = max(length(cum_pl), maximum(size(t[1], 2) for t in wtens))
    rot = make_rot_buffers(wtens)
    vbuf = zeros(ComplexF64, chimax + 1); wbuf = zeros(ComplexF64, 2chimax + 2)
    bits = zeros(Int, Nw)
    est_std = zeros(n_repeat, nobs); est_crm = zeros(n_repeat, nobs, np)
    xs = zeros(nobs)
    for rep in 1:n_repeat
        acc_s = zeros(nobs); acc_c = zeros(nobs, np)
        for _ in 1:nu
            basis = rand(1:3, Nw)
            rotate!(rot, wtens, basis)
            fill!(xs, 0.0)
            for _ in 1:nm
                sample_bits!(bits, rot, vbuf, wbuf, cum_pl)
                for k in 1:nobs
                    xs[k] += estimate_obs(obs_w[k], basis, bits)
                end
            end
            for k in 1:nobs
                mρ = xs[k] / nm
                acc_s[k] += mρ
                for p in 1:np
                    acc_c[k, p] += mρ - exact_prior_mean(obs_w[k], basis, Pσ_all[p][k])
                end
            end
        end
        est_std[rep, :] .= acc_s ./ nu
        for p in 1:np, k in 1:nobs
            est_crm[rep, k, p] = acc_c[k, p] / nu + trOσ_all[p][k]
        end
    end
    return est_std, est_crm
end

# ------------------------------------------------------------
# 5. 検証 (U=0): DMRG vs 厳密自由フェルミオン vs Wick
# ------------------------------------------------------------
# 検証は W=4 x Lx=2 (16 qubits) で行う: 中央カットのSchmidtランク上限が
# 2^8=256 なので chi=256 のMPSは厳密。DMRG収束性ではなく
# JW順序・シリンダー結合・C行列対応・Wick式の整合性を機械精度で検証する。
function validate_u0()
    t = 1.0
    Lxv = 2
    n = Lxv * W; Nup = n ÷ 2; Ndn = n ÷ 2
    T = hopping_matrix(Lxv, W, t)
    ε = sort(eigvals(Symmetric(T)))
    E_free = 2 * sum(ε[1:Nup])
    @printf("  [check] free-fermion E = %.8f  (gap at Fermi level: %.4f)\n",
            E_free, ε[Nup+1] - ε[Nup]); flush(stdout)

    BLAS.set_num_threads(4)
    E_dmrg, ψ = ground_state_cyl(Lxv, W, t, 0.0, 0.0; chi_max=256, nsweeps=10, nhole=0)
    BLAS.set_num_threads(1)
    @printf("  [check] E_DMRG(U=0) = %.8f  (diff %.2e)\n", E_dmrg, abs(E_dmrg - E_free))
    @assert abs(E_dmrg - E_free) < 1e-6

    # U=0ではUHF(=RHF)が厳密 → Wick期待値 vs DMRG期待値
    uhf = solve_uhf(Lxv, W, t, 0.0; Nup, Ndn)
    tens = extract_tensors(ψ)
    obs = build_observables_2d(Lxv, W)
    maxdiff = 0.0
    for o in obs
        v_dmrg = sum(tm.coeff * product_op_expect(tens, term_matrixdict(tm)) for tm in o.terms)
        v_wick = sum(tm.coeff * wick_term(tm.sup, uhf) for tm in o.terms)
        d = abs(v_dmrg - v_wick); maxdiff = max(maxdiff, d)
        @printf("  [check] %-12s DMRG=%9.5f  Wick=%9.5f  (diff %.1e)\n",
                o.name, v_dmrg, v_wick, d)
        @assert d < 1e-4 "$(o.name) mismatch"
    end
    println("  [check] all U=0 validations passed (maxdiff $(round(maxdiff, sigdigits=2)))")

    # 一般化Wick: θ=0 で共線Wickと一致 (有限UのNéel解で確認)
    uhf4 = solve_uhf(Lxv, W, t, 4.0; Nup, Ndn)
    C0 = full_C(uhf4, n)
    d1 = maximum(abs(wick_term_gen(tm.sup, C0) - wick_term(tm.sup, uhf4))
                 for o in obs for tm in o.terms)
    @printf("  [check] generalized Wick vs collinear Wick at theta=0: maxdiff %.1e\n", d1)
    @assert d1 < 1e-12

    # U=0 (RHF, スピン一重項) では回転平均が不変
    P0 = [[wick_term(tm.sup, uhf) for tm in o.terms] for o in obs]
    Pavg = symavg_P(obs, uhf, n)
    d2 = maximum(abs(P0[i][j] - Pavg[i][j])
                 for i in eachindex(obs) for j in eachindex(obs[i].terms))
    @printf("  [check] rotation-average invariance at U=0 (singlet): maxdiff %.1e\n", d2)
    @assert d2 < 1e-10
    println("  [check] symmetry-restoration checks passed")
    flush(stdout)
end

# ------------------------------------------------------------
# 6. メイン
# ------------------------------------------------------------
function main()
    t = 1.0
    chi_priors = haskey(ENV, "CHIP") ? parse.(Int, split(ENV["CHIP"], ",")) : [8, 16, 32, 64, 128]
    prior_labels = vcat("UHF", "UHF-sym", ["chi=$c" for c in chi_priors])
    nu = parse(Int, get(ENV, "NU", "40"))
    nm = parse(Int, get(ENV, "NM", "100"))
    n_repeat = parse(Int, get(ENV, "NREP", "25"))
    n = LX * W
    Nel = n - NHOLE                     # 全電子数 (half-filling n から NHOLE 個抜く)
    Nup = Nel ÷ 2; Ndn = Nel - Nup

    println("=== validation (U=0, $(W)x2 cylinder, exact MPS) ==="); flush(stdout)
    validate_u0()

    rows = []
    U_first = parse(Float64, get(ENV, "UVAL", "8.0"))
    for U in [U_first]
        println("\n", "="^70)
        @printf("W=%d x Lx=%d cylinder (N=%d qubits), U=%.1f, doping=%d/%d holes (Nel=%d)\n",
                W, LX, 2n, U, NHOLE, n, n - NHOLE)
        flush(stdout)

        tstart = time()
        BLAS.set_num_threads(4)
        E, ψ = ground_state_cyl(LX, W, t, U, 0.0; chi_max=CHI_EXP, nsweeps=14)
        BLAS.set_num_threads(1)
        @printf("  DMRG(chi=192): E0=%.6f  maxlinkdim=%d  (%.0fs)\n", E, maxlinkdim(ψ), time()-tstart)
        flush(stdout)

        uhf = solve_uhf_best(LX, W, t, U; Nup, Ndn)
        @printf("  UHF: E=%.6f  (E_UHF - E_DMRG = %.4f, staggered m=%.3f)\n",
                uhf.E, uhf.E - E, abs(uhf.nup[sidx(LX÷2,2)] - uhf.ndn[sidx(LX÷2,2)])/2)
        flush(stdout)

        exp_tens = extract_tensors(ψ)
        obs = build_observables_2d(LX, W)
        q1, q2 = obs_support(obs)
        obs_w = shift_obs(obs, q1 - 1)
        wtens, cum_pl = extract_window(ψ, q1, q2)
        @printf("  sampling window: qubits [%d, %d] (%d of %d)\n", q1, q2, q2-q1+1, 2n)

        Pρ = [[product_op_expect(exp_tens, term_matrixdict(tm)) for tm in o.terms] for o in obs]
        Otrue = [sum(tm.coeff * Pρ[k][ti] for (ti, tm) in enumerate(obs[k].terms)) for k in 1:length(obs)]

        # priors: UHF (Wick) + 対称性回復UHF (回転平均Wick) + 切断MPS
        Pσ_all = Vector{Vector{Vector{Float64}}}(); trOσ_all = Vector{Vector{Float64}}()
        fids = Float64[]
        push!(Pσ_all, [[wick_term(tm.sup, uhf) for tm in o.terms] for o in obs])
        push!(trOσ_all, [sum(tm.coeff * Pσ_all[1][k][ti] for (ti, tm) in enumerate(obs[k].terms)) for k in 1:length(obs)])
        push!(fids, NaN)   # UHFのグローバル忠実度は未計算 (局所Δで評価)
        push!(Pσ_all, symavg_P(obs, uhf, n))
        push!(trOσ_all, [sum(tm.coeff * Pσ_all[2][k][ti] for (ti, tm) in enumerate(obs[k].terms)) for k in 1:length(obs)])
        push!(fids, NaN)
        for chi_p in chi_priors
            σ = truncate(ψ; maxdim=chi_p); normalize!(σ)
            push!(fids, abs2(inner(σ, ψ)))
            st = extract_tensors(σ)
            Pσ = [[product_op_expect(st, term_matrixdict(tm)) for tm in o.terms] for o in obs]
            push!(Pσ_all, Pσ)
            push!(trOσ_all, [sum(tm.coeff * Pσ[k][ti] for (ti, tm) in enumerate(obs[k].terms)) for k in 1:length(obs)])
        end
        @printf("  MPS prior fidelities: %s\n",
                join([@sprintf("chi=%d: %.4f", c, f) for (c, f) in zip(chi_priors, fids[3:end])], ",  "))
        flush(stdout)

        tstart = time()
        est_std, est_crm = run_locals(wtens, cum_pl, obs_w, Pσ_all, trOσ_all;
                                      nu, nm, n_repeat, seed=2000 + round(Int, 10U))
        @printf("  local-observable run: %.0fs\n", time()-tstart); flush(stdout)

        @printf("\n  %-12s %9s |", "observable", "true")
        for lab in prior_labels; @printf(" %8s", "G($lab)"); end
        println()
        for (k, o) in enumerate(obs)
            v_std = var(est_std[:, k])
            @printf("  %-12s %9.4f |", o.name, Otrue[k])
            for (p, lab) in enumerate(prior_labels)
                G = v_std / var(est_crm[:, k, p])
                @printf(" %8.2f", G)
                Δ = Otrue[k] - trOσ_all[p][k]
                G_theo = NaN
                if o.pure
                    vs, vc = theory_var(length(o.terms[1].sup), Otrue[k], Δ, nu, nm)
                    G_theo = vs / vc
                end
                push!(rows, (U, lab, o.name, o.pure, Otrue[k], Δ, fids[p], G, G_theo))
            end
            println()
        end
        # Δの比較表 (UHF vs MPS priors)
        @printf("\n  %-12s |", "Delta")
        for lab in prior_labels; @printf(" %8s", lab); end
        println()
        for (k, o) in enumerate(obs)
            @printf("  %-12s |", o.name)
            for p in 1:length(prior_labels)
                @printf(" %8.4f", Otrue[k] - trOσ_all[p][k])
            end
            println()
        end
        flush(stdout)
    end

    utag = replace(string(U_first), "." => "p")
    out = joinpath(@__DIR__, "crm_2d_doped_results_U$(utag).tsv")
    open(out, "w") do io
        println(io, "U\tprior\tobservable\tpure\ttrue\tDelta\tprior_fid\tG_emp\tG_theo")
        for r in rows; println(io, join(r, "\t")); end
    end
    println("\nresults saved: $out"); flush(stdout)
    return rows, prior_labels
end

rows, prior_labels = main()

println("done (doped, NHOLE=$NHOLE, CHI_EXP=$CHI_EXP)")
