# ============================================================
# matchgate(フェルミオンGauss)シャドウ × CRM
#
# 動機: ランダムPauliシャドウでは、2DをJW変換した際のx方向ホッピング項に
#   長さ O(W) の演算子弦が付き、分散が 3^|A| で爆発して実質測定不能。
#   matchgateシャドウ (ランダムGauss回転 U_Q, Q∈SO(2n) + 占有数測定) では
#   ホッピングは2次Majorana演算子なので弦問題が消える。
#   さらに Slater行列式 (UHF) prior の厳密古典側が Pfaffian で計算できる
#   → CRMとの構造的相性が良い。
#
# 枠組み (Wan-Huggins-Lee-Babbush 2023 / Zhao-Rubin-Miyake PRL 2021):
#   チャネル固有値 λ_{2k} = C(n,k)/C(2n,2k) (次数2kのMajorana単項式)
#   スナップショット推定量:
#     X_S(b) = λ^{-1} Σ_{S'対角} det(Q[S,S']) <b|Γ_{S'}|b>,
#     <b|Γ_T|b> = i^{|T|} Π_{j∈T}(1-2b_j)
#   prior側の厳密条件付き平均 (CRM):
#     E[X_S|u] = λ^{-1} Σ_{S'} det(Q[S,S']) Pf(K'[S']),  K' = 回転後の2点関数
#
# 系: 4x2 ハバードシリンダー (8サイト, 16モード, 2^16次元の密ベクトルで厳密)
#   x方向ボンド (Pauliでは|A|=9) を含む。U=4, 8。
#   比較: Pauliシャドウ(標準/CRM) vs matchgateシャドウ(標準/CRM)、
#   個別観測量とエネルギー全体 (運動項 + U×二重占有)。
#
# 検証 (V1-V4, 8モード密行列で規約を機械精度確定):
#   V1: Givens分解ゲート列が U γ_μ U† = Σ_ν Q_{μν} γ_ν を実現
#   V2: 観測量のMajorana展開係数 = 密行列 (hop/docc/SzSz/長弦Pauli)
#   V3: シャドウ推定量の不偏性 (統計)
#   V4: K行列構成とPfaffian古典側 = 密行列総和
#
# 実行方法:
#  JULIA_LOAD_PATH="@:@v#.#:@stdlib" julia --project=Hubbard_MPS_Env_v2 crm_matchgate.jl
# ============================================================

using ITensors, ITensorMPS
using LinearAlgebra
using Statistics
using Printf
using Random

BLAS.set_num_threads(1)

# ------------------------------------------------------------
# 0. 格子 (4x2シリンダー) と規約
#    モード = JW qubit: qup(s)=2s-1, qdn(s)=2s / qubit 1 = MSB
#    Majorana: γ_{2m-1} = Z..Z X_m, γ_{2m} = Z..Z Y_m
# ------------------------------------------------------------
const W  = 4
const LX = 2
const NS = W * LX          # 8 sites
const NM = 2 * NS          # 16 modes
const NQ = NM              # 16 qubits
const DIM = 2^NQ

sidx(x, y) = (x - 1) * W + y
qup(s) = 2s - 1
qdn(s) = 2s
bitof(s::Int, q::Int, nq::Int) = (s >> (nq - q)) & 1

function cyl_edges(Lx, Wd)
    edges = Tuple{Int,Int}[]
    for x in 1:Lx, y in 1:Wd
        s = sidx(x, y)
        x < Lx && push!(edges, (s, sidx(x+1, y)))
        push!(edges, (s, sidx(x, y == Wd ? 1 : y + 1)))
    end
    return edges
end

# ------------------------------------------------------------
# 1. 密ベクトル演算 (任意のNQに対応)
# ------------------------------------------------------------
# Z回転 exp(-i θ/2 Z_q)
function apply_zrot!(ψ::Vector{ComplexF64}, θ::Float64, q::Int, nq::Int)
    p0, p1 = cis(-θ/2), cis(θ/2)
    mask = 1 << (nq - q)
    @inbounds for s in 0:length(ψ)-1
        ψ[s+1] *= (s & mask == 0) ? p0 : p1
    end
end

# XX回転 exp(-i θ/2 X_q X_{q+1})
function apply_xxrot!(ψ::Vector{ComplexF64}, θ::Float64, q::Int, nq::Int)
    c, s_ = cos(θ/2), sin(θ/2)
    m = (1 << (nq - q)) | (1 << (nq - q - 1))
    @inbounds for s in 0:length(ψ)-1
        t = s ⊻ m
        if s < t
            a, b = ψ[s+1], ψ[t+1]
            ψ[s+1] = c*a - im*s_*b
            ψ[t+1] = c*b - im*s_*a
        end
    end
end

# Majorana平面 (μ, μ+1) の回転角θのゲートを状態に適用
#   μ奇数 (同モード内): γ_{2j-1}γ_{2j} = iZ_j -> exp(-θ/2 γγ) = exp(-iθ/2 Z_j)
#   μ偶数 (モード跨ぎ):  γ_{2j}γ_{2j+1} = iX_jX_{j+1} -> exp(-iθ/2 XX)
function apply_plane_rot!(ψ, μ::Int, θ::Float64, nq::Int)
    if isodd(μ)
        apply_zrot!(ψ, θ, (μ + 1) ÷ 2, nq)
    else
        apply_xxrot!(ψ, θ, μ ÷ 2, nq)
    end
end

# ------------------------------------------------------------
# 2. Haar SO(2n) と隣接Givens分解
#    分解: G_m ... G_1 Q = I (atan2で対角を正に) -> Q = G_1' ... G_m'
#    ゲート適用順序の規約は V1 で機械的に確定する
# ------------------------------------------------------------
function haar_so(m::Int)
    A = randn(m, m)
    F = qr(A)
    Q = Matrix(F.Q) * Diagonal(sign.(diag(F.R)))
    if det(Q) < 0
        Q[:, end] .*= -1
    end
    return Q
end

# 隣接平面回転で Q を I に潰す。返り値: [(μ, θ), ...] (潰した順)
function givens_decompose(Q::Matrix{Float64})
    m = size(Q, 1)
    R = copy(Q)
    rots = Tuple{Int,Float64}[]
    for col in 1:m-1
        for i in m-1:-1:col
            a, b = R[i, col], R[i+1, col]
            (abs(b) < 1e-14 && a > 0) && continue
            θ = atan(b, a)
            c, s = cos(θ), sin(θ)
            for k in col:m
                ri, rj = R[i, k], R[i+1, k]
                R[i, k]   =  c*ri + s*rj
                R[i+1, k] = -s*ri + c*rj
            end
            push!(rots, (i, θ))
        end
    end
    @assert maximum(abs.(R - I)) < 1e-9 "Givens decomposition failed"
    return rots
end

# rots から状態への適用 (順序・符号の規約は V1 で決定した ORDER[] を使う)
const ORDER = Ref(:reverse_negate)   # V1で確定させる
function apply_Q!(ψ, rots, nq::Int)
    if ORDER[] == :reverse_negate
        for k in length(rots):-1:1
            μ, θ = rots[k]; apply_plane_rot!(ψ, μ, -θ, nq)
        end
    elseif ORDER[] == :forward_negate
        for k in 1:length(rots)
            μ, θ = rots[k]; apply_plane_rot!(ψ, μ, -θ, nq)
        end
    elseif ORDER[] == :reverse
        for k in length(rots):-1:1
            μ, θ = rots[k]; apply_plane_rot!(ψ, μ, θ, nq)
        end
    else
        for k in 1:length(rots)
            μ, θ = rots[k]; apply_plane_rot!(ψ, μ, θ, nq)
        end
    end
end

# ------------------------------------------------------------
# 3. Majorana観測量 (係数はV2で密行列と照合)
#    monomial = (係数::ComplexF64, S::Vector{Int})  Sはソート済みMajorana添字
#    <b|Γ_T|b> = i^|T| Π (1-2b_j)  (Tは対角ペア集合)
# ------------------------------------------------------------
struct MObs
    name::String
    mons::Vector{Tuple{ComplexF64,Vector{Int}}}   # 次数0項は S=[]
end

# c†_a c_b + c†_b c_a (モードa<b) = (i/2)γ_{2a-1}γ_{2b} - (i/2)γ_{2a}γ_{2b-1}
hop_mons(a, b) = [(0.5im, [2a-1, 2b]), (-0.5im, [2a, 2b-1])]
# n_m = 1/2 + (i/2)γ_{2m-1}γ_{2m}
dens_mons(m) = [(0.5+0im, Int[]), (0.5im, [2m-1, 2m])]
# n_a n_b (a<b, 異モード)
docc_mons(a, b) = [(0.25+0im, Int[]), (0.25im, [2a-1, 2a]), (0.25im, [2b-1, 2b]),
                   (-0.25+0im, [2a-1, 2a, 2b-1, 2b])]
# S^z_i S^z_j = (1/4)(n_iu - n_id)(n_ju - n_jd)
function szsz_mons(si, sj)
    mons = Tuple{ComplexF64,Vector{Int}}[]
    for (ma, sa) in ((qup(si), 1.0), (qdn(si), -1.0)),
        (mb, sb) in ((qup(sj), 1.0), (qdn(sj), -1.0))
        c = 0.25 * sa * sb
        for (cm, S) in docc_mons(min(ma,mb), max(ma,mb))
            push!(mons, (c*cm, S))
        end
    end
    return mons
end

mval(T_modes::Vector{Int}, b::Vector{Int}) =
    (im)^length(T_modes) * prod(1 - 2*b[j] for j in T_modes; init=1)

# 次数2k係数
lam(n, k) = binomial(n, k) / binomial(2n, 2k)

# 対角ペア集合の列挙 (次数2: 各モード / 次数4: モードペア)
function diag_subsets(n, k)
    if k == 1
        return [[j] for j in 1:n]
    else
        return [[j1, j2] for j1 in 1:n-1 for j2 in j1+1:n]
    end
end

majorana_of(T::Vector{Int}) = reduce(vcat, [[2j-1, 2j] for j in T]; init=Int[])

det2(M) = M[1,1]*M[2,2] - M[1,2]*M[2,1]
function subdet(Q, S::Vector{Int}, Sp::Vector{Int})
    k = length(S)
    if k == 2
        return Q[S[1],Sp[1]]*Q[S[2],Sp[2]] - Q[S[1],Sp[2]]*Q[S[2],Sp[1]]
    else
        return det(@view Q[S, Sp])
    end
end

pf2(A) = A[1,2]
pf4(A) = A[1,2]*A[3,4] - A[1,3]*A[2,4] + A[1,4]*A[2,3]

# u = (Q, rots) に対する1観測量の推定準備:
#   各単項式について w[T'] = λ^{-1} det(Q[S, majorana(T')]) を前計算
struct ObsEval
    consts::Float64                            # 次数0項の合計 (実部)
    packs::Vector{Tuple{ComplexF64,Int,Vector{Vector{Int}},Vector{Float64}}}
    # (係数, k, 対角モード集合リスト, w)
end

function prepare_obs(o::MObs, Q::Matrix{Float64}, n::Int, dsub::Dict{Int,Vector{Vector{Int}}})
    consts = 0.0
    packs = Tuple{ComplexF64,Int,Vector{Vector{Int}},Vector{Float64}}[]
    for (c, S) in o.mons
        if isempty(S)
            consts += real(c)
            continue
        end
        k = length(S) ÷ 2
        subs = dsub[k]
        w = Vector{Float64}(undef, length(subs))
        for (idx, T) in enumerate(subs)
            w[idx] = subdet(Q, S, majorana_of(T)) / lam(n, k)
        end
        push!(packs, (c, k, subs, w))
    end
    return ObsEval(consts, packs)
end

function eval_obs(oe::ObsEval, b::Vector{Int})
    x = oe.consts + 0im
    for (c, k, subs, w) in oe.packs
        acc = 0.0 + 0im
        for (idx, T) in enumerate(subs)
            acc += w[idx] * mval(T, b)
        end
        x += c * acc
    end
    return real(x)
end

# prior側の厳密条件付き平均: <Γ_T'>_{σ'} = Pf(K'[majorana(T')])
function exact_obs_mean(oe::ObsEval, Kp::Matrix{ComplexF64})
    x = oe.consts + 0im
    for (c, k, subs, w) in oe.packs
        acc = 0.0 + 0im
        for (idx, T) in enumerate(subs)
            Sm = majorana_of(T)
            pf = k == 1 ? pf2(view(Kp, Sm, Sm)) : pf4(Kp[Sm, Sm])
            acc += w[idx] * pf
        end
        x += c * acc
    end
    return real(x)
end

# Tr[Oσ] (回転なし)
function trace_obs(o::MObs, K::Matrix{ComplexF64})
    x = 0.0 + 0im
    for (c, S) in o.mons
        if isempty(S)
            x += c
        elseif length(S) == 2
            x += c * pf2(view(K, S, S))
        else
            x += c * pf4(K[S, S])
        end
    end
    return real(x)
end

# 数保存・実C行列 (フルモード) から K_{μν} = <γ_μ γ_ν> (μ≠ν, 対角0) を構成
function K_from_C(C::Matrix{Float64})
    n = size(C, 1)
    K = zeros(ComplexF64, 2n, 2n)
    for a in 1:n, b in 1:n
        AA = C[a,b] - C[b,a]                       # <A_a A_b> (a≠b)
        BB = C[a,b] - C[b,a]                       # <B_a B_b> (a≠b)
        AB = im * ((a == b ? 1.0 : 0.0) - C[a,b] - C[b,a])   # <A_a B_b>
        BA = -im * ((a == b ? 1.0 : 0.0) - C[a,b] - C[b,a])  # <B_a A_b> = -<A_b B_a>... V4で検証
        if a != b
            K[2a-1, 2b-1] = AA
            K[2a,   2b]   = BB
        end
        K[2a-1, 2b] = AB
        K[2a, 2b-1] = BA          # a==b でも <γ_{2a}γ_{2a-1}> = -i(1-2n) は非零
    end
    for μ in 1:2n; K[μ, μ] = 0; end
    return K
end

# ------------------------------------------------------------
# 4. 検証 (8モード密行列)
# ------------------------------------------------------------
function dense_paulis(nq)
    I2 = ComplexF64[1 0; 0 1]; X = ComplexF64[0 1; 1 0]
    Y = ComplexF64[0 -im; im 0]; Z = ComplexF64[1 0; 0 -1]
    return I2, X, Y, Z
end

function dense_gamma(nq)
    I2, X, Y, Z = dense_paulis(nq)
    γ = Vector{Matrix{ComplexF64}}(undef, 2nq)
    for j in 1:nq
        opsX = [q < j ? Z : (q == j ? X : I2) for q in 1:nq]
        opsY = [q < j ? Z : (q == j ? Y : I2) for q in 1:nq]
        γ[2j-1] = reduce(kron, opsX)
        γ[2j]   = reduce(kron, opsY)
    end
    return γ
end

function dense_mobs(o::MObs, γ, nq)
    M = zeros(ComplexF64, 2^nq, 2^nq)
    for (c, S) in o.mons
        term = Matrix{ComplexF64}(I, 2^nq, 2^nq)
        for μ in S; term = term * γ[μ]; end
        M .+= c .* term
    end
    return M
end

function dense_cdag(nq)
    cdag = [0.0 0; 1 0]; I2 = Matrix(1.0I, 2, 2); F = [1.0 0; 0 -1]
    return [reduce(kron, [q < m ? F : (q == m ? cdag : I2) for q in 1:nq]) for m in 1:nq]
end

function validate_matchgate()
    nq = 8; n = nq; dim = 2^nq
    Random.seed!(11)
    γ = dense_gamma(nq)
    I2, X, Y, Z = dense_paulis(nq)

    # --- V1: ゲート列の規約決定 ---
    Q = haar_so(2n)
    rots = givens_decompose(Q)
    best = :none
    for ord in (:reverse_negate, :forward_negate, :reverse, :forward)
        ORDER[] = ord
        # 密Uを構成: 基底ベクトルにapply_Q!を適用
        U = zeros(ComplexF64, dim, dim)
        for col in 1:dim
            e = zeros(ComplexF64, dim); e[col] = 1
            apply_Q!(e, rots, nq)
            U[:, col] = e
        end
        d = maximum(maximum(abs.(U * γ[μ] * U' .- sum(Q[μ, ν] .* γ[ν] for ν in 1:2n)))
                    for μ in 1:2n)
        if d < 1e-9
            best = ord
            @printf("  [V1] gate convention: %s  (U γ U† = Σ Q γ, maxdiff %.1e)\n", ord, d)
            break
        end
    end
    @assert best != :none "V1 failed: no gate ordering matches Q"
    ORDER[] = best

    # --- V2: 観測量のMajorana展開 = フェルミオン密行列 ---
    Cd = dense_cdag(nq); Cop = [Matrix(m') for m in Cd]
    tests = [
        ("hop(1,3)", MObs("h", hop_mons(1, 3)), Cd[1]*Cop[3] .+ Cd[3]*Cop[1]),
        ("hop(2,7)", MObs("h", hop_mons(2, 7)), Cd[2]*Cop[7] .+ Cd[7]*Cop[2]),
        ("dens(3)",  MObs("n", dens_mons(3)),   Cd[3]*Cop[3]),
        ("docc(1,2)",MObs("d", docc_mons(1, 2)), (Cd[1]*Cop[1])*(Cd[2]*Cop[2])),
        ("SzSz(1,2)",MObs("s", szsz_mons(1, 2)),
         0.25 .* (Cd[1]*Cop[1] .- Cd[2]*Cop[2]) * (Cd[3]*Cop[3] .- Cd[4]*Cop[4])),
    ]
    for (nm_, o, Mref) in tests
        d = maximum(abs.(dense_mobs(o, γ, nq) .- Mref))
        @assert d < 1e-10 "V2 failed: $nm_ ($d)"
    end
    println("  [V2] Majorana expansions match dense fermion operators")

    # Pauli長弦 = hop の確認 (x-bond相当): (XZ..ZX + YZ..ZY)/2 on qubits 2..7
    Zs = [q in 3:6 ? Z : I2 for q in 1:nq]
    XZX = reduce(kron, [q == 2 || q == 7 ? X : Zs[q] for q in 1:nq])
    YZY = reduce(kron, [q == 2 || q == 7 ? Y : Zs[q] for q in 1:nq])
    d = maximum(abs.(0.5 .* (XZX .+ YZY) .- (Cd[2]*Cop[7] .+ Cd[7]*Cop[2])))
    @assert d < 1e-10 "V2 failed: long-string Pauli hop"
    println("  [V2] long JW-string Pauli = fermionic hop confirmed")

    # --- V3: 不偏性 (ランダム状態, 統計テスト) ---
    ψ = randn(ComplexF64, dim); ψ ./= norm(ψ)
    obs3 = [MObs("hop13", hop_mons(1, 3)), MObs("hop27", hop_mons(2, 7)),
            MObs("docc12", docc_mons(1, 2)), MObs("szsz12", szsz_mons(1, 2))]
    truth = [real(dot(ψ, dense_mobs(o, γ, nq) * ψ)) for o in obs3]
    dsub = Dict(1 => diag_subsets(n, 1), 2 => diag_subsets(n, 2))
    nu3, nm3 = 4000, 20
    acc = zeros(length(obs3)); acc2 = zeros(length(obs3))
    b = zeros(Int, nq)
    for _ in 1:nu3
        Q3 = haar_so(2n); r3 = givens_decompose(Q3)
        φ = copy(ψ); apply_Q!(φ, r3, nq)
        probs = abs2.(φ); cp = cumsum(probs)
        oes = [prepare_obs(o, Q3, n, dsub) for o in obs3]
        m_o = zeros(length(obs3))
        for _ in 1:nm3
            so_ = searchsortedfirst(cp, rand()) - 1
            for q in 1:nq; b[q] = bitof(so_, q, nq); end
            for (k, oe) in enumerate(oes)
                m_o[k] += eval_obs(oe, b)
            end
        end
        acc .+= m_o ./ nm3
        acc2 .+= (m_o ./ nm3).^2
    end
    for (k, o) in enumerate(obs3)
        est = acc[k]/nu3
        se = sqrt(max(acc2[k]/nu3 - est^2, 1e-12) / nu3)
        z = abs(est - truth[k]) / se
        @printf("  [V3] %-8s truth=%8.4f  est=%8.4f ± %.4f  (z=%.1f)\n",
                o.name, truth[k], est, se, z)
        @assert z < 5 "V3 failed: $(o.name) biased"
    end
    println("  [V3] matchgate shadow estimators unbiased")

    # --- V4: K行列とPfaffian古典側 ---
    # ランダムSlater行列式 (4粒子)
    Φ = Matrix(qr(randn(n, 4)).Q)
    σv = zeros(ComplexF64, dim); σv[1] = 1
    for kk in 4:-1:1
        σv = sum(Φ[m, kk] .* (Cd[m] * σv) for m in 1:n)
    end
    σv ./= norm(σv)
    Cσ = Φ * Φ'
    K = K_from_C(Cσ)
    dK = maximum(abs(K[μ, ν] - (μ == ν ? 0 : dot(σv, γ[μ] * γ[ν] * σv)))
                 for μ in 1:2n, ν in 1:2n)
    @assert dK < 1e-9 "V4 failed: K matrix ($dK)"
    println("  [V4] K matrix from C matches dense <γγ>")
    # 数保存Slaterで <X Z..Z X> = <Y Z..Z Y> = 2C_ab (Pauli-CRM prior側の項別式)
    I2v, Xv, Yv, Zv = dense_paulis(nq)
    for (a, bq) in ((2, 7), (1, 3))
        mid = [q in a+1:bq-1 ? Zv : I2v for q in 1:nq]
        XZXv = reduce(kron, [q == a || q == bq ? Xv : mid[q] for q in 1:nq])
        YZYv = reduce(kron, [q == a || q == bq ? Yv : mid[q] for q in 1:nq])
        @assert abs(dot(σv, XZXv*σv) - 2Cσ[a,bq]) < 1e-9 "V4: <XZX> ≠ 2C ($a,$bq)"
        @assert abs(dot(σv, YZYv*σv) - 2Cσ[a,bq]) < 1e-9 "V4: <YZY> ≠ 2C ($a,$bq)"
    end
    println("  [V4] per-term <XZ..ZX> = <YZ..ZY> = 2C confirmed on Slater state")
    # 厳密条件付き平均 vs 全数和
    for trial in 1:3
        Q4 = haar_so(2n); r4 = givens_decompose(Q4)
        φσ = copy(σv); apply_Q!(φσ, r4, nq)
        probs = abs2.(φσ)
        Kp = transpose(Q4) * K * Q4      # K' = <γγ>_{UσU†}
        for o in obs3
            oe = prepare_obs(o, Q4, n, dsub)
            brute = 0.0
            for sidx_ in 0:dim-1
                probs[sidx_+1] < 1e-14 && continue
                bb = [bitof(sidx_, q, nq) for q in 1:nq]
                brute += probs[sidx_+1] * eval_obs(oe, bb)
            end
            ex = exact_obs_mean(oe, Kp)
            @assert abs(brute - ex) < 1e-8 "V4 failed: $(o.name) trial $trial ($brute vs $ex)"
        end
    end
    println("  [V4] Pfaffian exact prior mean = brute-force sum")
    println("  [validate] all matchgate validations passed"); flush(stdout)
end

# ------------------------------------------------------------
# 5. Pauliシャドウ側 (16量子ビット密ベクトル; 既存規約と同一)
# ------------------------------------------------------------
const SIGMA_P = (ComplexF64[0 1; 1 0], ComplexF64[0 -im; im 0], ComplexF64[1 0; 0 -1])
const UB = (ComplexF64[1 1; 1 -1]/sqrt(2), ComplexF64[1 -im; 1 im]/sqrt(2),
            ComplexF64[1 0; 0 1])
const DV = ntuple(β -> ntuple(a -> real.(diag(UB[β]*SIGMA_P[a]*UB[β]')), 3), 3)

function apply_1q!(ψ::Vector{ComplexF64}, u::Matrix{ComplexF64}, q::Int, nq::Int)
    mask = 1 << (nq - q)
    @inbounds for s in 0:length(ψ)-1
        if s & mask == 0
            a, b = ψ[s+1], ψ[s+mask+1]
            ψ[s+1]      = u[1,1]*a + u[1,2]*b
            ψ[s+mask+1] = u[2,1]*a + u[2,2]*b
        end
    end
end

struct PTerm
    coeff::Float64
    sup::Vector{Tuple{Int,Int}}
end
struct PObs
    name::String
    terms::Vector{PTerm}
end

# hop(モードa<b) のPauli表現: (X Z..Z X + Y Z..Z Y)/2 on qubits a..b
function hop_pauli(a, b)
    mid = [(q, 3) for q in a+1:b-1]
    return [PTerm(0.5, vcat([(a,1)], mid, [(b,1)])),
            PTerm(0.5, vcat([(a,2)], mid, [(b,2)]))]
end
docc_pauli(a, b) = [PTerm(0.25, Tuple{Int,Int}[]), PTerm(-0.25, [(a,3)]),
                    PTerm(-0.25, [(b,3)]), PTerm(0.25, [(a,3),(b,3)])]
function szsz_pauli(si, sj)
    t = PTerm[]
    for (qa, sa) in ((qup(si),1.0),(qdn(si),-1.0)), (qb, sb) in ((qup(sj),1.0),(qdn(sj),-1.0))
        push!(t, PTerm(0.25*sa*sb*0.25, Tuple{Int,Int}[]))
        push!(t, PTerm(-0.25*sa*sb*0.25, [(qa,3)]))
        push!(t, PTerm(-0.25*sa*sb*0.25, [(qb,3)]))
        push!(t, PTerm(0.25*sa*sb*0.25, [(qa,3),(qb,3)]))
    end
    return t
end

function eval_pauli(o::PObs, basis::Vector{Int}, b::Vector{Int})
    x = 0.0
    for t in o.terms
        v = t.coeff
        for (q, a) in t.sup
            v *= 3.0 * DV[basis[q]][a][b[q]+1]
            v == 0.0 && break
        end
        x += v
    end
    return x
end

# prior側 (UHF, Wick): 項ごとの <P>_σ を事前計算しておき基底一致で流す
function pauli_prior_means(o::PObs, C::Matrix{Float64})
    # 各項のσ期待値 (数保存Slater): Z列はK_from_C経由でなく直接Wick
    n = size(C, 1)
    vals = Float64[]
    for t in o.terms
        if isempty(t.sup)
            push!(vals, 1.0)
        elseif all(a == 3 for (_, a) in t.sup)
            qs = [q for (q, _) in t.sup]
            if length(qs) == 1
                push!(vals, 1 - 2C[qs[1], qs[1]])
            else
                a, b = qs
                push!(vals, (1-2C[a,a])*(1-2C[b,b]) - 4*C[a,b]*C[b,a] + 0.0)
            end
        else
            # (X Z..Z X) or (Y Z..Z Y): = c†c + h.c. ± pairing -> 数保存で 2C_ab
            a = t.sup[1][1]; b = t.sup[end][1]
            push!(vals, 2 * C[a, b])
        end
    end
    return vals
end

function exact_pauli_mean(o::PObs, basis::Vector{Int}, vals::Vector{Float64})
    m = 0.0
    for (ti, t) in enumerate(o.terms)
        ok = all(basis[q] == a for (q, a) in t.sup)
        ok || continue
        m += t.coeff * 3.0^length(t.sup) * vals[ti]
    end
    return m
end

# ------------------------------------------------------------
# 6. 状態構築 (DMRG -> 密ベクトル) と UHF
# ------------------------------------------------------------
function hubbard_cyl_mpo(sites, t, U, mu)
    os = OpSum()
    for (s, s2) in cyl_edges(LX, W)
        for (qa, qb) in ((qup(s), qup(s2)), (qdn(s), qdn(s2)))
            os += -t, "Cdag", qa, "C", qb
            os += -t, "Cdag", qb, "C", qa
        end
    end
    for s in 1:NS
        os += U, "N", qup(s), "N", qdn(s)
        os += -mu, "N", qup(s)
        os += -mu, "N", qdn(s)
    end
    return MPO(os, sites)
end

function ground_state_dense(t, U, mu)
    sites = siteinds("Fermion", NQ; conserve_qns=true)
    H = hubbard_cyl_mpo(sites, t, U, mu)
    init = fill("Emp", NQ)
    for x in 1:LX, y in 1:W
        s = sidx(x, y)
        init[iseven(x + y) ? qup(s) : qdn(s)] = "Occ"
    end
    ψ0 = productMPS(sites, init)
    E, ψ = dmrg(H, ψ0; nsweeps=10, maxdim=[20,50,100,256,256,256,256,256,256,256],
                cutoff=1e-12, noise=[1e-6,1e-7,1e-8,0,0,0,0,0,0,0], outputlevel=0)
    # MPS -> 密ベクトル (qubit1 = MSB)
    ψo = orthogonalize(ψ, 1)
    ss = siteinds(ψo); ls = linkinds(ψo)
    V = ones(ComplexF64, 1, 1)
    for j in 1:NQ
        A = if j == 1
            reshape(Array(ψo[j], ss[j], ls[j]), 2, 1, :)
        elseif j == NQ
            reshape(Array(ψo[j], ls[j-1], ss[j]), :, 2, 1) |> x -> permutedims(x, (2,1,3))
        else
            permutedims(Array(ψo[j], ls[j-1], ss[j], ls[j]), (2,1,3))
        end
        # A[b, l, r]; V[prefix, l] -> V'[(prefix,b), r]
        np, nl = size(V); nr = size(A, 3)
        Vn = zeros(ComplexF64, np*2, nr)
        for bb in 0:1
            Vn[(bb*np+1):(bb+1)*np, :] = V * A[bb+1, :, :]
        end
        # ここで prefix の桁: 新しいビットが下位に付く (qubit jはMSB側からj番目)
        V = Vn
    end
    vec_ = vec(V)
    # 上のスタッキングは bit_j が「ブロック上位」なので index = b1 + 2 b2 + ... の順になっている
    # 並べ替え: 実際は index = Σ b_j 2^{j-1} (qubit1が最下位) になっているため反転する
    ψd = zeros(ComplexF64, DIM)
    for s in 0:DIM-1
        s2 = 0
        for q in 1:NQ
            s2 |= ((s >> (q-1)) & 1) << (NQ - q)
        end
        ψd[s2+1] = vec_[s+1]
    end
    return E, ψd
end

function solve_uhf_cyl(t, U; Nup, Ndn, iters=800, tol=1e-10)
    T = zeros(NS, NS)
    for (s, s2) in cyl_edges(LX, W)
        T[s, s2] = T[s2, s] = -t
    end
    nup = zeros(NS); ndn = zeros(NS)
    for x in 1:LX, y in 1:W
        s = sidx(x, y)
        nup[s] = 0.5 + 0.4*(-1)^(x+y); ndn[s] = 0.5 - 0.4*(-1)^(x+y)
    end
    nup .*= Nup/sum(nup); ndn .*= Ndn/sum(ndn)
    local Fu, Fd
    for _ in 1:iters
        Fu = eigen(Symmetric(T + diagm(U .* ndn)))
        Fd = eigen(Symmetric(T + diagm(U .* nup)))
        nu_ = vec(sum(abs2, Fu.vectors[:, 1:Nup]; dims=2))
        nd_ = vec(sum(abs2, Fd.vectors[:, 1:Ndn]; dims=2))
        if max(maximum(abs.(nu_ .- nup)), maximum(abs.(nd_ .- ndn))) < tol
            nup, ndn = nu_, nd_; break
        end
        nup = 0.5 .* nu_ .+ 0.5 .* nup; ndn = 0.5 .* nd_ .+ 0.5 .* ndn
    end
    Cu = Fu.vectors[:, 1:Nup] * Fu.vectors[:, 1:Nup]'
    Cd_ = Fd.vectors[:, 1:Ndn] * Fd.vectors[:, 1:Ndn]'
    # フルモードC (mode qup(s), qdn(s))
    C = zeros(NM, NM)
    for a in 1:NS, b in 1:NS
        C[qup(a), qup(b)] = Cu[a, b]
        C[qdn(a), qdn(b)] = Cd_[a, b]
    end
    return C
end

# ------------------------------------------------------------
# 7. メイン実験
# ------------------------------------------------------------
function main()
    println("=== validation (8 modes, dense) ==="); flush(stdout)
    validate_matchgate()

    t = 1.0
    nu, nm, n_repeat = 50, 100, 40
    dsub = Dict(1 => diag_subsets(NM, 1), 2 => diag_subsets(NM, 2))
    rows = []

    for U in [4.0, 8.0]
        println("\n", "="^70)
        @printf("4x2 cylinder (16 qubits), U=%.1f, half-filling\n", U); flush(stdout)
        BLAS.set_num_threads(4)
        E0, ψd = ground_state_dense(t, U, U/2)
        BLAS.set_num_threads(1)
        C = solve_uhf_cyl(t, U; Nup=NS÷2, Ndn=NS÷2)
        K = K_from_C(C)
        @printf("  DMRG E0=%.6f  (dense norm %.6f)\n", E0, norm(ψd)); flush(stdout)

        # --- 観測量セット ---
        s0 = sidx(1, 2)
        sy = s0 + 1               # y方向隣
        sx = s0 + W               # x方向隣 (次列)
        mob = [MObs("hop y-bond", hop_mons(qup(s0), qup(sy))),
               MObs("hop x-bond", hop_mons(qup(s0), qup(sx))),
               MObs("DoubleOcc",  docc_mons(qup(s0), qdn(s0))),
               MObs("SzSz y",     szsz_mons(s0, sy))]
        pob = [PObs("hop y-bond", hop_pauli(qup(s0), qup(sy))),
               PObs("hop x-bond", hop_pauli(qup(s0), qup(sx))),
               PObs("DoubleOcc",  docc_pauli(qup(s0), qdn(s0))),
               PObs("SzSz y",     szsz_pauli(s0, sy))]
        # エネルギー (運動項 + U*docc; μ項は粒子数固定なので定数)
        emons = Tuple{ComplexF64,Vector{Int}}[]
        eterms = PTerm[]
        for (s, s2) in cyl_edges(LX, W), off in (qup, qdn)
            a, b = min(off(s), off(s2)), max(off(s), off(s2))
            for (c, S) in hop_mons(a, b); push!(emons, (-t*c, S)); end
            for pt in hop_pauli(a, b); push!(eterms, PTerm(-t*pt.coeff, pt.sup)); end
        end
        for s in 1:NS
            for (c, S) in docc_mons(qup(s), qdn(s)); push!(emons, (U*c, S)); end
            for pt in docc_pauli(qup(s), qdn(s)); push!(eterms, PTerm(U*pt.coeff, pt.sup)); end
        end
        push!(mob, MObs("Energy", emons))
        push!(pob, PObs("Energy", eterms))
        nobs = length(mob)

        # 真値 (密ベクトル): Pauli側の項評価で計算 (MajoranaとPauliは同一演算子)
        truth = zeros(nobs)
        for (k, o) in enumerate(pob)
            v = 0.0
            for t_ in o.terms
                if isempty(t_.sup)
                    v += t_.coeff
                else
                    # <ψ| Pauli列 |ψ> を密ベクトルで
                    φ = copy(ψd)
                    for (q, a) in t_.sup
                        mask = 1 << (NQ - q)
                        if a == 3
                            @inbounds for s_ in 0:DIM-1
                                (s_ & mask != 0) && (φ[s_+1] *= -1)
                            end
                        else
                            tmp = zeros(ComplexF64, DIM)
                            @inbounds for s_ in 0:DIM-1
                                s2 = s_ ⊻ mask
                                # X: 1 / Y: <1|Y|0>=+i, <0|Y|1>=-i
                                ph = a == 1 ? 1.0+0im : ((s_ & mask == 0) ? im : -im)
                                tmp[s2+1] += ph * φ[s_+1]
                            end
                            φ = tmp
                        end
                    end
                    v += t_.coeff * real(dot(ψd, φ))
                end
            end
            truth[k] = v
        end
        # エンドツーエンド検証: エネルギー真値は E0 + μN (μ=U/2, N=8) に一致するはず
        @assert abs(truth[end] - (E0 + U/2 * NS)) < 1e-6 "energy truth check failed: $(truth[end]) vs $(E0 + U/2*NS)"
        trOσ_m = [trace_obs(o, K) for o in mob]
        pvals  = [pauli_prior_means(o, C) for o in pob]
        trOσ_p = [sum(o.terms[ti].coeff * pvals[k][ti] for ti in 1:length(o.terms))
                  for (k, o) in enumerate(pob)]
        # PfaffianルートとWickルートの一致 (K_from_C の交差検証)
        for k in 1:nobs
            @assert abs(trOσ_m[k] - trOσ_p[k]) < 1e-8 "Tr[Oσ] mismatch: $(mob[k].name)"
        end
        @printf("  %-12s %9s %9s %9s\n", "observable", "true", "Tr[Oσ]_MG", "Tr[Oσ]_P")
        for k in 1:nobs
            @printf("  %-12s %9.4f %9.4f %9.4f\n", mob[k].name, truth[k], trOσ_m[k], trOσ_p[k])
        end
        flush(stdout)

        # --- 実験ループ ---
        est = zeros(n_repeat, nobs, 4)   # 1=P-std, 2=P-CRM, 3=MG-std, 4=MG-CRM
        b = zeros(Int, NQ)
        tstart = time()
        for rep in 1:n_repeat
            accP = zeros(nobs, 2); accM = zeros(nobs, 2)
            for _ in 1:nu
                # --- Pauli ---
                basis = rand(1:3, NQ)
                φ = copy(ψd)
                for q in 1:NQ; apply_1q!(φ, UB[basis[q]], q, NQ); end
                cp = cumsum(abs2.(φ))
                mP = zeros(nobs)
                for _ in 1:nm
                    s_ = searchsortedfirst(cp, rand()) - 1
                    for q in 1:NQ; b[q] = bitof(s_, q, NQ); end
                    for k in 1:nobs; mP[k] += eval_pauli(pob[k], basis, b); end
                end
                mP ./= nm
                for k in 1:nobs
                    accP[k,1] += mP[k]
                    accP[k,2] += mP[k] - exact_pauli_mean(pob[k], basis, pvals[k])
                end
                # --- matchgate ---
                Q = haar_so(2NM)
                rots = givens_decompose(Q)
                φ = copy(ψd); apply_Q!(φ, rots, NQ)
                cp = cumsum(abs2.(φ))
                oes = [prepare_obs(o, Q, NM, dsub) for o in mob]
                Kp = transpose(Q) * K * Q
                mM = zeros(nobs)
                for _ in 1:nm
                    s_ = searchsortedfirst(cp, rand()) - 1
                    for q in 1:NQ; b[q] = bitof(s_, q, NQ); end
                    for k in 1:nobs; mM[k] += eval_obs(oes[k], b); end
                end
                mM ./= nm
                for k in 1:nobs
                    accM[k,1] += mM[k]
                    accM[k,2] += mM[k] - exact_obs_mean(oes[k], Kp)
                end
            end
            for k in 1:nobs
                est[rep,k,1] = accP[k,1]/nu
                est[rep,k,2] = accP[k,2]/nu + trOσ_p[k]
                est[rep,k,3] = accM[k,1]/nu
                est[rep,k,4] = accM[k,2]/nu + trOσ_m[k]
            end
            if rep % 10 == 0
                @printf("  rep %d/%d (%.0fs)\n", rep, n_repeat, time()-tstart); flush(stdout)
            end
        end

        @printf("\n  %-12s %9s | %-19s %-19s | %-19s %-19s | %7s %7s\n",
                "observable", "true", "Pauli std", "Pauli CRM(UHF)",
                "MG std", "MG CRM(UHF)", "G_P", "G_MG")
        for k in 1:nobs
            m = [mean(est[:,k,e]) for e in 1:4]; s = [std(est[:,k,e]) for e in 1:4]
            @printf("  %-12s %9.4f | %8.3f ± %-8.3f %8.3f ± %-8.3f | %8.4f ± %-8.4f %8.4f ± %-8.4f | %7.2f %7.2f\n",
                    mob[k].name, truth[k], m[1], s[1], m[2], s[2], m[3], s[3], m[4], s[4],
                    (s[2]>0 ? (s[1]/s[2])^2 : NaN), (s[4]>0 ? (s[3]/s[4])^2 : NaN))
            push!(rows, (U, mob[k].name, truth[k], s[1], s[2], s[3], s[4]))
        end
        flush(stdout)
    end

    out = joinpath(@__DIR__, "crm_matchgate_results.tsv")
    open(out, "w") do io
        println(io, "U\tobservable\ttrue\terr_P_std\terr_P_crm\terr_MG_std\terr_MG_crm")
        for r in rows; println(io, join(r, "\t")); end
    end
    println("\nresults saved: $out"); flush(stdout)
    return rows
end

rows = main()

# ------------------------------------------------------------
# 8. プロット
# ------------------------------------------------------------
using Plots
const OI = ["#0072B2", "#D55E00", "#009E73", "#CC79A7"]
default(fontfamily="Helvetica", framestyle=:box, grid=true, gridalpha=0.12,
        linewidth=1.8, markersize=6, markerstrokewidth=0.6, markerstrokecolor=:white,
        guidefontsize=11, tickfontsize=9, legendfontsize=8, titlefontsize=11, dpi=200)

panels = []
for (iu, U) in enumerate((4.0, 8.0))
    sel = [r for r in rows if r[1] == U]
    names = [r[2] for r in sel]
    p = plot(xticks=(1:length(names), names), xrotation=25, yscale=:log10,
             ylabel=(iu == 1 ? "std error of estimate" : ""),
             legend=(iu == 1 ? :topleft : :none),
             title=@sprintf("4x2 cylinder,  U = %.0f", U))
    labels = ["Pauli shadow", "Pauli + CRM(UHF)", "matchgate shadow", "matchgate + CRM(UHF)"]
    marks = [:circle, :utriangle, :square, :diamond]
    for (e, lab) in enumerate(labels)
        errs = [sel[k][3+e] for k in 1:length(sel)]
        errs = [x <= 0 ? NaN : x for x in errs]   # std=0 (一度も測定されない項) は非表示
        plot!(p, 1:length(names), errs, color=OI[e], marker=marks[e], label=lab)
    end
    annotate!(p, 2, minimum(skipmissing([r for r in (sel[k][3+e] for k in 1:length(sel), e in 1:4) if r > 0])),
              text("Pauli: x-bond never sampled", 7, :gray40, :center))
    push!(panels, p)
end
final = plot(panels..., layout=(1,2), size=(1240, 470), margin=7Plots.mm)
outpng = joinpath(@__DIR__, "crm_matchgate.png")
savefig(final, outpng)
println("figure saved: $outpng")
