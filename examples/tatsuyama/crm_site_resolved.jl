# ============================================================
# CRM利得のサイト分解検証 (1Dハバード鎖, L=8/16/32/64)
#
# 動機:
#   README §2「局所性: グローバル忠実度が崩壊しても利得は生き残る」は
#   これまで鎖中央 i0=L/2 の観測量だけで示していた。中央は開放端鎖で
#   最もエンタングルメントが大きい＝切断priorが最も苦しい場所なので、
#   「中央で生き残るなら全域で生き残る」という主張は妥当に見えるが、
#   実際にサイト依存性を測っていない以上それは仮説にすぎない。
#
# 本スクリプトが答える問い:
#   Q1. 利得 G は鎖に沿って一様か? 端と中央でどう違うか?
#   Q2. サイト依存性は何で説明できるか?
#       予測: 切断priorの局所誤差 Δ(i) はボンドエンタングルメント
#             エントロピー S(i) が大きい場所ほど大きい → 中央で G 最小。
#   Q3. 2レジーム診断 (G/G_max) はサイトごとにどう分布するか?
#       G_max は「完全prior (σ=ρ, Δ=0)」を同じサンプルに適用して
#       経験的に測る。これなら多項観測量でも定義できる。
#
# 設計 (計算量):
#   - サイト i ごとに窓 [qup(i), qdn(i+1)] = 4 qubit のみをサンプル。
#     直交中心を窓左端に置く窓サンプリング（crm_mps_scaling.jl と同じ）
#     によりコストは L に依存しない。
#   - 実験側サンプルは全priorで共有（1回サンプル → 後処理のみprior毎）。
#   - 分散比には bootstrap 信頼区間を付ける（n_repeat 個の推定値を再標本化）。
#
# 実行:
#   JULIA_LOAD_PATH="@:@v#.#:@stdlib" julia --project=Hubbard_MPS_Env_v2 crm_site_resolved.jl
# 環境変数: L_LIST (既定 "8,16,32"), N_REPEAT (既定 300), CHI_EXP (既定 128)
# ============================================================

using ITensors, ITensorMPS
using LinearAlgebra
using Statistics
using Printf
using Random

BLAS.set_num_threads(parse(Int, get(ENV, "DMRG_BLAS", "1")))

qup(i) = 2i - 1
qdn(i) = 2i

# ------------------------------------------------------------
# 1. モデルと基底状態
# ------------------------------------------------------------
function hubbard_chain_mpo(sites, L, t, U, mu)
    os = OpSum()
    for i in 1:L-1
        os += -t, "Cdag", qup(i), "C", qup(i+1)
        os += -t, "Cdag", qup(i+1), "C", qup(i)
        os += -t, "Cdag", qdn(i), "C", qdn(i+1)
        os += -t, "Cdag", qdn(i+1), "C", qdn(i)
    end
    for i in 1:L
        os += U, "N", qup(i), "N", qdn(i)
        os += -mu, "N", qup(i)
        os += -mu, "N", qdn(i)
    end
    return MPO(os, sites)
end

"""δ=0 なら Néel 積状態（各サイト1電子＝ハーフフィリング）。
   δ>0 なら電子を round(δL) 個だけ等間隔に抜いてホールを入れる。
   QN保存DMRGでは初期状態が粒子数セクターを決めるので、これが doping の指定になる。"""
function initial_config(L, δ)
    occ = [isodd(div(q-1, 2) + 1) == isodd(q) for q in 1:2L]
    nh = round(Int, δ * L)
    if nh > 0
        occupied = findall(occ)
        # 端に寄せず鎖全体に等間隔で抜く
        pick = [occupied[round(Int, (k - 0.5) * length(occupied) / nh)] for k in 1:nh]
        occ[unique(pick)] .= false
    end
    return [o ? "Occ" : "Emp" for o in occ], count(occ)
end

function ground_state(L, t, U, mu; chi_max=128, nsweeps=12, δ=0.0)
    N = 2L
    sites = siteinds("Fermion", N; conserve_qns=true)
    H = hubbard_chain_mpo(sites, L, t, U, mu)
    init, nel = initial_config(L, δ)
    ψ0 = productMPS(sites, init)
    ramp = [20, 50, 100, 200, chi_max, chi_max, chi_max, chi_max,
            chi_max, chi_max, chi_max, chi_max]
    maxdims = min.(chi_max, ramp)[1:nsweeps]
    noise = [1e-6, 1e-6, 1e-7, 1e-8, 1e-9, 0, 0, 0, 0, 0, 0, 0][1:nsweeps]
    E, ψ = dmrg(H, ψ0; nsweeps, maxdim=maxdims, cutoff=1e-11, noise, outputlevel=0)
    return E, ψ, sites, H, nel
end

# 量子ビットボンド b (b と b+1 の間) のエンタングルメントエントロピー
function bond_entropies(ψ::MPS)
    N = length(ψ)
    ψo = copy(ψ); orthogonalize!(ψo, 1)
    S = zeros(N - 1)
    for b in 1:N-1
        orthogonalize!(ψo, b)
        inds_left = b == 1 ? (siteinds(ψo)[b],) : (linkinds(ψo)[b-1], siteinds(ψo)[b])
        _, Sv, _ = svd(ψo[b], inds_left)
        p = diag(Array(Sv, inds(Sv)...)).^2
        p = p[p .> 1e-14]
        S[b] = -sum(p .* log.(p))
    end
    return S
end

# ------------------------------------------------------------
# 2. MPS -> 密テンソル / 窓抽出
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

"""窓 [q1,q2] の縮約密度行列 ρ_w (2^Nw × 2^Nw, ビット s=(b_{q1}...b_{q2}), b_{q1} が最上位)。

直交中心を q1 に置くと、q1 より左は左正準・q1 より右は右正準になるので
左右の環境はともに恒等になる。よって窓のMPSブロック A[l,s,r] に対し
    ρ_w[s,s'] = Σ_{l,r} A[l,s,r] conj(A[l,s',r])
が厳密な縮約密度行列である。窓幅は 4 量子ビットなので ρ_w は 16×16 にすぎず、
以降のサンプリングはボンド次元 χ に一切依存しない（これが L=64 まで回せる理由）。"""
function rdm_from_window(tens)
    chain = [tens[1][1], tens[1][2]]                      # ビット列ごとの χl×χr 行列
    for j in 2:length(tens)
        A0, A1 = tens[j]
        chain = vcat([M * A0 for M in chain], [M * A1 for M in chain])
        # vcat の順序が「新ビットが最下位」になるよう並べ替える
        n = length(chain) ÷ 2
        chain = vec(permutedims(reshape(chain, n, 2), (2, 1)))
    end
    d = length(chain)
    ρ = Matrix{ComplexF64}(undef, d, d)
    for a in 1:d, b in 1:d
        ρ[a, b] = dot(chain[b], chain[a])                 # Σ_{l,r} M_a conj(M_b)
    end
    return Hermitian((ρ + ρ') / 2)
end

window_rdm(ψ::MPS, q1::Int, q2::Int) =
    (ψo = orthogonalize(ψ, q1); rdm_from_window([site_matrices(ψo, j) for j in q1:q2]))

"""窓の列（q1 昇順）に対する縮約密度行列を一度の前進掃引でまとめて作る。

各窓ごとに orthogonalize を呼び直すと直交中心の移動が O(N) 回起きて
L に対して二次のコストになる。中心を左から右へ一方向に動かせば全体で
O(N) 回の QR で済み、L=128 (256 量子ビット) でも軽い。"""
function window_rdms(ψ::MPS, windows)
    ψo = copy(ψ); orthogonalize!(ψo, first(windows)[1])
    out = Vector{Hermitian{ComplexF64,Matrix{ComplexF64}}}(undef, length(windows))
    for (k, (q1, q2)) in enumerate(windows)
        orthogonalize!(ψo, q1)
        out[k] = rdm_from_window([site_matrices(ψo, j) for j in q1:q2])
    end
    return out
end

expect_rdm(ρw, M) = real(tr(ρw * M))

"""基底 basis (1=X,2=Y,3=Z) での測定結果分布（窓の 2^Nw 個の確率）。"""
function basis_probs(ρw, basis::Vector{Int})
    Uw = reduce(kron, (UBASIS[b] for b in basis))
    return max.(real.(diag(Uw * ρw * Uw')), 0.0)
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

# ------------------------------------------------------------
# 3. シャドウ測定機構
# ------------------------------------------------------------
const SIGMA  = (ComplexF64[0 1; 1 0], ComplexF64[0 -im; im 0], ComplexF64[1 0; 0 -1])
const UBASIS = (ComplexF64[1 1; 1 -1]/sqrt(2), ComplexF64[1 -im; 1 im]/sqrt(2),
                ComplexF64[1 0; 0 1])
const DVAL = ntuple(β -> ntuple(a -> real.(diag(UBASIS[β]*SIGMA[a]*UBASIS[β]')), 3), 3)

"""確率ベクトルの累積和からビット列（b_{1} が最上位）を1つ引く。"""
function draw_bits!(bits::Vector{Int}, cum::Vector{Float64})
    idx = searchsortedfirst(cum, rand() * cum[end]) - 1
    Nw = length(bits)
    @inbounds for j in Nw:-1:1
        bits[j] = idx & 1
        idx >>= 1
    end
    return bits
end

# ------------------------------------------------------------
# 4. 観測量: サイト i を中心とする局所量一式
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

const I2C = Matrix{ComplexF64}(I, 2, 2)

"""窓 (1..Nw に添字を移した) 項 t に対応する 2^Nw × 2^Nw のPauli積行列。"""
term_window_matrix(t::Term, Nw::Int) = reduce(kron,
    ((a = findfirst(x -> x[1] == q, t.sup); a === nothing ? I2C : SIGMA[t.sup[a][2]])
     for q in 1:Nw))

"""サイト i (と、bond=true なら i+1 との結合) に載る観測量の一式。
   n_q = (1 - Z_q)/2 の JW 表現を用いる。"""
function site_observables(i::Int; bond::Bool)
    u, d = qup(i), qdn(i)
    obs = Obs[]
    # 電子数 n(i) = n_up + n_dn = 1 - (Z_u + Z_d)/2
    push!(obs, Obs("n", [Term(1.0, Tuple{Int,Int}[]),
                         Term(-0.5, [(u,3)]), Term(-0.5, [(d,3)])], false))
    # 局所磁化 Sz(i) = (n_up - n_dn)/2 = (Z_d - Z_u)/4
    push!(obs, Obs("Sz", [Term(-0.25, [(u,3)]), Term(0.25, [(d,3)])], false))
    # 二重占有 D(i) = n_up n_dn
    push!(obs, Obs("DoubleOcc", [Term(0.25, Tuple{Int,Int}[]), Term(-0.25, [(u,3)]),
                                 Term(-0.25, [(d,3)]), Term(0.25, [(u,3),(d,3)])], false))
    # 単一Pauli列 (理論式が閉形式で使える対照群)
    push!(obs, Obs("ZZ onsite", [Term(1.0, [(u,3),(d,3)])], true))
    if bond
        u2, d2 = qup(i+1), qdn(i+1)
        push!(obs, Obs("ZZ up-up r=1", [Term(1.0, [(u,3),(u2,3)])], true))
        push!(obs, Obs("SzSz r=1",
            [Term( 1/16, [(u,3),(u2,3)]), Term(-1/16, [(u,3),(d2,3)]),
             Term(-1/16, [(d,3),(u2,3)]), Term( 1/16, [(d,3),(d2,3)])], false))
        # ホッピング (up スピン): c†_u c_{u2} + h.c. = (X Z X + Y Z Y)/2
        push!(obs, Obs("hop up r=1",
            [Term(0.5, [(u,1),(d,3),(u2,1)]), Term(0.5, [(u,2),(d,3),(u2,2)])], false))
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
        all(basis[q] == a for (q, a) in t.sup) || continue
        m += t.coeff * 3.0^length(t.sup) * Pσ[ti]
    end
    return m
end

theory_var(absA, P, Δ, nu, nm) = (((3.0^absA - 1)*P^2 + 3.0^absA*(1-P^2)/nm)/nu,
                                  ((3.0^absA - 1)*Δ^2 + 3.0^absA*(1-P^2)/nm)/nu)

# ------------------------------------------------------------
# 5. サンプリング本体 (実験サンプルを全priorで共有)
# ------------------------------------------------------------
function run_locals(ρw, Nw, obs_w, Pσ_all, trOσ_all; nu, nm, n_repeat, seed)
    Random.seed!(seed)
    nobs = length(obs_w); np = length(Pσ_all)
    bits = zeros(Int, Nw)
    # 3^Nw 通りの基底は高々 3^4=81 通り。累積確率を事前計算しておけば
    # サンプリングは 1 回あたり二分探索のみになる。
    cum_cache = Dict{Vector{Int},Vector{Float64}}()
    est_std = zeros(n_repeat, nobs); est_crm = zeros(n_repeat, nobs, np)
    xs = zeros(nobs)
    for rep in 1:n_repeat
        acc_s = zeros(nobs); acc_c = zeros(nobs, np)
        for _ in 1:nu
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

"""分散比 var(a)/var(b) の bootstrap 信頼区間 (同一サンプルで対応あり)。"""
const NBOOT = parse(Int, get(ENV, "NBOOT", "600"))

function boot_ratio(a::Vector{Float64}, b::Vector{Float64}; nboot=NBOOT, q=(0.16, 0.84), rng=MersenneTwister(1))
    n = length(a)
    r = Vector{Float64}(undef, nboot)
    idx = Vector{Int}(undef, n)
    for t in 1:nboot
        rand!(rng, idx, 1:n)
        r[t] = var(@view a[idx]) / var(@view b[idx])
    end
    sort!(r)
    return (quantile(r, q[1]), quantile(r, q[2]))
end

# ------------------------------------------------------------
# 6. 検証 (L=4 密行列ED + 窓サンプラー検定)
# ------------------------------------------------------------
function dense_chain_check()
    L, t, U = 4, 1.0, 4.0; mu = U/2; N = 2L
    cdag = [0.0 0.0; 1.0 0.0]; I2 = Matrix(1.0I, 2, 2); Fm = [1.0 0.0; 0.0 -1.0]
    Cd = Vector{Matrix{Float64}}(undef, N)
    for k in 1:N
        Cd[k] = reduce(kron, [j < k ? Fm : (j == k ? cdag : I2) for j in 1:N])
    end
    C = [Matrix(m') for m in Cd]; Nop = [Cd[k]*C[k] for k in 1:N]
    H = zeros(2^N, 2^N)
    for i in 1:L-1, (a, b) in ((qup(i), qup(i+1)), (qdn(i), qdn(i+1)))
        H .+= -t .* (Cd[a]*C[b] .+ Cd[b]*C[a])
    end
    for i in 1:L
        H .+= U .* (Nop[qup(i)]*Nop[qdn(i)]) .- mu .* (Nop[qup(i)] .+ Nop[qdn(i)])
    end
    Fed = eigen(Symmetric(H)); E_ed = Fed.values[1]; ψ_ed = Fed.vectors[:, 1]

    E_dmrg, ψ, _, _, _ = ground_state(L, t, U, mu; chi_max=64, nsweeps=8)
    @printf("  [check] E_ED=%.8f  E_DMRG=%.8f  (diff %.2e)\n", E_ed, E_dmrg, abs(E_ed-E_dmrg))
    @assert abs(E_ed - E_dmrg) < 1e-6

    tens = extract_tensors(ψ)
    # 全サイトの観測量を ED と突き合わせる (サイト分解が主題なので全サイト検定)
    I2c = Matrix{ComplexF64}(I, 2, 2)
    dense_of(t::Term) = t.coeff * reduce(kron,
        [(a = findfirst(x -> x[1] == q, t.sup); a === nothing ? I2c : SIGMA[t.sup[a][2]]) for q in 1:N])
    maxerr = 0.0
    for i in 1:L
        for o in site_observables(i; bond = i < L)
            Mden = sum(dense_of(tm) for tm in o.terms)
            v_ed = real(dot(ψ_ed, Mden * ψ_ed))
            v_mps = sum(tm.coeff * product_op_expect(tens, term_matrixdict(tm)) for tm in o.terms)
            maxerr = max(maxerr, abs(v_ed - v_mps))
        end
    end
    @printf("  [check] all-site observables: max |ED - MPS| = %.2e\n", maxerr)
    @assert maxerr < 1e-6
    # ホッピングの符号規約が ED の c†c + h.c. と一致するか (別立てで明示)
    hop_ed = dot(ψ_ed, (Cd[qup(2)]*C[qup(3)] .+ Cd[qup(3)]*C[qup(2)])*ψ_ed)
    o_hop = site_observables(2; bond=true)[end]
    hop_mps = sum(tm.coeff * product_op_expect(tens, term_matrixdict(tm)) for tm in o_hop.terms)
    @printf("  [check] hop: ED=%.6f  (XZX+YZY)/2=%.6f\n", hop_ed, hop_mps)
    @assert abs(hop_ed - hop_mps) < 1e-6

    # 窓の縮約密度行列: 全窓・全基底で ED の周辺分布と厳密一致するか
    Random.seed!(7)
    maxpe = 0.0
    for q1 in 1:N-3
        q2 = q1 + 3; Nw = 4
        ρw = window_rdm(ψ, q1, q2)
        for _ in 1:5
            basis = rand(1:3, N)
            Uglob = reduce(kron, [Matrix(UBASIS[basis[q]]) for q in 1:N])
            probs = abs2.(Uglob * ψ_ed)
            marg = zeros(2^Nw)
            for idx in 0:2^N-1
                sub = 0
                for q in q1:q2; sub = 2sub + ((idx >> (N - q)) & 1); end
                marg[sub+1] += probs[idx+1]
            end
            maxpe = max(maxpe, maximum(abs.(basis_probs(ρw, basis[q1:q2]) .- marg)))
        end
    end
    @printf("  [check] window RDM vs ED marginals: max |Δp| = %.2e\n", maxpe)
    @assert maxpe < 1e-10

    # サンプラー自体の検定 (経験分布 vs ρ_w の Born 確率)
    ρw = window_rdm(ψ, 3, 6)
    basis = rand(1:3, 4)
    pex = basis_probs(ρw, basis); cum = cumsum(pex)
    counts = zeros(16); bits2 = zeros(Int, 4); nsamp = 200_000
    for _ in 1:nsamp
        draw_bits!(bits2, cum)
        sub = 0
        for j in 1:4; sub = 2sub + bits2[j]; end
        counts[sub+1] += 1
    end
    tv = 0.5 * sum(abs.(counts ./ nsamp .- pex))
    @printf("  [check] sampler TV = %.4f (expect ~ %.4f)\n", tv, 0.4*sqrt(16/nsamp))
    @assert tv < 0.01

    # ドープ系 (δ=1/8, L=8): 粒子数セクターが意図通りか、および
    # 窓RDM から計算した観測量が転送行列の値と一致するか (ED を使わない独立検証)
    Ld = 8
    Ed, ψd, _, _, nel = ground_state(Ld, t, U, U/2; chi_max=64, nsweeps=10, δ=0.125)
    tensd = extract_tensors(ψd)
    ntot = sum(0.5*(1 - product_op_expect(tensd, Dict{Int,Matrix{ComplexF64}}(q => SIGMA[3])))
               for q in 1:2Ld)
    @printf("  [check] doped δ=1/8: N_el(target)=%d  N_el(measured)=%.6f\n", nel, ntot)
    @assert abs(ntot - nel) < 1e-8 && nel == Ld - 1
    maxrd = 0.0
    for i in 1:Ld
        obs = site_observables(i; bond = i < Ld)
        q1, q2 = obs_support(obs)
        ρwd = window_rdm(ψd, q1, q2); Nw = q2 - q1 + 1
        for (o, ow) in zip(obs, shift_obs(obs, q1 - 1))
            v_tm = sum(tm.coeff * product_op_expect(tensd, term_matrixdict(tm)) for tm in o.terms)
            v_rd = 0.0
            for tm in ow.terms
                M = reduce(kron, ((a = findfirst(x -> x[1] == q, tm.sup);
                                   a === nothing ? I2c : SIGMA[tm.sup[a][2]]) for q in 1:Nw))
                v_rd += tm.coeff * real(tr(ρwd * M))
            end
            maxrd = max(maxrd, abs(v_tm - v_rd))
        end
    end
    @printf("  [check] doped window RDM vs transfer matrix: max diff = %.2e\n", maxrd)
    @assert maxrd < 1e-9
    println("  [check] all validations passed"); flush(stdout)
end

# ------------------------------------------------------------
# 7. メイン
# ------------------------------------------------------------
function main()
    t, U = 1.0, 4.0; mu = U/2
    L_list = parse.(Int, split(get(ENV, "L_LIST", "8,16,32"), ","))
    chi_priors = [2, 4, 8, 16, 32]
    nu, nm = 50, 100
    n_repeat = parse(Int, get(ENV, "N_REPEAT", "300"))
    chi_exp = parse(Int, get(ENV, "CHI_EXP", "128"))
    dope_list = parse.(Float64, split(get(ENV, "DOPE_LIST", "0.0,0.125"), ","))

    println("=== validation (L=4) ==="); flush(stdout)
    dense_chain_check()

    rows = []
    for L in L_list, δ in dope_list
        println("\n", "="^72)
        @printf("L=%d (N=%d qubits), U=%.1f, doping δ=%.3f, chi_exp=%d, n_repeat=%d\n",
                L, 2L, U, δ, chi_exp, n_repeat); flush(stdout)
        tstart = time()
        E, ψ, _, _, nel = ground_state(L, t, U, mu; chi_max=chi_exp, δ)
        @printf("  DMRG: E0=%.6f  N_el=%d  maxlinkdim=%d  (%.0fs)\n",
                E, nel, maxlinkdim(ψ), time()-tstart); flush(stdout)

        Sbond = bond_entropies(ψ)

        # 各サイトの観測量と、その台を覆う窓
        obs_by_site = [site_observables(i; bond = i < L) for i in 1:L]
        windows = [obs_support(o) for o in obs_by_site]

        # prior 群: chi 切断 MPS + 「完全prior (σ=ρ)」= 経験的な利得天井 G_max
        priors = MPS[]
        for chi_p in chi_priors
            σ = truncate(ψ; maxdim=chi_p); normalize!(σ)
            push!(priors, σ)
        end
        push!(priors, ψ)                       # 完全prior
        prior_fid = [abs2(inner(σ, ψ)) for σ in priors]

        # 期待値も窓RDMから取る (全鎖の転送行列は χ=256, N=256 では非現実的)
        tstart = time()
        rdm_ρ = window_rdms(ψ, windows)
        rdm_p = [window_rdms(σ, windows) for σ in priors]
        @printf("  window RDMs (rho + %d priors): %.0fs\n", length(priors), time()-tstart)
        prior_lab = vcat(string.(chi_priors), ["exact"])
        @printf("  prior fidelities |<σ|ρ>|²: %s\n",
                join([@sprintf("chi=%s: %.2e", l, f) for (l, f) in zip(prior_lab, prior_fid)], ",  "))
        flush(stdout)

        tstart = time()
        for i in 1:L
            obs = obs_by_site[i]
            q1, q2 = windows[i]
            Nw = q2 - q1 + 1
            obs_w = shift_obs(obs, q1 - 1)
            ρw = rdm_ρ[i]
            Mats = [[term_window_matrix(tm, Nw) for tm in o.terms] for o in obs_w]

            Pρ = [[expect_rdm(ρw, M) for M in Ms] for Ms in Mats]
            Otrue = [sum(tm.coeff * Pρ[k][ti] for (ti, tm) in enumerate(obs[k].terms))
                     for k in 1:length(obs)]
            Pσ_all = Vector{Vector{Vector{Float64}}}(); trOσ_all = Vector{Vector{Float64}}()
            for p in 1:length(priors)
                Pσ = [[expect_rdm(rdm_p[p][i], M) for M in Ms] for Ms in Mats]
                push!(Pσ_all, Pσ)
                push!(trOσ_all, [sum(tm.coeff * Pσ[k][ti] for (ti, tm) in enumerate(obs[k].terms))
                                 for k in 1:length(obs)])
            end

            est_std, est_crm = run_locals(ρw, Nw, obs_w, Pσ_all, trOσ_all;
                                          nu, nm, n_repeat, seed = 7000 + 100L + i)

            for (k, o) in enumerate(obs)
                vstd = var(est_std[:, k])
                a = est_std[:, k]
                Gmax = vstd / var(est_crm[:, k, end])          # 完全prior = ショットノイズ床
                for p in 1:length(chi_priors)
                    b = est_crm[:, k, p]
                    G = vstd / var(b)
                    glo, ghi = boot_ratio(a, b; rng = MersenneTwister(31 + i))
                    Δ = Otrue[k] - trOσ_all[p][k]
                    Gth = NaN
                    if o.pure
                        vs, vc = theory_var(length(o.terms[1].sup), Otrue[k], Δ, nu, nm)
                        Gth = vs / vc
                    end
                    # 観測量の台の中央にあたるボンドエントロピー
                    sb = Sbond[clamp(q2 - 1, 1, length(Sbond))]
                    push!(rows, (L, δ, i, chi_priors[p], o.name, o.pure, Otrue[k], Δ,
                                 prior_fid[p], G, glo, ghi, Gmax, Gth, sb))
                end
            end
            if i % max(1, L ÷ 8) == 0 || i == L
                @printf("    site %2d/%d done (%.0fs elapsed)\n", i, L, time()-tstart); flush(stdout)
            end
        end
        @printf("  site-resolved run: %.0fs\n", time()-tstart); flush(stdout)

        # 画面向けサマリ: chi=8 prior, 代表観測量のサイト依存
        for nm_obs in ["ZZ onsite", "DoubleOcc", "SzSz r=1"]
            sel = [r for r in rows if r[1]==L && r[2]==δ && r[4]==8 && r[5]==nm_obs]
            isempty(sel) && continue
            @printf("\n  [%s, chi_p=8]  site: G (G/Gmax)\n    ", nm_obs)
            for r in sel
                @printf("%d:%.1f(%.2f)  ", r[3], r[10], r[10]/r[13])
            end
            println()
        end
        flush(stdout)
    end

    out = joinpath(@__DIR__, "crm_site_resolved_results.tsv")
    open(out, "w") do io
        println(io, "L\tdoping\tsite\tchi_prior\tobservable\tpure\ttrue\tDelta\tprior_fid\tG_emp\tG_lo\tG_hi\tG_max\tG_theo\tS_bond")
        for r in rows; println(io, join(r, "\t")); end
    end
    println("\nresults saved: $out"); flush(stdout)
    return rows, L_list, chi_priors
end

rows, L_list, chi_priors = main()
