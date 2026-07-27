# ============================================================
# 1Dハバード鎖のCRM実験で共通に使う機構
#   - モデル / QN保存DMRG (doping δ 指定可)
#   - 窓の縮約密度行列 (ボンド次元に依存しないサンプリングと期待値)
#   - 観測量の定義とスナップショット推定量
#   - 分散比の bootstrap CI
# crm_site_resolved.jl から抽出したもので、数値的挙動は同一。
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
    # χ を段階的に上げ、最初の数スイープだけノイズを入れる（nsweeps 任意）
    ramp = vcat([20, 50, 100, 200], fill(chi_max, max(0, nsweeps - 4)))
    maxdims = min.(chi_max, ramp)[1:nsweeps]
    noise = vcat([1e-6, 1e-6, 1e-7, 1e-8, 1e-9], zeros(max(0, nsweeps - 5)))[1:nsweeps]
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

