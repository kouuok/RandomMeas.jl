# ============================================================
# CRMシャドウのスケーリング検証 (MPS版, 1Dハバード鎖 L=8/16/32)
#
# 仮説: 局所観測量のCRM利得は prior の「局所誤差 Δ」だけで決まり、
#       グローバル忠実度が L とともに崩壊しても生き残る。
#
# 設計（計算量に注意した点）:
#  - 実験状態 ρ: DMRG (chi_exp<=96) の基底状態。参照実験としては
#    厳密基底状態である必要はない（固定された「実験」であればよい）
#  - prior σ: ρ を chi_p = 2..32 に切断したMPS（prior品質のダイヤル）
#  - 実験側サンプルは全priorで共有（1回サンプル → 後処理のみprior毎）
#  - 局所Pauliのprior側厳密平均は「基底一致 × 3^|A| <P>_σ」の閉形式
#    → u毎の追加コストゼロ
#  - ★窓サンプリング: 観測量の台が載る連続窓 [q1,q2] のみをサンプルする。
#    直交中心を q1 に置くと左環境は恒等になるので、左Schmidt指標 l を
#    確率 p(l)=diag(Σ_b C[b]C[b]†) （回転に依存しない）で先に引き、
#    v0 = e_l から通常の逐次サンプリングをすればよい。
#    コストが窓幅のみに依存し、L に依存しなくなる。
#  - BLASは1スレッド（小さいgemvではスレッド化が逆効果）
#  - 忠実度推定は局所Pauliシャドウでは分散が qubit数に対し指数的に悪化
#    するため L=8 のみで実施（それ自体が「局所観測量に集中すべき」根拠）
#
# 実行方法:
#  JULIA_LOAD_PATH="@:@v#.#:@stdlib" julia --project=Hubbard_MPS_Env_v2 crm_mps_scaling.jl
# ============================================================

using ITensors, ITensorMPS
using LinearAlgebra
using Statistics
using Printf
using Random

BLAS.set_num_threads(1)

# ------------------------------------------------------------
# 1. モデル: 1Dハバード鎖, JW: qubit(2i-1)=site i up, qubit(2i)=site i dn
# ------------------------------------------------------------
qup(i) = 2i - 1
qdn(i) = 2i

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

function ground_state(L, t, U, mu; chi_max=96, nsweeps=10)
    N = 2L
    sites = siteinds("Fermion", N; conserve_qns=true)
    H = hubbard_chain_mpo(sites, L, t, U, mu)
    init = [isodd(div(q-1, 2) + 1) == isodd(q) ? "Occ" : "Emp" for q in 1:N]
    ψ0 = productMPS(sites, init)
    maxdims = min.(chi_max, [20, 50, 100, chi_max, chi_max, chi_max, chi_max, chi_max, chi_max, chi_max])[1:nsweeps]
    noise = [1e-6, 1e-6, 1e-7, 1e-8, 1e-9, 0, 0, 0, 0, 0][1:nsweeps]
    E, ψ = dmrg(H, ψ0; nsweeps, maxdim=maxdims, cutoff=1e-10, noise, outputlevel=0)
    return E, ψ, sites, H
end

# ------------------------------------------------------------
# 2. MPS -> 密テンソル抽出
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

# 期待値計算用 (任意ゲージでよい)
function extract_tensors(ψ::MPS)
    ψo = orthogonalize(ψ, 1)
    return [site_matrices(ψo, j) for j in 1:length(ψo)]
end

# ★窓サンプリング用: 直交中心を q1 に置き、窓テンソルと左指標分布を返す
function extract_window(ψ::MPS, q1::Int, q2::Int)
    ψo = orthogonalize(ψ, q1)
    tens = [site_matrices(ψo, j) for j in q1:q2]
    C1, C2 = tens[1]
    pl = vec(sum(abs2, C1; dims=2) .+ sum(abs2, C2; dims=2))
    pl ./= sum(pl)
    return tens, cumsum(pl)
end

# 積演算子 <ψ| ⊗_q O_q |ψ> の厳密期待値 (転送行列)
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
# 3. シャドウ測定機構 (基底 1=X, 2=Y, 3=Z)
# ------------------------------------------------------------
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

# 完全サンプリング。左指標 l を cum_pl から引き v0=e_l で開始
# (全鎖サンプリングは q1=1 で cum_pl=[1.0] の特殊例)
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
# 4. 観測量 (Pauli項のリスト; JW後のqubit Pauli列)
# ------------------------------------------------------------
struct Term
    coeff::Float64
    sup::Vector{Tuple{Int,Int}}   # (qubit, pauli 1=X,2=Y,3=Z)
end
struct Obs
    name::String
    terms::Vector{Term}
    pure::Bool
end

term_matrixdict(t::Term) = Dict{Int,Matrix{ComplexF64}}(q => SIGMA[a] for (q,a) in t.sup)

function build_observables(L)
    i0 = L ÷ 2
    obs = Obs[]
    push!(obs, Obs("ZZ onsite(i0)", [Term(1.0, [(qup(i0),3),(qdn(i0),3)])], true))
    for r in [1, 2, 4, 8]
        r <= L - i0 || continue
        push!(obs, Obs("ZZ up-up r=$r", [Term(1.0, [(qup(i0),3),(qup(i0+r),3)])], true))
    end
    push!(obs, Obs("DoubleOcc(i0)",
        [Term(0.25, Tuple{Int,Int}[]), Term(-0.25, [(qup(i0),3)]),
         Term(-0.25, [(qdn(i0),3)]),  Term(0.25, [(qup(i0),3),(qdn(i0),3)])], false))
    push!(obs, Obs("SzSz r=1",
        [Term( 1/16, [(qup(i0),3),(qup(i0+1),3)]), Term(-1/16, [(qup(i0),3),(qdn(i0+1),3)]),
         Term(-1/16, [(qdn(i0),3),(qup(i0+1),3)]), Term( 1/16, [(qdn(i0),3),(qdn(i0+1),3)])], false))
    q1, qm, q2 = qup(i0), qdn(i0), qup(i0+1)
    push!(obs, Obs("hop bond(i0)",
        [Term(0.5, [(q1,1),(qm,3),(q2,1)]), Term(0.5, [(q1,2),(qm,3),(q2,2)])], false))
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

# ------------------------------------------------------------
# 5. 局所観測量の実験 (窓サンプリング; 実験サンプルを全priorで共有)
# ------------------------------------------------------------
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
# 6. 忠実度の実験 (L=8のみ, 全鎖; prior側は安価な古典サンプル N_cl >> nm)
# ------------------------------------------------------------
function fid_snapshot(rotσ, bits)
    Lenv = ones(ComplexF64, 1, 1)
    for j in 1:length(rotσ)
        A1, A2 = rotσ[j]
        k1 = bits[j] == 0 ? 2.0 : -1.0
        k2 = bits[j] == 1 ? 2.0 : -1.0
        Lenv = k1 .* (A1' * Lenv * A1) .+ k2 .* (A2' * Lenv * A2)
    end
    return real(Lenv[1,1])
end

function run_fidelity(ρw, ρcpl, σw, σcpl; nu, nm, n_cl, n_repeat, seed)
    Random.seed!(seed)
    N = length(ρw)
    chiρ = max(length(ρcpl), maximum(size(t[1], 2) for t in ρw))
    chiσ = max(length(σcpl), maximum(size(t[1], 2) for t in σw))
    rotρ = make_rot_buffers(ρw); rotσ = make_rot_buffers(σw)
    vb = zeros(ComplexF64, chiρ + 1); wb = zeros(ComplexF64, 2chiρ + 2)
    vbσ = zeros(ComplexF64, chiσ + 1); wbσ = zeros(ComplexF64, 2chiσ + 2)
    bits = zeros(Int, N)
    est = zeros(n_repeat, 2)
    for rep in 1:n_repeat
        acc_s, acc_c = 0.0, 0.0
        for _ in 1:nu
            basis = rand(1:3, N)
            rotate!(rotρ, ρw, basis); rotate!(rotσ, σw, basis)
            mρ = 0.0
            for _ in 1:nm
                sample_bits!(bits, rotρ, vb, wb, ρcpl)
                mρ += fid_snapshot(rotσ, bits)
            end
            mρ /= nm
            mσ = 0.0
            for _ in 1:n_cl
                sample_bits!(bits, rotσ, vbσ, wbσ, σcpl)
                mσ += fid_snapshot(rotσ, bits)
            end
            mσ /= n_cl
            acc_s += mρ; acc_c += mρ - mσ
        end
        est[rep, 1] = acc_s / nu
        est[rep, 2] = acc_c / nu + 1.0
    end
    return est
end

# ------------------------------------------------------------
# 7. 検証 (L=4): 密行列EDとの突き合わせ + サンプラー検定
# ------------------------------------------------------------
function dense_chain_check()
    L, t, U = 4, 1.0, 4.0; mu = U/2; N = 2L
    cdag = [0.0 0.0; 1.0 0.0]; I2 = Matrix(1.0I, 2, 2); F = [1.0 0.0; 0.0 -1.0]
    Cd = Vector{Matrix{Float64}}(undef, N)
    for k in 1:N
        ops = [j < k ? F : (j == k ? cdag : I2) for j in 1:N]
        Cd[k] = reduce(kron, ops)
    end
    C = [Matrix(m') for m in Cd]; Nop = [Cd[k]*C[k] for k in 1:N]
    H = zeros(2^N, 2^N)
    for i in 1:L-1, (a, b) in ((qup(i), qup(i+1)), (qdn(i), qdn(i+1)))
        H .+= -t .* (Cd[a]*C[b] .+ Cd[b]*C[a])
    end
    for i in 1:L
        H .+= U .* (Nop[qup(i)]*Nop[qdn(i)]) .- mu .* (Nop[qup(i)] .+ Nop[qdn(i)])
    end
    F_ed = eigen(Symmetric(H)); E_ed = F_ed.values[1]; ψ_ed = F_ed.vectors[:, 1]

    E_dmrg, ψ, sites, _ = ground_state(L, t, U, mu; chi_max=64, nsweeps=8)
    @printf("  [check] E_ED=%.8f  E_DMRG=%.8f  (diff %.2e)\n", E_ed, E_dmrg, abs(E_ed-E_dmrg))
    @assert abs(E_ed - E_dmrg) < 1e-6

    tens = extract_tensors(ψ)
    pauliZ = [1.0 0.0; 0.0 -1.0]
    for (qs, name) in [([qup(2), qdn(2)], "ZZ onsite"), ([qup(2), qup(3)], "ZZ up r=1")]
        Mdense = reduce(kron, [q in qs ? pauliZ : I2 for q in 1:N])
        v_ed = dot(ψ_ed, Mdense*ψ_ed)
        v_mps = product_op_expect(tens, Dict{Int,Matrix{ComplexF64}}(q => SIGMA[3] for q in qs))
        @assert abs(v_ed - v_mps) < 1e-6 "$name mismatch: $v_ed vs $v_mps"
    end
    hop_ed = dot(ψ_ed, (Cd[qup(2)]*C[qup(3)] .+ Cd[qup(3)]*C[qup(2)])*ψ_ed)
    hop_obs = [(Dict{Int,Matrix{ComplexF64}}(qup(2)=>SIGMA[1], qdn(2)=>SIGMA[3], qup(3)=>SIGMA[1]), 0.5),
               (Dict{Int,Matrix{ComplexF64}}(qup(2)=>SIGMA[2], qdn(2)=>SIGMA[3], qup(3)=>SIGMA[2]), 0.5)]
    hop_mps = sum(c * product_op_expect(tens, d) for (d, c) in hop_obs)
    @printf("  [check] hop: ED=%.6f  MPS(XZX+YZY)/2=%.6f\n", hop_ed, hop_mps)
    @assert abs(hop_ed - hop_mps) < 1e-6

    # 全鎖サンプラー検定: 固定基底で経験分布 vs 厳密Born確率
    Random.seed!(7)
    basis = rand(1:3, N)
    wt, cpl = extract_window(ψ, 1, N)
    rot = make_rot_buffers(wt); rotate!(rot, wt, basis)
    Uglob = reduce(kron, [Matrix(UBASIS[basis[q]]) for q in 1:N])
    probs = abs2.(Uglob * ψ_ed)
    counts = zeros(2^N); bits = zeros(Int, N)
    vb = zeros(ComplexF64, 100); wb = zeros(ComplexF64, 200)
    nsamp = 200_000
    for _ in 1:nsamp
        sample_bits!(bits, rot, vb, wb, cpl)
        idx = 0
        for q in 1:N; idx = 2idx + bits[q]; end
        counts[idx+1] += 1
    end
    tv = 0.5 * sum(abs.(counts ./ nsamp .- probs))
    @printf("  [check] full sampler TV = %.4f (expect ~ %.4f)\n", tv, 0.4*sqrt(2^N/nsamp))
    @assert tv < 0.05

    # ★窓サンプラー検定: 窓 [3,6] の周辺分布 vs 密行列周辺分布
    q1, q2 = 3, 6; Nw = q2 - q1 + 1
    wt2, cpl2 = extract_window(ψ, q1, q2)
    rot2 = make_rot_buffers(wt2); rotate!(rot2, wt2, basis[q1:q2])
    marg = zeros(2^Nw)
    for idx in 0:2^N-1
        sub = 0
        for q in q1:q2
            sub = 2sub + ((idx >> (N - q)) & 1)
        end
        marg[sub+1] += probs[idx+1]
    end
    counts2 = zeros(2^Nw); bits2 = zeros(Int, Nw)
    for _ in 1:nsamp
        sample_bits!(bits2, rot2, vb, wb, cpl2)
        sub = 0
        for j in 1:Nw; sub = 2sub + bits2[j]; end
        counts2[sub+1] += 1
    end
    tv2 = 0.5 * sum(abs.(counts2 ./ nsamp .- marg))
    @printf("  [check] window sampler TV = %.4f (expect ~ %.4f)\n", tv2, 0.4*sqrt(2^Nw/nsamp))
    @assert tv2 < 0.02
    println("  [check] all validations passed"); flush(stdout)
end

# ------------------------------------------------------------
# 8. メイン
# ------------------------------------------------------------
function main()
    t, U = 1.0, 4.0; mu = U/2
    L_list = [8, 16, 32]
    chi_priors = [2, 4, 8, 16, 32]
    nu, nm, n_repeat = 50, 100, 100

    println("=== validation (L=4) ==="); flush(stdout)
    dense_chain_check()

    rows = []
    fid_table = Dict()

    for L in L_list
        println("\n", "="^70)
        @printf("L=%d (N=%d qubits), U=%.1f, half-filling\n", L, 2L, U); flush(stdout)
        chi_exp = L <= 8 ? 64 : 96
        tstart = time()
        E, ψ, sites, H = ground_state(L, t, U, mu; chi_max=chi_exp)
        @printf("  DMRG(chi=%d): E0=%.6f  maxlinkdim=%d  (%.0fs)\n",
                chi_exp, E, maxlinkdim(ψ), time()-tstart); flush(stdout)

        exp_tens = extract_tensors(ψ)
        obs = build_observables(L)
        q1, q2 = obs_support(obs)
        obs_w = shift_obs(obs, q1 - 1)
        wtens, cum_pl = extract_window(ψ, q1, q2)
        @printf("  sampling window: qubits [%d, %d] (%d of %d)\n", q1, q2, q2-q1+1, 2L)

        Pρ = [[product_op_expect(exp_tens, term_matrixdict(tm)) for tm in o.terms] for o in obs]
        Otrue = [sum(tm.coeff * Pρ[k][ti] for (ti, tm) in enumerate(obs[k].terms)) for k in 1:length(obs)]

        Pσ_all = Vector{Vector{Vector{Float64}}}(); trOσ_all = Vector{Vector{Float64}}()
        fids = Float64[]; σ_windows = []
        for chi_p in chi_priors
            σ = truncate(ψ; maxdim=chi_p); normalize!(σ)
            push!(fids, abs2(inner(σ, ψ)))
            st = extract_tensors(σ)
            push!(σ_windows, extract_window(σ, 1, 2L))
            Pσ = [[product_op_expect(st, term_matrixdict(tm)) for tm in o.terms] for o in obs]
            push!(Pσ_all, Pσ)
            push!(trOσ_all, [sum(tm.coeff * Pσ[k][ti] for (ti, tm) in enumerate(obs[k].terms)) for k in 1:length(obs)])
        end
        @printf("  prior fidelities |<σ|ρ>|²: %s\n",
                join([@sprintf("chi=%d: %.4f", c, f) for (c, f) in zip(chi_priors, fids)], ",  "))
        flush(stdout)

        tstart = time()
        est_std, est_crm = run_locals(wtens, cum_pl, obs_w, Pσ_all, trOσ_all;
                                      nu, nm, n_repeat, seed=100 + L)
        @printf("  local-observable run: %.0fs\n", time()-tstart); flush(stdout)

        @printf("\n  %-16s %9s %9s |", "observable", "true", "std_err")
        for c in chi_priors; @printf(" G(chi=%d)", c); end
        println()
        for (k, o) in enumerate(obs)
            v_std = var(est_std[:, k])
            @printf("  %-16s %9.4f %9.4f |", o.name, Otrue[k], sqrt(v_std))
            for p in 1:length(chi_priors)
                G = v_std / var(est_crm[:, k, p])
                @printf(" %8.2f", G)
                Δ = Otrue[k] - trOσ_all[p][k]
                G_theo = NaN
                if o.pure
                    vs, vc = theory_var(length(o.terms[1].sup), Otrue[k], Δ, nu, nm)
                    G_theo = vs / vc
                end
                push!(rows, (L, chi_priors[p], o.name, o.pure, Otrue[k], Δ, fids[p], G, G_theo))
            end
            println()
        end
        flush(stdout)

        if L == 8
            println("\n  fidelity estimation (L=8 only; local shadows scale exponentially):")
            ρw, ρcpl = extract_window(ψ, 1, 2L)
            for (p, chi_p) in enumerate(chi_priors)
                chi_p in (4, 16) || continue
                σw, σcpl = σ_windows[p]
                est = run_fidelity(ρw, ρcpl, σw, σcpl; nu, nm, n_cl=300,
                                   n_repeat=50, seed=999 + chi_p)
                G = var(est[:, 1]) / var(est[:, 2])
                fid_table[(L, chi_p)] = (fids[p], mean(est[:,1]), std(est[:,1]),
                                         mean(est[:,2]), std(est[:,2]), G)
                @printf("    chi=%2d (F=%.3f): std= %.3f ± %.3f   CRM= %.3f ± %.3f   G=%.2f\n",
                        chi_p, fids[p], mean(est[:,1]), std(est[:,1]),
                        mean(est[:,2]), std(est[:,2]), G)
                flush(stdout)
            end
        end
    end

    out = joinpath(@__DIR__, "crm_mps_scaling_results.tsv")
    open(out, "w") do io
        println(io, "L\tchi_prior\tobservable\tpure\ttrue\tDelta\tprior_fid\tG_emp\tG_theo")
        for r in rows
            println(io, join(r, "\t"))
        end
    end
    println("\nresults saved: $out"); flush(stdout)
    return rows, fid_table, L_list, chi_priors
end

rows, fid_table, L_list, chi_priors = main()

# ------------------------------------------------------------
# 9. プロット
# ------------------------------------------------------------
using Plots

p1 = plot(xlabel="L (chain length)", ylabel="G = Var_std / Var_CRM",
          yscale=:log10, legend=:topleft, title="Local gain vs L  (U=4)")
markers = Dict(4 => :circle, 16 => :square)
for chi_p in (4, 16)
    for (name, col) in [("ZZ onsite(i0)", :firebrick), ("ZZ up-up r=2", :steelblue)]
        g = [only([r[8] for r in rows if r[1]==L && r[2]==chi_p && r[3]==name]) for L in L_list]
        plot!(p1, L_list, g, marker=markers[chi_p], color=col, lw=2,
              label="$name chi=$chi_p")
    end
end
fidref = [only(unique([r[7] for r in rows if r[1]==L && r[2]==16])) for L in L_list]
plot!(p1, L_list, max.(fidref, 1e-3), color=:gray, ls=:dash, marker=:x, lw=1.5,
      label="prior fidelity (chi=16)")
hline!(p1, [1.0], color=:black, ls=:dot, label="")

gx = [r[9] for r in rows if r[4] && !isnan(r[9])]
gy = [r[8] for r in rows if r[4] && !isnan(r[9])]
lims = (0.5*min(minimum(gx), minimum(gy)), 2*max(maximum(gx), maximum(gy)))
p2 = scatter(gx, gy, xscale=:log10, yscale=:log10, xlims=lims, ylims=lims,
             color=:firebrick, ms=5, alpha=0.7, label="pure Pauli strings",
             xlabel="G theory", ylabel="G empirical",
             title="Scaling law: all (L, chi_p, obs)")
plot!(p2, [lims...], [lims...], color=:black, ls=:dash, label="y = x")

p3 = plot(xlabel="chi_prior", ylabel="G", xscale=:log10, yscale=:log10,
          legend=:bottomright, title="Gain vs prior quality (L=32)")
for name in unique(r[3] for r in rows if r[1]==32)
    g = [only([r[8] for r in rows if r[1]==32 && r[2]==c && r[3]==name]) for c in chi_priors]
    plot!(p3, chi_priors, g, marker=:circle, lw=2, label=name)
end
hline!(p3, [1.0], color=:black, ls=:dot, label="")

final = plot(p1, p2, p3, layout=(1,3), size=(1500, 430), margin=6Plots.mm)
outpng = joinpath(@__DIR__, "crm_mps_scaling.png")
savefig(final, outpng)
println("figure saved: $outpng")
