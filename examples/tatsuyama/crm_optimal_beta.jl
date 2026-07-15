# ============================================================
# 最適係数CRM と 多prior回帰CRM (方向1-2)
#
# CRMは制御変量法のβ=1特殊例:
#   o(β) = mean_u[mρ(u)] - β (mean_u[Y(u)] - Tr[Oσ]),  Y(u)=E[X_σ|u] (厳密)
# 最適係数 β* = Cov(mρ,Y)/Var(Y) を同一データからプラグイン推定すると
#   G_opt = 1/(1-corr²) >= 1
# となり「priorがどんなに悪くても(漸近的に)損をしない」。
# 複数prior {Y_i} は回帰 β = Σ_YY^+ c_Yρ に一般化され G = 1/(1-R²)。
#
# 実証: W=4x8シリンダー (crm_2d_symuhf.jl と同一系・同一seed=同一測定データ)
#   prior: UHF / UHF-sym / chi=8 / chi=32
#   比較: 標準 / CRM(β=1) / 最適β / 多prior回帰 {UHF,UHF-sym} と {全4prior}
#   検証: 達成G vs 理論 1/(1-R²) の一致、推定量の不偏性(真値との比較)
#
# 実行方法:
#  JULIA_LOAD_PATH="@:@v#.#:@stdlib" julia --project=Hubbard_MPS_Env_v2 crm_optimal_beta.jl
# ============================================================

using ITensors, ITensorMPS
using LinearAlgebra
using Statistics
using Printf
using Random

BLAS.set_num_threads(1)

# ------------------------------------------------------------
# 0. 格子
# ------------------------------------------------------------
const W  = 4
const LX = 8
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

function ground_state_cyl(Lx, W, t, U, mu; chi_max=192, nsweeps=10)
    N = 2 * Lx * W
    sites = siteinds("Fermion", N; conserve_qns=true)
    H = hubbard_cyl_mpo(sites, Lx, W, t, U, mu)
    # Néel初期状態: (x+y)偶 -> up占有, 奇 -> dn占有
    init = fill("Emp", N)
    for x in 1:Lx, y in 1:W
        s = sidx(x, y)
        init[iseven(x + y) ? qup(s) : qdn(s)] = "Occ"
    end
    ψ0 = productMPS(sites, init)
    maxdims = min.(chi_max, [20, 50, 100, 150, chi_max, chi_max, chi_max, chi_max, chi_max, chi_max])[1:nsweeps]
    noise = [1e-6, 1e-6, 1e-7, 1e-8, 1e-9, 0, 0, 0, 0, 0][1:nsweeps]
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

function solve_uhf(Lx, W, t, U; Nup, Ndn, iters=800, tol=1e-10)
    n = Lx * W
    T = hopping_matrix(Lx, W, t)
    nup = zeros(n); ndn = zeros(n)
    for x in 1:Lx, y in 1:W
        s = sidx(x, y)
        nup[s] = 0.5 + 0.4 * (-1)^(x + y)
        ndn[s] = 0.5 - 0.4 * (-1)^(x + y)
    end
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
    E_dmrg, ψ = ground_state_cyl(Lxv, W, t, 0.0, 0.0; chi_max=256, nsweeps=10)
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
# 5b. 実験 (per-uデータを保持し、β=1 / 最適β / 多prior回帰を同時評価)
# ------------------------------------------------------------
function run_locals_beta(wtens, cum_pl, obs_w, Pσ_all, trOσ_all, prior_sets;
                         nu, nm, n_repeat, seed)
    Random.seed!(seed)
    Nw = length(wtens); nobs = length(obs_w); np = length(Pσ_all)
    nsets = length(prior_sets)
    nmeth = 1 + 2np + nsets          # std | CRM(β=1)×np | optβ×np | multi×nsets
    chimax = max(length(cum_pl), maximum(size(t[1], 2) for t in wtens))
    rot = make_rot_buffers(wtens)
    vbuf = zeros(ComplexF64, chimax + 1); wbuf = zeros(ComplexF64, 2chimax + 2)
    bits = zeros(Int, Nw)
    est = zeros(n_repeat, nobs, nmeth)
    # プールした相関 (理論 1/(1-R²) 比較用)
    Sx = zeros(nobs); Sxx = zeros(nobs)
    Sy = zeros(nobs, np); Syy = zeros(nobs, np); Sxy = zeros(nobs, np)
    ntot = 0
    mρu = zeros(nu, nobs); Yu = zeros(nu, nobs, np)
    xs = zeros(nobs)
    for rep in 1:n_repeat
        for iu in 1:nu
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
                mρu[iu, k] = xs[k] / nm
                for p in 1:np
                    Yu[iu, k, p] = exact_prior_mean(obs_w[k], basis, Pσ_all[p][k])
                end
            end
        end
        ntot += nu
        for k in 1:nobs
            x = view(mρu, :, k)
            mx = mean(x)
            est[rep, k, 1] = mx
            Sx[k] += sum(x); Sxx[k] += sum(abs2, x)
            for p in 1:np
                y = view(Yu, :, k, p)
                my = mean(y)
                trv = trOσ_all[p][k]
                est[rep, k, 1+p] = mx - my + trv
                vy = var(y); cxy = cov(x, y)
                β = vy > 1e-12 ? cxy / vy : 0.0
                est[rep, k, 1+np+p] = mx - β * (my - trv)
                Sy[k,p] += sum(y); Syy[k,p] += sum(abs2, y); Sxy[k,p] += dot(x, y)
            end
            for (si, ps) in enumerate(prior_sets)
                q = length(ps)
                Ym = Yu[:, k, ps]
                my = vec(mean(Ym; dims=1))
                Σm = cov(Ym)
                cv = [cov(x, view(Ym, :, j)) for j in 1:q]
                β = pinv(Σm, rtol=1e-8) * cv
                trvs = [trOσ_all[p][k] for p in ps]
                est[rep, k, 1+2np+si] = mx - dot(β, my .- trvs)
            end
        end
    end
    # プール相関
    corr = zeros(nobs, np)
    for k in 1:nobs, p in 1:np
        vx = Sxx[k]/ntot - (Sx[k]/ntot)^2
        vy = Syy[k,p]/ntot - (Sy[k,p]/ntot)^2
        cxy = Sxy[k,p]/ntot - (Sx[k]/ntot)*(Sy[k,p]/ntot)
        corr[k,p] = (vx > 1e-14 && vy > 1e-14) ? cxy/sqrt(vx*vy) : 0.0
    end
    return est, corr
end

# ------------------------------------------------------------
# 6. メイン
# ------------------------------------------------------------
function main()
    t = 1.0
    chi_priors = [8, 32]
    prior_labels = ["UHF", "UHF-sym", "chi=8", "chi=32"]
    prior_sets = [[1, 2], [1, 2, 3, 4]]           # multi(free) / multi(all)
    set_labels = ["multi(free)", "multi(all)"]
    nu, nm, n_repeat = 50, 100, 50
    n = LX * W; Nup = n ÷ 2; Ndn = n ÷ 2

    println("=== validation (U=0, $(W)x2 cylinder, exact MPS) ==="); flush(stdout)
    validate_u0()

    rows = []
    for U in [4.0, 8.0]
        println("\n", "="^70)
        @printf("W=%d x Lx=%d cylinder (N=%d qubits), U=%.1f, half-filling\n", W, LX, 2n, U)
        flush(stdout)
        tstart = time()
        BLAS.set_num_threads(4)
        E, ψ = ground_state_cyl(LX, W, t, U, U/2; chi_max=192, nsweeps=10)
        BLAS.set_num_threads(1)
        @printf("  DMRG(chi=192): E0=%.6f  (%.0fs)\n", E, time()-tstart); flush(stdout)

        uhf = solve_uhf(LX, W, t, U; Nup, Ndn)
        exp_tens = extract_tensors(ψ)
        obs = build_observables_2d(LX, W)
        q1, q2 = obs_support(obs)
        obs_w = shift_obs(obs, q1 - 1)
        wtens, cum_pl = extract_window(ψ, q1, q2)

        Pρ = [[product_op_expect(exp_tens, term_matrixdict(tm)) for tm in o.terms] for o in obs]
        Otrue = [sum(tm.coeff * Pρ[k][ti] for (ti, tm) in enumerate(obs[k].terms)) for k in 1:length(obs)]

        Pσ_all = Vector{Vector{Vector{Float64}}}(); trOσ_all = Vector{Vector{Float64}}()
        push!(Pσ_all, [[wick_term(tm.sup, uhf) for tm in o.terms] for o in obs])
        push!(Pσ_all, symavg_P(obs, uhf, n))
        for chi_p in chi_priors
            σ = truncate(ψ; maxdim=chi_p); normalize!(σ)
            st = extract_tensors(σ)
            push!(Pσ_all, [[product_op_expect(st, term_matrixdict(tm)) for tm in o.terms] for o in obs])
        end
        for p in 1:length(Pσ_all)
            push!(trOσ_all, [sum(tm.coeff * Pσ_all[p][k][ti] for (ti, tm) in enumerate(obs[k].terms)) for k in 1:length(obs)])
        end

        tstart = time()
        est, corr = run_locals_beta(wtens, cum_pl, obs_w, Pσ_all, trOσ_all, prior_sets;
                                    nu, nm, n_repeat, seed=2000 + round(Int, 10U))
        @printf("  run: %.0fs\n", time()-tstart); flush(stdout)

        np = length(prior_labels)
        @printf("\n  ===== U=%.0f: G = Var_std / Var_method =====\n", U)
        @printf("  %-12s |", "observable")
        for lab in prior_labels; @printf(" %8s", "b1:$lab"); end
        for lab in prior_labels; @printf(" %8s", "op:$lab"); end
        for lab in set_labels; @printf(" %11s", lab); end
        println()
        for (k, o) in enumerate(obs)
            v_std = var(est[:, k, 1])
            @printf("  %-12s |", o.name)
            vals = Float64[]
            for m in 2:size(est, 3)
                G = v_std / var(est[:, k, m])
                push!(vals, G)
                if m <= 1 + 2np
                    @printf(" %8.2f", G)
                else
                    @printf(" %11.2f", G)
                end
            end
            println()
            # バイアス確認: 各推定法の平均が真値と整合するか (最悪ケースのz値)
            zmax = maximum(abs(mean(est[:,k,m]) - Otrue[k]) / (std(est[:,k,m])/sqrt(n_repeat))
                           for m in 1:size(est,3))
            push!(rows, (U, o.name, Otrue[k], sqrt(v_std), vals..., zmax,
                         [corr[k,p] for p in 1:np]...))
        end
        @printf("  (bias check: max |mean-true|/SE over methods, per obs -> tsvのzmax列)\n")
        flush(stdout)
    end

    out = joinpath(@__DIR__, "crm_optbeta_results.tsv")
    open(out, "w") do io
        hdr = ["U","observable","true","err_std",
               ["G_b1_$(l)" for l in ["UHF","UHFsym","chi8","chi32"]]...,
               ["G_opt_$(l)" for l in ["UHF","UHFsym","chi8","chi32"]]...,
               "G_multifree","G_multiall","zmax",
               ["corr_$(l)" for l in ["UHF","UHFsym","chi8","chi32"]]...]
        println(io, join(hdr, "\t"))
        for r in rows; println(io, join(r, "\t")); end
    end
    println("\nresults saved: $out"); flush(stdout)
    return rows
end

rows = main()

# ------------------------------------------------------------
# 7. プロット
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
             ylabel=(iu == 1 ? "G" : ""), legend=(iu == 1 ? :bottomright : :none),
             title=@sprintf("W=4 cylinder,  U=%.0f", U))
    # 列: r = (U,name,true,err, G_b1(4), G_opt(4), G_multi(2), zmax, corr(4))
    series = [(5,  "CRM β=1 (UHF)", :dash),
              (9,  "opt-β (UHF)", :solid),
              (13, "multi(free: UHF+UHF-sym)", :solid),
              (14, "multi(all: +MPS)", :solid)]
    for (e, (col, lab, ls)) in enumerate(series)
        plot!(p, 1:length(names), [r[col] for r in sel], color=OI[e], ls=ls,
              marker=:circle, label=lab)
    end
    hline!(p, [1.0], color=:black, ls=:dot, label="")
    push!(panels, p)
end

# 達成G(opt) vs 理論 1/(1-corr²) : 全 (U, obs, prior)
gx = Float64[]; gy = Float64[]
for r in rows, p in 1:4
    c = r[15+p]                # corr列
    push!(gx, 1 / max(1 - c^2, 1e-6))
    push!(gy, r[8+p])          # G_opt列
end
lims = (0.7, 2 * max(maximum(gx), maximum(gy)))
p3 = scatter(gx, gy, xscale=:log10, yscale=:log10, xlims=lims, ylims=lims,
             color=:firebrick, ms=5, alpha=0.7, label="(obs, prior, U)",
             xlabel="1 / (1 - corr²)", ylabel="achieved G (opt-β)",
             title="Optimal-coefficient theory check")
plot!(p3, [lims...], [lims...], color=:black, ls=:dash, label="y = x")

final = plot(panels..., p3, layout=(1, 3), size=(1650, 440), margin=7Plots.mm)
outpng = joinpath(@__DIR__, "crm_optbeta.png")
savefig(final, outpng)
println("figure saved: $outpng")
