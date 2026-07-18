# ============================================================
# 多観測量同時推定ベンチマーク: ランダムPauliシャドウ(±CRM) vs 貪欲derandomization
#
# 査読対策の中心実験。「単一観測量なら直接測定が最強」への回答として、
# 同一総ショット数 (n_settings × nm) で多数の観測量を同時推定する場合の
# 誤差を比較する:
#   (1) ランダムPauliシャドウ (標準)
#   (2)     + CRM (切断MPS χ=32 prior, 厳密古典側, β=1)
#   (3) 貪欲derandomization (HKP型: 指数減衰重みのカバレッジ最適化で
#       測定基底を決定論的に選ぶ。逆チャネル係数 3^|A| は不要)
#
# セットA (構造化): 窓内の物理観測量 (二重占有/SzSz/ホッピング/距離別ZZ)
# セットB (非構造): ランダムな2-3局所Pauli列 150個
# 期待: 構造化少数 → derand有利 / 多数非構造 → シャドウ系が拮抗し、
#       <P>≠0 の物理観測量には CRM が追加利得
#
# 系: 1Dハバード鎖 L=16, U=4, half-filling (参照 DMRG χ=96, 窓18量子ビット)
# 実行: JULIA_LOAD_PATH="@:@v#.#:@stdlib" julia --project=Hubbard_MPS_Env_v2 crm_derand_benchmark.jl
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
# ------------------------------------------------------------
# 5. 観測量セット
# ------------------------------------------------------------
docc_terms(s) = [Term(0.25, Tuple{Int,Int}[]), Term(-0.25, [(qup(s),3)]),
                 Term(-0.25, [(qdn(s),3)]), Term(0.25, [(qup(s),3),(qdn(s),3)])]
szsz_terms(si, sj) = [Term( 1/16, [(qup(si),3),(qup(sj),3)]), Term(-1/16, [(qup(si),3),(qdn(sj),3)]),
                      Term(-1/16, [(qdn(si),3),(qup(sj),3)]), Term( 1/16, [(qdn(si),3),(qdn(sj),3)])]
hop_terms(s) = [Term(0.5, [(qup(s),1),(qdn(s),3),(qup(s+1),1)]),
                Term(0.5, [(qup(s),2),(qdn(s),3),(qup(s+1),2)])]

function build_setA(i0)
    obs = Obs[]
    for s in i0-3:i0+3
        push!(obs, Obs("docc(s$s)", docc_terms(s), false))
    end
    for s in i0-3:i0+2
        push!(obs, Obs("SzSz($s)", szsz_terms(s, s+1), false))
        push!(obs, Obs("hop($s)", hop_terms(s), false))
    end
    for r in 1:4
        push!(obs, Obs("ZZ r=$r", [Term(1.0, [(qup(i0),3),(qup(i0+r),3)])], true))
    end
    return obs
end

function build_setB(q1, q2; M=150, seed=42)
    rng = MersenneTwister(seed)
    obs = Obs[]
    for i in 1:M
        k = rand(rng, (2, 3))
        qs = sort(collect(q1:q2)[randperm(rng, q2-q1+1)][1:k])
        push!(obs, Obs("rnd$i", [Term(1.0, [(q, rand(rng, 1:3)) for q in qs])], true))
    end
    return obs
end

# ------------------------------------------------------------
# 6. 貪欲derandomization (HKP型)
#    ユニークPauli列単位で命中回数を追跡し、指数減衰重み exp(-α·hits) の
#    命中期待値を最大化する基底を qubit 順に貪欲選択する。
# ------------------------------------------------------------
function unique_sups(obs_list)
    sups = Vector{Vector{Tuple{Int,Int}}}()
    wts = Float64[]
    idx = Dict{Vector{Tuple{Int,Int}},Int}()
    for o in obs_list, t in o.terms
        isempty(t.sup) && continue
        if haskey(idx, t.sup)
            wts[idx[t.sup]] += abs(t.coeff)
        else
            push!(sups, t.sup); push!(wts, abs(t.coeff)); idx[t.sup] = length(sups)
        end
    end
    return sups, wts, idx
end

function greedy_settings(sups, wts, Nw, nsettings; alpha=1.0)
    hits = zeros(Int, length(sups))
    settings = Vector{Vector{Int}}()
    for _ in 1:nsettings
        basis = zeros(Int, Nw)
        for q in 1:Nw
            bestβ, bestsc = 1, -Inf
            for β in 1:3
                basis[q] = β
                sc = 0.0
                for (k, sup) in enumerate(sups)
                    p = 1.0; ok = true
                    for (qq, a) in sup
                        if qq <= q
                            (basis[qq] == a) || (ok = false; break)
                        else
                            p /= 3
                        end
                    end
                    ok || continue
                    sc += wts[k] * exp(-alpha * hits[k]) * p
                end
                sc > bestsc && (bestsc = sc; bestβ = β)
            end
            basis[q] = bestβ
        end
        for (k, sup) in enumerate(sups)
            all(basis[q] == a for (q, a) in sup) && (hits[k] += 1)
        end
        push!(settings, copy(basis))
    end
    return settings, hits
end

# derandomized推定: 命中設定のショットをプールした符号積平均 (3^|A|係数なし)
sign_prod(sup, bits) = prod(1 - 2*bits[q] for (q, _) in sup; init=1)

# ------------------------------------------------------------
# 7. ベンチマーク本体
# ------------------------------------------------------------
function run_benchmark(setname, obs_w, wtens, cum_pl, Pσ, trOσ; nsettings, nm, n_repeat, seed)
    Nw = length(wtens)
    nobs = length(obs_w)
    sups, wts, supidx = unique_sups(obs_w)
    println("  [$setname] observables=$nobs, unique Pauli strings=$(length(sups))")

    # derandomization の設定列 (決定論的; 一度だけ構築)
    dsettings, hits = greedy_settings(sups, wts, Nw, nsettings)
    ncov = count(>(0), hits)
    @printf("  [%s] greedy coverage: %d/%d strings hit (min hits=%d, median=%.0f)\n",
            setname, ncov, length(sups), minimum(hits), median(hits))

    rot = make_rot_buffers(wtens)
    chimax = max(length(cum_pl), maximum(size(t[1], 2) for t in wtens))
    vbuf = zeros(ComplexF64, chimax + 1); wbuf = zeros(ComplexF64, 2chimax + 2)
    bits = zeros(Int, Nw)

    est = zeros(n_repeat, nobs, 3)    # 1=shadow std, 2=shadow CRM, 3=derand
    Random.seed!(seed)
    for rep in 1:n_repeat
        # --- (1)(2) ランダムシャドウ ---
        accS = zeros(nobs); accC = zeros(nobs)
        for _ in 1:nsettings
            basis = rand(1:3, Nw)
            rotate!(rot, wtens, basis)
            xs = zeros(nobs)
            for _ in 1:nm
                sample_bits!(bits, rot, vbuf, wbuf, cum_pl)
                for k in 1:nobs
                    xs[k] += estimate_obs(obs_w[k], basis, bits)
                end
            end
            for k in 1:nobs
                m = xs[k] / nm
                accS[k] += m
                accC[k] += m - exact_prior_mean(obs_w[k], basis, Pσ[k])
            end
        end
        est[rep, :, 1] .= accS ./ nsettings
        est[rep, :, 2] .= accC ./ nsettings .+ trOσ
        # --- (3) derandomization (決定論的設定列で新データ取得) ---
        ssum = zeros(length(sups)); scnt = zeros(Int, length(sups))
        for basis in dsettings
            rotate!(rot, wtens, basis)
            live = [k for k in 1:length(sups) if all(basis[q] == a for (q, a) in sups[k])]
            for _ in 1:nm
                sample_bits!(bits, rot, vbuf, wbuf, cum_pl)
                for k in live
                    ssum[k] += sign_prod(sups[k], bits); scnt[k] += 1
                end
            end
        end
        for k in 1:nobs
            v = 0.0
            for t in obs_w[k].terms
                if isempty(t.sup)
                    v += t.coeff
                else
                    ki = supidx[t.sup]
                    v += scnt[ki] > 0 ? t.coeff * ssum[ki] / scnt[ki] : 0.0
                end
            end
            est[rep, k, 3] = v
        end
    end
    return est
end

function summarize(setname, est, obs_w, truth)
    nobs = length(obs_w)
    err = [std(est[:, k, e]) for k in 1:nobs, e in 1:3]
    bias = [abs(mean(est[:, k, e]) - truth[k]) for k in 1:nobs, e in 1:3]
    @printf("\n  ===== %s: 誤差まとめ (中央値 / 最悪値 over %d observables) =====\n", setname, nobs)
    @printf("  %-24s %10s %10s %10s\n", "method", "median", "max", "max|bias|")
    for (e, lab) in enumerate(("shadow std", "shadow CRM", "derandomized"))
        @printf("  %-24s %10.4f %10.4f %10.4f\n",
                lab, median(err[:, e]), maximum(err[:, e]), maximum(bias[:, e]))
    end
    return err
end

function main()
    t, U = 1.0, 4.0
    L = 16; i0 = L ÷ 2
    nsettings, nm, n_repeat = 50, 100, 30

    E, ψ = ground_state(L, t, U, U/2; chi_max=96)
    @printf("DMRG: E0=%.6f\n", E)
    exp_tens = extract_tensors(ψ)
    σ = truncate(ψ; maxdim=32); normalize!(σ)
    σ_tens = extract_tensors(σ)
    @printf("prior fidelity (chi=32): %.4f\n", abs2(inner(σ, ψ))); flush(stdout)

    rows = []
    for (setname, obs) in [("SetA(structured)", build_setA(i0)),
                           ("SetB(random)", build_setB(qup(i0-4), qdn(i0+4); M=150))]
        q1, q2 = obs_support(obs)
        obs_w = shift_obs(obs, q1 - 1)
        wtens, cum_pl = extract_window(ψ, q1, q2)
        Pρ = [[product_op_expect(exp_tens, term_matrixdict(tm)) for tm in o.terms] for o in obs]
        truth = [sum(tm.coeff * Pρ[k][ti] for (ti, tm) in enumerate(obs[k].terms)) for k in 1:length(obs)]
        Pσ = [[product_op_expect(σ_tens, term_matrixdict(tm)) for tm in o.terms] for o in obs]
        trOσ = [sum(tm.coeff * Pσ[k][ti] for (ti, tm) in enumerate(obs[k].terms)) for k in 1:length(obs)]
        tstart = time()
        est = run_benchmark(setname, obs_w, wtens, cum_pl, Pσ, trOσ;
                            nsettings, nm, n_repeat, seed=777)
        @printf("  [%s] run: %.0fs\n", setname, time() - tstart)
        err = summarize(setname, est, obs, truth)
        for k in 1:length(obs)
            push!(rows, (setname, obs[k].name, truth[k], err[k,1], err[k,2], err[k,3],
                         abs(mean(est[:,k,3]) - truth[k])))
        end
        flush(stdout)
    end

    out = joinpath(@__DIR__, "crm_derand_results.tsv")
    open(out, "w") do io
        println(io, "set\tobservable\ttrue\terr_shadow\terr_crm\terr_derand\tbias_derand")
        for r in rows; println(io, join(r, "\t")); end
    end
    println("\nresults saved: $out")
    return rows
end

rows = main()
println("done")
