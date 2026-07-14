# ============================================================
# CRMシャドウのパラメータ掃引 (nm依存 / U・ドーピング依存)
#
# 実験A (nm掃引, L=8):
#   利得の飽和則 G_max(nm) = (V_u + v_s(nm)) / ((3^|A|-1)Δ² + v_s(nm))
#   を nm = 10..1000 で検証する。CRMは「少ない設定×多ショット」配分で
#   最大化される、という測定設計指針の定量化。
#
# 実験B (U・ドーピング掃引, L=16):
#   U = 1, 4, 8 (half-filling) と U=4 ドープ (n=7/8) で、
#   prior品質(χ切断)ごとの局所観測量CRM利得を比較する。
#
# 共通機構は crm_mps_scaling.jl と同一（検証済みコードの複製）。
# 実行方法:
#  JULIA_LOAD_PATH="@:@v#.#:@stdlib" julia --project=Hubbard_MPS_Env_v2 crm_param_sweep.jl
# ============================================================

using ITensors, ITensorMPS
using LinearAlgebra
using Statistics
using Printf
using Random

BLAS.set_num_threads(1)

# ------------------------------------------------------------
# 1. モデルとDMRG (crm_mps_scaling.jl と同一 + 粒子数指定)
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

# n_elec 個の電子 (Sz=0付近) を最初の n_elec サイトに交互スピンで配置
function ground_state(L, t, U, mu; chi_max=96, nsweeps=10, n_elec=L)
    N = 2L
    sites = siteinds("Fermion", N; conserve_qns=true)
    H = hubbard_chain_mpo(sites, L, t, U, mu)
    init = fill("Emp", N)
    for i in 1:n_elec
        init[isodd(i) ? qup(i) : qdn(i)] = "Occ"
    end
    ψ0 = productMPS(sites, init)
    maxdims = min.(chi_max, [20, 50, 100, chi_max, chi_max, chi_max, chi_max, chi_max, chi_max, chi_max])[1:nsweeps]
    noise = [1e-6, 1e-6, 1e-7, 1e-8, 1e-9, 0, 0, 0, 0, 0][1:nsweeps]
    E, ψ = dmrg(H, ψ0; nsweeps, maxdim=maxdims, cutoff=1e-10, noise, outputlevel=0)
    return E, ψ
end

# ------------------------------------------------------------
# 2. テンソル抽出・期待値・サンプラー (crm_mps_scaling.jl と同一)
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
# 3. 観測量 (crm_mps_scaling.jl と同一)
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
# 4. 1パラメータ点の実験一式 (状態構築 -> prior -> run_locals)
# ------------------------------------------------------------
function prepare_point(L, t, U, mu; chi_exp, chi_priors, n_elec=L)
    E, ψ = ground_state(L, t, U, mu; chi_max=chi_exp, n_elec)
    exp_tens = extract_tensors(ψ)
    obs = build_observables(L)
    q1, q2 = obs_support(obs)
    obs_w = shift_obs(obs, q1 - 1)
    wtens, cum_pl = extract_window(ψ, q1, q2)
    Pρ = [[product_op_expect(exp_tens, term_matrixdict(tm)) for tm in o.terms] for o in obs]
    Otrue = [sum(tm.coeff * Pρ[k][ti] for (ti, tm) in enumerate(obs[k].terms)) for k in 1:length(obs)]
    Pσ_all = Vector{Vector{Vector{Float64}}}(); trOσ_all = Vector{Vector{Float64}}()
    fids = Float64[]
    for chi_p in chi_priors
        σ = truncate(ψ; maxdim=chi_p); normalize!(σ)
        push!(fids, abs2(inner(σ, ψ)))
        st = extract_tensors(σ)
        Pσ = [[product_op_expect(st, term_matrixdict(tm)) for tm in o.terms] for o in obs]
        push!(Pσ_all, Pσ)
        push!(trOσ_all, [sum(tm.coeff * Pσ[k][ti] for (ti, tm) in enumerate(obs[k].terms)) for k in 1:length(obs)])
    end
    return (; E, obs, obs_w, wtens, cum_pl, Otrue, Pσ_all, trOσ_all, fids)
end

# ------------------------------------------------------------
# 5. メイン
# ------------------------------------------------------------
function main()
    t = 1.0
    chi_priors = [2, 4, 8, 16, 32]
    nu = 50
    rows_A = []; rows_B = []

    # ===== 実験A: nm掃引 (L=8, U=4, half-filling) =====
    println("="^70)
    println("Experiment A: nm sweep (L=8, U=4, half-filling)"); flush(stdout)
    L, U = 8, 4.0
    pt = prepare_point(L, t, U, U/2; chi_exp=64, chi_priors)
    @printf("  E0=%.6f  prior fids: %s\n", pt.E,
            join([@sprintf("chi=%d: %.4f", c, f) for (c,f) in zip(chi_priors, pt.fids)], ", "))
    flush(stdout)
    nm_list = [10, 30, 100, 300, 1000]
    for nm in nm_list
        n_repeat = nm >= 1000 ? 50 : 100
        tstart = time()
        est_std, est_crm = run_locals(pt.wtens, pt.cum_pl, pt.obs_w, pt.Pσ_all, pt.trOσ_all;
                                      nu, nm, n_repeat, seed=7000 + nm)
        @printf("  nm=%4d (n_rep=%d): %.0fs\n", nm, n_repeat, time()-tstart); flush(stdout)
        for (k, o) in enumerate(pt.obs)
            v_std = var(est_std[:, k])
            for (p, chi_p) in enumerate(chi_priors)
                G = v_std / var(est_crm[:, k, p])
                Δ = pt.Otrue[k] - pt.trOσ_all[p][k]
                G_theo = NaN
                if o.pure
                    vs, vc = theory_var(length(o.terms[1].sup), pt.Otrue[k], Δ, nu, nm)
                    G_theo = vs / vc
                end
                push!(rows_A, (nm, chi_p, o.name, o.pure, pt.Otrue[k], Δ, G, G_theo))
            end
        end
    end

    # ===== 実験B: U・ドーピング掃引 (L=16) =====
    println("\n", "="^70)
    println("Experiment B: U & doping sweep (L=16, chi_exp=96)"); flush(stdout)
    L = 16
    nm, n_repeat = 100, 60
    points = [(U=1.0, n_elec=16, label="U=1 half"),
              (U=4.0, n_elec=16, label="U=4 half"),
              (U=8.0, n_elec=16, label="U=8 half"),
              (U=4.0, n_elec=14, label="U=4 doped(7/8)")]
    for pset in points
        tstart = time()
        pt = prepare_point(L, t, pset.U, pset.U/2; chi_exp=96, chi_priors, n_elec=pset.n_elec)
        @printf("\n  --- %s ---  E0=%.6f  (prep %.0fs)\n", pset.label, pt.E, time()-tstart)
        @printf("  prior fids: %s\n",
                join([@sprintf("chi=%d: %.4f", c, f) for (c,f) in zip(chi_priors, pt.fids)], ", "))
        flush(stdout)
        tstart = time()
        est_std, est_crm = run_locals(pt.wtens, pt.cum_pl, pt.obs_w, pt.Pσ_all, pt.trOσ_all;
                                      nu, nm, n_repeat, seed=8000 + round(Int, 10pset.U) + pset.n_elec)
        @printf("  run: %.0fs\n", time()-tstart); flush(stdout)
        @printf("  %-16s %9s |", "observable", "true")
        for c in chi_priors; @printf(" G(chi=%d)", c); end
        println()
        for (k, o) in enumerate(pt.obs)
            v_std = var(est_std[:, k])
            @printf("  %-16s %9.4f |", o.name, pt.Otrue[k])
            for (p, chi_p) in enumerate(chi_priors)
                G = v_std / var(est_crm[:, k, p])
                @printf(" %8.2f", G)
                Δ = pt.Otrue[k] - pt.trOσ_all[p][k]
                push!(rows_B, (pset.label, pset.U, pset.n_elec, chi_p, o.name,
                               pt.Otrue[k], Δ, pt.fids[p], G))
            end
            println()
        end
        flush(stdout)
    end

    # 保存
    outA = joinpath(@__DIR__, "crm_sweep_nm.tsv")
    open(outA, "w") do io
        println(io, "nm\tchi_prior\tobservable\tpure\ttrue\tDelta\tG_emp\tG_theo")
        for r in rows_A; println(io, join(r, "\t")); end
    end
    outB = joinpath(@__DIR__, "crm_sweep_U.tsv")
    open(outB, "w") do io
        println(io, "label\tU\tn_elec\tchi_prior\tobservable\ttrue\tDelta\tprior_fid\tG_emp")
        for r in rows_B; println(io, join(r, "\t")); end
    end
    println("\nresults saved: $outA, $outB"); flush(stdout)
    return rows_A, rows_B, nm_list, chi_priors, points
end

rows_A, rows_B, nm_list, chi_priors, points = main()

# ------------------------------------------------------------
# 6. プロット
# ------------------------------------------------------------
using Plots

# P1: G vs nm (理論曲線を重ねる; 純Pauli列, chi=4と32)
p1 = plot(xlabel="nm (shots per setting)", ylabel="G", xscale=:log10, yscale=:log10,
          legend=:topleft, title="Gain vs shot allocation (L=8, U=4)")
colors = Dict("ZZ onsite(i0)" => :firebrick, "ZZ up-up r=1" => :steelblue)
styles = Dict(4 => :dash, 32 => :solid)
for chi_p in (4, 32), name in keys(colors)
    sel = [(r[1], r[7], r[8]) for r in rows_A if r[2]==chi_p && r[3]==name]
    nms = [s[1] for s in sel]; ge = [s[2] for s in sel]
    plot!(p1, nms, ge, marker=:circle, color=colors[name], ls=styles[chi_p], lw=2,
          label="$name chi=$chi_p")
    # 理論曲線 (連続nm)
    r0 = first(r for r in rows_A if r[2]==chi_p && r[3]==name)
    P = r0[5]; Δ = r0[6]; absA = 2
    nmc = 10 .^ range(1, 3.2; length=50)
    gt = [(vs = theory_var(absA, P, Δ, 50, m); vs[1]/vs[2]) for m in nmc]
    plot!(p1, nmc, gt, color=colors[name], ls=:dot, lw=1, alpha=0.7, label="")
end
hline!(p1, [1.0], color=:black, ls=:dot, label="")

# P2: G経験 vs G理論 (実験A全点)
gx = [r[8] for r in rows_A if r[4] && !isnan(r[8])]
gy = [r[7] for r in rows_A if r[4] && !isnan(r[8])]
lims = (0.5*min(minimum(gx), minimum(gy)), 2*max(maximum(gx), maximum(gy)))
p2 = scatter(gx, gy, xscale=:log10, yscale=:log10, xlims=lims, ylims=lims,
             color=:firebrick, ms=5, alpha=0.7, label="pure Pauli strings",
             xlabel="G theory", ylabel="G empirical", title="nm sweep: theory check")
plot!(p2, [lims...], [lims...], color=:black, ls=:dash, label="y = x")

# P3: G vs U (chi=8; half-filling 3点 + doped)
p3 = plot(xlabel="U", ylabel="G (chi_prior=8)", yscale=:log10, legend=:topleft,
          title="Gain vs interaction & doping (L=16)")
sel_names = ["ZZ onsite(i0)", "DoubleOcc(i0)", "SzSz r=1", "hop bond(i0)"]
for name in sel_names
    Us = Float64[]; gs = Float64[]
    for r in rows_B
        (r[3] == 16 && r[4] == 8 && r[5] == name) || continue
        push!(Us, r[2]); push!(gs, r[9])
    end
    perm = sortperm(Us)
    plot!(p3, Us[perm], gs[perm], marker=:circle, lw=2, label=name)
    # doped点 (U=4) を白抜きマーカーで
    for r in rows_B
        (r[3] == 14 && r[4] == 8 && r[5] == name) || continue
        scatter!(p3, [r[2] + 0.15], [r[9]], marker=:diamond, ms=7, mc=:white,
                 msc=:auto, label="")
    end
end
hline!(p3, [1.0], color=:black, ls=:dot, label="")
annotate!(p3, 4.6, 1.35, text("◇ = doped n=7/8", 8, :gray))

final = plot(p1, p2, p3, layout=(1,3), size=(1500, 430), margin=6Plots.mm)
outpng = joinpath(@__DIR__, "crm_param_sweep.png")
savefig(final, outpng)
println("figure saved: $outpng")
