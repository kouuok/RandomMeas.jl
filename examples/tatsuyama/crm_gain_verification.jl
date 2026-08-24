# ============================================================
# CRMシャドウの分散利得検証（古典側厳密化版）
#
# 目的:
#   1. 従来実装（prior側もnmショットでサンプル）と修正版
#      （prior側を厳密計算）のCRM推定器を比較する
#   2. 局所Pauli観測量について、分散の閉形式理論
#        Var_std  = [(3^|A|-1)<P>_ρ² + 3^|A|(1-<P>_ρ²)/nm] / nu
#        Var_CRM  = [(3^|A|-1)Δ²     + 3^|A|(1-<P>_ρ²)/nm] / nu,  Δ = <P>_ρ - <P>_σ
#      と経験分散を突き合わせ、「利得は prior の局所誤差 Δ だけで決まる」
#      ことを数値的に検証する
#
# 系: 4サイト・ハバードプラケット (8量子ビット, JW)
#     実験状態 ρ = ED基底状態 / prior σ = HF Slater行列式
# ============================================================

using LinearAlgebra
using Statistics
using Printf
using Random

const NQ  = 8
const DIM = 256

# ------------------------------------------------------------
# 1. 演算子・ハミルトニアン（既存ノートブックと同じ構成）
#    make_matrix(ops...): ops[1] が最上位ビット(MSB)に作用
#    JW順序: "1u"がLSB(位置8), "4d"がMSB(位置1)
# ------------------------------------------------------------
make_matrix(ops...) = reduce(kron, ops)

function get_creation_ops()
    cdag = [0.0 0.0; 1.0 0.0]
    I2   = [1.0 0.0; 0.0 1.0]
    F    = [1.0 0.0; 0.0 -1.0]
    Cdag = Dict{String,Matrix{Float64}}()
    labels = ["1u","1d","2u","2d","3u","3d","4u","4d"]  # k番目がLSB側からk番目
    for (k, lab) in enumerate(labels)
        ops = Matrix{Float64}[]
        for pos in 1:NQ                       # pos=1 が MSB
            q = NQ - pos + 1                  # LSBから数えた番号
            if q < k;      push!(ops, I2)
            elseif q == k; push!(ops, cdag)
            else;          push!(ops, F)
            end
        end
        Cdag[lab] = make_matrix(ops...)
    end
    return Cdag
end

# ラベル -> make_matrix引数位置 (1=MSB) : "1u"->8, "1d"->7, ..., "4d"->1
qpos(lab::String) = NQ - findfirst(==(lab), ["1u","1d","2u","2d","3u","3d","4u","4d"]) + 1

function build_hamiltonian(t_L, t_S, U, mu)
    Cdag = get_creation_ops()
    C    = Dict(k => Matrix(v') for (k, v) in Cdag)
    Nop  = Dict(k => Cdag[k] * C[k] for k in keys(Cdag))
    H = zeros(DIM, DIM)
    edges_L = [(1,2),(3,4),(4,1)]
    edges_S = [(2,3)]
    for s in ("u","d")
        for (i,j) in edges_L
            H .+= -t_L .* (Cdag["$i$s"]*C["$j$s"] .+ Cdag["$j$s"]*C["$i$s"])
        end
        for (i,j) in edges_S
            H .+= -t_S .* (Cdag["$i$s"]*C["$j$s"] .+ Cdag["$j$s"]*C["$i$s"])
        end
    end
    for i in 1:4
        H .+= U .* (Nop["$(i)u"]*Nop["$(i)d"])
    end
    for op in values(Nop)
        H .+= -mu .* op
    end
    return H
end

# ------------------------------------------------------------
# 2. HF prior（既存ノートブックの動的HFを再利用）
# ------------------------------------------------------------
function solve_hf_orbitals(t_L, t_S, U, mu, Nup, Ndn)
    T = zeros(4,4)
    for (i,j) in [(1,2),(3,4),(4,1)]; T[i,j] = T[j,i] = -t_L; end
    for (i,j) in [(2,3)];             T[i,j] = T[j,i] = -t_S; end
    n_up, n_dn = zeros(4), zeros(4)
    for i in 1:Nup; n_up[i] = 0.8; end
    for i in 1:Ndn; n_dn[i] = 0.8; end
    if Nup > 0; n_up[1] += 0.1; end
    if Ndn > 0; n_dn[1] -= 0.1; end
    ev_up, ev_dn = zeros(4,4), zeros(4,4)
    for _ in 1:500
        F_up = eigen(Symmetric(T + diagm(U .* n_dn) - mu*I))
        F_dn = eigen(Symmetric(T + diagm(U .* n_up) - mu*I))
        ev_up, ev_dn = F_up.vectors, F_dn.vectors
        new_up = sum(ev_up[:,i].^2 for i in 1:Nup; init=zeros(4))
        new_dn = sum(ev_dn[:,i].^2 for i in 1:Ndn; init=zeros(4))
        if maximum(abs.(new_up .- n_up)) < 1e-10; break; end
        n_up = 0.5 .* new_up .+ 0.5 .* n_up
        n_dn = 0.5 .* new_dn .+ 0.5 .* n_dn
    end
    return ev_up, ev_dn
end

function build_hf_state(ev_up, ev_dn, Nup, Ndn)
    Cdag = get_creation_ops()
    D_up = [sum(ev_up[s,k] * Cdag["$(s)u"] for s in 1:4) for k in 1:4]
    D_dn = [sum(ev_dn[s,k] * Cdag["$(s)d"] for s in 1:4) for k in 1:4]
    v = zeros(DIM); v[1] = 1.0
    for k in reverse(1:Ndn); v = D_dn[k] * v; end
    for k in reverse(1:Nup); v = D_up[k] * v; end
    return v ./ norm(v)
end

# ------------------------------------------------------------
# 3. シャドウ測定の機構
#    基底 1=X, 2=Y, 3=Z を測るユニタリ u_β（u†|b><b|u が測定基底）
# ------------------------------------------------------------
const SIGMA = (ComplexF64[0 1; 1 0], ComplexF64[0 -im; im 0], ComplexF64[1 0; 0 -1])
const UBASIS = (ComplexF64[1 1; 1 -1]/sqrt(2),      # X基底
                ComplexF64[1 -im; 1 im]/sqrt(2),    # Y基底
                ComplexF64[1 0; 0 1])               # Z基底
# DVAL[β][a][b+1] = Re <b| u_β σ_a u_β† |b>  （基底βでPauli aを測った際の値）
const DVAL = ntuple(β -> ntuple(a -> real.(diag(UBASIS[β]*SIGMA[a]*UBASIS[β]')), 3), 3)

# 状態ベクトルの qubit pos (1=MSB) に 2x2 行列を作用
function apply_1q!(out::Vector{ComplexF64}, ψ::Vector{ComplexF64}, u::Matrix{ComplexF64}, pos::Int)
    shift = NQ - pos                    # そのqubitのビット位置
    mask  = 1 << shift
    @inbounds for s in 0:DIM-1
        if s & mask == 0
            a, b = ψ[s+1], ψ[s+mask+1]
            out[s+1]      = u[1,1]*a + u[1,2]*b
            out[s+mask+1] = u[2,1]*a + u[2,2]*b
        end
    end
    return out
end

function rotated_probs(ψ::Vector{ComplexF64}, basis::Vector{Int})
    φ = copy(ψ); tmp = similar(φ)
    for pos in 1:NQ
        apply_1q!(tmp, φ, UBASIS[basis[pos]], pos)
        φ, tmp = tmp, φ
    end
    return abs2.(φ)
end

sample_outcome(cq::Vector{Float64}) = searchsortedfirst(cq, rand()) - 1  # 0-based

# W[b+1,s+1] = Π_i (3δ_{b_i s_i} - 1) : 忠実度推定量の基礎行列（uに依存しない）
function build_W()
    W = Matrix{Float64}(undef, DIM, DIM)
    for b in 0:DIM-1, s in 0:DIM-1
        h = count_ones(b ⊻ s)
        W[b+1, s+1] = 2.0^(NQ - h) * (-1.0)^h
    end
    return W
end

# ------------------------------------------------------------
# 4. 観測量の定義（Pauli和）
#    term = (係数, Dict(qubit位置 => Pauli番号 1=X,2=Y,3=Z))
# ------------------------------------------------------------
struct PauliObs
    name::String
    terms::Vector{Tuple{Float64,Dict{Int,Int}}}
end

function dense_matrix(obs::PauliObs)
    M = zeros(ComplexF64, DIM, DIM)
    for (c, sup) in obs.terms
        ops = [haskey(sup, p) ? SIGMA[sup[p]] : ComplexF64[1 0; 0 1] for p in 1:NQ]
        M .+= c .* make_matrix(ops...)
    end
    return M
end

# アウトカム s (0-based) に対する1スナップショット推定値 X(s)
function pauli_estimate(obs::PauliObs, basis::Vector{Int}, s::Int)
    x = 0.0
    for (c, sup) in obs.terms
        t = c
        for (p, a) in sup
            bit = (s >> (NQ - p)) & 1
            t *= 3.0 * DVAL[basis[p]][a][bit+1]
            t == 0.0 && break
        end
        x += t
    end
    return x
end

# 分布 q の下での厳密な条件付き期待値 E[X|u]（prior側の厳密計算に使用）
function pauli_exact_mean(obs::PauliObs, basis::Vector{Int}, q::Vector{Float64})
    m = 0.0
    @inbounds for s in 0:DIM-1
        q[s+1] == 0.0 && continue
        m += q[s+1] * pauli_estimate(obs, basis, s)
    end
    return m
end

# ------------------------------------------------------------
# 5. メイン: 3推定器の同時シミュレーション
#    返り値: est[repeat, obs, estimator]  (1=標準, 2=CRM旧(サンプル), 3=CRM新(厳密))
# ------------------------------------------------------------
function run_experiment(ψρ, ψσ, obs_list, W, tr_Oσ, fid_idx;
                        nu::Int, nm::Int, n_repeat::Int, rng_seed::Int)
    Random.seed!(rng_seed)
    nobs = length(obs_list) + 1                     # +1 = 忠実度
    est = zeros(n_repeat, nobs, 3)
    ψρc = ComplexF64.(ψρ); ψσc = ComplexF64.(ψσ)
    for rep in 1:n_repeat
        acc = zeros(nobs, 3)
        for _ in 1:nu
            basis = rand(1:3, NQ)
            qρ = rotated_probs(ψρc, basis)
            qσ = rotated_probs(ψσc, basis)
            cqρ = cumsum(qρ); cqσ = cumsum(qσ)
            v = W * qσ                              # 忠実度: X(s) = v[s+1]
            sρ = [sample_outcome(cqρ) for _ in 1:nm]
            sσ = [sample_outcome(cqσ) for _ in 1:nm]
            for (k, obs) in enumerate(obs_list)
                mρ  = mean(pauli_estimate(obs, basis, s) for s in sρ)
                mσs = mean(pauli_estimate(obs, basis, s) for s in sσ)
                mσe = pauli_exact_mean(obs, basis, qσ)
                acc[k,1] += mρ
                acc[k,2] += mρ - mσs
                acc[k,3] += mρ - mσe
            end
            # 忠実度
            mρ  = mean(v[s+1] for s in sρ)
            mσs = mean(v[s+1] for s in sσ)
            mσe = dot(qσ, v)
            acc[fid_idx,1] += mρ
            acc[fid_idx,2] += mρ - mσs
            acc[fid_idx,3] += mρ - mσe
        end
        for k in 1:nobs
            est[rep,k,1] = acc[k,1]/nu
            est[rep,k,2] = acc[k,2]/nu + tr_Oσ[k]   # CRMは Tr[Oσ] を加算
            est[rep,k,3] = acc[k,3]/nu + tr_Oσ[k]
        end
    end
    return est
end

# 局所Pauli(単一列)の理論分散
function theory_var(absA::Int, Pρ::Float64, Δ::Float64, nu::Int, nm::Int)
    shot = 3.0^absA * (1 - Pρ^2) / nm
    v_std = ((3.0^absA - 1) * Pρ^2 + shot) / nu
    v_crm = ((3.0^absA - 1) * Δ^2  + shot) / nu
    return v_std, v_crm
end

# ------------------------------------------------------------
# 6. 実行
# ------------------------------------------------------------
"""1つのホッピング配置について全 (U, 設定) を回す。

t_L=t_S=1 が本編(一様鎖)。t_S=0.5 は「結合を1本弱めた場合でも同じ結論か」を
確かめるための補足で、以前はこちらが本編だった(弱化の理由づけが誤っていた
経緯は README §1 を参照)。乱数種は配置に依存させないので、両者の違いは
ハミルトニアンの違いだけから来る。"""
function run_config(t_L, t_S, tag, obs_list, W, fid_idx, names, U_list, settings, n_repeat)
    results = Dict()
    for U in U_list
        mu = U / 2               # half-filling
        println("="^70)
        @printf("[%s] U = %.1f  (t_L=%.1f, t_S=%.1f, mu=U/2)\n", tag, U, t_L, t_S)

        H = build_hamiltonian(t_L, t_S, U, mu)
        F = eigen(Symmetric(H))
        ψρ = F.vectors[:, 1]
        gap = F.values[2]-F.values[1]
        @printf("  ED: E0=%.6f, gap=%.4f\n", F.values[1], gap)
        @assert gap > 1e-6 "基底状態が縮退している (tag=$tag, U=$U)"

        Cdag = get_creation_ops()
        Nup_op = sum(Cdag["$(i)u"]*Cdag["$(i)u"]' for i in 1:4)
        Ndn_op = sum(Cdag["$(i)d"]*Cdag["$(i)d"]' for i in 1:4)
        Nup = round(Int, dot(ψρ, Nup_op*ψρ)); Ndn = round(Int, dot(ψρ, Ndn_op*ψρ))
        @printf("  particle sector: (N_up, N_dn) = (%d, %d)\n", Nup, Ndn)

        ev_up, ev_dn = solve_hf_orbitals(t_L, t_S, U, mu, Nup, Ndn)
        ψσ = build_hf_state(ev_up, ev_dn, Nup, Ndn)
        F_true = abs(dot(ψσ, ψρ))^2
        @printf("  HF prior fidelity: %.4f\n", F_true)

        Oρ = Float64[]; Oσ = Float64[]
        for obs in obs_list
            M = dense_matrix(obs)
            push!(Oρ, real(dot(ψρ, M*ψρ))); push!(Oσ, real(dot(ψσ, M*ψσ)))
        end
        push!(Oρ, F_true); push!(Oσ, 1.0)

        # 検証: 高速Pauli推定器 vs 密行列スナップショット
        Random.seed!(1)
        basis_t = rand(1:3, NQ); s_t = rand(0:DIM-1)
        ρ̂ = ones(ComplexF64,1,1)
        for pos in 1:NQ
            b = (s_t >> (NQ-pos)) & 1
            u = UBASIS[basis_t[pos]]
            sb = b == 0 ? ComplexF64[1,0] : ComplexF64[0,1]
            ρ̂ = kron(ρ̂, 3 .* (u'*(sb*sb')*u) .- ComplexF64[1 0; 0 1])
        end
        for (k, obs) in enumerate(obs_list)
            @assert isapprox(pauli_estimate(obs, basis_t, s_t),
                             real(tr(dense_matrix(obs)*ρ̂)); atol=1e-9) "fast-path mismatch"
        end

        for (si, st) in enumerate(settings)
            est = run_experiment(ψρ, ψσ, obs_list, W, Oσ, fid_idx;
                                 nu=st.nu, nm=st.nm, n_repeat=n_repeat,
                                 rng_seed=1000*si + round(Int, 10U))
            results[(U, st.nu, st.nm)] = (est=est, Oρ=Oρ, Oσ=Oσ, F_true=F_true, gap=gap)

            @printf("\n  --- nu=%d, nm=%d (total shots %d), n_repeat=%d ---\n",
                    st.nu, st.nm, st.nu*st.nm, n_repeat)
            @printf("  %-18s %9s | %-19s %-19s %-19s | %7s %7s %8s\n",
                    "observable", "true", "std shadow", "CRM old(sampled)", "CRM new(exact)",
                    "G_emp", "G_theo", "Delta")
            for k in 1:fid_idx
                m = [mean(est[:,k,e]) for e in 1:3]
                sd = [std(est[:,k,e])  for e in 1:3]
                G_emp = (sd[3] > 0) ? (sd[1]/sd[3])^2 : NaN
                G_theo = NaN; Δ = Oρ[k] - Oσ[k]
                if k <= length(obs_list) && length(obs_list[k].terms) == 1
                    absA = length(obs_list[k].terms[1][2])
                    vs, vc = theory_var(absA, Oρ[k], Δ, st.nu, st.nm)
                    G_theo = vs / vc
                end
                @printf("  %-18s %9.4f | %8.4f ± %-8.4f %8.4f ± %-8.4f %8.4f ± %-8.4f | %7.2f %7.2f %8.4f\n",
                        names[k], Oρ[k], m[1], sd[1], m[2], sd[2], m[3], sd[3],
                        G_emp, G_theo, Δ)
            end
        end
    end
    return results
end

function main()
    # ホッピング配置。**本編は一様 (t_L=t_S=1)**。
    # t_S=0.5(結合 2--3 を1本だけ弱める)は補足で、以前はこちらを本編にしていた。
    # 当初のコメントは「二量体化で縮退を回避」としていたが誤りで、4サイトリングは
    # 二部格子(|A|=|B|=2)なので Lieb の定理から半充填・U>0 の基底状態は一意な
    # singlet であり、一様でも縮退しない(E1-E0 = 4.8e-2 @U=1, 3.3e-1 @U=8)。
    # 縮退するのは U=0 の開殻だけで、ここでは使わない。詳細は
    # crm_plaquette_degeneracy.jl と README §1。
    configs  = [("uniform", 1.0, 1.0), ("weak", 1.0, 0.5)]
    U_list   = [1.0, 8.0]        # 弱相関 / 強相関 (HF priorの質が変わる)
    settings = [(nu=50, nm=10), (nu=50, nm=100), (nu=50, nm=1000)]
    n_repeat = 200

    W = build_W()

    p1u, p1d = qpos("1u"), qpos("1d")
    p3u, p3d = qpos("3u"), qpos("3d")
    obs_list = PauliObs[
        PauliObs("Z(1u)  |A|=1", [(1.0, Dict(p1u=>3))]),
        PauliObs("ZZ(1u,1d) |A|=2", [(1.0, Dict(p1u=>3, p1d=>3))]),
        PauliObs("ZZZZ(1,3) |A|=4", [(1.0, Dict(p1u=>3, p1d=>3, p3u=>3, p3d=>3))]),
        PauliObs("DoubleOcc site1", [(0.25, Dict{Int,Int}()), (-0.25, Dict(p1u=>3)),
                                     (-0.25, Dict(p1d=>3)), (0.25, Dict(p1u=>3, p1d=>3))]),
        PauliObs("Sz1*Sz3", [( 1/16, Dict(p1d=>3, p3d=>3)), (-1/16, Dict(p1d=>3, p3u=>3)),
                             (-1/16, Dict(p1u=>3, p3d=>3)), ( 1/16, Dict(p1u=>3, p3u=>3))]),
    ]
    fid_idx = length(obs_list) + 1
    names = vcat([o.name for o in obs_list], "Fidelity")

    all_results = Dict{String,Any}()
    for (tag, tL, tS) in configs
        all_results[tag] = run_config(tL, tS, tag, obs_list, W, fid_idx, names,
                                      U_list, settings, n_repeat)
    end
    return all_results, names, obs_list, settings, n_repeat, configs
end

all_results, names, obs_list, settings, n_repeat, configs = main()

# ------------------------------------------------------------
# 7. 結果の書き出し (README/論文の表はここから作る)
# ------------------------------------------------------------
open(joinpath(@__DIR__, "crm_gain_verification_results.tsv"), "w") do io
    println(io, "config\tt_L\tt_S\tU\tnu\tnm\tobservable\ttrue\tprior\tDelta\t" *
                "sd_std\tsd_crm_old\tsd_crm_new\tG_emp\tG_theo\tHF_fid\tgap")
    for (tag, tL, tS) in configs, U in [1.0, 8.0], st in settings
        r = all_results[tag][(U, st.nu, st.nm)]
        for k in 1:length(names)
            sd = [std(r.est[:,k,e]) for e in 1:3]
            G_emp = sd[3] > 0 ? (sd[1]/sd[3])^2 : NaN
            G_theo = NaN
            if k <= length(obs_list) && length(obs_list[k].terms) == 1
                absA = length(obs_list[k].terms[1][2])
                vs, vc = theory_var(absA, r.Oρ[k], r.Oρ[k]-r.Oσ[k], st.nu, st.nm)
                G_theo = vs/vc
            end
            println(io, join((tag, tL, tS, U, st.nu, st.nm, names[k],
                              r.Oρ[k], r.Oσ[k], r.Oρ[k]-r.Oσ[k],
                              sd[1], sd[2], sd[3], G_emp, G_theo, r.F_true, r.gap), "\t"))
        end
    end
end
println("\nresults saved: crm_gain_verification_results.tsv")

# 一様 vs 弱化 の比較サマリ
println("\n" * "="^70)
println("一様 (t=1) と 1本弱化 (t_S=0.5) の比較")
@printf("%-6s %-5s %10s %10s | %-18s %10s %10s\n",
        "config","U","ED gap","HF忠実度","observable","G_emp","G_theo")
for (tag, _, _) in configs, U in [1.0, 8.0]
    r = all_results[tag][(U, 50, 100)]
    for (k, o) in enumerate(obs_list)
        length(o.terms) == 1 || continue
        sd = [std(r.est[:,k,e]) for e in 1:3]
        absA = length(o.terms[1][2])
        vs, vc = theory_var(absA, r.Oρ[k], r.Oρ[k]-r.Oσ[k], 50, 100)
        @printf("%-6s %-5.0f %10.4f %10.4f | %-18s %10.2f %10.2f\n",
                tag, U, r.gap, r.F_true, o.name, (sd[1]/sd[3])^2, vs/vc)
    end
end

# ------------------------------------------------------------
# 8. プロット (配置ごとに1枚)
# ------------------------------------------------------------
using Plots

for (tag, tL, tS) in configs
    results = all_results[tag]
    plots_fid = []
    for U in [1.0, 8.0]
        shots = [st.nu*st.nm for st in settings]
        s_std = Float64[]; s_old = Float64[]; s_new = Float64[]
        local F_true = 0.0
        for st in settings
            r = results[(U, st.nu, st.nm)]
            k = length(names)
            push!(s_std, std(r.est[:,k,1])); push!(s_old, std(r.est[:,k,2]))
            push!(s_new, std(r.est[:,k,3])); F_true = r.F_true
        end
        p = plot(shots, s_std, marker=:circle, label="Standard shadow",
                 xscale=:log10, yscale=:log10, lw=2, color=:steelblue)
        plot!(p, shots, s_old, marker=:utriangle, label="CRM old (sampled prior)", lw=2, color=:gray)
        plot!(p, shots, s_new, marker=:square, label="CRM new (exact prior)", lw=2, color=:firebrick)
        title!(p, @sprintf("Fidelity est.  U=%.0f  (HF fid=%.3f)", U, F_true))
        xlabel!(p, "total shots (nu*nm, nu=50)"); ylabel!(p, "std error")
        push!(plots_fid, p)
    end

    gx = Float64[]; gy = Float64[]
    for U in [1.0, 8.0], st in settings
        r = results[(U, st.nu, st.nm)]
        for (k, obs) in enumerate(obs_list)
            length(obs.terms) == 1 || continue
            absA = length(obs.terms[1][2])
            vs, vc = theory_var(absA, r.Oρ[k], r.Oρ[k]-r.Oσ[k], st.nu, st.nm)
            push!(gx, vs/vc); push!(gy, var(r.est[:,k,1]) / var(r.est[:,k,3]))
        end
    end
    lims = (0.5*min(minimum(gx), minimum(gy)), 2*max(maximum(gx), maximum(gy)))
    p3 = scatter(gx, gy, xscale=:log10, yscale=:log10, xlims=lims, ylims=lims,
                 label="single Pauli strings", color=:firebrick, ms=6, alpha=0.8)
    plot!(p3, [lims...], [lims...], label="y = x", color=:black, ls=:dash)
    ttl = tag == "uniform" ? "CRM gain: empirical vs theory (uniform t=1)" :
                             "CRM gain: empirical vs theory (one bond t=0.5)"
    title!(p3, ttl)
    xlabel!(p3, "G theory = Var_std/Var_CRM"); ylabel!(p3, "G empirical")

    final = plot(plots_fid..., p3, layout=(1,3), size=(1500, 420),
                 margin=6Plots.mm, legend=:bottomleft)
    fn = tag == "uniform" ? "crm_gain_verification.png" : "crm_gain_verification_weak.png"
    out = joinpath(@__DIR__, fn)
    savefig(final, out)
    println("Figure saved: $out")
end
