# ============================================================
# 利得は「窓(部分系)の忠実度」で決まるか
#
# 動機:
#   §2 の主張は「大域忠実度 F は利得を決めない、観測量ごとの局所誤差 Δ が決める」
#   だった。その中間に **窓の縮約密度行列どうしの近さ** という量がある。
#   これが利得を決めるなら、観測量を選ばずに prior を評価できる1つの数になる。
#
# 理論的な関係:
#   Pauli 列 P は ||P||_∞ = 1 なので
#       |Δ| = |Tr[P(ρ_w - σ_w)]| ≤ ||ρ_w - σ_w||_1 = 2 D_tr
#   つまり窓のトレース距離は |Δ| の上界を与え、したがって G の下界を与える。
#   問題はこの不等式がどれだけタイトか。緩ければ「窓の忠実度では利得は決まらない」
#   ことになり、§2 の主張が一段強くなる。
#
#   平均場 prior の窓RDMは、窓内の全 Pauli 列(4^Nw 個)の期待値から
#       σ_w = 2^{-Nw} Σ_P ⟨P⟩_σ P
#   と再構成する。crm_wick_pauli.jl の一般Wickがあるので任意の P を評価できる。
#
# 実行: julia --project=Hubbard_MPS_Env_v2 crm_window_fidelity.jl
# 環境変数: L_LIST(既定 "16,32"), U(既定 4.0), DOPING(既定 "0.0,0.125")
# ============================================================
include(joinpath(@__DIR__, "crm_chain_common.jl"))
include(joinpath(@__DIR__, "crm_wick_pauli.jl"))

"""窓内の全 Pauli 列の期待値から縮約密度行列を再構成する。"""
function rdm_from_paulis(expect_fn, Nw::Int)
    d = 1 << Nw
    ρ = zeros(ComplexF64, d, d)
    for idx in 0:(4^Nw - 1)
        labs = Vector{Int}(undef, Nw); x = idx
        for j in Nw:-1:1; labs[j] = x % 4; x ÷= 4; end     # 0=I,1=X,2=Y,3=Z
        sup = [(q, labs[q]) for q in 1:Nw if labs[q] != 0]
        v = expect_fn(sup)
        v == 0 && continue
        M = reduce(kron, [labs[q] == 0 ? Matrix{ComplexF64}(I,2,2) : SIGMA[labs[q]] for q in 1:Nw])
        ρ .+= v .* M
    end
    Hermitian(ρ ./ d)
end

"""Uhlmann 忠実度 F(ρ,σ) = (Tr√(√ρ σ √ρ))^2 と トレース距離。"""
function fid_trdist(ρ, σ)
    A = Matrix(Hermitian((ρ+ρ')/2)); B = Matrix(Hermitian((σ+σ')/2))
    Fa = eigen(Hermitian(A)); sq = Fa.vectors * Diagonal(sqrt.(max.(Fa.values,0))) * Fa.vectors'
    M = Hermitian(sq*B*sq)
    F = (sum(sqrt.(max.(eigen(M).values, 0))))^2
    D = 0.5*sum(abs.(eigen(Hermitian(A-B)).values))
    return real(F), real(D)
end

function main()
    t = 1.0
    U = parse(Float64, get(ENV, "U", "4.0"))
    L_list = parse.(Int, split(get(ENV, "L_LIST", "16,32"), ","))
    dops   = parse.(Float64, split(get(ENV, "DOPING", "0.0,0.125"), ","))
    chi_priors = [2, 4, 8, 16, 32]
    nm = 100      # 天井の計算に使う(測定はしない)

    rows = []
    for L in L_list, δ in dops
        E, ψ, _, _, nel = ground_state(L, t, U, U/2; chi_max=256, nsweeps=20, δ)
        Nup = (nel+1)÷2; Ndn = nel÷2
        @printf("\n%s\nL=%d δ=%.3f  E0=%.6f chi=%d N_el=%d\n", "="^92, L, δ, E, maxlinkdim(ψ), nel)
        flush(stdout)

        # UHF(共線)と対称性回復UHF
        T = zeros(L,L); for i in 1:L-1; T[i,i+1]=T[i+1,i]=-t; end
        nup = [ (Nup+Ndn)/(2L) + 0.4*(-1)^i for i in 1:L]
        ndn = [ (Nup+Ndn)/(2L) - 0.4*(-1)^i for i in 1:L]
        clamp!(nup,0.02,0.98); clamp!(ndn,0.02,0.98)
        nup .*= Nup/sum(nup); ndn .*= Ndn/sum(ndn)
        local Fu, Fd
        for _ in 1:4000
            Fu = eigen(Symmetric(T + diagm(U .* ndn))); Fd = eigen(Symmetric(T + diagm(U .* nup)))
            nu2 = vec(sum(abs2, Fu.vectors[:,1:Nup]; dims=2)); nd2 = vec(sum(abs2, Fd.vectors[:,1:Ndn]; dims=2))
            if max(maximum(abs.(nu2.-nup)), maximum(abs.(nd2.-ndn))) < 1e-12
                nup, ndn = nu2, nd2; break
            end
            nup = 0.5.*nu2 .+ 0.5.*nup; ndn = 0.5.*nd2 .+ 0.5.*ndn
        end
        Φu = Fu.vectors[:,1:Nup]; Φd = Fd.vectors[:,1:Ndn]
        C0 = zeros(ComplexF64, 2L, 2L)
        for a in 1:L, b in 1:L
            C0[qup(a),qup(b)] = (Φu*Φu')[a,b]; C0[qdn(a),qdn(b)] = (Φd*Φd')[a,b]
        end
        rot(θ,φ) = begin
            c,s = cos(θ/2), sin(θ/2); Um = zeros(ComplexF64,2L,2L)
            for i in 1:L
                a,b = qup(i), qdn(i)
                Um[a,a]=c; Um[a,b]=-exp(-im*φ)*s; Um[b,a]=exp(im*φ)*s; Um[b,b]=c
            end
            Um*C0*Um'
        end
        Muhf = majorana_M(C0)
        θs = [acos(-1 + (k-0.5)*2/24) for k in 1:24]; φs = [2π*(l-0.5)/8 for l in 1:8]
        Msym = [majorana_M(rot(θ,φ)) for θ in θs, φ in φs]

        sites = 2:L-2
        for i in sites
            obs = site_observables(i; bond=true)
            q1,q2 = obs_support(obs); Nw = q2-q1+1
            ow = shift_obs(obs, q1-1)
            Mats = [[term_window_matrix(tm,Nw) for tm in o.terms] for o in ow]
            ρw = window_rdms(ψ, [(q1,q2)])[1]
            Otrue = [sum(tm.coeff*expect_rdm(ρw, Mats[k][ti]) for (ti,tm) in enumerate(ow[k].terms))
                     for k in 1:length(obs)]

            cands = Tuple{String,Any}[]
            for cp in chi_priors
                σ = truncate(ψ; maxdim=cp); normalize!(σ)
                push!(cands, ("chi$cp", window_rdms(σ, [(q1,q2)])[1]))
            end
            # 平均場 prior の窓RDMを Pauli 期待値から再構成
            shift = q1 - 1
            push!(cands, ("UHF", rdm_from_paulis(
                sup -> gauss_pauli_expect([(q+shift,a) for (q,a) in sup], Muhf), Nw)))
            push!(cands, ("UHFsym", rdm_from_paulis(
                sup -> mean(gauss_pauli_expect([(q+shift,a) for (q,a) in sup], M) for M in Msym), Nw)))

            for (lb, σw) in cands
                Fw, Dtr = fid_trdist(Matrix(ρw), Matrix(σw))
                for (k,o) in enumerate(obs)
                    Pσ = [expect_rdm(σw, Mats[k][ti]) for ti in 1:length(o.terms)]
                    tr_o = sum(tm.coeff*Pσ[ti] for (ti,tm) in enumerate(o.terms))
                    Δ = Otrue[k] - tr_o
                    absA = o.pure ? length(o.terms[1].sup) : 0
                    gmax = o.pure && abs(Otrue[k])<1 ?
                        1 + nm*(3.0^absA-1)*Otrue[k]^2/(3.0^absA*(1-Otrue[k]^2)) : NaN
                    gth = o.pure && abs(Otrue[k])<1 ?
                        ((3.0^absA-1)*Otrue[k]^2 + 3.0^absA*(1-Otrue[k]^2)/nm) /
                        ((3.0^absA-1)*Δ^2      + 3.0^absA*(1-Otrue[k]^2)/nm) : NaN
                    push!(rows, (L, δ, i, o.name, o.pure, absA, lb, Otrue[k], Δ, abs(Δ),
                                 Fw, 1-Fw, Dtr, gth, gmax))
                end
            end
        end
        @printf("  サイト %d〜%d を処理\n", first(sites), last(sites)); flush(stdout)
    end

    out = joinpath(@__DIR__, "crm_window_fidelity_results.tsv")
    open(out,"w") do io
        println(io, "L\tdoping\tsite\tobservable\tpure\tabsA\tprior\ttrue\tDelta\tabsDelta\t" *
                    "F_window\tinfid\tD_trace\tG_theory\tG_max")
        for r in rows; println(io, join(r,"\t")); end
    end
    println("\nresults saved: $out")
end
main()
