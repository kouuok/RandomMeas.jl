# ============================================================
# 1つのシャドウを全サイト・全観測量で共有する(本来のプロトコル)
#
# これまでの実験(run_locals)は **観測量ごと・サイトごとに独立な乱数種**で
# 窓の4量子ビットだけをサンプルしていた。個々の観測量の分散はそれで正しく
# 再現される(基底は量子ビットごとに独立なので、窓に制限した周辺分布は
# 全系から引いた場合と同じ)。しかし古典シャドウの本来のプロトコルは
#   「系全体で n_u 個の基底を引き、n_m ショット測り、**同じデータ**を
#    後処理して任意の観測量を推定する」
# であり、**観測量どうしが同じ設定・同じショットを共有する**。
# したがって観測量の和(全エネルギー、構造因子)の分散や、
# 「全観測量を同時に賄うのに n_u がいくつ要るか」は再現できていなかった。
#
# ここでは L=8(16量子ビット)を密ベクトルで扱い、全系のビット列を厳密に
# サンプルして本来のプロトコルをそのまま模擬する。測るのは
#   (a) サイトごとの局所観測量 — 独立サンプルの結果と一致するか
#   (b) 全エネルギー — 多数の項の和。相関が効く
#   (c) スタッガード構造因子 S(π) — L^2 個の項の和
#
# 実行: julia --project=Hubbard_MPS_Env_v2 crm_global_shadow.jl
# 環境変数: L(既定 8), U(既定 4.0), NU_LIST(既定 "50,100,200,400,800"),
#           NM(既定 100), N_REPEAT(既定 200)
# ============================================================
include(joinpath(@__DIR__, "crm_chain_common.jl"))

"""MPS を密ベクトルに展開する(量子ビット n が bit n-1)。"""
function mps_dense(ψ::MPS)
    tens = extract_tensors(ψ)
    Ψ = ComplexF64[1.0;;]
    for (A0, A1) in tens
        Ψ = vcat(Ψ*ComplexF64.(A0), Ψ*ComplexF64.(A1))
    end
    vec(Ψ)
end

"""密ベクトルに1量子ビットずつ基底回転をかける(全系 N 量子ビット)。"""
function rotate_all(v::Vector{ComplexF64}, basis::Vector{Int}, N::Int)
    w = copy(v)
    for q in 1:N
        Uq = UBASIS[basis[q]]
        stride = 1 << (q-1)
        @inbounds for blk in 0:(1 << (N-q))-1, off in 0:stride-1
            i0 = (blk << q) | off; i1 = i0 | stride
            a = w[i0+1]; b = w[i1+1]
            w[i0+1] = Uq[1,1]*a + Uq[1,2]*b
            w[i1+1] = Uq[2,1]*a + Uq[2,2]*b
        end
    end
    w
end

"""Pauli 項の集合。coeff と台 (量子ビット, 1=X/2=Y/3=Z) を持つ。"""
struct GTerm
    coeff::Float64
    qs::Vector{Int}
    as::Vector{Int}
end

"""観測量 = 名前 + 項のリスト。"""
struct GObs
    name::String
    terms::Vector{GTerm}
end

function build_observables(L::Int, t, U, mu)
    obs = GObs[]
    # --- サイトごとの局所観測量(鎖中央付近の代表サイト) ---
    for i in (L÷2,)
        u, d = qup(i), qdn(i); u2, d2 = qup(i+1), qdn(i+1)
        push!(obs, GObs("ZZ onsite", [GTerm(1.0,[u,d],[3,3])]))
        push!(obs, GObs("DoubleOcc", [GTerm(0.25,Int[],Int[]), GTerm(-0.25,[u],[3]),
                                      GTerm(-0.25,[d],[3]), GTerm(0.25,[u,d],[3,3])]))
        push!(obs, GObs("SzSz r=1", [GTerm(1/16,[u,u2],[3,3]), GTerm(-1/16,[u,d2],[3,3]),
                                     GTerm(-1/16,[d,u2],[3,3]), GTerm(1/16,[d,d2],[3,3])]))
        push!(obs, GObs("hop r=1", [GTerm(0.5,[u,d,u2],[1,3,1]), GTerm(0.5,[u,d,u2],[2,3,2])]))
    end
    # --- 全エネルギー(全ボンド + 全サイト) ---
    et = GTerm[]
    for i in 1:L-1
        for (a,b) in ((qup(i),qup(i+1)), (qdn(i),qdn(i+1)))
            mid = a+1                                  # 間の1モード
            push!(et, GTerm(-t, [a,mid,b], [1,3,1]))
            push!(et, GTerm(-t, [a,mid,b], [2,3,2]))
        end
    end
    for i in 1:L
        u,d = qup(i), qdn(i)
        push!(et, GTerm(U*0.25, Int[], Int[])); push!(et, GTerm(-U*0.25,[u],[3]))
        push!(et, GTerm(-U*0.25,[d],[3]));      push!(et, GTerm(U*0.25,[u,d],[3,3]))
        push!(et, GTerm(-mu*1.0, Int[], Int[]))
        push!(et, GTerm(mu*0.5,[u],[3]));       push!(et, GTerm(mu*0.5,[d],[3]))
    end
    push!(obs, GObs("全エネルギー", et))
    # --- スタッガード構造因子 S(π) = (1/L) Σ_ij (-1)^{i-j} <Sz_i Sz_j> ---
    st = GTerm[]
    for i in 1:L, j in 1:L
        sgn = (-1)^(i-j) / L
        u,d = qup(i), qdn(i); u2,d2 = qup(j), qdn(j)
        if i == j
            push!(st, GTerm(sgn*0.25, Int[], Int[]))          # (Sz_i)^2 = 1/4 - (1/2)n_up n_dn ...
            push!(st, GTerm(-sgn*0.5, [u,d], [3,3]))
            push!(st, GTerm(sgn*0.25, Int[], Int[]))
        else
            push!(st, GTerm( sgn/16,[u,u2],[3,3])); push!(st, GTerm(-sgn/16,[u,d2],[3,3]))
            push!(st, GTerm(-sgn/16,[d,u2],[3,3])); push!(st, GTerm( sgn/16,[d,d2],[3,3]))
        end
    end
    push!(obs, GObs("構造因子 S(π)", st))
    obs
end

"""1設定 u とビット列から、項ごとのスナップショット推定値を足し込む。"""
@inline function term_estimate(tm::GTerm, basis::Vector{Int}, bits::Vector{Int})
    isempty(tm.qs) && return tm.coeff
    v = tm.coeff
    @inbounds for k in 1:length(tm.qs)
        q = tm.qs[k]
        basis[q] == tm.as[k] || return 0.0
        v *= 3.0 * (bits[q] == 0 ? 1.0 : -1.0)
    end
    v
end

@inline function term_prior(tm::GTerm, basis::Vector{Int}, Pσ::Float64)
    isempty(tm.qs) && return tm.coeff
    @inbounds for k in 1:length(tm.qs)
        basis[tm.qs[k]] == tm.as[k] || return 0.0
    end
    tm.coeff * 3.0^length(tm.qs) * Pσ
end

function main()
    t = 1.0
    U  = parse(Float64, get(ENV, "U", "4.0"))
    L  = parse(Int, get(ENV, "L", "8"))
    nm = parse(Int, get(ENV, "NM", "100"))
    n_repeat = parse(Int, get(ENV, "N_REPEAT", "200"))
    nu_list = parse.(Int, split(get(ENV, "NU_LIST", "50,100,200,400,800"), ","))
    chi_priors = [2, 4, 8, 32]
    N = 2L; mu = U/2

    E, ψ, _, _, _ = ground_state(L, t, U, mu; chi_max=128, nsweeps=20)
    @printf("L=%d U=%.1f E0=%.6f chi=%d (%d量子ビット, 密ベクトル %d 次元)\n",
            L, U, E, maxlinkdim(ψ), N, 1<<N); flush(stdout)

    vρ = mps_dense(ψ); vρ ./= norm(vρ)
    priors = MPS[]
    for cp in chi_priors
        σ = truncate(ψ; maxdim=cp); normalize!(σ); push!(priors, σ)
    end
    push!(priors, ψ)
    labels = vcat(["chi$c" for c in chi_priors], ["exact"])
    vσ = [ (w = mps_dense(σ); w ./ norm(w)) for σ in priors ]

    obs = build_observables(L, t, U, mu)
    # 項ごとの厳密期待値(ρ と 各 prior)
    function pauli_expect(v::Vector{ComplexF64}, qs, as)
        isempty(qs) && return 1.0
        acc = 0.0
        @inbounds for idx in 0:(1<<N)-1
            a = abs2(v[idx+1]); a == 0 && continue
            s = 1.0; ok = true
            for k in 1:length(qs)
                as[k] == 3 || (ok = false; break)
                s *= ((idx >> (qs[k]-1)) & 1) == 0 ? 1.0 : -1.0
            end
            ok || return NaN
            acc += a*s
        end
        acc
    end
    # X/Y を含む項は回転してから対角期待値を取る
    function pexp(v, qs, as)
        isempty(qs) && return 1.0
        b = fill(3, N); for k in 1:length(qs); b[qs[k]] = as[k]; end
        w = rotate_all(v, b, N)
        acc = 0.0
        @inbounds for idx in 0:(1<<N)-1
            p = abs2(w[idx+1]); p == 0 && continue
            s = 1.0
            for q in qs; s *= ((idx >> (q-1)) & 1) == 0 ? 1.0 : -1.0; end
            acc += p*s
        end
        acc
    end
    Pρ = [[pexp(vρ, tm.qs, tm.as) for tm in o.terms] for o in obs]
    Pσs = [[[pexp(v, tm.qs, tm.as) for tm in o.terms] for o in obs] for v in vσ]
    Otrue = [sum(o.terms[j].coeff == 0 ? 0.0 :
                 (isempty(o.terms[j].qs) ? o.terms[j].coeff : o.terms[j].coeff*Pρ[k][j])
                 for j in 1:length(o.terms)) for (k,o) in enumerate(obs)]
    trOσ = [[sum(isempty(o.terms[j].qs) ? o.terms[j].coeff : o.terms[j].coeff*Pσs[p][k][j]
                 for j in 1:length(o.terms)) for (k,o) in enumerate(obs)]
            for p in 1:length(priors)]
    @printf("\n%-18s %14s %12s\n", "観測量", "真値", "項数")
    for (k,o) in enumerate(obs); @printf("%-18s %14.6f %12d\n", o.name, Otrue[k], length(o.terms)); end
    flush(stdout)

    rows = []
    for nu in nu_list
        Random.seed!(20260826 + nu)
        nobs = length(obs); npri = length(priors)
        est_s = zeros(n_repeat, nobs); est_c = zeros(n_repeat, nobs, npri)
        basis = zeros(Int, N); bits = zeros(Int, N)
        t0 = time()
        for rep in 1:n_repeat
            acc_s = zeros(nobs); acc_c = zeros(nobs, npri)
            for _ in 1:nu
                rand!(basis, 1:3)
                w = rotate_all(vρ, basis, N)
                cum = cumsum(abs2.(w))
                xs = zeros(nobs)
                for _ in 1:nm
                    idx = searchsortedfirst(cum, rand()*cum[end]) - 1
                    @inbounds for q in 1:N; bits[q] = (idx >> (q-1)) & 1; end
                    for (k,o) in enumerate(obs)
                        v = 0.0
                        for tm in o.terms; v += term_estimate(tm, basis, bits); end
                        xs[k] += v
                    end
                end
                for k in 1:nobs
                    mρ = xs[k]/nm; acc_s[k] += mρ
                    for p in 1:npri
                        y = 0.0
                        for (j,tm) in enumerate(obs[k].terms)
                            y += term_prior(tm, basis, Pσs[p][k][j])
                        end
                        acc_c[k,p] += mρ - y
                    end
                end
            end
            est_s[rep,:] .= acc_s ./ nu
            for p in 1:npri, k in 1:nobs; est_c[rep,k,p] = acc_c[k,p]/nu + trOσ[p][k]; end
        end
        @printf("\n-- n_u=%-5d (n_repeat=%d, %.0fs) 1つのシャドウを全観測量で共有 --\n",
                nu, n_repeat, time()-t0)
        @printf("  %-18s %10s", "観測量", "sd(標準)")
        for lb in labels; @printf("%11s", "G($lb)"); end; println()
        for (k,o) in enumerate(obs)
            v = var(est_s[:,k])
            @printf("  %-18s %10.5f", o.name, sqrt(v))
            for p in 1:npri; @printf("%11.2f", v/var(est_c[:,k,p])); end
            println()
            for (p,lb) in enumerate(labels)
                push!(rows, (L, nu, nm, n_repeat, o.name, length(o.terms), lb, Otrue[k],
                             Otrue[k]-trOσ[p][k], sqrt(v), sqrt(var(est_c[:,k,p])),
                             v/var(est_c[:,k,p])))
            end
        end
        flush(stdout)
    end

    out = joinpath(@__DIR__, "crm_global_shadow_results.tsv")
    open(out,"w") do io
        println(io, "L\tn_u\tn_m\tn_repeat\tobservable\tnterms\tprior\ttrue\tDelta\tsd_std\tsd_crm\tG")
        for r in rows; println(io, join(r,"\t")); end
    end
    println("\nresults saved: $out")
end
main()
