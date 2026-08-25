# ============================================================
# Gauss 状態(Slater行列式・その回転平均)に対する任意 Pauli 列の期待値
#
# 動機:
#   crm_1d_hf.jl / crm_2d_site.jl の wick_gen は「全Z列(長さ1–2)」と
#   ホッピングの1パターンしか扱えず、それ以外は error を投げる。
#   そのため横スピン相関 <S^x_i S^x_j> のような XXXX / XXYY 型の項を
#   UHF prior で評価できず、**測れる観測量が Z 型に縛られていた**。
#   CRM は prior の「項ごとの期待値 <P_t>_σ」しか要求しないので、
#   ここを一般化すれば prior 側の制約が外れる。
#
# 方法:
#   Jordan-Wigner で Majorana を
#     γ_{2q-1} = (∏_{k<q} Z_k) X_q,   γ_{2q} = (∏_{k<q} Z_k) Y_q
#   と定めると
#     X_q = (∏_{k<q} Z_k) γ_{2q-1},  Y_q = (∏_{k<q} Z_k) γ_{2q},
#     Z_q = -i γ_{2q-1} γ_{2q},      ∏_{k<q} Z_k = ∏_{k<q} (-i γ_{2k-1} γ_{2k})
#   なので任意の Pauli 列は Majorana の単項式(位相つき)になる。
#   粒子数保存の Gauss 状態では ⟨γ_j γ_k⟩ = δ_jk + i M_jk で
#     M[2p-1, 2q] = δ_pq - 2C[p,q],  奇-奇 と 偶-偶 ブロックは 0
#   (C[p,q] = ⟨c†_p c_q⟩、実対称)であり、
#     ⟨γ_{j1}…γ_{j2m}⟩ = i^m Pf(M[j,j])
#   で閉じる。
#
# 検証: crm_wick_pauli_check() が
#   (1) 既存 wick_gen の対応範囲と一致するか
#   (2) 密行列で作った Slater 行列式の期待値と一致するか
#   を確かめる。
# ============================================================
using LinearAlgebra
using Random

"""Majorana 単項式の積。添字は昇順・重複なしに正規化し、符号を追跡する。"""
function mmul(c1::ComplexF64, v1::Vector{Int}, c2::ComplexF64, v2::Vector{Int})
    v = vcat(v1, v2); c = c1 * c2
    for i in 2:length(v)                      # 挿入ソート(交換のたびに符号反転)
        j = i
        while j > 1 && v[j-1] > v[j]
            v[j-1], v[j] = v[j], v[j-1]; c = -c; j -= 1
        end
    end
    out = Int[]; i = 1                        # γ^2 = 1 で隣接する同一添字を消す
    while i <= length(v)
        if i < length(v) && v[i] == v[i+1]
            i += 2
        else
            push!(out, v[i]); i += 1
        end
    end
    return c, out
end

"""Pauli 列 sup = [(量子ビット, 1=X/2=Y/3=Z), …] を Majorana 単項式に変換。"""
function pauli_majorana(sup)
    c = ComplexF64(1); v = Int[]
    for (q, a) in sort(collect(sup), by = first)
        if a == 3
            c2, v2 = ComplexF64(-im), [2q-1, 2q]
        else
            c2 = ComplexF64(1); v2 = Int[]
            for k in 1:q-1                     # JW 弦 ∏_{k<q} Z_k
                c2, v2 = mmul(c2, v2, ComplexF64(-im), [2k-1, 2k])
            end
            c2, v2 = mmul(c2, v2, ComplexF64(1), [a == 1 ? 2q-1 : 2q])
        end
        c, v = mmul(c, v, c2, v2)
    end
    return c, v
end

"""実反対称行列の Pfaffian(小さい行列専用の再帰展開)。"""
function pfaffian(A::AbstractMatrix)
    n = size(A, 1)
    n == 0 && return 1.0
    isodd(n) && return 0.0
    n == 2 && return A[1, 2]
    s = 0.0
    for j in 2:n
        a = A[1, j]
        a == 0 && continue
        idx = [k for k in 2:n if k != j]
        s += (isodd(j) ? -1 : 1) * a * pfaffian(@view A[idx, idx])
    end
    return s
end

"""相関行列 C[p,q] = ⟨c†_p c_q⟩(エルミート、複素可)から Majorana 共分散 M を作る。

⟨γ_j γ_k⟩ = δ_jk + i M_jk と置くと
    M[2p-1, 2q-1] = M[2p, 2q]   = 2 Im C[p,q]
    M[2p-1, 2q]                  = δ_pq - 2 Re C[p,q]
となる(M は実反対称)。C が実対称なら奇-奇・偶-偶ブロックが消えて
従来の式に戻る。スピン回転を方位角 φ まで含めると C は複素になるので
この一般形が要る。"""
function majorana_M(C::AbstractMatrix)
    n = size(C, 1); M = zeros(2n, 2n)
    for p in 1:n, q in 1:n
        re = (p == q ? 1.0 : 0.0) - 2real(C[p, q])
        im_ = 2imag(C[p, q])
        M[2p-1, 2q]   = re;  M[2q,   2p-1] = -re
        M[2p-1, 2q-1] = im_
        M[2p,   2q]   = im_
    end
    for j in 1:2n; M[j, j] = 0.0; end
    return M
end

"""Gauss 状態での任意 Pauli 列の期待値。"""
function gauss_pauli_expect(sup, M::AbstractMatrix)
    isempty(sup) && return 1.0
    c, v = pauli_majorana(sup)
    isempty(v) && return real(c)
    isodd(length(v)) && return 0.0
    m = length(v) ÷ 2
    return real(c * (im^m) * pfaffian(M[v, v]))
end

# ------------------------------------------------------------
# 検証
# ------------------------------------------------------------
"""密行列で Slater 行列式を作り、Pauli 期待値を直接計算して突き合わせる。"""
function crm_wick_pauli_check(; n = 4, seed = 7)
    Random.seed!(seed)
    # ランダムな実直交軌道から粒子数 Nocc の Slater 行列式を作る
    Nocc = 2
    A = randn(n, n); Q = Matrix(qr(A).Q); Φ = Q[:, 1:Nocc]
    C = Φ * Φ'                                    # C[p,q] = ⟨c†_p c_q⟩
    M = majorana_M(C)

    # 密行列: JW で c_p を作り、Slater 行列式ベクトルを構成
    sx = ComplexF64[0 1; 1 0]; sy = ComplexF64[0 -im; im 0]; sz = ComplexF64[1 0; 0 -1]
    id = Matrix{ComplexF64}(I, 2, 2)
    cdag_op(p) = reduce(kron, [k < p ? sz : (k == p ? ComplexF64[0 0; 1 0] : id) for k in 1:n])
    Cd = [cdag_op(p) for p in 1:n]
    v = zeros(ComplexF64, 2^n); v[1] = 1
    for k in 1:Nocc
        d = sum(Φ[p, k] * Cd[p] for p in 1:n)
        v = d * v
    end
    v ./= norm(v)

    pauli_dense(sup) = reduce(kron,
        [(j = findfirst(x -> x[1] == q, sup);
          j === nothing ? id : (sup[j][2] == 1 ? sx : sup[j][2] == 2 ? sy : sz)) for q in 1:n])

    maxerr = 0.0; ntest = 0
    for a in 1:3, b in 1:3, p in 1:n-1, q in p+1:n
        for sup in ([(p, a)], [(p, a), (q, b)])
            ref = real(dot(v, pauli_dense(sup) * v))
            got = gauss_pauli_expect(sup, M)
            maxerr = max(maxerr, abs(ref - got)); ntest += 1
        end
    end
    # 4体項(横スピン相関に現れる型)も確認
    for (a, b, c2, d) in ((1,1,1,1), (1,1,2,2), (2,2,1,1), (2,2,2,2), (3,3,3,3), (1,2,1,2))
        sup = [(1, a), (2, b), (3, c2), (4, d)]
        ref = real(dot(v, pauli_dense(sup) * v))
        got = gauss_pauli_expect(sup, M)
        maxerr = max(maxerr, abs(ref - got)); ntest += 1
    end
    return maxerr, ntest, C
end
