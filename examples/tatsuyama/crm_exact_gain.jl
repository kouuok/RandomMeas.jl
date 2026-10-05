# ============================================================
# パウリ和の観測量に対する、標準シャドウと CRM シャドウの厳密な分散比
#
#   O = Σ_k c_k P_k,  x_k = <P_k>_ρ,  y_k = <P_k>_σ (prior),  Δ_k = x_k - y_k
#
# 1基底・n_m ショットの推定量について (n_u は比で約分される)
#   Var_std = Σ_{k,l 両立} c_k c_l 3^{|A_k ∩ A_l|} [<P_k P_l>/n_m + (1-1/n_m) x_k x_l] - (Σ c_k x_k)^2
#   Var_CRM = Σ_{k,l 両立} c_k c_l 3^{|A_k ∩ A_l|} [(<P_k P_l> - x_k x_l)/n_m + Δ_k Δ_l] - (Σ c_k Δ_k)^2
# 「両立」は共通の量子ビットで同じパウリ文字であること(同じ基底で同時に測れる)。
# 単一パウリでは元論文の式(25) [PRX: (C7)] の比 (theory_var) に一致する。
# crm_pauli_sum_variance.py と同じ式。モンテカルロの n_repeat が小さい値を置き換えるために使う。
#
# 使い方: Term / Obs / SIGMA / product_op_expect を定義したスクリプトから include する。
# ============================================================

function pauli_product_sup(sk::Vector{Tuple{Int,Int}}, sl::Vector{Tuple{Int,Int}})
    dk = Dict(sk); dl = Dict(sl)
    compatible = all(get(dl, q, a) == a for (q, a) in dk)
    compatible || return (false, 0, Tuple{Int,Int}[])
    ov = count(q -> haskey(dl, q), keys(dk))
    # 同じ量子ビットで同じ文字の積は恒等になるので、重ならない部分だけが残る
    sup = vcat([(q, a) for (q, a) in sk if !haskey(dl, q)],
               [(q, a) for (q, a) in sl if !haskey(dk, q)])
    return (true, ov, sup)
end

# tens: 真の状態のテンソル列 (product_op_expect に渡す形), Pσ: prior の項ごとの値
function exact_gain(o::Obs, tens, Pσ::Vector{Float64}, nm)
    T = o.terms; K = length(T)
    x = [product_op_expect(tens, term_matrixdict(t)) for t in T]
    Δ = x .- Pσ
    vs = 0.0; vc = 0.0
    for k in 1:K, l in 1:K
        ok, ov, sup = pauli_product_sup(T[k].sup, T[l].sup)
        ok || continue
        pkl = product_op_expect(tens, Dict{Int,Matrix{ComplexF64}}(q => SIGMA[a] for (q, a) in sup))
        w = T[k].coeff * T[l].coeff * 3.0^ov
        vs += w * (pkl / nm + (1 - 1 / nm) * x[k] * x[l])
        vc += w * ((pkl - x[k] * x[l]) / nm + Δ[k] * Δ[l])
    end
    sx = sum(T[k].coeff * x[k] for k in 1:K)
    sd = sum(T[k].coeff * Δ[k] for k in 1:K)
    return (vs - sx^2) / (vc - sd^2)
end
