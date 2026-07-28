# U=4 と U=8 の中央 Schmidt スペクトルを比べる。
# 仮説: 強結合ではスピンセクターが SU(2) 対称に近く、Schmidt 値が多重項で
# 縮退するため、χ_p=8 のような「多重項の途中で切る」切断が相対的に大きな
# 誤差を生む。エントロピーは小さいのに prior 誤差が大きい理由の候補。
include("crm_chain_common.jl")

function central_spectrum(ψ::MPS; ntop=24)
    N=length(ψ); ψo=copy(ψ); orthogonalize!(ψo, N÷2)
    _,Sv,_ = svd(ψo[N÷2], (linkinds(ψo)[N÷2-1], siteinds(ψo)[N÷2]))
    s = diag(Array(Sv, inds(Sv)...))
    p = sort(s.^2; rev=true)
    return p[1:min(ntop,length(p))], p
end

for U in (4.0, 8.0)
    E, ψ, _, _, _ = ground_state(64, 1.0, U, U/2; chi_max=256, nsweeps=30)
    top, p = central_spectrum(ψ)
    S = -sum(p[p.>1e-14] .* log.(p[p.>1e-14]))
    @printf("\n=== U=%.1f, L=64 : S_center=%.4f ===\n", U, S)
    @printf("上位24個の Schmidt 重み (p_k) と、その累積・χ切断で捨てる重み\n")
    @printf("%4s %14s %12s %14s\n","k","p_k","p_k/p_{k-1}","1 - Σ_{j<=k} p_j")
    c = 0.0
    for k in 1:length(top)
        c += top[k]
        @printf("%4d %14.6e %12.4f %14.3e\n", k, top[k], k>1 ? top[k]/top[k-1] : NaN, 1-c)
    end
    for x in (4,8,16,32,64)
        @printf("  χ_p=%3d で捨てる重み: %.4e\n", x, 1-sum(p[1:min(x,length(p))]))
    end
end
