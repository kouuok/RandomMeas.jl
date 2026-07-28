# 仮説: 強結合では二重占有 <D> 自体が小さく、しかも仮想的な電荷ゆらぎ
# （＝Schmidt 重みの小さい成分）が担っている。したがって「重みの 99.3% を
# 保つ切断」でも、その 0.7% に住む量は大きく壊れる。
# 検証: <D> の絶対誤差と相対誤差を U=4 / U=8 で比べる。
include("crm_chain_common.jl")
L = 64
@printf("%5s %10s %12s %12s %12s %12s %12s\n",
        "U","<D>_ρ","<D>_σ(χ8)","|Δ_D|","相対誤差","捨てる重み","S_center")
for U in (2.0, 4.0, 8.0)
    E, ψ, _, _, _ = ground_state(L, 1.0, U, U/2; chi_max=256, nsweeps=30)
    N=length(ψ); ψo=copy(ψ); orthogonalize!(ψo, N÷2)
    _,Sv,_ = svd(ψo[N÷2], (linkinds(ψo)[N÷2-1], siteinds(ψo)[N÷2]))
    p = sort(diag(Array(Sv,inds(Sv)...)).^2; rev=true)
    S = -sum(p[p.>1e-14].*log.(p[p.>1e-14]))
    disc = 1 - sum(p[1:min(8,length(p))])
    σ = truncate(ψ; maxdim=8); normalize!(σ)
    obs = site_observables(L÷2; bond=true); q1,q2 = obs_support(obs); Nw=q2-q1+1
    ow = shift_obs(obs, q1-1)
    k = findfirst(o->o.name=="DoubleOcc", obs)
    rρ = window_rdms(ψ, [(q1,q2)])[1]; rσ = window_rdms(σ, [(q1,q2)])[1]
    Ms = [term_window_matrix(tm,Nw) for tm in ow[k].terms]
    dρ = sum(tm.coeff*expect_rdm(rρ,Ms[t]) for (t,tm) in enumerate(ow[k].terms))
    dσ = sum(tm.coeff*expect_rdm(rσ,Ms[t]) for (t,tm) in enumerate(ow[k].terms))
    @printf("%5.1f %10.5f %12.5f %12.3e %12.3e %12.3e %12.4f\n",
            U, dρ, dσ, abs(dρ-dσ), abs((dρ-dσ)/dρ), disc, S)
end
