# U=8 半充填で ε が L とともに増える件が、DMRG の残留収束誤差かどうかの診断。
# ρ 自体が未収束だと、その χ_p 切断との差 Δ は「切断誤差」ではなく
# 「収束誤差の差」を測ってしまう。スイープ数と χ_exp を変えて ε が動くかを見る。
include(joinpath(@__DIR__, "crm_chain_common.jl"))

function eps_median(ψ, L, chi_p)
    obs_by_site = [site_observables(i; bond = i < L) for i in 1:L]
    windows = [obs_support(o) for o in obs_by_site]
    rρ = window_rdms(ψ, windows)
    σ = truncate(ψ; maxdim=chi_p); normalize!(σ)
    rσ = window_rdms(σ, windows)
    acc = Dict{String,Vector{Float64}}()
    for i in 1:L
        obs = obs_by_site[i]; q1,q2 = windows[i]; Nw = q2-q1+1
        for (o, ow) in zip(obs, shift_obs(obs, q1-1))
            Ms = [term_window_matrix(tm, Nw) for tm in ow.terms]
            vρ = sum(tm.coeff*expect_rdm(rρ[i], Ms[t]) for (t,tm) in enumerate(ow.terms))
            vσ = sum(tm.coeff*expect_rdm(rσ[i], Ms[t]) for (t,tm) in enumerate(ow.terms))
            abs(vρ) > 1e-3 && push!(get!(acc, o.name, Float64[]), abs((vρ-vσ)/vρ))
        end
    end
    return acc
end

@printf("%5s %5s %8s %6s %14s %8s %9s %11s %11s\n",
        "U","L","nsweeps","chi","E0","maxlink","S_center","ε(ZZ onsite)","ε(DoubleOcc)")
for U in (4.0, 8.0), L in (32, 128), (ns, cx) in ((16,256),(30,256),(45,256),(45,512))
    (L == 32 && cx == 512) && continue
    t0=time()
    E, ψ, _, _, _ = ground_state(L, 1.0, U, U/2; chi_max=cx, nsweeps=ns)
    N=length(ψ); ψo=copy(ψ); orthogonalize!(ψo, N÷2)
    _,Sv,_ = svd(ψo[N÷2], (linkinds(ψo)[N÷2-1], siteinds(ψo)[N÷2]))
    p=diag(Array(Sv,inds(Sv)...)).^2; p=p[p.>1e-14]
    a = eps_median(ψ, L, 8)
    @printf("%5.1f %5d %8d %6d %14.6f %8d %9.4f %11.3e %11.3e   (%.0fs)\n",
            U, L, ns, cx, E, maxlinkdim(ψ), -sum(p.*log.(p)),
            median(a["ZZ onsite"]), median(a["DoubleOcc"]), time()-t0)
    flush(stdout)
end
