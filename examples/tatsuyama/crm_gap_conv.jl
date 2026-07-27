# U=8 half-filling の S_center が L=64→128 で減少した件の収束チェック。
# c=1 の gapless スピンセクターなら S_c は (1/6)ln2 ≈ 0.116 だけ増えるはず。
# スイープ数と χ を上げて、増加に転じるかを見る。
include(joinpath(@__DIR__, "crm_chain_common.jl"))

function center_entropy(ψ::MPS)
    N = length(ψ); ψo = copy(ψ); orthogonalize!(ψo, 1); vals = Float64[]
    for b in (N ÷ 2 - 1):(N ÷ 2 + 2)
        (2 <= b <= N - 1) || continue
        orthogonalize!(ψo, b)
        _, Sv, _ = svd(ψo[b], (linkinds(ψo)[b-1], siteinds(ψo)[b]))
        p = diag(Array(Sv, inds(Sv)...)).^2; p = p[p .> 1e-14]
        push!(vals, -sum(p .* log.(p)))
    end
    return mean(vals)
end

@printf("%6s %6s %8s %6s %14s %10s %10s\n", "U", "L", "nsweeps", "chi", "E0", "maxlink", "S_center")
for (U, δ) in ((8.0, 0.0), (8.0, 0.125))
    for L in (64, 128), (ns, cx) in ((16, 256), (30, 256), (30, 512))
        t0 = time()
        E, ψ, _, _, _ = ground_state(L, 1.0, U, U/2; chi_max=cx, nsweeps=ns, δ)
        @printf("%6.1f %6d %8d %6d %14.6f %10d %10.4f   (%.0fs)\n",
                U, L, ns, cx, E, maxlinkdim(ψ), center_entropy(ψ), time()-t0)
        flush(stdout)
    end
end
