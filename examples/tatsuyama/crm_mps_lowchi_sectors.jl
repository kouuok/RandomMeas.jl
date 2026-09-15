# 低い結合次元の MPS が、中央の結合の少ない状態を何に使っているかを見る。
# 中央の結合で Schmidt 分解し、左半分の電子数 N_L とスピン 2S^z_L ごとに、残した状態の数と重み(特異値の2乗の和)を出す。
#   厳密な状態(states/eigcorr_*.jls)、それを χ に切り詰めた切断 MPS、変分 MPS(states/mpslow_*.jls、crm_mps_lowchi.jl)
#   TAGS="W1L64_cylinder_U12.0" julia --project=. crm_mps_lowchi_sectors.jl
using ITensors, ITensorMPS, LinearAlgebra, Printf, Serialization

"""結合 (c, c+1) の Schmidt 係数を量子数のブロックごとにまとめる。"""
function sectors(ψ0::MPS, c::Int)
    ψ = orthogonalize(ψ0, c)
    li = c > 1 ? (linkind(ψ, c-1), siteind(ψ, c)) : (siteind(ψ, c),)
    U, S, V = svd(ψ[c], li...)
    u = commonind(U, S)
    s = diag(Array(dense(S), inds(S)...))
    out = Dict{Tuple{Int,Int},Tuple{Int,Float64}}(); k = 0
    for (q, d) in space(u)
        w = sum(abs2, s[k+1:k+d]); k += d
        key = (val(q, "Nf"), val(q, "Sz"))
        n0, w0 = get(out, key, (0, 0.0))
        w > 1e-14 && (out[key] = (n0 + count(>(1e-7), s[k-d+1:k]), w0 + w))
    end
    out
end

function show_sectors(label, sec, nL)
    tot = sum(last, values(sec))
    ch = sum((w for ((N, Sz), (_, w)) in sec if abs(N) != nL); init=0.0)
    @printf("  %-12s 電荷がゆらいだセクター(N_L≠%d)の重み %.4f\n", label, nL, ch/tot)
    for (key, (nst, w)) in sort(collect(sec), by = t -> -t[2][2])
        w/tot > 1e-6 || continue
        @printf("      N_L=%3d 2S^z_L=%+d : 状態 %3d  重み %.6f\n", key[1], key[2], nst, w/tot)
    end
end

for tag in split(get(ENV, "TAGS", "W1L64_cylinder_U12.0"), ",")
    d = deserialize(joinpath(@__DIR__, "states", "eigcorr_$tag.jls")); ψ = d.ψ
    N = length(ψ); c = N ÷ 2
    println("== $tag  中央の結合 ($c, $(c+1))  厳密な状態の χ=$(maxlinkdim(ψ))")
    show_sectors("厳密", sectors(ψ, c), c)
    for χ in (4, 8)
        σ = copy(ψ); orthogonalize!(σ, 1); truncate!(σ; maxdim=χ, cutoff=0.0); normalize!(σ)
        show_sectors("切断 χ=$χ", sectors(σ, c), c)
    end
    for χ in (4, 8, 16)
        f = joinpath(@__DIR__, "states", "mpslow_$(tag)_chi$χ.jls")
        isfile(f) && show_sectors("変分 χ=$χ", sectors(deserialize(f).ψ, c), c)
    end
end
