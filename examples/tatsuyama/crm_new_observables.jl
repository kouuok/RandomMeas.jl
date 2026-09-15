# 物理的に自然な「和の観測量」とその成分を、保存済みの厳密な ρ・UHF・切断 MPS で計算する。
#   中央の結合 (i, i+1) の2サイト縮約密度行列から:
#     hopU = c†_{i↑}c_{j↑}+h.c.   hopD = c†_{i↓}c_{j↓}+h.c.
#     hopU·(1-2n_{i↓}), hopU·(1-2n_{j↓}), hopD·(1-2n_{j↑}), hopD·(1-2n_{i↑})(和の分散の交差項に要る)
#     Z_{i↑}, Z_{i↓}, Z_{j↑}, Z_{j↓} のあらゆる積(対角)
#   距離 r: 電荷の偶奇の2点相関 <F_c F_{c+r}>、ホッピング <c†_{c↑} c_{c+r,↑}>
# 2サイト演算子は同じ局所基底の2サイト系で OpSum から作るので、フェルミオンの符号は ITensor に任せる
# (隣接サイトの c†c では、それより左のサイトの JW の弦は2つの演算子で打ち消し合う)。
#   W=1 LX=32 PBC=0 CRM_U=8.0 julia --project=. examples/tatsuyama/crm_new_observables.jl
include("crm_2d_ed_fid_table.jl")
using Serialization

const S2 = siteinds("Electron", 2)
function op2(terms)   # terms: [(係数, 演算子名..., サイト)...] を OpSum に → 16×16 行列(|a b⟩ の順)
    os = OpSum()
    for t in terms; os += t; end
    M = MPO(os, S2); T = M[1] * M[2]
    A = Array(T, prime(S2[1]), prime(S2[2]), S2[1], S2[2])
    out = zeros(ComplexF64, 16, 16)
    for ia in 1:4, ib in 1:4, ja in 1:4, jb in 1:4
        out[(ia-1)*4+ib, (ja-1)*4+jb] = A[ia, ib, ja, jb]
    end
    out
end
hop(σ) = [(1.0, "Cdag$σ", 1, "C$σ", 2), (1.0, "Cdag$σ", 2, "C$σ", 1)]
withn(ts, nname, site) = vcat(ts, [(-2.0, t[2], t[3], t[4], t[5], nname, site) for t in ts])
const OPS = Dict(
    "hopU" => op2(hop("up")), "hopD" => op2(hop("dn")),
    "hopU_Zb" => op2(withn(hop("up"), "Ndn", 1)), "hopU_Zd" => op2(withn(hop("up"), "Ndn", 2)),
    "hopD_Zc" => op2(withn(hop("dn"), "Nup", 2)), "hopD_Za" => op2(withn(hop("dn"), "Nup", 1)))

function bond_values(ψ::MPS, i::Int)
    ρ = site_rdm(ψ, [i, i+1]); ρ ./= real(tr(ρ))
    vals = Dict{String,Float64}()
    for (k, O) in OPS; vals[k] = real(tr(ρ * O)); end
    # Z の積は占有数基底で対角: 状態 (ia, ib) の (n_{i↑}, n_{i↓}, n_{j↑}, n_{j↓})
    for mask in 1:15
        v = 0.0
        for ia in 1:4, ib in 1:4
            (nu1, nd1) = locbits(ia); (nu2, nd2) = locbits(ib)
            bits = (nu1, nd1, nu2, nd2)
            sgn = prod((mask >> (q-1)) & 1 == 1 ? (1 - 2bits[q]) : 1 for q in 1:4)
            v += sgn * real(ρ[(ia-1)*4+ib, (ia-1)*4+ib])
        end
        vals["Z" * join(((mask >> (q-1)) & 1 == 1 ? "abcd"[q] : "") for q in 1:4)] = v
    end
    vals
end

function run_new()
    geo = PBC ? "torus" : "cylinder"; tag = @sprintf("W%dL%d_%s_U%.1f", W, LX, geo, U)
    d = deserialize(joinpath(@__DIR__, "states", "eigcorr_$tag.jls")); ψ, ψu = d.ψ, d.ψu
    @printf("%s  分散=%.3e  χ=%d\n", tag, d.var, maxlinkdim(ψ))
    states = Any[ψ, ψu]; labels = ["true", "UHF"]
    for c in (2,4,8,16,32,64,128,256,512)
        c < maxlinkdim(ψ) || continue
        σ = copy(ψ); orthogonalize!(σ,1); truncate!(σ; maxdim=c, cutoff=0.0); normalize!(σ)
        push!(states, σ); push!(labels, "chi$c")
    end
    N = length(ψ); c0 = PBC ? 1 : N÷2
    bonds = PBC ? [1] : [c0, c0+1]                     # 開放端は強い結合と弱い結合の両方
    R = PBC ? N÷2 : N÷2 - 1
    fn = joinpath(@__DIR__, "crm_newobs_$tag.tsv")
    open(fn, "w") do io
        println(io, "W\tLX\tgeometry\tU\tHvar\tstate\tkind\ti\tj\tname\tvalue")
        for (lab, φ) in zip(labels, states)
            for i in bonds
                for (k, v) in sort(collect(bond_values(φ, i)))
                    @printf(io, "%d\t%d\t%s\t%.1f\t%.3e\t%s\tbond\t%d\t%d\t%s\t%.12f\n", W, LX, geo, U, d.var, lab, i, i+1, k, v)
                end
            end
            rng_ = c0:min(N, c0+R)
            CF = correlation_matrix(φ, "F", "F"; sites=rng_)
            CH = correlation_matrix(φ, "Cdagup", "Cup"; sites=rng_)
            for (n_, j) in enumerate(rng_)
                n_ == 1 && continue
                @printf(io, "%d\t%d\t%s\t%.1f\t%.3e\t%s\tdist\t%d\t%d\tFF\t%.12f\n", W, LX, geo, U, d.var, lab, c0, j, real(CF[1, n_]))
                @printf(io, "%d\t%d\t%s\t%.1f\t%.3e\t%s\tdist\t%d\t%d\thopU\t%.12f\n", W, LX, geo, U, d.var, lab, c0, j, 2real(CH[1, n_]))
            end
            @printf("  %-6s 完了\n", lab)
        end
    end
    println("書き出し: ", fn)
end

if get(ENV, "NEWOBS_SELFTEST", "0") == "1"
    d = deserialize(ENV["NEWOBS_STATE"]); ψ = d.ψ
    v = bond_values(ψ, 4)
    CH = correlation_matrix(ψ, "Cdagup", "Cup"); CD = correlation_matrix(ψ, "Cdagdn", "Cdn")
    @printf("  hopU: RDM %.10f  相関行列 %.10f\n", v["hopU"], 2real(CH[4,5]))
    @printf("  hopD: RDM %.10f  相関行列 %.10f\n", v["hopD"], 2real(CD[4,5]))
    @printf("  Zab(オンサイト): %.10f   Zac(Z↑Z↑ 隣): %.10f   Zabcd(偶奇 P_2): %.10f\n", v["Zab"], v["Zac"], v["Zabcd"])
    # 密度付きホッピングを MPO で直接計算して照合
    os = OpSum(); os += 1.0,"Cdagup",4,"Cup",5; os += 1.0,"Cdagup",5,"Cup",4
    os += -2.0,"Cdagup",4,"Ndn",4,"Cup",5; os += -2.0,"Cdagup",5,"Cup",4,"Ndn",4
    s = siteinds(ψ); @printf("  hopU·(1-2n_{i↓}): RDM %.10f  MPO %.10f\n", v["hopU_Zb"], real(inner(ψ', MPO(os, s), ψ)))
    os = OpSum(); os += 1.0,"Cdagdn",4,"Cdn",5; os += 1.0,"Cdagdn",5,"Cdn",4
    os += -2.0,"Nup",4,"Cdagdn",4,"Cdn",5; os += -2.0,"Cdagdn",5,"Nup",4,"Cdn",4
    @printf("  hopD·(1-2n_{i↑}): RDM %.10f  MPO %.10f\n", v["hopD_Za"], real(inner(ψ', MPO(os, s), ψ)))
else
    run_new()
end
