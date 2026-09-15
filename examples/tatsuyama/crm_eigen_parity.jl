# ブロックの偶奇の文字列で「どういう観測量なら HF を使った CRM の利得が大きいか」を確かめる。
#   P_ℓ = Π_{i∈B} Z_{i↑}Z_{i↓} = (-1)^{N_B}   ブロック内の電荷の偶奇(台 2ℓ 量子ビット、スピン回転で不変)
#   Q_ℓ = Π_{i∈B} Z_{i↑}       = (-1)^{N↑_B}  ブロック内の上向き電子の偶奇(台 ℓ 量子ビット、スピン回転で変わる)
# crm_eigen_corr.jl が states/ に保存した厳密な ρ と UHF の MPS を読むだけで、DMRG はしない。
#   W=1 LX=32 PBC=0 CRM_U=8.0 julia --project=. examples/tatsuyama/crm_eigen_parity.jl
include("crm_2d_ed_fid_table.jl")    # 定数(W, LX, PBC, U)
using Serialization

"""直交中心を s0 に置き、s0 から右へ1サイトずつ演算子を掛けて、ℓ=1..ℓmax の Π_{j=s0}^{s0+ℓ-1} O_j を一度に求める
   (ITensorMPS の correlation_matrix と同じ転送行列の組み方)。"""
function string_expect(psi0::MPS, opname::String, s0::Int, lmax::Int)
    psi = orthogonalize(psi0, s0); s = siteinds(psi); N = length(psi)
    nrm = norm(psi[s0])^2
    L = ITensor(1.0)
    if s0 > 1
        li = commonind(psi[s0], psi[s0-1]); L = delta(dag(li), li')
    end
    vals = Float64[]
    for j in s0:min(N, s0+lmax-1)
        T = L * psi[j] * op(opname, s[j])
        closer = j > 1 ? prime(dag(psi[j]), (s[j], commonind(psi[j], psi[j-1]))) : prime(dag(psi[j]), s[j])
        push!(vals, real(scalar(T * closer)) / nrm)
        L = T * dag(psi[j])'
    end
    vals
end

"""検証用: ブロックに演算子を掛けた MPS との内積で愚直に計算する。"""
function string_bruteforce(psi::MPS, opname::String, s0::Int, l::Int)
    φ = copy(psi); s = siteinds(psi)
    for j in s0:s0+l-1
        φ[j] = noprime(φ[j] * op(opname, s[j]))
    end
    real(inner(psi, φ)) / real(inner(psi, psi))
end

function run_parity()
    geo = PBC ? "torus" : "cylinder"
    tag = @sprintf("W%dL%d_%s_U%.1f", W, LX, geo, U)
    d = deserialize(joinpath(@__DIR__, "states", "eigcorr_$tag.jls"))
    ψ, ψu = d.ψ, d.ψu
    @printf("%s  E0=%.12f  分散=%.3e  χ=%d\n", tag, d.E, d.var, maxlinkdim(ψ))
    priors = Any[ψu]; labels = ["UHF"]
    for c in (2,4,8,16,32,64,128,256,512)
        c < maxlinkdim(ψ) || continue
        σ = copy(ψ); orthogonalize!(σ,1); truncate!(σ; maxdim=c, cutoff=0.0); normalize!(σ)
        push!(priors, σ); push!(labels, "chi$c")
    end
    N = length(ψ)
    # 開放端: 中央の半分から始まるブロック(始点の偶奇で境界の結合の強弱が変わるので2通り)。周期: 1通り
    starts = PBC ? [1] : [N÷4 + 1, N÷4 + 2]
    lmax = PBC ? N÷2 : N÷2
    fn = joinpath(@__DIR__, "crm_eigparity_$tag.tsv")
    open(fn, "w") do io
        println(io, "W\tLX\tgeometry\tU\tHvar\tprior\tkind\tstart\tl\ttrue\tprior_val")
        for (kind, opn) in (("P", "F"), ("Q", "Fup")), s0 in starts
            tv = string_expect(ψ, opn, s0, lmax)
            for (lab, σ) in zip(labels, priors)
                pv = string_expect(σ, opn, s0, lmax)
                for l in eachindex(tv)
                    @printf(io, "%d\t%d\t%s\t%.1f\t%.3e\t%s\t%s\t%d\t%d\t%.12f\t%.12f\n", W, LX, geo, U, d.var, lab, kind, s0, l, tv[l], pv[l])
                end
            end
        end
    end
    println("書き出し: ", fn)
end

if get(ENV, "PARITY_SELFTEST", "0") == "1"
    d = deserialize(ENV["PARITY_STATE"]); ψ, ψu = d.ψ, d.ψu
    for (opn, s0) in (("F", 2), ("Fup", 3), ("F", 1))
        fast = string_expect(ψ, opn, s0, 6)
        slow = [string_bruteforce(ψ, opn, s0, l) for l in eachindex(fast)]
        @printf("  %-3s 始点 %d: 最大差 %.1e   値 %s\n", opn, s0, maximum(abs.(fast .- slow)), join([@sprintf("%+.5f", v) for v in fast], " "))
    end
    # ℓ=1 の F はオンサイト ZZ、ℓ=2 の Fup は隣接 Z↑Z↑ に一致するはず
    println("  L=8 U=8 の中央サイト 4: F(ℓ=1) = ", round(string_expect(ψ, "F", 4, 1)[1], digits=10), "   Fup(ℓ=2, 4-5) = ", round(string_expect(ψ, "Fup", 4, 2)[2], digits=10))
else
    run_parity()
end
