# 大域忠実度と局所利得の対応を、同じ系で UHF と MPS prior の両方について測る。
#
# これまで UHF の大域忠実度はどのデータセットにも入っていなかった(prior_fid が NaN)。
# ここでは UHF の Slater 行列式を MPS として組み、DMRG 基底状態との重なりを直接取る。
# 出力: crm_fid_vs_gain_L{L}_U{U}.tsv  (§2m/図(a) 用)
using ITensors, ITensorMPS, LinearAlgebra, Printf, Random

# --- UHF 自己無撞着解(1次元鎖、Néel 初期条件) -------------------------
function solve_uhf_chain(L, t, U; Nup, Ndn, iters=6000, tol=1e-13)
    nup = [isodd(i) ? 0.9 : 0.1 for i in 1:L]
    ndn = [isodd(i) ? 0.1 : 0.9 for i in 1:L]
    evu = zeros(L,L); evd = zeros(L,L)
    for _ in 1:iters
        Hu = zeros(L,L); Hd = zeros(L,L)
        for i in 1:L-1; Hu[i,i+1]=Hu[i+1,i]=-t; Hd[i,i+1]=Hd[i+1,i]=-t; end
        for i in 1:L; Hu[i,i]=U*ndn[i]; Hd[i,i]=U*nup[i]; end
        Fu = eigen(Symmetric(Hu)); Fd = eigen(Symmetric(Hd))
        evu = Fu.vectors; evd = Fd.vectors
        nu = [sum(abs2, evu[i,1:Nup]) for i in 1:L]
        nd = [sum(abs2, evd[i,1:Ndn]) for i in 1:L]
        δ = maximum(abs.(nu.-nup)) + maximum(abs.(nd.-ndn))
        nup = 0.5*nup + 0.5*nu; ndn = 0.5*ndn + 0.5*nd
        δ < tol && break
    end
    return evu, evd, nup, ndn
end

# --- Slater 行列式を MPS として組む -------------------------------------
# 生成演算子1つはパリティ奇で MPO にできないので、2つずつ積にして作用させる。
# 並べ替えは全体の位相を変えるだけで、忠実度 |<σ|ψ>|^2 には影響しない。
function slater_mps(sites, evu, evd, Nup, Ndn; maxdim=600, cutoff=1e-15)
    L = length(sites)
    ops = vcat([(:dn,k) for k in Ndn:-1:1], [(:up,k) for k in Nup:-1:1])
    @assert iseven(length(ops)) "生成演算子の総数が奇数"
    ψ = MPS(sites, ["Emp" for _ in 1:L])
    for i in 1:2:length(ops)
        (σa,ka) = ops[i]; (σb,kb) = ops[i+1]
        va = σa === :up ? evu[:,ka] : evd[:,ka]
        vb = σb === :up ? evu[:,kb] : evd[:,kb]
        na = σa === :up ? "Cdagup" : "Cdagdn"
        nb = σb === :up ? "Cdagup" : "Cdagdn"
        os = OpSum()
        for s in 1:L, s2 in 1:L
            c = va[s]*vb[s2]
            if abs(c) > 1e-14
                os += c, na, s, nb, s2
            end
        end
        ψ = apply(MPO(os, sites), ψ; maxdim, cutoff)
        normalize!(ψ)
    end
    normalize!(ψ); return ψ
end

# --- 観測量(Pauli 列としての台の大きさ nA も返す) ---------------------
# Jordan-Wigner: n_{i,σ} = (1 - Z)/2。単一 Pauli 列のみ利得法則が厳密。
function observables(sites, L)
    c = L ÷ 2
    obs = Tuple{String,OpSum,Int,Bool}[]           # 名前, 演算子, 台の大きさ, 単一Pauli列か
    o = OpSum(); o += 4.0,"Nup",c; o += -2.0,"Ntot",c; o += 1.0,"Id",c
    push!(obs, ("ZZ onsite", let q=OpSum()
        q += 4.0,"Nupdn",c; q += -2.0,"Nup",c; q += -2.0,"Ndn",c; q += 1.0,"Id",c; q end, 2, true))
    push!(obs, ("ZZ up-up r=1", let q=OpSum()
        q += 4.0,"Nup",c,"Nup",c+1; q += -2.0,"Nup",c; q += -2.0,"Nup",c+1; q += 1.0,"Id",c; q end, 2, true))
    push!(obs, ("DoubleOcc", let q=OpSum(); q += 1.0,"Nupdn",c; q end, 2, false))
    push!(obs, ("SzSz r=1", let q=OpSum(); q += 1.0,"Sz",c,"Sz",c+1; q end, 4, false))
    push!(obs, ("Sz", let q=OpSum(); q += 1.0,"Sz",c; q end, 2, false))
    return obs
end

gain(P, Δ, nA, nm) = begin
    a = 3.0^nA - 1; vs = 3.0^nA*(1-P^2)/nm
    (a*P^2 + vs) / (a*Δ^2 + vs)
end

function main()
    L  = parse(Int, get(ENV,"CRM_L","8"))
    U  = parse(Float64, get(ENV,"CRM_U","8.0"))
    t  = 1.0
    nm = parse(Int, get(ENV,"CRM_NM","100"))
    Nup = L÷2; Ndn = L÷2
    sites = siteinds("Electron", L; conserve_qns=true)

    os = OpSum()
    for i in 1:L-1
        os += -t,"Cdagup",i,"Cup",i+1;  os += -t,"Cdagup",i+1,"Cup",i
        os += -t,"Cdagdn",i,"Cdn",i+1;  os += -t,"Cdagdn",i+1,"Cdn",i
    end
    for i in 1:L; os += U,"Nupdn",i; end
    H = MPO(os, sites)
    st = [isodd(i) ? "Up" : "Dn" for i in 1:L]
    ψ0 = random_mps(sites, st; linkdims=20)
    E, ψ = dmrg(H, ψ0; nsweeps=24,
                maxdim=[20,40,80,160,320,640,800], cutoff=1e-13, outputlevel=0)
    normalize!(ψ)
    @printf("L=%d U=%.1f  E0=%.10f  最大結合次元 %d\n", L, U, E, maxlinkdim(ψ))

    evu, evd, nup, ndn = solve_uhf_chain(L, t, U; Nup, Ndn)
    ψu = slater_mps(sites, evu, evd, Nup, Ndn)
    @printf("UHF: サイト磁化 m=%.4f  結合次元 %d\n", nup[1]-ndn[1], maxlinkdim(ψu))

    chis = [2,4,8,16,32,64]
    priors = Any[]; labels = String[]; fids = Float64[]
    for c in chis
        σ = copy(ψ); orthogonalize!(σ, 1)
        truncate!(σ; maxdim=c, cutoff=0.0); normalize!(σ)
        push!(priors, σ); push!(labels, "chi$c"); push!(fids, abs2(inner(σ, ψ)))
    end
    push!(priors, ψu); push!(labels, "UHF"); push!(fids, abs2(inner(ψu, ψ)))

    for (l,f) in zip(labels,fids); @printf("  %-8s 大域忠実度 F = %.6e\n", l, f); end

    obs = observables(sites, L)
    open(joinpath(@__DIR__, @sprintf("crm_fid_vs_gain_L%d_U%.1f.tsv", L, U)), "w") do io
        println(io, "L\tU\tprior\tglobal_fid\tobservable\tnA\tsingle_pauli\ttrue\tprior_val\tDelta\teps\tG\tG_max")
        for (name, op, nA, single) in obs
            Op = MPO(op, sites)
            P  = real(inner(ψ', Op, ψ))
            for (l,σ,F) in zip(labels, priors, fids)
                Pσ = real(inner(σ', Op, σ))
                Δ  = P - Pσ
                G  = gain(P, Δ, nA, nm); Gm = gain(P, 0.0, nA, nm)
                eps = abs(P) > 1e-12 ? abs(Δ)/abs(P) : NaN
                @printf(io, "%d\t%.1f\t%s\t%.10e\t%s\t%d\t%s\t%.10f\t%.10f\t%.6e\t%.6e\t%.6f\t%.6f\n",
                        L, U, l, F, name, nA, single, P, Pσ, Δ, eps, G, Gm)
            end
        end
    end
    println("書き出し: ", @sprintf("crm_fid_vs_gain_L%d_U%.1f.tsv", L, U))
end
main()
