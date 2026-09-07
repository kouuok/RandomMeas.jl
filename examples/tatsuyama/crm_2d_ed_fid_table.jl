# 2次元(シリンダー / トーラス)で、ρ を厳密基底状態とし、σ に HF と各ボンド次元の
# MPS を使ったときの「利得・大域忠実度・台の忠実度」を1つの表にまとめる。
#
# ρ の厳密性: cutoff=0、maxdim ≥ 4^(n/2) で DMRG を回す。この結合次元では MPS は
# セクター内の任意の状態を誤差なく表せるので切断が起こりえない。到達結合次元と
# エネルギー分散 <H^2>-<H>^2 を出力して、実際に固有状態であることを確認する。
#   W=2 LX=4 PBC=0 julia --project=. examples/tatsuyama/crm_2d_ed_fid_table.jl
using ITensors, ITensorMPS, LinearAlgebra, Statistics, Printf

const W   = parse(Int, get(ENV,"W","2"))
const LX  = parse(Int, get(ENV,"LX","4"))
const PBC = get(ENV,"PBC","0") == "1"
const U   = parse(Float64, get(ENV,"CRM_U","8.0"))
const NM  = parse(Int, get(ENV,"CRM_NM","100"))
const CHIMAX = parse(Int, get(ENV,"CHIMAX","0"))   # 0 なら 4^(n/2) を使う
sidx(x,y) = (x-1)*W + y

function lat_edges(Lx, Wd; pbc_x::Bool)
    e = Tuple{Int,Int}[]
    for x in 1:Lx, y in 1:Wd
        s = sidx(x,y)
        if x < Lx;                push!(e, minmax(s, sidx(x+1,y)))
        elseif pbc_x && Lx > 2;   push!(e, minmax(s, sidx(1,y))) end
        if Wd > 2;                push!(e, minmax(s, sidx(x, y == Wd ? 1 : y+1)))
        elseif Wd == 2 && y == 1; push!(e, minmax(s, sidx(x,2))) end
    end
    unique(e)
end

function bipartite(edges, n)
    col = fill(0,n); col[1] = 1; q = [1]; ok = true
    adj = [Int[] for _ in 1:n]; for (a,b) in edges; push!(adj[a],b); push!(adj[b],a) end
    while !isempty(q)
        s = popfirst!(q)
        for t in adj[s]
            col[t] == 0 ? (col[t] = -col[s]; push!(q,t)) : (col[t] == col[s] && (ok = false))
        end
    end
    ok
end

function solve_uhf(edges, n, t, U; Nup, Ndn, iters=8000, tol=1e-13)
    T = zeros(n,n); for (a,b) in edges; T[a,b] = T[b,a] = -t end
    nup = zeros(n); ndn = zeros(n); fa = (Nup+Ndn)/(2n)
    for x in 1:LX, y in 1:W
        s = sidx(x,y); g = isodd(x+y) ? 1 : -1
        nup[s] = fa + 0.4g; ndn[s] = fa - 0.4g
    end
    clamp!(nup,0.02,0.98); clamp!(ndn,0.02,0.98)
    nup .*= Nup/sum(nup); ndn .*= Ndn/sum(ndn)
    local Fu, Fd
    for _ in 1:iters
        Fu = eigen(Symmetric(T + diagm(U .* ndn))); Fd = eigen(Symmetric(T + diagm(U .* nup)))
        nu = vec(sum(abs2, Fu.vectors[:,1:Nup]; dims=2))
        nd = vec(sum(abs2, Fd.vectors[:,1:Ndn]; dims=2))
        if max(maximum(abs.(nu.-nup)), maximum(abs.(nd.-ndn))) < tol
            nup, ndn = nu, nd; break
        end
        nup = 0.5.*nu .+ 0.5.*nup; ndn = 0.5.*nd .+ 0.5.*ndn
    end
    (Φu=Fu.vectors[:,1:Nup], Φd=Fd.vectors[:,1:Ndn], m=mean(abs.(nup.-ndn))/2)
end

"""Slater 行列式を MPS に(生成演算子はパリティ奇なので2つずつ積にする)。"""
function slater_mps(sites, Φu, Φd, Nup, Ndn; maxdim=2048, cutoff=1e-15)
    L = length(sites)
    ops = vcat([(:dn,k) for k in Ndn:-1:1], [(:up,k) for k in Nup:-1:1])
    ψ = MPS(sites, ["Emp" for _ in 1:L])
    for i in 1:2:length(ops)
        (σa,ka) = ops[i]; (σb,kb) = ops[i+1]
        va = σa === :up ? Φu[:,ka] : Φd[:,ka]; vb = σb === :up ? Φu[:,kb] : Φd[:,kb]
        na = σa === :up ? "Cdagup" : "Cdagdn"; nb = σb === :up ? "Cdagup" : "Cdagdn"
        os = OpSum()
        for s in 1:L, s2 in 1:L
            c = va[s]*vb[s2]
            if abs(c) > 1e-13
                os += c, na, s, nb, s2
            end
        end
        ψ = apply(MPO(os, sites), ψ; maxdim, cutoff); normalize!(ψ)
    end
    normalize!(ψ); ψ
end

"""サイト集合 ss(1〜2サイト)の縮約密度行列。ITensor の Electron 基底
   |0>,|up>,|dn>,|updn> は そのまま (upキュービット ⊗ dnキュービット) である。"""
function site_rdm(ψ::MPS, ss::Vector{Int})
    ϕ = orthogonalize(ψ, first(ss))
    s = siteinds(ψ)
    if length(ss) == 1
        A = ϕ[ss[1]]
        ρ = prime(A, s[ss[1]]) * dag(A)
        return Array(ρ, prime(s[ss[1]]), s[ss[1]])
    else
        a, b = ss
        A = ϕ[a]; for k in a+1:b; A *= ϕ[k] end
        ρ = prime(A, s[a], s[b]) * dag(A)
        M = Array(ρ, prime(s[a]), prime(s[b]), s[a], s[b])
        return reshape(permutedims(M, (1,2,3,4)), 16, 16)
    end
end

"""4^k の RDM から、残すキュービット集合 keep(サイトごとに (:up,:dn) を並べた
   2k キュービットのうち)に縮約する。"""
function qubit_rdm(ρ::Matrix, nq::Int, keep::Vector{Int})
    d = 2^nq; @assert size(ρ,1) == d
    out = zeros(ComplexF64, 2^length(keep), 2^length(keep))
    tr_ = setdiff(1:nq, keep)
    bit(i,q) = (i >> (nq-q)) & 1
    idx(i) = sum(bit(i,q) << (length(keep)-k) for (k,q) in enumerate(keep); init=0)
    for i in 0:d-1, j in 0:d-1
        all(bit(i,q) == bit(j,q) for q in tr_) || continue
        out[idx(i)+1, idx(j)+1] += ρ[i+1, j+1]
    end
    out
end

fid(ρ::Matrix, σ::Matrix) = begin
    w, V = eigen(Hermitian(ρ)); w = clamp.(real(w), 0, Inf)
    sr = V * Diagonal(sqrt.(w)) * V'
    m = clamp.(real(eigvals(Hermitian(sr*σ*sr))), 0, Inf)
    real(sum(sqrt.(m))^2)
end
gain(P, Δ, nA) = begin
    a = 3.0^nA - 1; vs = 3.0^nA*(1-P^2)/NM
    (a*P^2 + vs)/(a*Δ^2 + vs)
end

function main()
    n = LX*W; Nup = n÷2; Ndn = n÷2
    edges = lat_edges(LX, W; pbc_x=PBC)
    geo = PBC ? "torus" : "cylinder"
    bip = bipartite(edges, n)
    chimax = CHIMAX > 0 ? CHIMAX : min(4^(n÷2), 8192)
    @printf("%s %dx%d  n=%d  ボンド %d  二部格子:%s  厳密化に要る χ=4^%d=%d (使用 %d)\n",
            geo, LX, W, n, length(edges), bip ? "はい" : "いいえ(奇環)", n÷2, 4^(n÷2), chimax)

    sites = siteinds("Electron", n; conserve_qns=true)
    os = OpSum()
    for (a,b) in edges
        os += -1.0,"Cdagup",a,"Cup",b; os += -1.0,"Cdagup",b,"Cup",a
        os += -1.0,"Cdagdn",a,"Cdn",b; os += -1.0,"Cdagdn",b,"Cdn",a
    end
    for s in 1:n; os += U,"Nupdn",s end
    H = MPO(os, sites)
    st = [isodd(sum(divrem(s-1,W))) ? "Up" : "Dn" for s in 1:n]
    E, ψ = dmrg(H, random_mps(sites, st; linkdims=32);
                nsweeps=60, maxdim=chimax, cutoff=0.0, outputlevel=0)
    normalize!(ψ)
    Hψ = apply(H, ψ; cutoff=1e-16)
    var = real(inner(Hψ,Hψ) - E^2)
    @printf("E0=%.12f  到達χ=%d  <H^2>-<H>^2=%.3e  (厳密なら 0)\n", E, maxlinkdim(ψ), var)

    uhf = solve_uhf(edges, n, 1.0, U; Nup, Ndn)
    ψu = slater_mps(sites, uhf.Φu, uhf.Φd, Nup, Ndn)
    @printf("UHF: 磁化 m=%.4f  χ=%d\n", uhf.m, maxlinkdim(ψu))

    chis = [2,4,8,16,32,64,128,256]
    filter!(c -> c < maxlinkdim(ψ), chis)
    priors = Any[]; labels = String[]; fids = Float64[]
    for c in chis
        σ = copy(ψ); orthogonalize!(σ,1); truncate!(σ; maxdim=c, cutoff=0.0); normalize!(σ)
        push!(priors,σ); push!(labels,"chi$c"); push!(fids, abs2(inner(σ,ψ)))
    end
    push!(priors, ψu); push!(labels,"UHF"); push!(fids, abs2(inner(ψu,ψ)))

    c0 = sidx(max(1,cld(LX,2)), max(1,cld(W,2)))
    nb = first(b for (a,b) in edges if a == c0)
    # (名前, 演算子, キュービット台の大きさ, 台のサイト, サイト内で残すキュービット)
    obs = Any[
      ("ZZ onsite", let q=OpSum(); q += 4.0,"Nupdn",c0; q += -2.0,"Nup",c0; q += -2.0,"Ndn",c0; q += 1.0,"Id",c0; q end, 2, [c0], [1,2]),
      ("ZZ up-up nb", let q=OpSum(); q += 4.0,"Nup",c0,"Nup",nb; q += -2.0,"Nup",c0; q += -2.0,"Nup",nb; q += 1.0,"Id",c0; q end, 2, sort([c0,nb]), [1,3]),
      ("SzSz nb", let q=OpSum(); q += 1.0,"Sz",c0,"Sz",nb; q end, 4, sort([c0,nb]), [1,2,3,4]),
      ("DoubleOcc", let q=OpSum(); q += 1.0,"Nupdn",c0; q end, 2, [c0], [1,2]),
      ("Sz", let q=OpSum(); q += 1.0,"Sz",c0; q end, 2, [c0], [1,2]) ]

    fn = joinpath(@__DIR__, @sprintf("crm_2d_edfid_W%dL%d_%s_U%.1f.tsv", W, LX, geo, U))
    open(fn,"w") do io
        println(io, "W\tLX\tnsites\tgeometry\tbipartite\tU\tEref\tHvar\tprior\tglobal_fid\t"*
                    "observable\tnA\tF_supp\tD_supp\ttrue\tprior_val\tDelta\teps\tG\tG_max")
        for (name, op, nA, ss, keep) in obs
            Op = MPO(op, sites); P = real(inner(ψ', Op, ψ))
            nq = 2*length(ss)
            ρq = qubit_rdm(site_rdm(ψ, ss), nq, keep)
            for (l,σ,F) in zip(labels, priors, fids)
                Pσ = real(inner(σ', Op, σ)); Δ = P - Pσ
                σq = qubit_rdm(site_rdm(σ, ss), nq, keep)
                Fs = fid(ρq, σq); Ds = 0.5*sum(abs, eigvals(Hermitian(ρq - σq)))
                eps = abs(P) > 1e-12 ? abs(Δ)/abs(P) : NaN
                @printf(io, "%d\t%d\t%d\t%s\t%s\t%.1f\t%.10f\t%.3e\t%s\t%.10e\t%s\t%d\t%.8f\t%.8e\t%.10f\t%.10f\t%.6e\t%.6e\t%.6f\t%.6f\n",
                        W, LX, n, geo, bip, U, E, var, l, F, name, nA, Fs, Ds, P, Pσ, Δ, eps, gain(P,Δ,nA), gain(P,0.0,nA))
            end
        end
    end
    println("書き出し: ", fn)
end
main()
