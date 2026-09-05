# 2次元(シリンダー / トーラス)で、大域忠実度と局所利得を同じ系で測る。
# crm_fid_vs_gain.jl の2次元版。UHF の Slater 行列式を MPS として組み、
# DMRG 参照状態との重なりから大域忠実度を出す。
#   W=4 LX=4 PBC=1 julia --project=. examples/tatsuyama/crm_2d_pbc_fid_vs_gain.jl
using ITensors, ITensorMPS, LinearAlgebra, Statistics, Printf

const W   = parse(Int, get(ENV,"W","4"))
const LX  = parse(Int, get(ENV,"LX","4"))
const PBC = get(ENV,"PBC","1") == "1"
const U   = parse(Float64, get(ENV,"CRM_U","8.0"))
const NM  = parse(Int, get(ENV,"CRM_NM","100"))
const CHIREF = parse(Int, get(ENV,"CHIREF","1024"))
sidx(x, y) = (x - 1) * W + y

function lat_edges(Lx, Wd; pbc_x::Bool)
    e = Tuple{Int,Int}[]
    for x in 1:Lx, y in 1:Wd
        s = sidx(x, y)
        if x < Lx;                push!(e, minmax(s, sidx(x+1, y)))
        elseif pbc_x && Lx > 2;   push!(e, minmax(s, sidx(1, y))) end
        if Wd > 2;                push!(e, minmax(s, sidx(x, y == Wd ? 1 : y+1)))
        elseif Wd == 2 && y == 1; push!(e, minmax(s, sidx(x, 2))) end
    end
    unique(e)
end

"""UHF 自己無撞着解。軌道係数 Φ も返す(Slater MPS を組むのに要る)。"""
function solve_uhf(edges, n, t, U; Nup, Ndn, iters=6000, tol=1e-13)
    T = zeros(n,n); for (a,b) in edges; T[a,b] = T[b,a] = -t; end
    coord = Dict(sidx(x,y) => (x,y) for x in 1:LX, y in 1:W)
    fill_avg = (Nup+Ndn)/(2n)
    nup = zeros(n); ndn = zeros(n)
    for s in 1:n
        (x,y) = coord[s]; g = isodd(x+y) ? 1 : -1
        nup[s] = fill_avg + 0.4g; ndn[s] = fill_avg - 0.4g
    end
    clamp!(nup,0.02,0.98); clamp!(ndn,0.02,0.98)
    nup .*= Nup/sum(nup); ndn .*= Ndn/sum(ndn)
    local Fu, Fd
    for _ in 1:iters
        Fu = eigen(Symmetric(T + diagm(U .* ndn))); Fd = eigen(Symmetric(T + diagm(U .* nup)))
        nu2 = vec(sum(abs2, Fu.vectors[:,1:Nup]; dims=2))
        nd2 = vec(sum(abs2, Fd.vectors[:,1:Ndn]; dims=2))
        if max(maximum(abs.(nu2 .- nup)), maximum(abs.(nd2 .- ndn))) < tol
            nup, ndn = nu2, nd2; break
        end
        nup = 0.5.*nu2 .+ 0.5.*nup; ndn = 0.5.*nd2 .+ 0.5.*ndn
    end
    (Φu=Fu.vectors[:,1:Nup], Φd=Fd.vectors[:,1:Ndn],
     m=mean(abs(nup[s]-ndn[s])/2 for s in 1:n))
end

"""Slater 行列式を MPS に。生成演算子はパリティ奇で MPO にできないので2つずつ。"""
function slater_mps(sites, Φu, Φd, Nup, Ndn; maxdim=1024, cutoff=1e-14)
    L = length(sites)
    ops = vcat([(:dn,k) for k in Ndn:-1:1], [(:up,k) for k in Nup:-1:1])
    @assert iseven(length(ops))
    ψ = MPS(sites, ["Emp" for _ in 1:L])
    for i in 1:2:length(ops)
        (σa,ka) = ops[i]; (σb,kb) = ops[i+1]
        va = σa === :up ? Φu[:,ka] : Φd[:,ka]
        vb = σb === :up ? Φu[:,kb] : Φd[:,kb]
        na = σa === :up ? "Cdagup" : "Cdagdn"
        nb = σb === :up ? "Cdagup" : "Cdagdn"
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

gain(P, Δ, nA) = begin
    a = 3.0^nA - 1; vs = 3.0^nA*(1-P^2)/NM
    (a*P^2 + vs) / (a*Δ^2 + vs)
end

function main()
    n = LX*W; Nup = n÷2; Ndn = n÷2
    edges = lat_edges(LX, W; pbc_x=PBC)
    geo = PBC ? "torus" : "cylinder"
    # 二部格子か(奇環があると反強磁性がフラストレートする)
    col = fill(0, n); col[1] = 1; q = [1]; bip = true
    adj = [Int[] for _ in 1:n]; for (a,b) in edges; push!(adj[a],b); push!(adj[b],a); end
    while !isempty(q)
        s = popfirst!(q)
        for t2 in adj[s]
            col[t2] == 0 ? (col[t2] = -col[s]; push!(q,t2)) : (col[t2] == col[s] && (bip = false))
        end
    end
    @printf("%s %dx%d  ボンド %d  二部格子: %s\n", geo, LX, W, length(edges), bip ? "はい" : "いいえ(奇環)")

    sites = siteinds("Electron", n; conserve_qns=true)
    os = OpSum()
    for (a,b) in edges
        os += -1.0,"Cdagup",a,"Cup",b;  os += -1.0,"Cdagup",b,"Cup",a
        os += -1.0,"Cdagdn",a,"Cdn",b;  os += -1.0,"Cdagdn",b,"Cdn",a
    end
    for s in 1:n; os += U,"Nupdn",s; end
    H = MPO(os, sites)
    st = [isodd(sum(divrem(s-1,W))) ? "Up" : "Dn" for s in 1:n]
    ψ0 = random_mps(sites, st; linkdims=32)
    E, ψ = dmrg(H, ψ0; nsweeps=40,
        maxdim=[32,64,128,256,512,CHIREF], cutoff=1e-12, outputlevel=0)
    normalize!(ψ)
    @printf("E0=%.10f  参照の結合次元 %d\n", E, maxlinkdim(ψ))

    uhf = solve_uhf(edges, n, 1.0, U; Nup, Ndn)
    ψu = slater_mps(sites, uhf.Φu, uhf.Φd, Nup, Ndn)
    @printf("UHF: 磁化 m=%.4f  結合次元 %d\n", uhf.m, maxlinkdim(ψu))

    chis = [2,4,8,16,32,64,128]
    priors = Any[]; labels = String[]; fids = Float64[]
    for c in chis
        σ = copy(ψ); orthogonalize!(σ,1); truncate!(σ; maxdim=c, cutoff=0.0); normalize!(σ)
        push!(priors,σ); push!(labels,"chi$c"); push!(fids, abs2(inner(σ,ψ)))
    end
    push!(priors, ψu); push!(labels,"UHF"); push!(fids, abs2(inner(ψu,ψ)))
    for (l,f) in zip(labels,fids); @printf("  %-8s F = %.6e\n", l, f); end

    c0 = sidx(max(1,LX÷2), max(1,W÷2))
    nb = first(b for (a,b) in edges if a == c0)          # c0 に接続する隣
    obs = Tuple{String,OpSum,Int}[]
    push!(obs, ("ZZ onsite", let q=OpSum()
        q += 4.0,"Nupdn",c0; q += -2.0,"Nup",c0; q += -2.0,"Ndn",c0; q += 1.0,"Id",c0; q end, 2))
    push!(obs, ("ZZ up-up nb", let q=OpSum()
        q += 4.0,"Nup",c0,"Nup",nb; q += -2.0,"Nup",c0; q += -2.0,"Nup",nb; q += 1.0,"Id",c0; q end, 2))
    push!(obs, ("SzSz nb", let q=OpSum(); q += 1.0,"Sz",c0,"Sz",nb; q end, 4))
    push!(obs, ("DoubleOcc", let q=OpSum(); q += 1.0,"Nupdn",c0; q end, 2))
    push!(obs, ("Sz", let q=OpSum(); q += 1.0,"Sz",c0; q end, 2))

    fn = joinpath(@__DIR__, @sprintf("crm_2d_pbc_fid_W%dL%d_%s_U%.1f.tsv", W, LX, geo, U))
    open(fn,"w") do io
        println(io, "W\tLX\tgeometry\tbipartite\tU\tprior\tglobal_fid\tobservable\tnA\ttrue\tprior_val\tDelta\teps\tG\tG_max")
        for (name, op, nA) in obs
            Op = MPO(op, sites); P = real(inner(ψ', Op, ψ))
            for (l,σ,F) in zip(labels, priors, fids)
                Pσ = real(inner(σ', Op, σ)); Δ = P - Pσ
                eps = abs(P) > 1e-12 ? abs(Δ)/abs(P) : NaN
                @printf(io, "%d\t%d\t%s\t%s\t%.1f\t%s\t%.10e\t%s\t%d\t%.10f\t%.10f\t%.6e\t%.6e\t%.6f\t%.6f\n",
                        W, LX, geo, bip, U, l, F, name, nA, P, Pσ, Δ, eps, gain(P,Δ,nA), gain(P,0.0,nA))
            end
        end
    end
    println("書き出し: ", fn)
end
main()
