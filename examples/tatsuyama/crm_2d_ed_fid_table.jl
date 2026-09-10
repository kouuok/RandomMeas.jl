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
const CUTOFF = parse(Float64, get(ENV,"CUTOFF","0.0"))  # 0 なら切断なし(小さい系)
const SWPER  = parse(Int, get(ENV,"SWEEPS_PER","6"))
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

"""サイト集合 ss(1〜2サイト)の縮約密度行列を、ITensor の Electron 局所基底
   |0>,|up>,|dn>,|updn> の順序のまま返す。"""
function site_rdm(ψ::MPS, ss::Vector{Int})
    ϕ = orthogonalize(ψ, first(ss)); s = siteinds(ψ)
    if length(ss) == 1
        A = ϕ[ss[1]]; ρ = prime(A, s[ss[1]]) * dag(A)
        return Array(ρ, prime(s[ss[1]]), s[ss[1]])
    else
        a, b = ss
        A = ϕ[a]; for k in a+1:b; A *= ϕ[k] end
        ρ = prime(A, s[a], s[b]) * dag(A)
        M = Array(ρ, prime(s[a]), prime(s[b]), s[a], s[b])   # (ia', ib', ia, ib)
        out = zeros(ComplexF64, 16, 16)
        for ia in 1:4, ib in 1:4, ja in 1:4, jb in 1:4
            out[(ia-1)*4+ib, (ja-1)*4+jb] = M[ia, ib, ja, jb]   # |a>⊗|b> の順に並べ直す
        end
        return out
    end
end

# Electron 局所基底の番号 i(1..4)を (上向き占有, 下向き占有) に分解する。
#   1=|0>→(0,0)  2=|up>→(1,0)  3=|dn>→(0,1)  4=|updn>→(1,1)
locbits(i) = ((i-1) & 1, ((i-1) >> 1) & 1)

"""k サイト(k≤2)の RDM から、指定したモードだけに縮約する。
   modes は (サイトの何番目, :up または :dn) の並び。残りのモードはトレースする。"""
function mode_rdm(ρ::Matrix, k::Int, modes::Vector{Tuple{Int,Symbol}})
    nk = length(modes); d = 4^k
    @assert size(ρ,1) == d
    sel(i_vec, (p, m)) = m === :up ? locbits(i_vec[p])[1] : locbits(i_vec[p])[2]
    allm = [(p,m) for p in 1:k for m in (:up,:dn)]
    tr_  = setdiff(allm, modes)
    unpack(I) = k == 1 ? (I,) : (div(I-1,4)+1, mod(I-1,4)+1)
    idx(v) = 1 + sum((sel(v,mm) << (nk-t)) for (t,mm) in enumerate(modes); init=0)
    out = zeros(ComplexF64, 2^nk, 2^nk)
    for I in 1:d, J in 1:d
        vi = unpack(I); vj = unpack(J)
        all(sel(vi,mm) == sel(vj,mm) for mm in tr_) || continue
        out[idx(vi), idx(vj)] += ρ[I, J]
    end
    out
end

"""自己検証: 縮約した RDM から求めた期待値が MPO の値と一致するか。"""
function check_rdm(name, ρq, modes, P_mpo)
    nk = length(modes)
    diagv = [prod((( (i >> (nk-t)) & 1) == 1 ? -1.0 : 1.0) for t in 1:nk) for i in 0:2^nk-1]
    val = real(sum(diagv[i+1]*ρq[i+1,i+1] for i in 0:2^nk-1))   # Z...Z の期待値
    (name, val, P_mpo)
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
    chimax = CHIMAX > 0 ? CHIMAX : min(n÷2 >= 32 ? typemax(Int) : 4^(n÷2), 8192)
    @printf("%s %dx%d  n=%d  ボンド %d  二部格子:%s  最悪値 χ=4^%d (使用 %d)\n",
            geo, LX, W, n, length(edges), bip ? "はい" : "いいえ(奇環)", n÷2, chimax)

    sites = siteinds("Electron", n; conserve_qns=true)
    os = OpSum()
    for (a,b) in edges
        os += -1.0,"Cdagup",a,"Cup",b; os += -1.0,"Cdagup",b,"Cup",a
        os += -1.0,"Cdagdn",a,"Cdn",b; os += -1.0,"Cdagdn",b,"Cdn",a
    end
    for s in 1:n; os += U,"Nupdn",s end
    H = MPO(os, sites)
    st = [isodd(sum(divrem(s-1,W))) ? "Up" : "Dn" for s in 1:n]
    # 結合次元をランプさせる。必要な χ が小さい系で上限から始めると初期スイープが
    # 無駄に重くなる(L=28 は χ≈380 で足りるのに上限 8192 だと5.5時間、ランプ後は数分)。
    # ノイズ項は Néel 初期状態からの局所解を抜けるために要る(L≥40 で収束不足になった)。
    sched = Int[]; cc = 64
    while cc < chimax; push!(sched, cc); cc *= 2 end
    push!(sched, chimax)
    md = vcat([fill(d, SWPER) for d in sched]...)
    E, ψ = dmrg(H, random_mps(sites, st; linkdims=32);
                nsweeps=length(md)+30, maxdim=vcat(md, fill(chimax,30)),
                cutoff=CUTOFF, noise=[1e-6,1e-7,1e-8,1e-9,0.0], outputlevel=0)
    normalize!(ψ)
    Hψ = apply(H, ψ; cutoff=1e-16)
    var = real(inner(Hψ,Hψ) - E^2)
    @printf("E0=%.12f  到達χ=%d  <H^2>-<H>^2=%.3e  (厳密なら 0)\n", E, maxlinkdim(ψ), var)

    uhf = solve_uhf(edges, n, 1.0, U; Nup, Ndn)
    ψu = slater_mps(sites, uhf.Φu, uhf.Φd, Nup, Ndn)
    @printf("UHF: 磁化 m=%.4f  χ=%d\n", uhf.m, maxlinkdim(ψu))

    chis = [2,4,8,16,32,64,128,256,512]
    filter!(c -> c < maxlinkdim(ψ), chis)
    priors = Any[]; labels = String[]; fids = Float64[]
    for c in chis
        σ = copy(ψ); orthogonalize!(σ,1); truncate!(σ; maxdim=c, cutoff=0.0); normalize!(σ)
        push!(priors,σ); push!(labels,"chi$c"); push!(fids, abs2(inner(σ,ψ)))
    end
    push!(priors, ψu); push!(labels,"UHF"); push!(fids, abs2(inner(ψu,ψ)))

    c0 = sidx(max(1,cld(LX,2)), max(1,cld(W,2)))
    nb = first(b for (a,b) in edges if a == c0)
    # (名前, 演算子, キュービット台の大きさ, 台のサイト, 残すモード)
    ss1 = [c0]; ss2 = sort([c0,nb]); pa = c0 < nb ? 1 : 2; pb = 3 - pa
    obs = Any[
      ("ZZ onsite", let q=OpSum(); q += 4.0,"Nupdn",c0; q += -2.0,"Nup",c0; q += -2.0,"Ndn",c0; q += 1.0,"Id",c0; q end,
        2, ss1, [(1,:up),(1,:dn)]),
      ("ZZ up-up nb", let q=OpSum(); q += 4.0,"Nup",c0,"Nup",nb; q += -2.0,"Nup",c0; q += -2.0,"Nup",nb; q += 1.0,"Id",c0; q end,
        2, ss2, [(pa,:up),(pb,:up)]),
      ("SzSz nb", let q=OpSum(); q += 1.0,"Sz",c0,"Sz",nb; q end,
        4, ss2, [(pa,:up),(pa,:dn),(pb,:up),(pb,:dn)]),
      ("SxSx nb", let q=OpSum()
          q += 0.25,"S+",c0,"S+",nb; q += 0.25,"S+",c0,"S-",nb
          q += 0.25,"S-",c0,"S+",nb; q += 0.25,"S-",c0,"S-",nb; q end,
        4, ss2, [(pa,:up),(pa,:dn),(pb,:up),(pb,:dn)]),
      ("DoubleOcc", let q=OpSum(); q += 1.0,"Nupdn",c0; q end, 2, ss1, [(1,:up),(1,:dn)]),
      ("n", let q=OpSum(); q += 1.0,"Nup",c0; q += 1.0,"Ndn",c0; q end, 2, ss1, [(1,:up),(1,:dn)]),
      ("Sz", let q=OpSum(); q += 1.0,"Sz",c0; q end, 2, ss1, [(1,:up),(1,:dn)]) ]

    # 自己検証: ZZ onsite / ZZ up-up は Z...Z 型なので RDM から出した値が MPO と一致するはず
    for (name, op, nA, ss, modes) in obs
        name in ("ZZ onsite","ZZ up-up nb") || continue
        Pm = real(inner(ψ', MPO(op, sites), ψ))
        ρq = mode_rdm(site_rdm(ψ, ss), length(ss), modes)
        (_, v, _) = check_rdm(name, ρq, modes, Pm)
        @printf("  検証 %-14s RDMから %+.10f  MPOから %+.10f  差 %.2e\n", name, v, Pm, abs(v-Pm))
        abs(v-Pm) < 1e-8 || error("台の縮約が MPO と一致しない: $name")
    end

    fn = joinpath(@__DIR__, @sprintf("crm_2d_edfid_W%dL%d_%s_U%.1f.tsv", W, LX, geo, U))
    open(fn,"w") do io
        println(io, "W\tLX\tnsites\tgeometry\tbipartite\tU\tEref\tHvar\tprior\tglobal_fid\t"*
                    "observable\tnA\tF_supp\tD_supp\ttrue\tprior_val\tDelta\teps\tG\tG_max")
        for (name, op, nA, ss, modes) in obs
            Op = MPO(op, sites); P = real(inner(ψ', Op, ψ))
            ρq = mode_rdm(site_rdm(ψ, ss), length(ss), modes)
            for (l,σ,F) in zip(labels, priors, fids)
                Pσ = real(inner(σ', Op, σ)); Δ = P - Pσ
                σq = mode_rdm(site_rdm(σ, ss), length(ss), modes)
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
