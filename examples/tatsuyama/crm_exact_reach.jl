# 「厳密な ρ」をどこまで大きい系で作れるかを測る。
#
# 現状の crm_2d_ed_fid_table.jl は cutoff=0 で走らせているので、結合次元は必ず
# 上限まで伸び、「本当はいくつ必要だったか」が分からない。ここでは打ち切り誤差を
# 与えて結合次元を適応的にし、到達 χ とエネルギー分散を測る。
#   W=1 LX=20 PBC=0 julia --project=. examples/tatsuyama/crm_exact_reach.jl
using ITensors, ITensorMPS, LinearAlgebra, Printf

const W    = parse(Int, get(ENV,"W","1"))
const LX   = parse(Int, get(ENV,"LX","16"))
const PBC  = get(ENV,"PBC","0") == "1"
const U    = parse(Float64, get(ENV,"CRM_U","8.0"))
const CUT  = parse(Float64, get(ENV,"CUTOFF","1e-13"))
const CHI  = parse(Int, get(ENV,"CHIMAX","8192"))
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

function main()
    n = LX*W; Nup = n÷2; Ndn = n÷2
    edges = lat_edges(LX, W; pbc_x=PBC)
    geo = PBC ? "torus" : "cylinder"
    sites = siteinds("Electron", n; conserve_qns=true)
    os = OpSum()
    for (a,b) in edges
        os += -1.0,"Cdagup",a,"Cup",b; os += -1.0,"Cdagup",b,"Cup",a
        os += -1.0,"Cdagdn",a,"Cdn",b; os += -1.0,"Cdagdn",b,"Cdn",a
    end
    for s in 1:n; os += U,"Nupdn",s end
    H = MPO(os, sites)
    st = [isodd(sum(divrem(s-1,W))) ? "Up" : "Dn" for s in 1:n]
    t0 = time()
    # 結合次元をランプさせる。必要な χ が小さい系で上限をそのまま使うと、
    # 初期スイープが無駄に重くなる(1D OBC は χ≈400 で足りるのに 8192 だと数時間)。
    sched = Int[]; c = 64
    while c < CHI; push!(sched, c); c *= 2 end
    push!(sched, CHI)
    nsw = parse(Int, get(ENV,"SWEEPS_PER","6"))
    md = vcat([fill(d, nsw) for d in sched]...)
    E, ψ = dmrg(H, random_mps(sites, st; linkdims=32);
                noise=[1e-6,1e-7,1e-8,1e-9,0.0],
                nsweeps=length(md)+30, maxdim=vcat(md, fill(CHI,30)),
                cutoff=CUT, outputlevel=0)
    normalize!(ψ)
    Hψ = apply(H, ψ; cutoff=1e-16, maxdim=4*CHI)
    var = real(inner(Hψ,Hψ) - E^2)
    dt = time() - t0
    @printf("%s %dx%d n=%2d  cutoff=%.0e  到達χ=%5d (上限 %d)  E0=%.10f  分散=%.2e  %.0f秒  %s\n",
            geo, LX, W, n, CUT, maxlinkdim(ψ), CHI, E, var, dt,
            abs(var) < 1e-8 ? "厳密" : "近似")
end
main()
