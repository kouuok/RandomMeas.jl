# 低い結合次元の MPS を prior にしたときの値を、スピン射影 HF と同じ表(crm_phf_compare.py)に並べるために求める。
#
#   変分 MPS: 結合次元を χ に制限した DMRG の基底状態。厳密な状態を知らなくても作れる prior。
#   切断 MPS: 厳密な状態を χ に切り詰めたもの(統合表 crm_edfid_all.tsv の chi*)。中央の値はもう統合表に
#             あるので、ここではエネルギーを足し、照合のために中央の値も出し直す。
#
# 系ごとに Julia を起動するとコンパイルに時間がかかるので、1つのプロセスで系を順に回す。
#   SYSTEMS="1:cylinder:8:t;1:torus:16:tp" US="4,8,12" CHIS="2,4,8,16" julia --project=. crm_mps_lowchi.jl
# フラグ  t: 切断 MPS(states/eigcorr_*.jls の厳密な状態。なければ 16 サイト以下に限り DMRG で作る)。
#            TRUNC_US に含まれる U だけで作る。
#         p: 全サイト対の Z 相関(crm_eigen_corr.jl と同じ量)
#         q: ブロックの偶奇の文字列(crm_eigen_parity.jl と同じ量)
using ITensors, ITensorMPS, LinearAlgebra, Printf, Random, Serialization

const CHIS   = parse.(Int, split(get(ENV, "CHIS", "2,4,8,16"), ","))
const US     = parse.(Float64, split(get(ENV, "US", "2,4,8,12"), ","))
const NSTART = parse(Int, get(ENV, "NSTART", "4"))
const NSWEEP = parse(Int, get(ENV, "NSWEEP", "60"))
const TRUNC_US = parse.(Float64, split(get(ENV, "TRUNC_US", get(ENV, "US", "2,4,8,12")), ","))   # 切断 MPS を作る U

"""crm_2d_ed_fid_table.jl の lat_edges と同じ辺の並び(幅 Wd を定数ではなく引数にした)。"""
function edges_of(Lx, Wd, pbc_x)
    sidx(x, y) = (x-1)*Wd + y
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

function hubbard_mpo(sites, edges, n, U)
    os = OpSum()
    for (a, b) in edges
        os += -1.0,"Cdagup",a,"Cup",b; os += -1.0,"Cdagup",b,"Cup",a
        os += -1.0,"Cdagdn",a,"Cdn",b; os += -1.0,"Cdagdn",b,"Cdn",a
    end
    for s in 1:n; os += U,"Nupdn",s end
    MPO(os, sites)
end

"""統合表と同じ中央のサイトと結合の観測量(演算子の定義も crm_2d_ed_fid_table.jl と同じ)。"""
function centre_ops(sites, c0, nb)
    ops = [
      ("ZZ onsite", let q=OpSum(); q += 4.0,"Nupdn",c0; q += -2.0,"Nup",c0; q += -2.0,"Ndn",c0; q += 1.0,"Id",c0; q end),
      ("ZZ up-up nb", let q=OpSum(); q += 4.0,"Nup",c0,"Nup",nb; q += -2.0,"Nup",c0; q += -2.0,"Nup",nb; q += 1.0,"Id",c0; q end),
      ("SzSz nb", let q=OpSum(); q += 1.0,"Sz",c0,"Sz",nb; q end),
      ("SxSx nb", let q=OpSum()
          q += 0.25,"S+",c0,"S+",nb; q += 0.25,"S+",c0,"S-",nb
          q += 0.25,"S-",c0,"S+",nb; q += 0.25,"S-",c0,"S-",nb; q end),
      ("DoubleOcc", let q=OpSum(); q += 1.0,"Nupdn",c0; q end),
      ("n", let q=OpSum(); q += 1.0,"Nup",c0; q += 1.0,"Ndn",c0; q end),
      ("Sz", let q=OpSum(); q += 1.0,"Sz",c0; q end) ]
    [(name, MPO(q, sites)) for (name, q) in ops]
end

"""crm_eigen_corr.jl と同じ: 全サイト対の Z↑、Z↑Z↑、Z↑Z↓。"""
function zobs(ψ::MPS)
    nu = expect(ψ, "Nup"); nd = expect(ψ, "Ndn")
    Cuu = correlation_matrix(ψ, "Nup", "Nup"); Cud = correlation_matrix(ψ, "Nup", "Ndn")
    n = length(ψ); out = Dict{Tuple{String,Int,Int},Float64}()
    for i in 1:n
        out[("Zup",i,i)] = 1 - 2nu[i]
        for j in 1:n
            i < j && (out[("ZupZup",i,j)] = 1 - 2nu[i] - 2nu[j] + 4real(Cuu[i,j]))
            out[("ZupZdn",i,j)] = 1 - 2nu[i] - 2nd[j] + 4real(Cud[i,j])
        end
    end
    out
end

"""crm_eigen_parity.jl と同じ: s0 から右へ ℓ=1..ℓmax の Π O_j を一度に求める。"""
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

"""結合次元を χ に制限した DMRG。局所解を避けるため、Néel の積状態と乱数の MPS から始めて、一番低いものを採る。"""
function lowchi_ground(H, sites, st, χ)
    nz = vcat(fill(1e-4, 4), fill(1e-5, 4), fill(1e-6, 4), fill(1e-7, 4), fill(1e-8, 4), 0.0)
    best = nothing; Es = Float64[]
    for k in 1:NSTART
        Random.seed!(7919k + χ)
        ψ0 = k == 1 ? MPS(sites, st) : random_mps(sites, st; linkdims=χ)
        _, ψ = dmrg(H, ψ0; nsweeps=NSWEEP, maxdim=χ, cutoff=0.0, noise=nz, outputlevel=0)
        normalize!(ψ); E = real(inner(ψ', H, ψ)); push!(Es, E)
        (best === nothing || E < best[1]) && (best = (E, ψ))
    end
    best[1], best[2], Es
end

"""保存済みの厳密な状態がない小さい系(16 サイト以下)で、統合表と同じ手順(ランプ + ノイズ項)の DMRG を回す。"""
function exact_dmrg(H, sites, st)
    chimax = 2048; sched = Int[]; cc = 64
    while cc < chimax; push!(sched, cc); cc *= 2 end
    push!(sched, chimax)
    md = vcat([fill(d, 6) for d in sched]...)
    E, ψ = dmrg(H, random_mps(sites, st; linkdims=32); nsweeps=length(md)+30, maxdim=vcat(md, fill(chimax, 30)),
                cutoff=1e-14, noise=[1e-6,1e-7,1e-8,1e-9,0.0], outputlevel=0)
    normalize!(ψ); Hψ = apply(H, ψ; cutoff=1e-16); var = real(inner(Hψ, Hψ) - E^2)
    ψ, E, var
end

function run_system(Wd, geo, Lx, flags, U)
    t0 = time()
    n = Wd*Lx; pbc = geo == "torus"; edges = edges_of(Lx, Wd, pbc)
    tag = @sprintf("W%dL%d_%s_U%.1f", Wd, Lx, geo, U)
    c0 = (max(1, cld(Lx, 2))-1)*Wd + max(1, cld(Wd, 2)); nb = first(b for (a, b) in edges if a == c0)
    st = [isodd(sum(divrem(s-1, Wd))) ? "Up" : "Dn" for s in 1:n]
    rows = Tuple{String,Int,String,Int,Int,Float64}[]            # (method, χ, quantity, i, j, value)

    # 切断 MPS を作る系では、厳密な状態のサイトの添字に合わせて H を作る(保存済みの状態は1回だけ読む)
    ψex = nothing
    if 't' in flags && U in TRUNC_US
        f = joinpath(@__DIR__, "states", "eigcorr_$tag.jls")
        if isfile(f)
            d = deserialize(f); ψex, Eex, var = d.ψ, d.E, d.var
            sites = siteinds(ψex); H = hubbard_mpo(sites, edges, n, U)
        else
            n <= 16 || error("厳密な状態がない: $tag")
            sites = siteinds("Electron", n; conserve_qns=true); H = hubbard_mpo(sites, edges, n, U)
            ψex, Eex, var = exact_dmrg(H, sites, st)
        end
        @printf("%s  厳密: E0=%.12f 分散=%.2e χ=%d (%.0f 秒)\n", tag, Eex, var, maxlinkdim(ψex), time()-t0)
        push!(rows, ("exact", 0, "E", 0, 0, Eex), ("exact", 0, "Hvar", 0, 0, var),
                    ("exact", 0, "chi_reached", 0, 0, maxlinkdim(ψex)))
    else
        sites = siteinds("Electron", n; conserve_qns=true)
        H = hubbard_mpo(sites, edges, n, U)
    end
    ops = centre_ops(sites, c0, nb)

    # 全サイト対とブロックの偶奇は変分 MPS だけ(厳密な状態と切断 MPS の値は crm_eigcorr_*.tsv / crm_eigparity_*.tsv にある)
    function record!(method, χ, σ; full=false)
        push!(rows, (method, χ, "E", 0, 0, real(inner(σ', H, σ))), (method, χ, "chi_reached", 0, 0, maxlinkdim(σ)))
        for (name, Op) in ops
            push!(rows, (method, χ, name, c0, nb, real(inner(σ', Op, σ))))
        end
        sz = expect(σ, "Sz"); push!(rows, (method, χ, "mz_max", 0, 0, maximum(abs, sz)))
        full || return
        if 'p' in flags
            for ((kind, i, j), v) in zobs(σ)
                push!(rows, (method, χ, kind, i, j, v))
            end
        end
        if 'q' in flags
            for (kind, opn) in (("P", "F"), ("Q", "Fup")), s0 in (pbc ? [1] : [n÷4 + 1, n÷4 + 2])
                for (l, v) in enumerate(string_expect(σ, opn, s0, n÷2))
                    push!(rows, (method, χ, kind, s0, l, v))
                end
            end
        end
    end

    if ψex !== nothing
        for χ in CHIS
            χ < maxlinkdim(ψex) || continue
            σ = copy(ψex); orthogonalize!(σ, 1); truncate!(σ; maxdim=χ, cutoff=0.0); normalize!(σ)
            record!("trunc", χ, σ)
        end
    end
    for χ in CHIS
        E, σ, Es = lowchi_ground(H, sites, st, χ)
        record!("var", χ, σ; full=true)
        for (k, e) in enumerate(Es); push!(rows, ("var", χ, "E_start", k, 0, e)) end
        mkpath(joinpath(@__DIR__, "states"))
        serialize(joinpath(@__DIR__, "states", "mpslow_$(tag)_chi$χ.jls"), (ψ=σ, E=E, Es=Es))
        @printf("%s  変分 χ=%-3d E=%.8f  出発点ごと %s  (%.0f 秒)\n", tag, χ, E,
                join([@sprintf("%.6f", e) for e in Es], " "), time()-t0)
    end

    fn = joinpath(@__DIR__, "crm_mpslow_$tag.tsv")
    open(fn, "w") do io
        println(io, "W\tLX\tnsites\tgeometry\tU\tmethod\tchi\tquantity\ti\tj\tvalue")
        for (m, χ, q, i, j, v) in rows
            @printf(io, "%d\t%d\t%d\t%s\t%.1f\t%s\t%d\t%s\t%d\t%d\t%.12f\n", Wd, Lx, n, geo, U, m, χ, q, i, j, v)
        end
    end
    @printf("書き出し: %s (%d 行, %.0f 秒)\n", fn, length(rows), time()-t0)
end

for spec in split(ENV["SYSTEMS"], ";"), U in US
    parts = split(spec, ":")
    Wd, geo, Lx = parse(Int, parts[1]), String(parts[2]), parse(Int, parts[3])
    flags = length(parts) >= 4 ? String(parts[4]) : ""
    run_system(Wd, geo, Lx, flags, U)
end
