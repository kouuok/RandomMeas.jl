# 仮説「prior が観測量の固有状態に近いほど利得が大きい」を、距離 r の Z 相関で連続的に検証する。
#
# 総当たり(crm_2d_ed_fid_table.jl)は中央の1結合しか測っていない。ここでは厳密な ρ と各 prior
# について全サイト対の <n_iσ n_jσ'> と <n_iσ> を相関行列で一度に取り、Z = 1-2n から
#   Zup     : Z_{i↑}                  (台 1 量子ビット)
#   ZupZup  : Z_{i↑} Z_{j↑}   (i<j)   (台 2 量子ビット)
#   ZupZdn  : Z_{i↑} Z_{j↓}   (全 i,j。i=j はオンサイト)
# を組み立てる。UHF は長距離秩序を持つので、どの距離でも固有状態の近くに留まる。一方で
# 厳密な ρ の相関は距離とともに減衰する。つまり「prior は固有状態に近いまま、ρ だけが
# 離れていく」系列が1つの系の中で作れる。
#
# 厳密な ρ と UHF の MPS は states/ に保存する(別の観測量を足すとき DMRG をやり直さなくてよい)。
#   W=1 LX=8 PBC=0 CRM_U=8.0 julia --project=. examples/tatsuyama/crm_eigen_corr.jl
include("crm_2d_ed_fid_table.jl")   # 定数・lat_edges・solve_uhf・slater_mps(main は走らない)
using Serialization
# REFINE=1: 保存済みの状態の分散が VARTOL を超えていたら、そこからノイズ付きでスイープを足す。
#   ランダム初期状態からの DMRG は L=64 で局所解に留まることがある(U=12 で分散 4.3e-6、χ=485)。
const REFINE = get(ENV, "REFINE", "0") == "1"
const VARTOL = parse(Float64, get(ENV, "VARTOL", "1e-8"))

function hubbard_mpo(sites, edges, n)
    os = OpSum()
    for (a,b) in edges
        os += -1.0,"Cdagup",a,"Cup",b; os += -1.0,"Cdagup",b,"Cup",a
        os += -1.0,"Cdagdn",a,"Cdn",b; os += -1.0,"Cdagdn",b,"Cdn",a
    end
    for s in 1:n; os += U,"Nupdn",s end
    MPO(os, sites)
end
energy_var(H, ψ, E) = (Hψ = apply(H, ψ; cutoff=1e-16); real(inner(Hψ,Hψ) - E^2))

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

function run()
    t0 = time()
    n = LX*W; Nup = n÷2; Ndn = n÷2
    edges = lat_edges(LX, W; pbc_x=PBC); geo = PBC ? "torus" : "cylinder"
    chimax = CHIMAX > 0 ? CHIMAX : min(n÷2 >= 32 ? typemax(Int) : 4^(n÷2), 8192)
    @printf("%s %dx%d  n=%d  U=%.1f  二部格子:%s  χ上限 %d\n", geo, LX, W, n, U,
            bipartite(edges, n) ? "はい" : "いいえ", chimax)
    tag = @sprintf("W%dL%d_%s_U%.1f", W, LX, geo, U)
    mkpath(joinpath(@__DIR__, "states"))
    fstate = joinpath(@__DIR__, "states", "eigcorr_$tag.jls")

    if isfile(fstate)
        d = deserialize(fstate); ψ, ψu, E, var = d.ψ, d.ψu, d.E, d.var
        @printf("保存済みの状態を読み込み: %s  (E0=%.12f 分散=%.3e χ=%d)\n", fstate, E, var, maxlinkdim(ψ))
        if REFINE && var > VARTOL
            H = hubbard_mpo(siteinds(ψ), edges, n)
            for round in 1:8
                E, ψ = dmrg(H, ψ; nsweeps=20, maxdim=chimax, cutoff=CUTOFF,
                            noise=[1e-5,1e-6,1e-7,1e-8,1e-9,1e-10,0.0], outputlevel=0)
                normalize!(ψ); var = energy_var(H, ψ, E)
                @printf("  追加スイープ %d: E0=%.12f  χ=%d  分散=%.3e  (%.0f 秒)\n", round, E, maxlinkdim(ψ), var, time()-t0)
                var < VARTOL && break
            end
            serialize(fstate, (ψ=ψ, ψu=ψu, E=E, var=var, m=d.m))
        end
    else
        sites = siteinds("Electron", n; conserve_qns=true)
        H = hubbard_mpo(sites, edges, n)
        # DMRG の設定は crm_2d_ed_fid_table.jl の main() と同一(ランプ + ノイズ項)
        st = [isodd(sum(divrem(s-1,W))) ? "Up" : "Dn" for s in 1:n]
        sched = Int[]; cc = 64
        while cc < chimax; push!(sched, cc); cc *= 2 end
        push!(sched, chimax)
        md = vcat([fill(d, SWPER) for d in sched]...)
        E, ψ = dmrg(H, random_mps(sites, st; linkdims=32);
                    nsweeps=length(md)+30, maxdim=vcat(md, fill(chimax,30)),
                    cutoff=CUTOFF, noise=[1e-6,1e-7,1e-8,1e-9,0.0], outputlevel=0)
        normalize!(ψ); var = energy_var(H, ψ, E)
        uhf = solve_uhf(edges, n, 1.0, U; Nup, Ndn)
        ψu = slater_mps(sites, uhf.Φu, uhf.Φd, Nup, Ndn)
        serialize(fstate, (ψ=ψ, ψu=ψu, E=E, var=var, m=uhf.m))
        @printf("UHF: 磁化 m=%.4f  χ=%d\n", uhf.m, maxlinkdim(ψu))
    end
    @printf("E0=%.12f  到達χ=%d  <H^2>-<H>^2=%.3e  (%.0f 秒)\n", E, maxlinkdim(ψ), var, time()-t0)

    priors = Any[ψu]; labels = ["UHF"]
    for c in (2,4,8,16,32,64,128,256,512)
        c < maxlinkdim(ψ) || continue
        σ = copy(ψ); orthogonalize!(σ,1); truncate!(σ; maxdim=c, cutoff=0.0); normalize!(σ)
        push!(priors, σ); push!(labels, "chi$c")
    end

    zr = zobs(ψ); @printf("ρ の相関行列 (%.0f 秒)\n", time()-t0)
    coord(s) = divrem(s-1, W) .+ (1, 1)          # サイト番号 → (x, y)
    fn = joinpath(@__DIR__, "crm_eigcorr_$tag.tsv")
    open(fn, "w") do io
        println(io, "W\tLX\tnsites\tgeometry\tU\tEref\tHvar\tchi_ref\tprior\tkind\ti\tj\txi\tyi\txj\tyj\ttrue\tprior_val")
        for (l, σ) in zip(labels, priors)
            zs = zobs(σ)
            for k in sort(collect(keys(zr)))
                (kind, i, j) = k; (xi, yi) = coord(i); (xj, yj) = coord(j)
                @printf(io, "%d\t%d\t%d\t%s\t%.1f\t%.12f\t%.3e\t%d\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%.12f\t%.12f\n",
                        W, LX, n, geo, U, E, var, maxlinkdim(ψ), l, kind, i, j, xi, yi, xj, yj, zr[k], zs[k])
            end
            @printf("  %-6s 完了 (%.0f 秒)\n", l, time()-t0)
        end
    end
    println("書き出し: ", fn)
end
run()
