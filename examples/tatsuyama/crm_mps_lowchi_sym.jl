# 変分 MPS をスピンの向きについて平均した prior(UHF-sym と同じ回転平均)の値を求める。
#
#   変分 MPS(crm_mps_lowchi.jl、states/mpslow_*.jls に保存済み)は電荷の量では良い prior だが、
#   UHF と同じくスピンの向きを選ぶ。全スピンを一様に回して平均した混合状態 σ̄ = ∫dΩ R σ R† を使うと、
#     <S^z_i>        -> 0
#     <S^z_i S^z_j>  -> <S_i·S_j>/3
#     <n_i>, <n_i n_j>, <n_i↑ n_i↓>, <S_i·S_j> は回転で変わらない
#   なので、Z_{i↑} = 1 - n_i - 2 S^z_i から
#     <Z_{i↑}>_sym          = 1 - <n_i>
#     <Z_{i↑} Z_{j↑}>_sym   = <(1-n_i)(1-n_j)> + (4/3) <S_i·S_j>
#     <Z_{i↑} Z_{i↓}>_sym   = <Z_{i↑} Z_{i↓}>   (オンサイトは回転不変)
#   となる。状態ごとに中央のサイト c0 と隣 nb の <n_i>、<n_j>、<n_i n_j>、<S_i·S_j> を出す。
#   (S_i·S_j は S+S- を含むが、S± は同じサイトの c†↑c↓ なので JW の符号は付かない)
#
#   julia --project=Hubbard_MPS_Env_v2 crm_mps_lowchi_sym.jl   (出力: crm_mpslow_sym.tsv)
using ITensors, ITensorMPS, LinearAlgebra, Printf, Serialization

"""crm_mps_lowchi.jl の edges_of と同じ。"""
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

expect_mpo(ψ, os) = real(inner(ψ', MPO(os, siteinds(ψ)), ψ))

dir = joinpath(@__DIR__, "states")
files = sort(filter(f -> startswith(f, "mpslow_") && endswith(f, ".jls"), readdir(dir)))
out = joinpath(@__DIR__, "crm_mpslow_sym.tsv")
open(out, "w") do io
    println(io, "W\tLX\tgeometry\tU\tchi\tc0\tnb\tn_i\tn_j\tnn\tSS\tSz_i\tZZ_onsite\tZupZup_nb")
    for f in files
        m = match(r"mpslow_W(\d+)L(\d+)_(\w+)_U([\d.]+)_chi(\d+)\.jls", f)
        m === nothing && continue
        Wd, Lx, geo, U, χ = parse(Int, m[1]), parse(Int, m[2]), m[3], parse(Float64, m[4]), parse(Int, m[5])
        ψ = deserialize(joinpath(dir, f)).ψ
        edges = edges_of(Lx, Wd, geo == "torus")
        c0 = (max(1, cld(Lx, 2))-1)*Wd + max(1, cld(Wd, 2)); nb = first(b for (a, b) in edges if a == c0)
        q(ops...) = (o = OpSum(); for t in ops; o += t; end; o)
        ni = expect_mpo(ψ, q((1.0, "Ntot", c0)))
        nj = expect_mpo(ψ, q((1.0, "Ntot", nb)))
        nn = expect_mpo(ψ, q((1.0, "Ntot", c0, "Ntot", nb)))
        ss = expect_mpo(ψ, q((1.0, "Sz", c0, "Sz", nb), (0.5, "S+", c0, "S-", nb), (0.5, "S-", c0, "S+", nb)))
        sz = expect_mpo(ψ, q((1.0, "Sz", c0)))
        zon = expect_mpo(ψ, q((4.0, "Nupdn", c0), (-2.0, "Nup", c0), (-2.0, "Ndn", c0), (1.0, "Id", c0)))
        zuu = expect_mpo(ψ, q((4.0, "Nup", c0, "Nup", nb), (-2.0, "Nup", c0), (-2.0, "Nup", nb), (1.0, "Id", c0)))
        @printf(io, "%d\t%d\t%s\t%.1f\t%d\t%d\t%d\t%.12f\t%.12f\t%.12f\t%.12f\t%.12f\t%.12f\t%.12f\n",
                Wd, Lx, geo, U, χ, c0, nb, ni, nj, nn, ss, sz, zon, zuu)
        flush(io)
    end
end
println("書き出し: $out")
