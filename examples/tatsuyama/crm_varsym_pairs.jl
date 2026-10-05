# 変分 MPS の全サイト対について、回転平均した prior と「連結相関の prior」に要る量を求める。
#
#   回転平均した prior(crm_mps_lowchi_sym.jl と同じ):
#     <Z_{i↑}Z_{j↑}>_sym = <(1-n_i)(1-n_j)> + (4/3)<S_i·S_j>
#   スピンの向きを選んだ状態の <S_i·S_j> には、距離によらない <S_i>·<S_j>(偽の長距離秩序)が含まれ、
#   S_i·S_j は回転で変わらないので回転平均でも残る。そこで、それを引いた「連結相関の prior」
#     <Z_{i↑}Z_{j↑}>_conn = <(1-n_i)(1-n_j)> + (4/3)(<S_i·S_j> - <S_i>·<S_j>)
#   も比べる。CRM はどんな数値を prior に使っても偏らないので、これが物理的な状態の期待値でなくてもよい。
#   変分 MPS は S^z の総和を保つので <S^x_i> = <S^y_i> = 0 で、<S_i>·<S_j> = <S^z_i><S^z_j>。
#
#   対象: 1 次元 L=64 開放端と L=16 周期境界(厳密な状態の全対は crm_eigcorr_*.tsv にある)
#   julia --project=Hubbard_MPS_Env_v2 crm_varsym_pairs.jl   (出力: crm_varsym_pairs.tsv)
using ITensors, ITensorMPS, LinearAlgebra, Printf, Serialization

dir = joinpath(@__DIR__, "states")
files = sort(filter(f -> occursin(r"^mpslow_W1L(64_cylinder|16_torus)_U[\d.]+_chi\d+\.jls$", f), readdir(dir)))
out = joinpath(@__DIR__, "crm_varsym_pairs.tsv")
open(out, "w") do io
    println(io, "LX\tgeometry\tU\tchi\ti\tj\tn_i\tn_j\tnn\tSS\tSz_i\tSz_j\tZupZup")
    for f in files
        m = match(r"mpslow_W1L(\d+)_(\w+)_U([\d.]+)_chi(\d+)\.jls", f)
        Lx, geo, U, χ = parse(Int, m[1]), m[2], parse(Float64, m[3]), parse(Int, m[4])
        ψ = deserialize(joinpath(dir, f)).ψ
        n  = expect(ψ, "Ntot"); sz = expect(ψ, "Sz"); nup = expect(ψ, "Nup")
        Cnn = correlation_matrix(ψ, "Ntot", "Ntot")
        Czz = correlation_matrix(ψ, "Sz", "Sz")
        Cpm = correlation_matrix(ψ, "S+", "S-")
        Cmp = correlation_matrix(ψ, "S-", "S+")
        Cuu = correlation_matrix(ψ, "Nup", "Nup")
        N = length(ψ)
        for i in 1:N, j in i+1:N
            ss = real(Czz[i, j]) + 0.5*real(Cpm[i, j] + Cmp[i, j])
            zuu = 1 - 2nup[i] - 2nup[j] + 4real(Cuu[i, j])        # 照合用(crm_mpslow_*.tsv の ZupZup と同じはず)
            @printf(io, "%d\t%s\t%.1f\t%d\t%d\t%d\t%.12f\t%.12f\t%.12f\t%.12f\t%.12f\t%.12f\t%.12f\n",
                    Lx, geo, U, χ, i, j, n[i], n[j], real(Cnn[i, j]), ss, sz[i], sz[j], zuu)
        end
        println("done: $f"); flush(stdout)
    end
end
println("書き出し: $out")
