# ============================================================
# 帰無検定: サイト依存性は本当に「開放端の効果」か
#
# §2b/§2c では、利得 G がサイトによって変わり(端で1.3-3.6倍、ボンド量では
# 隣り合うボンドで最大17倍)、その原因は prior 品質ではなく観測量自身の
# <P>(i) が開放端由来の Friedel 振動で変わるためだと結論した。
#
# この説明が正しければ、**周期境界条件にすると並進対称性から全サイトが
# 等価になり、サイト依存性は完全に消えるはず**である。消えなければ説明が誤り。
#
# G_max = 1 + n_m (3^|A|-1)<P>^2 / [3^|A|(1-<P>^2)] は <P> だけの関数なので、
# <P>(i) のサイト依存を見れば十分で、Monte Carlo は不要。
#
# 実行: julia --project=Hubbard_MPS_Env_v2 crm_pbc_check.jl
# 環境変数: L_LIST(既定 "16,32"), CHI_EXP(既定 512), NSWEEPS(既定 30)
# ============================================================
include(joinpath(@__DIR__, "crm_chain_common.jl"))

"""周期境界(pbc=true)にも対応した1Dハバード鎖のMPO。"""
function chain_mpo_bc(sites, L, t, U, mu; pbc::Bool)
    os = OpSum()
    bonds = pbc ? [(i, mod1(i+1, L)) for i in 1:L] : [(i, i+1) for i in 1:L-1]
    for (i, j) in bonds
        for q in (qup, qdn)
            os += -t, "Cdag", q(i), "C", q(j)
            os += -t, "Cdag", q(j), "C", q(i)
        end
    end
    for i in 1:L
        os += U, "N", qup(i), "N", qdn(i)
        os += -mu, "N", qup(i); os += -mu, "N", qdn(i)
    end
    return MPO(os, sites)
end

function gs_bc(L, t, U, mu; chi_max, nsweeps, pbc)
    sites = siteinds("Fermion", 2L; conserve_qns=true)
    H = chain_mpo_bc(sites, L, t, U, mu; pbc)
    init, nel = initial_config(L, 0.0)
    ψ0 = productMPS(sites, init)
    ramp = vcat([50, 100, 200, 400], fill(chi_max, max(0, nsweeps-4)))
    noise = vcat([1e-5,1e-6,1e-7,1e-8,1e-9], zeros(max(0, nsweeps-5)))[1:nsweeps]
    E, ψ = dmrg(H, ψ0; nsweeps, maxdim=min.(chi_max, ramp)[1:nsweeps],
                cutoff=1e-11, noise, outputlevel=0)
    return E, ψ, nel
end

# G_max（単一Pauli列の閉形式）
gmax(P, absA, nm) = 1 + nm*(3.0^absA - 1)*P^2 / (3.0^absA*(1-P^2))

function main()
    t, U, nm = 1.0, 4.0, 100
    L_list = parse.(Int, split(get(ENV, "L_LIST", "16,32"), ","))
    chi_exp = parse(Int, get(ENV, "CHI_EXP", "512"))
    nsweeps = parse(Int, get(ENV, "NSWEEPS", "30"))

    for L in L_list, pbc in (false, true)
        E, ψ, nel = gs_bc(L, t, U, U/2; chi_max=chi_exp, nsweeps, pbc)
        @printf("\n=== L=%d  %s : E0=%.6f  N_el=%d  chi=%d ===\n",
                L, pbc ? "PBC" : "OBC", E, nel, maxlinkdim(ψ)); flush(stdout)

        # 周期境界でも窓が連続になるサイトのみ（i=1..L-1 のボンド）
        obs_by_site = [site_observables(i; bond=true) for i in 1:L-1]
        windows = [obs_support(o) for o in obs_by_site]
        rdms = window_rdms(ψ, windows)
        res = Dict{String,Vector{Float64}}()
        for i in 1:L-1
            q1, q2 = windows[i]; Nw = q2-q1+1
            for (o, ow) in zip(obs_by_site[i], shift_obs(obs_by_site[i], q1-1))
                v = sum(tm.coeff * expect_rdm(rdms[i], term_window_matrix(tm, Nw))
                        for tm in ow.terms)
                push!(get!(res, o.name, Float64[]), v)
            end
        end

        for nmn in ("ZZ onsite", "ZZ up-up r=1", "DoubleOcc")
            v = res[nmn]
            absA = nmn == "DoubleOcc" ? 2 : 2
            gs = [gmax(x, absA, nm) for x in v]
            @printf("  %-14s <P>: min=%+.4f max=%+.4f  変動係数=%.4f | G_max: min=%.1f max=%.1f 比=%.2f\n",
                    nmn, minimum(v), maximum(v), std(v)/abs(mean(v)),
                    minimum(gs), maximum(gs), maximum(gs)/minimum(gs))
        end
        # 先頭8ボンドの <P> を並べる（交替が見えるか）
        @printf("  ZZ up-up r=1 の <P>(i), i=1..8: %s\n",
                join([@sprintf("%+.4f", x) for x in res["ZZ up-up r=1"][1:min(8,end)]], " "))
        flush(stdout)
    end
end

main()
