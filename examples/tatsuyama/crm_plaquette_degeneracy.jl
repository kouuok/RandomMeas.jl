# ============================================================
# §1 の4サイトリングは、ホッピングを一様にすると本当に縮退するのか
#
# 背景: crm_gain_verification.jl は 2--3 の結合だけを t=0.5 に弱めており、
#   コメントは「二量体化で縮退を回避」としていた。しかし4サイトリングは
#   二部格子(|A|=|B|=2)なので、Lieb の定理から半充填・U>0 の基底状態は
#   S=0 の一意な singlet であり、一様でも縮退しないはずである。直接確かめる。
#
# 結果: 縮退するのは U=0 の開殻(1粒子準位 -2t,0,0,+2t で各スピン2電子)だけ。
#   U=1, 8 では一様でもギャップは 4.8e-2, 3.3e-1 と有限。U=8 では一様の方が
#   むしろ広い(0.332 対 0.287)。UHF も一様リングで正常に収束する。
#   → 弱化は必要ではなく、実効果は U=1 でギャップを約8倍にすることのみ。
#
# 実行: julia --project=Hubbard_MPS_Env_v2 crm_plaquette_degeneracy.jl
# ============================================================
using LinearAlgebra, Printf
const NQ=8; const DIM=256
make_matrix(ops...) = reduce(kron, ops)
function get_creation_ops()
    cdag=[0.0 0.0; 1.0 0.0]; I2=[1.0 0.0;0.0 1.0]; F=[1.0 0.0;0.0 -1.0]
    Cdag=Dict{String,Matrix{Float64}}()
    labels=["1u","1d","2u","2d","3u","3d","4u","4d"]
    for (k,lab) in enumerate(labels)
        ops=Matrix{Float64}[]
        for pos in 1:NQ
            q=NQ-pos+1
            q<k ? push!(ops,I2) : (q==k ? push!(ops,cdag) : push!(ops,F))
        end
        Cdag[lab]=make_matrix(ops...)
    end
    Cdag
end
function build_H(t_L,t_S,U,mu)
    Cdag=get_creation_ops(); C=Dict(k=>Matrix(v') for (k,v) in Cdag)
    Nop=Dict(k=>Cdag[k]*C[k] for k in keys(Cdag))
    H=zeros(DIM,DIM)
    for s in ("u","d")
        for (i,j) in [(1,2),(3,4),(4,1)]
            H .+= -t_L .* (Cdag["$i$s"]*C["$j$s"] .+ Cdag["$j$s"]*C["$i$s"])
        end
        for (i,j) in [(2,3)]
            H .+= -t_S .* (Cdag["$i$s"]*C["$j$s"] .+ Cdag["$j$s"]*C["$i$s"])
        end
    end
    for i in 1:4; H .+= U .* (Nop["$(i)u"]*Nop["$(i)d"]); end
    for op in values(Nop); H .+= -mu .* op; end
    H, Nop
end

println("全ヒルベルト空間 (2^8=256次元、μ=U/2) での最低5準位")
@printf("%-8s %-8s %12s %12s %12s %12s %12s %10s %10s\n",
        "t_S","U","E0","E1","E2","E3","E4","E1-E0","<N>_gs")
for U in (0.0, 1.0, 8.0), tS in (1.0, 0.5)
    H, Nop = build_H(1.0, tS, U, U/2)
    F = eigen(Symmetric(H)); ev = F.values
    Ntot = sum(values(Nop))
    ngs = real(dot(F.vectors[:,1], Ntot*F.vectors[:,1]))
    @printf("%-8.1f %-8.1f %12.6f %12.6f %12.6f %12.6f %12.6f %10.2e %10.3f\n",
            tS, U, ev[1],ev[2],ev[3],ev[4],ev[5], ev[2]-ev[1], ngs)
end

println("\n半充填セクター (N=4, Nup=Ndn=2) に射影した最低4準位")
@printf("%-8s %-8s %12s %12s %12s %10s\n","t_S","U","E0","E1","E2","E1-E0")
for U in (0.0, 1.0, 8.0), tS in (1.0, 0.5)
    H, Nop = build_H(1.0, tS, U, U/2)
    Nup = Nop["1u"]+Nop["2u"]+Nop["3u"]+Nop["4u"]
    Ndn = Nop["1d"]+Nop["2d"]+Nop["3d"]+Nop["4d"]
    idx = [i for i in 1:DIM if abs(Nup[i,i]-2)<1e-9 && abs(Ndn[i,i]-2)<1e-9]
    Hs = Symmetric(H[idx,idx]); ev = eigen(Hs).values
    @printf("%-8.1f %-8.1f %12.6f %12.6f %12.6f %10.2e\n", tS,U,ev[1],ev[2],ev[3],ev[2]-ev[1])
end

# --- UHF 側: 一様リングだと自己無撞着解が縮退/不定になるか ---
function hf_run(t_L,t_S,U,mu,Nup,Ndn)
    T=zeros(4,4)
    for (i,j) in [(1,2),(3,4),(4,1)]; T[i,j]=T[j,i]=-t_L; end
    for (i,j) in [(2,3)]; T[i,j]=T[j,i]=-t_S; end
    n_up=zeros(4); n_dn=zeros(4)
    for i in 1:Nup; n_up[i]=0.8; end
    for i in 1:Ndn; n_dn[i]=0.8; end
    n_up[1]+=0.1; n_dn[1]-=0.1
    conv=false; it=0
    local Fu,Fd
    for k in 1:500
        Fu=eigen(Symmetric(T+diagm(U.*n_dn)-mu*I)); Fd=eigen(Symmetric(T+diagm(U.*n_up)-mu*I))
        nu2=sum(Fu.vectors[:,i].^2 for i in 1:Nup; init=zeros(4))
        nd2=sum(Fd.vectors[:,i].^2 for i in 1:Ndn; init=zeros(4))
        if maximum(abs.(nu2.-n_up))<1e-10; conv=true; it=k; n_up,n_dn=nu2,nd2; break; end
        n_up=0.5.*nu2.+0.5.*n_up; n_dn=0.5.*nd2.+0.5.*n_dn
    end
    # HOMO-LUMO ギャップ (占有軌道の選択が一意か)
    gu = Fu.values[Nup+1]-Fu.values[Nup]
    (conv=conv, it=it, gap=gu, m=maximum(abs.(n_up.-n_dn)))
end
println("\nUHF: 一様リング vs 1本弱化 (Nup=Ndn=2)")
@printf("%-8s %-6s %8s %6s %14s %12s\n","t_S","U","収束","反復","HOMO-LUMOギャップ","max|n↑-n↓|")
for U in (1.0,8.0), tS in (1.0,0.5)
    r=hf_run(1.0,tS,U,U/2,2,2)
    @printf("%-8.1f %-6.1f %8s %6d %14.3e %12.4f\n",tS,U,r.conv ? "した" : "しない",r.it,r.gap,r.m)
end
