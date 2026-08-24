# ============================================================
# §2 の参照状態 ρ を「厳密対角化の基底状態」にして計算し直す
#
# 動機:
#   §2 の ρ は DMRG(有限 χ)で作った状態だった。主張は「固定された ρ に
#   対する分散比」なので自己無撞着ではあるが、**ρ 自身が切断MPSであること**
#   が結論に効いていないかは別途確かめる価値がある。特に prior も切断MPSなので
#   「同じ MPS 族の中での比較」になっており、有利に働いている可能性がある。
#
# 手法:
#   - 粒子数セクター (N_up, N_dn) 内で、行列を保持しない Lanczos
#     (完全再直交化・thick restart)で基底状態を求める。
#   - 量子ビット占有基底そのもので張るので JW 符号は「間のモード」の占有で決まる:
#       up の i↔i+1: qup(i)=2i-1 と qup(i+1)=2i+1 の間は qdn(i)=2i のみ
#       dn の i↔i+1: qdn(i)=2i と qdn(i+1)=2i+2 の間は qup(i+1)=2i+1 のみ
#   - 全 2^(2L) ベクトルに展開し、窓の縮約密度行列は reshape で厳密に取る。
#   - prior は DMRG(χ切断)と UHF。ρ(ED)と ρ(DMRG)の両方で同じ量を測る。
#
# 規模: セクター次元 C(L,L/2)^2、全ベクトル 2^(2L)*16 バイト
#   L=12: 8.5e5 / 0.27 GB   L=14: 1.2e7 / 4.3 GB   L=16: 1.7e8 / 69 GB(不可)
# ============================================================
include(joinpath(@__DIR__, "crm_chain_common.jl"))

function combos(L, k)
    out = Int[]
    for x in 0:(1<<L - 1); count_ones(x) == k && push!(out, x); end
    return out
end

struct Sector
    L::Int; ups::Vector{Int}; dns::Vector{Int}
    upidx::Dict{Int,Int}; dnidx::Dict{Int,Int}
end
function Sector(L, Nup, Ndn)
    ups = combos(L, Nup); dns = combos(L, Ndn)
    Sector(L, ups, dns, Dict(a=>i for (i,a) in enumerate(ups)),
           Dict(b=>i for (i,b) in enumerate(dns)))
end
Base.length(s::Sector) = length(s.ups)*length(s.dns)
bit(x, i) = (x >> (i-1)) & 1

"""ホッピング先を前計算しておく(内側ループから Dict 参照を除く)。"""
struct HopTables
    up::Vector{Vector{Tuple{Int,Int}}}   # ia -> [(bond i, 行き先 ia2)]
    dn::Vector{Vector{Tuple{Int,Int}}}
end
function HopTables(S::Sector)
    L = S.L
    mk(cfgs, idx) = [begin
        h = Tuple{Int,Int}[]
        for i in 1:L-1
            bit(x,i) == bit(x,i+1) && continue
            push!(h, (i, idx[x ⊻ ((1<<(i-1)) | (1<<i))]))
        end
        h
    end for x in cfgs]
    HopTables(mk(S.ups, S.upidx), mk(S.dns, S.dnidx))
end

function apply_H!(w, v, S::Sector, HT::HopTables, t, U, mu)
    nu = length(S.ups); nd = length(S.dns)
    fill!(w, 0)
    @inbounds for ia in 1:nu
        a = S.ups[ia]; base = (ia-1)*nd
        na = count_ones(a)
        # 対角
        for ib in 1:nd
            b = S.dns[ib]
            w[base+ib] += (U*count_ones(a & b) - mu*(na + count_ones(b))) * v[base+ib]
        end
        # up ホッピング: 行き先は ib について連続
        for (i, ia2) in HT.up[ia]
            tb = (ia2-1)*nd
            for ib in 1:nd
                sgn = iseven(bit(S.dns[ib], i)) ? 1.0 : -1.0
                w[tb+ib] += -t*sgn*v[base+ib]
            end
        end
        # dn ホッピング: 符号は a_{i+1} で決まるので ia 固定なら定数
        for ib in 1:nd
            amp = v[base+ib]
            amp == 0 && continue
            for (i, ib2) in HT.dn[ib]
                sgn = iseven(bit(a, i+1)) ? 1.0 : -1.0
                w[base+ib2] += -t*sgn*amp
            end
        end
    end
    return w
end

function lanczos_gs(S::Sector, t, U, mu; m=50, maxrestart=60, tol=1e-11, seed=1234)
    D = length(S); HT = HopTables(S); Random.seed!(seed)
    v = randn(D); v ./= norm(v); Eprev = Inf; w = similar(v)
    for restart in 1:maxrestart
        Q = [copy(v)]; α = Float64[]; β = Float64[]
        for j in 1:m
            apply_H!(w, Q[j], S, HT, t, U, mu)
            a = dot(Q[j], w); push!(α, a)
            w .-= a .* Q[j]
            j > 1 && (w .-= β[j-1] .* Q[j-1])
            for q in Q; w .-= dot(q, w) .* q; end
            b = norm(w)
            if b < 1e-12 || j == m; push!(β, b); break; end
            push!(β, b); push!(Q, w ./ b)
        end
        k = length(α)
        F = eigen(SymTridiagonal(α, β[1:k-1]))
        E0 = F.values[1]; c = F.vectors[:, 1]
        v = zeros(D); for j in 1:k; v .+= c[j] .* Q[j]; end
        v ./= norm(v)
        abs(E0 - Eprev) < tol && return E0, v, restart
        Eprev = E0
    end
    return Eprev, v, maxrestart
end

function expand_full(v, S::Sector)
    L = S.L; nd = length(S.dns)
    full = zeros(ComplexF64, 1 << (2L))
    @inbounds for (ia, a) in enumerate(S.ups), (ib, b) in enumerate(S.dns)
        q = 0
        for i in 1:L
            q |= bit(a,i) << (2i-2)
            q |= bit(b,i) << (2i-1)
        end
        full[q+1] = v[(ia-1)*nd + ib]
    end
    return full
end

"""全ベクトルから窓 [q1,q2] の縮約密度行列。窓内は b_{q1} が最上位になるよう並べ替える。"""
function window_rdm_full(full, N, q1, q2)
    Nw = q2 - q1 + 1
    dl = 1 << (q1-1); dw = 1 << Nw; dr = 1 << (N - q2)
    A = reshape(full, dl, dw, dr)
    ρ = zeros(ComplexF64, dw, dw)
    @inbounds for r in 1:dr, s2 in 1:dw, s1 in 1:dw
        acc = zero(ComplexF64)
        for l in 1:dl; acc += A[l, s1, r] * conj(A[l, s2, r]); end
        ρ[s1, s2] += acc
    end
    perm = [begin
                x = s - 1; y = 0
                for k in 0:Nw-1; y |= ((x >> k) & 1) << (Nw-1-k); end
                y + 1
            end for s in 1:dw]
    R = ρ[perm, perm]
    return Hermitian((R + R')/2)
end

# ---------------- 1D鎖 UHF(crm_1d_hf.jl と同一) ----------------
chain_hop_matrix(L, t) = (T = zeros(L, L);
    for i in 1:L-1; T[i,i+1] = T[i+1,i] = -t; end; T)

function solve_uhf_chain(L, t, U; Nup, Ndn, iters=4000, tol=1e-12, init=:neel)
    T = chain_hop_matrix(L, t); fill_avg = (Nup+Ndn)/(2L)
    nup = zeros(L); ndn = zeros(L)
    for i in 1:L
        if init == :neel
            nup[i] = fill_avg + 0.4*(-1)^i; ndn[i] = fill_avg - 0.4*(-1)^i
        else
            nup[i] = fill_avg + 0.3*(rand()-0.5); ndn[i] = fill_avg + 0.3*(rand()-0.5)
        end
    end
    clamp!(nup, 0.02, 0.98); clamp!(ndn, 0.02, 0.98)
    nup .*= Nup/sum(nup); ndn .*= Ndn/sum(ndn)
    local Fu, Fd
    for _ in 1:iters
        Fu = eigen(Symmetric(T + diagm(U .* ndn))); Fd = eigen(Symmetric(T + diagm(U .* nup)))
        nu2 = vec(sum(abs2, Fu.vectors[:, 1:Nup]; dims=2))
        nd2 = vec(sum(abs2, Fd.vectors[:, 1:Ndn]; dims=2))
        if max(maximum(abs.(nu2 .- nup)), maximum(abs.(nd2 .- ndn))) < tol
            nup, ndn = nu2, nd2; break
        end
        nup = 0.5.*nu2 .+ 0.5.*nup; ndn = 0.5.*nd2 .+ 0.5.*ndn
    end
    Φu = Fu.vectors[:, 1:Nup]; Φd = Fd.vectors[:, 1:Ndn]
    Cu = Φu*Φu'; Cd = Φd*Φd'
    return (Cu=Cu, Cd=Cd, E=tr(T*Cu)+tr(T*Cd)+U*sum(diag(Cu).*diag(Cd)),
            m=mean((-1)^i*(Cu[i,i]-Cd[i,i])/2 for i in 1:L))
end

full_C(uhf, L) = (C = zeros(2L, 2L);
    for a in 1:L, b in 1:L
        C[qup(a), qup(b)] = uhf.Cu[a,b]; C[qdn(a), qdn(b)] = uhf.Cd[a,b]
    end; C)

function rotate_C(C0, L, θ)
    c, s = cos(θ/2), sin(θ/2); Um = zeros(2L, 2L)
    for i in 1:L
        a, b = qup(i), qdn(i)
        Um[a,a]=c; Um[a,b]=-s; Um[b,a]=s; Um[b,b]=c
    end
    return Um*C0*Um'
end

function wick_gen(sup, C)
    isempty(sup) && return 1.0
    ps = [a for (_,a) in sup]
    if all(==(3), ps)
        qs = [q for (q,_) in sup]
        length(qs)==1 && return 1 - 2C[qs[1],qs[1]]
        a,b = qs
        return 1 - 2C[a,a] - 2C[b,b] + 4*(C[a,a]*C[b,b] - C[a,b]*C[b,a])
    elseif length(sup)==3 && (ps==[1,3,1] || ps==[2,3,2])
        return 2*C[sup[1][1], sup[3][1]]
    end
    error("wick_gen: 未対応 $(sup)")
end

symavg(obs, C0, L; nθ=64) = begin
    P = [[0.0 for _ in o.terms] for o in obs]
    for k in 1:nθ
        θ = acos(-1 + (k-0.5)*2/nθ); C = rotate_C(C0, L, θ)
        for (i,o) in enumerate(obs), (j,tm) in enumerate(o.terms)
            P[i][j] += wick_gen(tm.sup, C)/nθ
        end
    end
    P
end

# ---------------- メイン ----------------
function main()
    t, U, nu, nm = 1.0, 4.0, 50, 100
    L_list = parse.(Int, split(get(ENV, "L_LIST", "4,6,8,10,12"), ","))
    n_repeat = parse(Int, get(ENV, "N_REPEAT", "500"))
    chi_exp = parse(Int, get(ENV, "CHI_EXP", "256"))
    chi_priors = [2, 4, 8, 16, 32]

    println("=== validation (L=4, DMRG系) ==="); flush(stdout)
    dense_chain_check()

    rows = []
    for L in L_list
        N = 2L; Nup = L ÷ 2; Ndn = L - Nup
        S = Sector(L, Nup, Ndn)
        @printf("\n%s\nL=%d (N=%d qubits)  セクター次元=%d  全ベクトル=%.2f GB\n",
                "="^72, L, N, length(S), (1<<N)*16/1e9); flush(stdout)

        t0 = time()
        E_ed, vsec, nres = lanczos_gs(S, t, U, U/2)
        @printf("  ED(Lanczos): E0=%.10f  restart=%d  (%.0fs)\n", E_ed, nres, time()-t0)
        flush(stdout)
        full = expand_full(vsec, S)

        E_dm, ψ, _, _, _ = ground_state(L, t, U, U/2; chi_max=chi_exp, nsweeps=20)
        @printf("  DMRG(χ=%d): E0=%.10f  maxlinkdim=%d  |ΔE|=%.2e\n",
                chi_exp, E_dm, maxlinkdim(ψ), abs(E_dm-E_ed)); flush(stdout)

        uhf = solve_uhf_chain(L, t, U; Nup, Ndn); C0 = full_C(uhf, L)

        sites = 1:L-1
        obs_by_site = [site_observables(i; bond=true) for i in sites]
        windows = [obs_support(o) for o in obs_by_site]
        rdm_dm = window_rdms(ψ, windows)
        priors = MPS[]
        for cp in chi_priors
            σ = truncate(ψ; maxdim=cp); normalize!(σ); push!(priors, σ)
        end
        push!(priors, ψ)
        rdm_p = [window_rdms(σ, windows) for σ in priors]
        labels = vcat(["chi$c" for c in chi_priors], ["exact", "UHF", "UHFsym"])

        maxdiff = 0.0
        for (si, i) in enumerate(sites)
            obs = obs_by_site[si]; q1, q2 = windows[si]; Nw = q2-q1+1
            ow = shift_obs(obs, q1-1)
            Mats = [[term_window_matrix(tm, Nw) for tm in o.terms] for o in ow]
            ρ_ed = window_rdm_full(full, N, q1, q2)
            for (name, ρref, rdmp) in (("ED", ρ_ed, nothing), ("DMRG", rdm_dm[si], nothing))
                Otrue = [sum(tm.coeff*expect_rdm(ρref, Mats[k][ti])
                             for (ti,tm) in enumerate(ow[k].terms)) for k in 1:length(obs)]
                if name == "ED"
                    Odm = [sum(tm.coeff*expect_rdm(rdm_dm[si], Mats[k][ti])
                               for (ti,tm) in enumerate(ow[k].terms)) for k in 1:length(obs)]
                    maxdiff = max(maxdiff, maximum(abs.(Otrue .- Odm)))
                end
                Pσ_all = Vector{Vector{Vector{Float64}}}(); trOσ_all = Vector{Vector{Float64}}()
                for p in 1:length(priors)
                    Pσ = [[expect_rdm(rdm_p[p][si], M) for M in Ms] for Ms in Mats]
                    push!(Pσ_all, Pσ)
                    push!(trOσ_all, [sum(tm.coeff*Pσ[k][ti] for (ti,tm) in enumerate(obs[k].terms))
                                     for k in 1:length(obs)])
                end
                Pu = [[wick_gen(tm.sup, C0) for tm in o.terms] for o in obs]
                Ps = symavg(obs, C0, L)
                for P in (Pu, Ps)
                    push!(Pσ_all, P)
                    push!(trOσ_all, [sum(tm.coeff*P[k][ti] for (ti,tm) in enumerate(obs[k].terms))
                                     for k in 1:length(obs)])
                end
                es, ec = run_locals(ρref, Nw, ow, Pσ_all, trOσ_all;
                                    nu, nm, n_repeat, seed = 9100 + 100L + i)
                for (k,o) in enumerate(obs)
                    v = var(es[:,k]); gmx = v/var(ec[:,k,length(chi_priors)+1])
                    for (p, lb) in enumerate(labels)
                        push!(rows, (L, name, i, lb, o.name, Otrue[k],
                                     Otrue[k]-trOσ_all[p][k], v/var(ec[:,k,p]), gmx))
                    end
                end
            end
        end
        @printf("  局所観測量の ED vs DMRG 最大差: %.3e\n", maxdiff); flush(stdout)

        @printf("\n  %-14s %-6s %8s %8s %8s %8s %8s | %8s\n",
                "observable","ρ","χ=2","χ=4","χ=8","χ=32","UHFsym","天井")
        for o in unique(r[5] for r in rows if r[1]==L)
            for nmr in ("ED","DMRG")
                f(lb)=(h=[r[8] for r in rows if r[1]==L&&r[2]==nmr&&r[4]==lb&&r[5]==o];
                       isempty(h) ? NaN : median(h))
                gm=median([r[9] for r in rows if r[1]==L&&r[2]==nmr&&r[5]==o])
                @printf("  %-14s %-6s %8.2f %8.2f %8.2f %8.2f %8.2f | %8.2f\n",
                        nmr=="ED" ? o : "", nmr, f("chi2"),f("chi4"),f("chi8"),
                        f("chi32"),f("UHFsym"),gm)
            end
        end
        flush(stdout)
    end

    out = joinpath(@__DIR__, "crm_ed_reference_results.tsv")
    open(out,"w") do io
        println(io, "L\trho\tsite\tprior\tobservable\ttrue\tDelta\tG\tG_max")
        for r in rows; println(io, join(r,"\t")); end
    end
    println("\nresults saved: $out")
end
main()
