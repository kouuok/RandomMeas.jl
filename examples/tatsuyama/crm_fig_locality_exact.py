"""修士論文 第4章の図(局所性)を、分散の厳密な式の値で描く。

  (a) prior(切断 MPS)の大域的な忠実度 F と系の長さ L
  (b) 利得 G と L(厳密な値)
  (c) 利得 G と 1 設定あたりのショット数 n_m(L=8、U=4)。線は厳密な式、点は反復の少ないモンテカルロ
  データ: crm_mps_scaling_results.tsv(G_exact 列)、crm_sweep_nm.tsv(G_exact、G_emp 列)
  出力: thesis/figures/fig_locality.png

  忠実度と利得は尺度の違う量なので、同じ図に 2 本の縦軸で重ねず、別の図にする。
  配色は dataviz の既定パレットの slot 1〜4 を固定順で使い、線種とマーカーでも系列を区別する。
"""
import csv, collections
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams.update({
    "font.family": ["Hiragino Sans", "Hiragino Kaku Gothic Pro", "sans-serif"],
    "font.size": 9, "axes.linewidth": 0.6, "axes.edgecolor": "#52514e",
    "axes.labelcolor": "#0b0b0b", "xtick.color": "#52514e", "ytick.color": "#52514e",
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.grid": True, "grid.color": "#e4e3df", "grid.linewidth": 0.5,
    "legend.frameon": False, "figure.facecolor": "white", "axes.facecolor": "white",
})
C = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"]          # 既定パレット slot 1〜4(固定順)
MK = ["o", "s", "^", "D"]; LS = ["-", "--", "-.", ":"]

rows = list(csv.DictReader(open("crm_mps_scaling_results.tsv"), delimiter="\t"))
def get(L, chi, obs, col):
    for r in rows:
        if int(r["L"]) == L and int(r["chi_prior"]) == chi and r["observable"] == obs:
            return float(r[col])
Ls = [8, 16, 32]

fig, ax = plt.subplots(1, 3, figsize=(7.2, 3.3), constrained_layout=True)

# (a) 忠実度
for k, chi in enumerate([2, 4, 8]):
    F = [get(L, chi, "ZZ onsite(i0)", "prior_fid") for L in Ls]
    ax[0].plot(Ls, F, LS[k], marker=MK[k], color=C[k], lw=1.5, ms=5, label=f"$\\chi_p={chi}$")
ax[0].set_xscale("log", base=2); ax[0].set_yscale("log")
ax[0].set_xticks(Ls, [str(L) for L in Ls]); ax[0].set_xlabel("系の長さ $L$")
ax[0].set_ylabel("大域的な忠実度 $F$"); ax[0].set_title("(a) prior の忠実度", loc="left", fontsize=9)
ax[0].legend(fontsize=7, loc="upper center", bbox_to_anchor=(0.5, -0.25), ncol=3)

# (b) 利得と L
series = [("ZZ onsite(i0)", 16, "オンサイト $ZZ$($\\chi_p=16$)"),
          ("ZZ up-up r=1", 2, "隣接 $Z_\\uparrow Z_\\uparrow$($\\chi_p=2$)"),
          ("SzSz r=1", 2, "隣接 $S^zS^z$($\\chi_p=2$)"),
          ("DoubleOcc(i0)", 16, "二重占有($\\chi_p=16$)")]
for k, (obs, chi, lab) in enumerate(series):
    G = [get(L, chi, obs, "G_exact") for L in Ls]
    ax[1].plot(Ls, G, LS[k], marker=MK[k], color=C[k], lw=1.5, ms=5, label=lab)
ax[1].set_xscale("log", base=2); ax[1].set_yscale("log")
ax[1].set_xticks(Ls, [str(L) for L in Ls]); ax[1].set_xlabel("系の長さ $L$")
ax[1].set_ylabel("利得 $G$($n_m=100$)"); ax[1].set_title("(b) 利得は $L$ によらない", loc="left", fontsize=9)
ax[1].set_ylim(4, 90); ax[1].legend(fontsize=7, loc="upper center", bbox_to_anchor=(0.5, -0.25), ncol=1)

# (c) 利得と n_m
nmr = list(csv.DictReader(open("crm_sweep_nm.tsv"), delimiter="\t"))
def nrep(nm): return 50 if nm >= 1000 else (400 if nm <= 3 else 100)
def gain(x, d, nA, nm):
    K = 3.0**nA - 1; vs = 3.0**nA*(1 - x*x)/nm
    return (K*x*x + vs)/(K*d*d + vs)
nmc = np.logspace(0, 3, 200)
for k, (obs, chi, lab) in enumerate([("ZZ onsite(i0)", 32, "オンサイト $ZZ$($\\chi_p=32$)"),
                                      ("ZZ onsite(i0)", 4, "オンサイト $ZZ$($\\chi_p=4$)"),
                                      ("ZZ up-up r=1", 32, "隣接 $Z_\\uparrow Z_\\uparrow$($\\chi_p=32$)")]):
    sel = [r for r in nmr if r["observable"] == obs and int(r["chi_prior"]) == chi]
    x, d = float(sel[0]["true"]), float(sel[0]["Delta"])
    ax[2].plot(nmc, gain(x, d, 2, nmc), LS[k], color=C[k], lw=1.5, label=lab)
    nm = np.array([int(r["nm"]) for r in sel]); ge = np.array([float(r["G_emp"]) for r in sel])
    err = np.array([np.sqrt(4/(nrep(m) - 1)) for m in nm])      # 分散の比の対数の標準偏差
    ax[2].errorbar(nm, ge, yerr=[ge*(1 - np.exp(-err)), ge*(np.exp(err) - 1)], fmt=MK[k], color=C[k],
                   ms=4, mfc="white", mew=1.2, elinewidth=0.8, capsize=0)
ax[2].set_xscale("log"); ax[2].set_yscale("log")
ax[2].set_xlabel("1 設定あたりのショット数 $n_m$"); ax[2].set_ylabel("利得 $G$($L=8$、$U=4$)")
ax[2].set_title("(c) 線: 厳密な式、点: 模擬", loc="left", fontsize=9)
ax[2].legend(fontsize=7, loc="upper center", bbox_to_anchor=(0.5, -0.25), ncol=1)

fig.savefig("thesis/figures/fig_locality.png", dpi=300)
print("saved thesis/figures/fig_locality.png")
