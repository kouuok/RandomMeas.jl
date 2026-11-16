#!/usr/bin/env python3
r"""README_prior_comparison.md の図を作る。

新しい計算はしない。数値はすべて既存の文書の表からの再掲である
(どれも n_m = 100 での厳密な分散の比)。

  図1  crm_new_fig_prior_compare.png
       U = 12、厳密な参照状態を持つ 20 系の利得の中央値
       出典: README_phf_mps.md 5 章、README_varsym.md 2.1 節、
             README_eigenstate.md 3.6 節 (RHF)、README_observables.md 5.2 節 (RHF の二重占有)
  図2  crm_new_fig_prior_distance.png
       1 次元 L = 64 開放端、U = 12、距離 r の Z↑Z↑ の利得の中央値
       出典: README_varsym_distance.md 2.1 節

実行: /opt/anaconda3/bin/python3 crm_prior_comparison_figure.py
"""
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams["font.family"] = "Hiragino Sans"
plt.rcParams["axes.unicode_minus"] = False

INK = "#1f2328"      # 文字
MUTED = "#656d76"    # 補助の文字・軸
GRID = "#d8dee4"
GAIN = "#0072B2"     # 得 (G > 1)
LOSS = "#D55E00"     # 損 (G < 1)
EVEN = "#8c959f"     # 得も損もしない (G ≈ 1)
BOX = dict(facecolor="white", edgecolor="none", pad=1.2)   # 値のラベルの下地

# ---------------------------------------------------------------- 図1
PRIORS = [
    "RHF", "UHF", "UHF-sym", "PHF",
    "切断 χ=4", "切断 χ=8",
    "変分 χ=4", "変分 χ=8", "変分 χ=16",
    "変分 χ=4-sym", "変分 χ=8-sym", "変分 χ=16-sym",
]
GROUP_BREAKS = [4, 6, 9]          # この番号の prior の前に区切り線
GROUP_LABELS = [(1.5, "平均場"), (4.5, "切断 MPS\n(厳密な状態が要る)"),
                (7.0, "変分 MPS"), (10.0, "変分 MPS\nの回転平均")]

PANELS = [
    ("オンサイト ZZ\n(電荷の量)",
     [1.00, 465, 465, 513, 151, 546, 553, 555, 555, 553, 555, 555], 555),
    ("隣接 Z↑Z↑\n(隣のスピン相関)",
     [8.67, 1.22, 6.95, 14.2, 29.2, 29.4, 15.3, 16.8, 28.8, 18.0, 19.8, 26.7], 31.6),
    ("単一サイト Z↑\n(向きを持つ量、真の値は 0)",
     [1.00, 0.02, 1.00, 1.00, 1.00, 1.00, 0.04, 0.11, 0.46, 1.00, 1.00, 1.00], 1.0),
    ("二重占有\n(4項の和として測る)",
     [1.00, 1.81, 123, 127, 77.8, 128, 4.75, 13.3, 53.2, 129, 129, 129], None),
]


def color_of(g):
    if g < 0.99:
        return LOSS
    if g <= 1.01:
        return EVEN
    return GAIN


def fmt(g):
    if g >= 100:
        return f"{g:.0f}"
    if g >= 10:
        return f"{g:.1f}"
    return f"{g:.2f}"


def figure1():
    n = len(PRIORS)
    fig, axes = plt.subplots(1, 4, figsize=(13.5, 6.4), sharey=True)
    ys = list(range(n))
    for ax, (title, vals, ceil) in zip(axes, PANELS):
        ax.set_xscale("log")
        ax.set_xlim(0.008, 4000)
        ax.set_ylim(n - 0.4, -0.6)
        for y, g in zip(ys, vals):
            c = color_of(g)
            if abs(g - 1) <= 0.01:
                ax.plot([1], [y], marker="o", ms=6, color=c, zorder=3)
                ax.text(1.5, y, fmt(g), va="center", ha="left", fontsize=9, color=INK, bbox=BOX, zorder=5)
            elif g > 1:
                ax.barh(y, g - 1, left=1, height=0.62, color=c, zorder=3)
                ax.text(g * 1.25, y, fmt(g), va="center", ha="left", fontsize=9, color=INK, bbox=BOX, zorder=5)
            else:
                ax.barh(y, 1 - g, left=g, height=0.62, color=c, zorder=3)
                ax.text(g / 1.25, y, fmt(g), va="center", ha="right", fontsize=9, color=INK, bbox=BOX, zorder=5)
        ax.axvline(1, color=INK, lw=1.0, zorder=2)
        if ceil is not None and ceil > 1:
            ax.axvline(ceil, color=MUTED, lw=1.0, ls=(0, (4, 3)), zorder=1)
            ax.text(ceil, -0.72, f"天井 {fmt(ceil)}", ha="center", va="bottom",
                    fontsize=8.5, color=MUTED)
        elif ceil is not None:
            ax.text(1, -0.72, "天井 1(得はできない)", ha="center", va="bottom",
                    fontsize=8.5, color=MUTED)
        else:
            ax.text(1, -0.72, "損得の境目", ha="center", va="bottom",
                    fontsize=8.5, color=MUTED)
        for b in GROUP_BREAKS:
            ax.axhline(b - 0.5, color=GRID, lw=1.0, zorder=0)
        ax.set_title(title, fontsize=11, color=INK, pad=22)
        ax.set_xticks([0.01, 0.1, 1, 10, 100, 1000])
        ax.set_xticklabels(["0.01", "0.1", "1", "10", "100", "1000"])
        ax.tick_params(axis="x", colors=MUTED, labelsize=8.5, length=3)
        ax.minorticks_off()
        ax.tick_params(axis="y", length=0)
        ax.xaxis.grid(True, color=GRID, lw=0.6, zorder=0)
        ax.set_axisbelow(True)
        for s in ("top", "right", "left"):
            ax.spines[s].set_visible(False)
        ax.spines["bottom"].set_color(GRID)
        ax.set_xlabel("利得 G(対数目盛)", fontsize=9, color=MUTED)
    axes[0].set_yticks(ys)
    axes[0].set_yticklabels(PRIORS, fontsize=10, color=INK)
    for y, lab in GROUP_LABELS:
        axes[0].text(-0.60, y, lab, transform=axes[0].get_yaxis_transform(),
                     ha="right", va="center", fontsize=9, color=MUTED)

    handles = [
        plt.Rectangle((0, 0), 1, 1, color=GAIN),
        plt.Rectangle((0, 0), 1, 1, color=LOSS),
        plt.Line2D([0], [0], marker="o", ls="", ms=6, color=EVEN),
        plt.Line2D([0], [0], color=MUTED, lw=1.0, ls=(0, (4, 3))),
    ]
    fig.legend(handles,
               ["得(G > 1、右へ伸びる)", "損(G < 1、左へ伸びる)",
                "得も損もしない(G = 1)", "天井(完全な prior の利得)"],
               loc="lower center", ncol=4, frameon=False, fontsize=9.5,
               bbox_to_anchor=(0.56, 0.0))
    fig.suptitle("prior ごとの利得 — U = 12、厳密な参照状態を持つ 20 系の中央値、"
                 "1設定あたり 100 ショット",
                 fontsize=12.5, color=INK, x=0.56, y=0.985)
    fig.subplots_adjust(left=0.215, right=0.985, top=0.83, bottom=0.14, wspace=0.10)
    fig.savefig("crm_new_fig_prior_compare.png", dpi=170, facecolor="white")
    plt.close(fig)


# ---------------------------------------------------------------- 図2
R_LABELS = ["1", "2", "3", "4", "5〜8", "9〜16", "17〜31"]
# (名前, 値, 色, 線種, 右端のラベルの縦位置)
SERIES = [
    ("変分 χ=16(対称な解)", [56.4, 5.7, 5.2, 2.5, 1.8, 1.21, 1.06], "#0072B2", "-", 1.95),
    ("変分 χ=8-sym", [43.4, 5.7, 4.9, 2.5, 1.7, 1.21, 1.06], "#56B4E9", "-", 1.30),
    ("UHF-conn(連結相関)", [1.09, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00], "#8c959f", "-", 0.86),
    ("変分 χ=4-sym", [31.2, 5.7, 3.9, 2.5, 1.7, 0.86, 0.63], "#009E73", "-", 0.56),
    ("UHF-sym", [4.5, 3.3, 2.7, 0.61, 0.32, 0.17, 0.13], "#E69F00", "-", 0.13),
    ("UHF", [3.3, 0.12, 0.10, 0.04, 0.03, 0.02, 0.01], "#D55E00", "-", 0.012),
]


def figure2():
    fig, ax = plt.subplots(figsize=(9.2, 5.4))
    xs = list(range(len(R_LABELS)))
    ax.set_yscale("log")
    ax.set_ylim(0.007, 110)
    ax.set_xlim(-0.25, len(xs) - 1 + 1.95)
    ax.axhspan(0.007, 1, color=LOSS, alpha=0.07, lw=0, zorder=0)
    ax.axhline(1, color=INK, lw=1.0, zorder=1)
    ax.text(4.3, 40, "得をする領域(G > 1)", fontsize=9.5, color=MUTED, va="center")
    ax.text(-0.15, 0.013, "損をする領域(G < 1)", fontsize=9.5, color=MUTED, va="center")
    for name, vals, c, ls, ylab in SERIES:
        ax.plot(xs, vals, color=c, lw=2.0, ls=ls, marker="o", ms=5.5,
                mec="white", mew=1.2, zorder=3)
        ax.plot([xs[-1] + 0.08, xs[-1] + 0.30], [vals[-1], ylab], color=c, lw=0.8, zorder=2)
        ax.text(xs[-1] + 0.36, ylab, name, color=INK, fontsize=9.5, va="center")
    ax.set_xticks(xs)
    ax.set_xticklabels(R_LABELS)
    ax.set_yticks([0.01, 0.1, 1, 10, 100])
    ax.set_yticklabels(["0.01", "0.1", "1", "10", "100"])
    ax.tick_params(colors=MUTED, labelsize=9.5, length=3)
    ax.minorticks_off()
    ax.yaxis.grid(True, color=GRID, lw=0.6, zorder=0)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.spines["left"].set_color(GRID)
    ax.spines["bottom"].set_color(GRID)
    ax.set_xlabel("サイト間の距離 r", fontsize=10, color=MUTED)
    ax.set_ylabel("利得 G の中央値(対数目盛)", fontsize=10, color=MUTED)
    ax.set_title("距離 r の Z↑Z↑ の利得 — 1次元 L = 64 開放端、U = 12、"
                 "1設定あたり 100 ショット",
                 fontsize=11.5, color=INK, pad=12)
    fig.subplots_adjust(left=0.09, right=0.985, top=0.90, bottom=0.12)
    fig.savefig("crm_new_fig_prior_distance.png", dpi=170, facecolor="white")
    plt.close(fig)


if __name__ == "__main__":
    figure1()
    figure2()
    print("wrote crm_new_fig_prior_compare.png, crm_new_fig_prior_distance.png")
