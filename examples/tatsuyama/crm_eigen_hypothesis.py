"""仮説「prior が観測量の固有状態に近いほど利得が大きい」を統合表で検証する。

単一パウリ文字列 P(固有値 ±1)について、状態 τ が固有状態からどれだけ離れているかを
    η_τ = (1 - |<P>_τ|) / 2      (最も近い固有空間の外にある確率の重み。0 なら固有状態)
で測る。x=<P>_ρ, y=<P>_σ とし、同符号なら |Δ| = 2|η_σ - η_ρ| が恒等的に成り立つ。
"""
import csv, collections, math

rows = [r for r in csv.DictReader(open("crm_edfid_all.tsv"), delimiter="\t")
        if r["single_pauli"] == "true" and r["exact"] == "yes"]
eta = lambda v: (1 - abs(v)) / 2
TOL = 1e-6

def classify(x, y):
    if abs(y) < TOL: return "零"          # prior が P について何も言わない (y=0)
    if x * y < 0:    return "逆符号"
    return "過信" if eta(y) < eta(x) - TOL else ("控えめ" if eta(y) > eta(x) + TOL else "一致")

for r in rows:
    r["x"], r["y"] = float(r["true"]), float(r["prior_val"])
    r["eta_r"], r["eta_s"] = eta(r["x"]), eta(r["y"])
    r["cls"] = classify(r["x"], r["y"])
    r["Gf"], r["Gm"] = float(r["G"]), float(r["G_max"])
    r["sys"] = (r["W"], r["geometry"], int(r["LX"]), float(r["U"]))
    r["fam"] = "UHF" if r["prior"] == "UHF" else ("MPS χ≤4" if int(r["prior"][3:]) <= 4 else "MPS χ≥8")

def spearman(a, b):
    def rank(v):
        o = sorted(range(len(v)), key=lambda i: v[i]); rk = [0.0]*len(v); i = 0
        while i < len(v):
            j = i
            while j+1 < len(v) and v[o[j+1]] == v[o[i]]: j += 1
            for k in range(i, j+1): rk[o[k]] = (i+j)/2
            i = j+1
        return rk
    ra, rb = rank(a), rank(b); n = len(a)
    ma, mb = sum(ra)/n, sum(rb)/n
    num = sum((p-ma)*(q-mb) for p, q in zip(ra, rb))
    den = math.sqrt(sum((p-ma)**2 for p in ra) * sum((q-mb)**2 for q in rb))
    return num/den if den > 0 else float("nan")

OBS = ("ZZ onsite", "ZZ up-up nb")
print(f"厳密な系の単一パウリ行: {len(rows)}\n")

# ---- 1. 恒等式の照合 -------------------------------------------------------
bad = sum((r["Gf"] > 1 + 1e-9) != (0 < r["y"]/r["x"] < 2) for r in rows if abs(r["Gf"]-1) > 1e-6)
print("[1] G>1 ⟺ 0 < y/x < 2 が破れた行:", bad)

# ---- 2. 実際の prior はどちら側にいるか -----------------------------------
print("\n[2] prior の分類 (過信 = ρ より固有状態に近い)")
for obs in OBS:
    for fam in ("UHF", "MPS χ≤4", "MPS χ≥8"):
        c = collections.Counter(r["cls"] for r in rows if r["observable"] == obs and r["fam"] == fam)
        n = sum(c.values())
        print(f"  {obs:12s} {fam:9s} n={n:4d}  " + "  ".join(f"{k} {c[k]:4d} ({100*c[k]/n:4.1f}%)"
              for k in ("過信", "控えめ", "一致", "逆符号", "零")))

# ---- 3. 仮説そのもの: UHF, 格子を固定して U を振る -------------------------
print("\n[3] UHF・格子固定で U を振ったとき、η_σ が小さい U ほど G が大きいか")
for obs in OBS:
    agree = rev = mixed = 0
    lat = collections.defaultdict(list)
    for r in rows:
        if r["observable"] == obs and r["prior"] == "UHF" and r["cls"] != "零":
            lat[r["sys"][:3]].append(r)
    for k, v in lat.items():
        if len(v) < 3: continue
        s = spearman([r["eta_s"] for r in v], [r["Gf"] for r in v])
        if s <= -0.99: agree += 1
        elif s >= 0.99: rev += 1
        else: mixed += 1
    print(f"  {obs:12s} 格子 {agree+rev+mixed:2d}:  仮説どおり(ρ_s=-1) {agree:2d}  逆向き(ρ_s=+1) {rev:2d}  混在 {mixed:2d}")

# ---- 4. ρ を固定して prior を振る ------------------------------------------
print("\n[4] 系を固定して prior を振る: 最も固有状態に近い prior が最良か")
for obs in OBS:
    by = collections.defaultdict(list)
    for r in rows:
        if r["observable"] == obs: by[r["sys"]].append(r)
    n = best_is_eig = 0; rs_eta = []; rs_dist = []
    for k, v in by.items():
        if len(v) < 4: continue
        n += 1
        eig = min(v, key=lambda r: r["eta_s"]); best = max(v, key=lambda r: r["Gf"])
        best_is_eig += (eig is best)
        rs_eta.append(spearman([r["eta_s"] for r in v], [r["Gf"] for r in v]))
        rs_dist.append(spearman([abs(r["eta_s"]-r["eta_r"]) for r in v], [r["Gf"] for r in v]))
    rs_eta.sort()
    q = lambda a, p: a[min(len(a)-1, int(p*len(a)))]
    print(f"  {obs:12s} 系 {n}: 最も固有状態に近い prior が最良なのは {best_is_eig} 系 ({100*best_is_eig/n:.0f}%)")
    print(f"     Spearman(η_σ, G)        中央値 {q(rs_eta,.5):+.2f}  [10%,90%] = [{q(rs_eta,.1):+.2f}, {q(rs_eta,.9):+.2f}]")
    print(f"     Spearman(|η_σ-η_ρ|, G)  最大値 {max(rs_dist):+.3f}  (恒等的に -1 のはず)")

# ---- 5. 自然実験: OBC 鎖の強い結合と弱い結合 -------------------------------
print("\n[5] OBC 鎖・UHF・ZZ up-up nb: η_σ はほぼ同じで、中央の結合の強弱だけが違う")
print(f"  {'L':>3s} {'結合':4s} {'U':>3s} {'η_σ':>7s} {'η_ρ':>7s} {'G':>7s}")
for r in sorted((r for r in rows if r["observable"] == "ZZ up-up nb" and r["prior"] == "UHF"
                 and r["W"] == "1" and r["geometry"] == "cylinder" and float(r["U"]) in (8.0, 12.0)),
                key=lambda r: (float(r["U"]), int(r["LX"]))):
    L = int(r["LX"]); strong = (L//2) % 2 == 1        # 結合 (L/2, L/2+1): 左端が奇数なら強い結合
    print(f"  {L:3d} {'強' if strong else '弱':4s} {float(r['U']):3.0f} {r['eta_s']:7.4f} {r['eta_r']:7.4f} {r['Gf']:7.2f}")

# ---- 6. ほぼ固有状態の prior だけを取り出す --------------------------------
print("\n[6] ほぼ固有状態の prior (η_σ<0.02, 同符号): G>1 は |x|>|y|/2 で決まる")
sub = [r for r in rows if r["eta_s"] < 0.02 and r["x"]*r["y"] > 0]
win = [r for r in sub if r["Gf"] > 1]; lose = [r for r in sub if r["Gf"] < 1]
print(f"  該当 {len(sub)} 行: 勝ち {len(win)}, 負け {len(lose)}")
if win:  print(f"  勝ちの |x| 最小値 {min(abs(r['x']) for r in win):.4f}")
if lose: print(f"  負けの |x| 最大値 {max(abs(r['x']) for r in lose):.4f}")
for obs in OBS:
    s = [r for r in sub if r["observable"] == obs]
    if s: print(f"  {obs:12s} {len(s):3d} 行  負け {sum(r['Gf']<1 for r in s):3d}  prior 内訳 {dict(collections.Counter(r['fam'] for r in s))}")
