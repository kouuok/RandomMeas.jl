#!/usr/bin/env python3
"""GitHub 上で数式が正しく表示されるように markdown を整形・検査する。

README 付録「GitHub で数式を正しく表示させるための書式規則」を実装したもの。
規則は gh api /markdown(GitHub 実物のレンダラ)に通した実測から決めた。

  python3 mdmath.py fix    README.md   # 書式を直す
  python3 mdmath.py verify README.md   # 実際にレンダリングして未変換の $ を探す
"""
import json, re, subprocess, sys

FENCE = re.compile(r'^(\s*)(```|~~~)')
OK_BEFORE = set(' \t([*_>|~"\'')
ESC = {'{': r'\lbrace ', '}': r'\rbrace ', '|': r'\vert '}   # 後ろの空白は、次の文字と命令名がつながらないようにするため

def code_lines(text):
    out, inf, fc = [], False, None
    for ln in text.split('\n'):
        m = FENCE.match(ln)
        if m:
            if not inf: inf, fc = True, m.group(2)
            elif ln.strip().startswith(fc): inf = False
            out.append(True); continue
        out.append(inf)
    return out

def inline_code(line):
    return [(m.start(), m.end()) for m in re.finditer(r'`[^`\n]*`', line)]

def fix_line(line, st):
    prot = inline_code(line)
    def repl(m):
        body = m.group(0)
        new = re.sub(r'\\([{}|])', lambda k: ESC[k.group(1)], body)
        if not body.startswith('$$'):
            new = new.replace('*', r'\ast ')         # 生の * は強調に食われる(空白で次の文字とつながらないようにする)
        # 1行の $$…$$ も中身が Markdown として処理されるので、以下は両方に当てる
        new = re.sub(r'\\[,;:!]', r'\\thinspace ', new)  # \, などは Markdown にバックスラッシュを剥がされるので、剥がされない \thinspace に置き換える。
        # 以前はここで空文字に置き換えていたが、それだと \gamma\,U が未定義の \gammaU になって表示が壊れた
        new = re.sub(r'(?<!\\)\\%', r'\\\\%', new)   # \% は % になり TeX の注釈記号として式の残りを消すので \\% にする
        if not body.startswith('$$'):
            new = re.sub(r'\s+\$$', '$', new)           # 閉じ $ の直前に空白を残さない
        if new != body: st['esc'] += 1
        return new
    out, last = [], 0
    for m in re.finditer(r'\$\$.+?\$\$|\$[^$\n]+?\$', line):
        if any(a <= m.start() < b for a, b in prot): continue
        out.append(line[last:m.start()]); out.append(repl(m)); last = m.end()
    out.append(line[last:]); line = ''.join(out)

    prot = inline_code(line); res, i, n = [], 0, len(line)
    while i < n:
        if line[i] == '$' and not any(a <= i < b for a, b in prot):
            d = 2 if line.startswith('$$', i) else 1
            close = line.find('$$' if d == 2 else '$', i + d)
            if close == -1 or (d == 1 and close == i + 1):
                res.append(line[i]); i += 1; continue
            if d == 1 and line[i+d:close].strip() == '':
                res.append(line[i]); i += 1; continue
            prev = res[-1] if res else ''
            if d == 1 and prev != '' and prev not in OK_BEFORE:
                res.append(' '); st['open'] += 1          # 開き $ の直前は空白が要る
            res.append(line[i:close+d]); j = close + d
            if d == 1 and j < n and re.match(r'[0-9A-Za-z]', line[j]):
                res.append(' '); st['close'] += 1         # 閉じ $ の直後に英数字は不可
            i = j
        else:
            res.append(line[i]); i += 1
    return ''.join(res)

def do_fix(path):
    text = open(path, encoding='utf-8').read()
    inc, lines = code_lines(text), text.split('\n')
    st = {'open': 0, 'close': 0, 'esc': 0}
    for i, ln in enumerate(lines):
        if inc[i] or '$' not in ln: continue
        lines[i] = fix_line(ln, st)
    open(path, 'w', encoding='utf-8').write('\n'.join(lines))
    print(f'{path}: 開き空白 {st["open"]} / 閉じ空白 {st["close"]} / エスケープ {st["esc"]}')

def _render(text):
    p = subprocess.run(['gh', 'api', '-X', 'POST', '/markdown', '--input', '-'],
                       input=json.dumps({'text': text, 'mode': 'gfm'}),
                       capture_output=True, text=True)
    if p.returncode: sys.exit('gh api 失敗: ' + p.stderr[:300])
    return p.stdout

def _payload(text):
    """API に送る JSON の実バイト数。日本語はユニコードエスケープに展開されて
    約2倍になるので、生テキストの長さで測ると上限判定を誤る。"""
    return len(json.dumps({'text': text, 'mode': 'gfm'}))

def _chunks(text, limit=350_000):
    """API の 400 KB 上限を避けて分割する。見出し境界で切るので
    コードフェンスや <details> を途中で割らない。"""
    parts, cur = [], []
    for block in re.split(r'(?m)(?=^## )', text):
        if cur and _payload(''.join(cur) + block) > limit:
            parts.append(''.join(cur)); cur = []
        cur.append(block)
    if cur: parts.append(''.join(cur))
    return parts

# この研究の文書で使っている LaTeX 命令。ここにない命令名は、命令と次の文字が
# つながってしまった(例: \gamma\,U を整形して \gammaU になった)可能性が高い。
# GitHub の API は数式を描画しないので、この確認をしないと見逃す。
KNOWN_CMDS = set(['Big', 'Bigl', 'Bigr', 'Delta', 'Gamma', 'Longleftrightarrow', 'Longrightarrow', 'Omega', 'Phi', 'Pi', 'Pr', 'Rightarrow', 'Sigma', 'Tr', 'alpha', 'approx', 'arg', 'ast', 'bar', 'begin', 'beta', 'bigl', 'bigotimes', 'bigr', 'binom', 'boxed', 'cap', 'cdot', 'cdots', 'chi', 'cos', 'cup', 'dagger', 'delta', 'det', 'dfrac', 'dots', 'downarrow', 'ell', 'end', 'epsilon', 'equiv', 'eta', 'exp', 'frac', 'gamma', 'ge', 'gg', 'gtrsim', 'hat', 'iff', 'in', 'infty', 'int', 'kappa', 'lVert', 'lambda', 'langle', 'lbrace', 'lceil', 'le', 'left', 'leftrightarrow', 'lesssim', 'll', 'ln', 'log', 'longrightarrow', 'lvert', 'mathbb', 'mathbf', 'mathcal', 'mathrm', 'mathsf', 'max', 'mid', 'min', 'mp', 'mu', 'ne', 'neq', 'nu', 'otimes', 'overline', 'partial', 'phi', 'pi', 'pm', 'pmod', 'prod', 'propto', 'psi', 'qquad', 'quad', 'rVert', 'rangle', 'rbrace', 'rceil', 'rho', 'right', 'rm', 'rvert', 'sigma', 'sim', 'simeq', 'sin', 'sqrt', 'subseteq', 'sum', 'sup', 'tau', 'text', 'tfrac', 'theta', 'thinspace', 'tilde', 'times', 'to', 'underbrace', 'uparrow', 'varepsilon', 'varphi', 'vert', 'wedge', 'widehat', 'xi', 'xrightarrow'])

def unknown_commands(text):
    t = re.sub(r'```math\n(.*?)```', lambda m: '$$' + m.group(1) + '$$', text, flags=re.S)
    t = re.sub(r'```.*?```', '', t, flags=re.S)
    bad = []
    for m in re.finditer(r'\$\$.+?\$\$|\$[^$\n]+?\$', t, flags=re.S):
        body = m.group(0).replace('\\\\', ' ')      # 行替えの \\ を除く
        for c in re.findall(r'\\([A-Za-z]+)', body):
            if c not in KNOWN_CMDS:
                bad.append((c, m.group(0)[:100]))
    return bad

def compare_inline(text, html_out):
    r"""文中の数式(1行の $$…$$ を含む)を1つずつ、GitHub が数式描画に渡す文字列と突き合わせる。
    \, が , に、\% が % になるなど、Markdown が式を書き換えた箇所を見つける。
    < と > が実体参照になる違いは描画に影響しないので無視する。
    規則の説明のために $…$ をわざと文字として書いた箇所(… を含むもの)は除く。"""
    import difflib, html as H
    lines = text.split('\n'); inc = code_lines(text); src = []
    for i, ln in enumerate(lines):
        if inc[i]: continue
        prot = inline_code(ln)
        for mm in re.finditer(r'\$\$.+?\$\$|\$`.+?`\$|\$[^$\n]+?\$', ln):
            if any(a <= mm.start() < b for a, b in prot): continue
            s = mm.group(0).replace('$`', '$').replace('`$', '$')
            s = s.replace('\\\\%', '\\%')   # \\% が \% になるのは意図どおり
            if '…' in s: continue
            src.append((i + 1, s))
    got = [H.unescape(H.unescape(x)) for x in
           re.findall(r'<math-renderer[^>]*js-inline-math[^>]*>(.*?)</math-renderer>', html_out, flags=re.S)]
    sm = difflib.SequenceMatcher(a=[s for _, s in src], b=got, autojunk=False)
    out = []
    for op, a1, a2, b1, b2 in sm.get_opcodes():
        if op == 'equal': continue
        for k in range(max(a2 - a1, b2 - b1)):
            ln, s = src[a1 + k] if a1 + k < a2 else (None, '(なし)')
            g = got[b1 + k] if b1 + k < b2 else None
            out.append((ln, s, g))
    return out

def do_verify(path):
    text = open(path, encoding='utf-8').read()
    html = ''.join(_render(c) for c in _chunks(text))
    n = len(re.findall(r'<math-renderer', html))
    # ```math が数式にならなかった場合 GitHub は <pre lang="math"> のまま返す。
    # 下で <pre> を落とすので、先に数えておかないと見逃す。
    dead = re.findall(r'<pre lang="math"[^>]*><code[^>]*>(.*?)</code></pre>', html, re.S)
    h = re.sub(r'<pre.*?</pre>', '', html, flags=re.S)
    h = re.sub(r'<code.*?</code>', '', h, flags=re.S)
    h = re.sub(r'<math-renderer[^>]*>.*?</math-renderer>', '', h, flags=re.S)
    h = re.sub(r'<[^>]+>', '', h)
    bad = [l.strip() for l in h.split('\n') if '$' in l]
    print(f'{path}: 数式 {n}個がレンダリング / 未変換の $ を含む行 {len(bad)}個'
          f' / 表示されない math ブロック {len(dead)}個')
    for l in bad[:10]: print('   ', l[:120])
    mism = compare_inline(text, html)
    print(f'{path}: 元の式と描画に渡る式の食い違い {len(mism)}個')
    for ln, s, g in mism[:10]: print(f'    行{ln}\n      元  : {s[:120]}\n      描画: {str(g)[:120]}')
    unk = unknown_commands(text)
    print(f'{path}: 知らない命令名 {len(unk)}個(命令と文字がつながっていないか確認)')
    for c, s in unk[:10]: print(f'    [\\{c}]', s.replace(chr(10), ' '))
    for d in dead[:10]: print('    [math未描画]', d.strip()[:100])
    return len(bad) + len(dead) + len(unk) + len(mism)

if __name__ == '__main__':
    if len(sys.argv) < 3: sys.exit(__doc__)
    cmd, paths = sys.argv[1], sys.argv[2:]
    bad = 0
    for p in paths:
        if cmd == 'fix': do_fix(p)
        elif cmd == 'verify': bad += do_verify(p)
        else: sys.exit(f'不明なコマンド: {cmd}')
    sys.exit(1 if bad else 0)
