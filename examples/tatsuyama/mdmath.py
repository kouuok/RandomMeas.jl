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
ESC = {'{': r'\lbrace', '}': r'\rbrace', '|': r'\vert'}

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
            new = new.replace('*', r'\ast')          # 生の * は強調に食われる
            new = re.sub(r'\\[,;:!]', '', new)       # 間隔命令はバックスラッシュを剥がされる
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

def do_verify(path):
    text = open(path, encoding='utf-8').read()
    html = ''.join(_render(c) for c in _chunks(text))
    n = len(re.findall(r'<math-renderer', html))
    h = re.sub(r'<pre.*?</pre>', '', html, flags=re.S)
    h = re.sub(r'<code.*?</code>', '', h, flags=re.S)
    h = re.sub(r'<math-renderer[^>]*>.*?</math-renderer>', '', h, flags=re.S)
    h = re.sub(r'<[^>]+>', '', h)
    bad = [l.strip() for l in h.split('\n') if '$' in l]
    print(f'{path}: 数式 {n}個がレンダリング / 未変換の $ を含む行 {len(bad)}個')
    for l in bad[:10]: print('   ', l[:120])
    return len(bad)

if __name__ == '__main__':
    if len(sys.argv) < 3: sys.exit(__doc__)
    cmd, paths = sys.argv[1], sys.argv[2:]
    bad = 0
    for p in paths:
        if cmd == 'fix': do_fix(p)
        elif cmd == 'verify': bad += do_verify(p)
        else: sys.exit(f'不明なコマンド: {cmd}')
    sys.exit(1 if bad else 0)
