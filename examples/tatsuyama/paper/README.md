# 論文原稿 (LaTeX)

`crm_hubbard.tex` — 本研究(examples/tatsuyama/)を論文形式にまとめた作業原稿。図は `../crm_*.png` を参照する。

## ビルド

```bash
cd examples/tatsuyama/paper
lualatex crm_hubbard.tex   # 目次・参照のため2回実行
lualatex crm_hubbard.tex
```

LuaLaTeX + luatexja(Harano Aji フォント同梱、TeX Live に標準)で日本語をコンパイルする。生成物 `crm_hubbard.pdf`(12ページ)もコミット済み。

`crm_tutorial.tex` — 前提知識を最小限にした自習用の詳細解説。物性の用語(第二量子化・ハバード模型・Mott絶縁体・平均場・Slater行列式・Wickの定理・DMRG/MPS・JW変換)を数式で導入し、古典シャドウの分散〜CRMの利得法則・制御変量・Wick prior・matchgateまでを一から導出する。同じく `lualatex crm_tutorial.tex` を2回でPDF化(12ページ)。

`crm_observables.tex` — 本研究で推定している各物理量(観測量)の物理的意味を初歩から解説。占有数・オンサイトZZ・二重占有率・スピン相関・電荷相関・ホッピング(運動エネルギー)・全エネルギー・忠実度について、定義、量子ビット表現、物理的な読み方、本研究での実測値、そして「なぜ観測量ごとに測定効率(CRM利得)が違うのか」を整理する。同じく `lualatex crm_observables.tex` を2回でPDF化(9ページ)。

内容の対応: 詳細な実験記録は `../README.md`、論文の主張構成・査読対応は `../CRM_RESEARCH_NOTES.md`。
