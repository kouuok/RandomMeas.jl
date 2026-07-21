# 論文原稿 (LaTeX)

`crm_hubbard.tex` — 本研究(examples/tatsuyama/)を論文形式にまとめた作業原稿。図は `../crm_*.png` を参照する。

## ビルド

```bash
cd examples/tatsuyama/paper
lualatex crm_hubbard.tex   # 目次・参照のため2回実行
lualatex crm_hubbard.tex
```

LuaLaTeX + luatexja(Harano Aji フォント同梱、TeX Live に標準)で日本語をコンパイルする。生成物 `crm_hubbard.pdf`(12ページ)もコミット済み。

内容の対応: 詳細な実験記録は `../README.md`、論文の主張構成・査読対応は `../CRM_RESEARCH_NOTES.md`。
