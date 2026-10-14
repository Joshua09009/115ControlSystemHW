# Homework 01 Problem 3（電容麥克風）

編譯（需 XeLaTeX）：

```bash
cd hw01_p3
latexmk -xelatex main.tex
```

或在 Overleaf：Compiler 選 **XeLaTeX**。若 Overleaf 找不到「Microsoft YaHei」，把 `main.tex` 裡的

```latex
\usepackage{fontspec}
\usepackage{xeCJK}
\setmainfont{Times New Roman}
\setCJKmainfont{DFKai-SB}[AutoFakeBold,AutoFakeSlant]
```

改成

```latex
\usepackage{ctex}
```

圖片在 `figures/`。正文請自行填入姓名學號。
