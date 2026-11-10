# AGENTS.md —— 给 AI 的作图规范

> **任何 AI 进入本仓库，写或改绘图代码之前先读这份文件。**
> 适用于 Claude Code / Cursor / GitHub Copilot / ChatGPT / DeepSeek 等。

---

## 这个仓库是什么

GWAS 与群体遗传学论文配图的统一作图库。`R/` 下按分析类型分三个子文件夹：

```
R/
├── theme_popgen.R                          # 统一样式（必须最先加载）
├── load_all.R                              # 递归载入 R/ 下全部函数
├── Population Genetic Structure Analysis/  # 群体结构
│   ├── PCA_2D_Plot.R
│   ├── PCA_3D_Plot.R
│   ├── Admixture_Plot.R
│   └── Tree_Pie_Plot.R
├── Population Genetic Diversity Analysis/  # 遗传多样性
│   ├── LD_Decay_Plot.R
│   └── SNP_Density_Plot.R
├── Selection Signature Analysis/           # 选择信号
│   ├── FST_Manhattan_Plot.R
│   ├── XPCLR_Manhattan_Plot.R
│   ├── XPEHH_Manhattan_Plot.R
│   └── legacy/                             # 改造前的原始手写脚本，仅作参考，
│                                           # load_all.R 会跳过这个目录
└── General Plots/                          # 通用
    └── Venn_Plot.R                          # 韦恩图（2~4 集合，复刻 Venny）
```

### 命名约定

**函数名 = 词之间用下划线分隔、以 `_Plot` 结尾，缩写保持全大写**，
与文件名完全一致（`LD_Decay_Plot.R` 里就是 `LD_Decay_Plot()`）。

| 缩写 | 写法 | 例 |
|---|---|---|
| PCA / LD / SNP / FST | 保持全大写 | `PCA_2D_Plot` `LD_Decay_Plot` `SNP_Density_Plot` |
| XP-CLR / XP-EHH | 去掉连字符，全大写 | `XPCLR_Manhattan_Plot` `XPEHH_Manhattan_Plot` |

不要再用 `plot_xxx` 前缀，也不要中文文件名。

| 场景 | 函数 | 文件 |
|---|---|---|
| PCA 二维散点（含 95% 椭圆） | `PCA_2D_Plot()` | `Population Genetic Structure Analysis/PCA_2D_Plot.R` |
| PCA 碎石图 | `PCA_Scree_Plot()` | `Population Genetic Structure Analysis/PCA_2D_Plot.R` |
| PCA 三维散点 | `PCA_3D_Plot()` | `Population Genetic Structure Analysis/PCA_3D_Plot.R` |
| ADMIXTURE 堆叠柱（多 K） | `Admixture_Plot()` | `Population Genetic Structure Analysis/Admixture_Plot.R` |
| 系统发育树 + 祖先成分饼图 | `Tree_Pie_Plot()` | `Population Genetic Structure Analysis/Tree_Pie_Plot.R` |
| LD 衰减曲线 | `LD_Decay_Plot()` | `Population Genetic Diversity Analysis/LD_Decay_Plot.R` |
| SNP 密度图 | `SNP_Density_Plot()` | `Population Genetic Diversity Analysis/SNP_Density_Plot.R` |
| Fst 曼哈顿图 | `FST_Manhattan_Plot()` | `Selection Signature Analysis/FST_Manhattan_Plot.R` |
| XP-CLR 曼哈顿图 | `XPCLR_Manhattan_Plot()` | `Selection Signature Analysis/XPCLR_Manhattan_Plot.R` |
| XP-EHH 曼哈顿图 | `XPEHH_Manhattan_Plot()` | `Selection Signature Analysis/XPEHH_Manhattan_Plot.R` |
| 韦恩图（2~4 集合） | `Venn_Plot()` | `General Plots/Venn_Plot.R` |

**唯一目标：不管哪个 AI、哪台机器、哪次生成，图必须长得一样。**

---

## 铁律（违反即视为改坏）

### 1. 样式只能从 `R/theme_popgen.R` 取

禁止在任何作图函数里写死颜色、字号、点大小、导出宽高、dpi。

```r
# 错
geom_point(size = 2, colour = "red")
ggsave(f, plot, width = 10, height = 4.5, dpi = 500)

# 对
geom_point(size = POPGEN_POINT_SIZE)
scale_colour_manual(values = popgen_colors(3))
save_popgen(plot, f)          # 内部用 POPGEN_WIDTH / HEIGHT / DPI
```

要改风格 → **只改 `theme_popgen.R`**，全仓库同步生效。不要在单个函数里打补丁。

### 2. 禁止硬编码路径

所有输入输出路径必须走函数参数。代码里不许出现 `C:/Users/...`、`H:/...`、`Desktop`、`D:/`。

### 3. 改动画法必须说明理由

这是论文用图。几何对象（`geom_*` / `barplot` / `CMplot`）、aes 映射、图层顺序、坐标轴含义，
一旦改动会直接影响已发表的图，**改之前必须说明为什么**。

### 4. 不新增依赖

现有依赖已经够用：
`ggplot2`、`readr`、`data.table`、`scatterplot3d`、`CMplot`。
确实需要新包，先问，不要直接加。

### 5. 改完必须跑验证

```bash
Rscript tools/run_demo.R
```

然后检查 `demo_out/` 里的图是否正常。图跑不出来或明显变形 = 改坏了。

---

## 怎么加载

在仓库根目录：

```r
source("R/load_all.R")
```

`load_all.R` 会先切 UTF-8 locale 再依次载入各文件（原因见文末「已知坑 1」）。

---

## 函数签名

### PCA_2D_Plot()

```r
PCA_2D_Plot(
  eigenvec_path,            # PLINK eigenvec（FID IID PC1..PCn）
  eigenval_path,            # PLINK eigenval
  group_path,               # 分组表，需含 SampleID + 分组列
  group_col   = "Breed",    # 分组列的列名
  sample_col  = "SampleID",
  out_file    = NULL,       # NULL = 不落盘，只返回 ggplot 对象
  colors      = NULL,       # NULL = 用统一色板 popgen_colors()
  point_size  = POPGEN_POINT_SIZE,
  ellipse     = TRUE,       # 是否画 95% 置信椭圆
  ellipse_level = 0.95,
  ellipse_alpha = 0.2,
  n_pc_axis   = 10,         # 碎石图/方差解释取前几个 PC
  width, height, dpi
)
```

返回 `ggplot` 对象；方差解释数据在 `attr(p, "pve")`。

### PCA_Scree_Plot()

```r
PCA_Scree_Plot(eigenval_path, out_file = NULL, n_pc = 10, width, height, dpi)
```

### PCA_3D_Plot()

```r
PCA_3D_Plot(
  eigenvec_path, eigenval_path, group_path,
  group_col = "Breed", sample_col = "SampleID",
  out_file = NULL,
  colors = NULL,
  drop_first_row = FALSE,   # eigenvec 若带表头行才设 TRUE
  legend_pos = "topright",
  legend_cex = 1,
  width = 8, height = 6, dpi = 300
)
```

### Admixture_Plot()

```r
Admixture_Plot(
  sname,                    # 含 breed + sampleID 的 csv
  file_prefix,              # Q 文件前缀，如 "admixture/bestqc."
  out,                      # 输出前缀（不含后缀）
  k_range   = 2:9,
  height = 18, width = 14,
  cexlab = 3,
  out_type  = c("pdf", "png", "none"),
  liner = 3, adj = -0.25,
  mfrow = 6,                # 与原脚本一致；见「已知坑 3」
  font_family = "",         # PDF 中文字体，Windows 可试 "GB1"
  breed_col = "breed", sample_col = "sampleID",
  dpi = 900, png_width, png_height
)
```

Q 文件命名规则：`paste0(file_prefix, K, ".Q")`，即 `bestqc.3.Q` 这种形式。

### FST_Manhattan_Plot()

```r
FST_Manhattan_Plot(
  fst_file,                 # 4 列以上：SNP, CHR, BP, 值...
  out_dir = NULL,           # 输出目录；NULL = 当前目录
  memo    = "FST",          # 输出文件名标记
  plot_type = "m",
  ylab = "FST", LOG10 = FALSE,
  cex = 0.5, signal_cex = 0.5,
  threshold = 0.3, threshold_col = "red",
  threshold_lty = 1, threshold_lwd = 1,
  width = 14, height = 6,
  file = "pdf",             # "pdf" / "png" / "jpg" / "tiff"
  col = PALETTE_CHR
)
```

### XPCLR_Manhattan_Plot()

```r
XPCLR_Manhattan_Plot(
  xpclr_file,               # 制表符分隔，带表头
  out_file   = NULL,
  chr_col    = "chrom",
  start_col  = "start",
  end_col    = "stop",
  value_col  = "xpclr_norm",
  chr_levels = NULL,        # 染色体排序，默认 1..33 + "Z"
  keep_axis  = c(...),      # 只在 x 轴标出这些染色体
  colors     = NULL,        # NULL = PALETTE_CHR8
  highlight_chr = "Z",      # 强制用第一色的染色体
  threshold_q   = 0.99,     # 分位虚线
  point_size = 1, point_alpha = 0.8,
  title, xlab, ylab,
  top_pct = NULL,           # 如 0.05 = 同时导出前 5% 位点
  top_out = NULL,
  width = 12, height = 6, dpi = 300
)
```

### XPEHH_Manhattan_Plot()

```r
XPEHH_Manhattan_Plot(
  xpehh_file,               # 8 列无表头，第 1 列是 "chr:pos"
  out_file    = NULL,
  col_names   = c("pos","gpos","p1","ihh1","p2","ihh2","xpehh_raw","xpehh_norm"),
  value_col   = "xpehh_norm",
  show_chr_labels = c(...), # 只在 x 轴标出这些染色体
  colors      = NULL,       # NULL = PALETTE_CHR8
  top_pct     = 0.05,       # 上、下各取 5% 作双尾阈值
  add_threshold = TRUE,
  point_size  = 0.5, point_alpha = 0.7,
  xlab = "Chr", ylab = "XP-EHH",
  width = 14, height = 7, dpi = 300
)
```

### Venn_Plot()

```r
Venn_Plot(
  sets,                     # 带名字的 list，每个元素是字符向量
  out_file    = NULL,       # 输出 png / pdf；NULL = 画到当前设备
  colors      = NULL,       # NULL = PALETTE_VENN（Venny 的红黄绿蓝）
  alpha       = 0.4,        # 填充透明度
  title       = NULL,
  cex         = 1.6,        # 区域数字字号
  cat_cex     = 1.6,        # 集合名字号
  line_width  = 2,
  font_family = "sans",     # "sans" / "serif" / "mono"
  detail_out  = NULL,       # 可选：每个区域的元素明细写成 tsv
  width = 8, height = 8, dpi = 300
)
```

```r
Venn_Plot(list(`List 1` = genesA, `List 2` = genesB, `List 3` = genesC),
         out_file = "out/venn.png")
```

返回 `invisible(list(plot, counts, regions))` —— `counts` 是各区域元素数，
`regions` 是各区域的元素名。

**注意**：返回的是 grid grob，不是 ggplot 对象，**用不了 `save_popgen()`**，
函数内部自己开关设备。

### LD_Decay_Plot()

```r
LD_Decay_Plot(
  files,                    # 命名向量 c("25R" = "25R.stat.gz", ...)，名字就是图例名
  pop_names  = NULL,        # files 没 names 时用它指定
  out_file   = NULL,
  bin_size   = NULL,        # NULL = 不平滑；给数值 = 按该距离分箱取均值
  xlim       = c(0, 500),
  ylim       = c(0, 1),
  breaks_x   = seq(0, 500, 100),
  title      = NULL,
  xlab       = "Distance (bp)",
  ylab       = expression(LD~(r^2)),
  colors     = NULL,
  linewidth  = 1,
  add_smooth = FALSE,       # 是否再叠一层 geom_smooth
  smooth_span = 0.3,
  skip       = 1,           # 跳过文件开头几行（PLINK 表头）
  width, height, dpi
)
```

输入是 PLINK `--r2` 输出的 `*.stat.gz`：第一行表头，之后每行两列 `Distance  R2`。

### SNP_Density_Plot()

```r
SNP_Density_Plot(
  snp_file,                 # 两列：CHR  BP（不含表头）
  out_dir  = NULL,
  memo     = "SNP_density",
  bin_size = 1e6,           # 密度窗口，1e6 = 1Mb
  col      = c("darkgreen", "yellow", "red"),   # 低 -> 中 -> 高
  file     = "pdf",
  width    = 14, height = 6, dpi = 300,
  chr_labels = NULL,        # 自定义染色体标签
  verbose  = TRUE
)
```

本质是 `CMplot(plot.type = "d")`。输出文件名形如 `Marker_Density.<memo>.pdf`。

### Tree_Pie_Plot()

```r
Tree_Pie_Plot(
  tree_file,                # newick 树
  pie_file,                 # csv：Tip + 各祖先成分列，行和约为 1
  out_file    = NULL,
  tip_col     = "Tip",      # 与树叶节点名对应的列
  pie_cols    = NULL,       # NULL = 除 tip_col 外的全部列
  layout      = "rectangular",
  tiplab      = TRUE,
  tiplab_size = 2.5,
  pie_radius  = 0.45,       # 饼半径，0.5 = 正好填满样本间距
  pie_nseg    = 36,         # 扇形分段数
  colors      = NULL,
  legend      = TRUE,
  width = 7, height = 12, dpi = 300
)
```

依赖 `ggtree` / `treeio`（Bioconductor）+ `patchwork`（CRAN），未安装会直接报错并给出安装命令。

**实现方式**：不是常见的 `ggtreeExtra::geom_fruit`，而是 ggtree 画树 + 手工算扇形多边形
（`geom_polygon`）+ patchwork 横向拼接。原因是 ggtreeExtra 1.16 的 `is.layer()` 会用空参数
试调用 geom，`geom_scatterpie` 的 `cols`/`data` 是必需参数直接崩；且 scatterpie 的 geom 在
ggplot2 3.5.2 下还会因 `from_theme()` 报错。详见「已知坑 8」。

---

## 样式常量（全部定义在 `R/theme_popgen.R`）

| 常量 | 默认值 | 含义 |
|---|---|---|
| `PALETTE_POPGEN` | 8 色（红蓝绿橙紫棕青粉） | 群体配色，前三位刻意贴近红/蓝/绿 |
| `PALETTE_CHR` | 6 色 | Fst 曼哈顿图染色体交替色 |
| `PALETTE_CHR8` | 8 色 | XP-CLR / XP-EHH 多染色体循环色 |
| `PALETTE_VENN` | 红黄绿蓝 | 韦恩图配色（Venny 默认） |
| `POPGEN_POINT_SIZE` | 2 | 散点大小 |
| `POPGEN_BASE_SIZE` | 11 | theme 基准字号 |
| `POPGEN_FONT` | `""` | 字体族，留空 = 系统默认 |
| `POPGEN_WIDTH` / `HEIGHT` | 10 / 4.5 | 默认导出尺寸（英寸） |
| `POPGEN_DPI` | 500 | 默认导出分辨率 |

辅助函数：

- `popgen_colors(n)` —— 取 n 个群体色（超过 8 个循环取用）
- `popgen_colors_ramp(n)` —— 取 n 个可区分色（n ≤ 8 用离散色，n > 8 才插值）
- `theme_popgen(base_size, base_family)` —— 统一 theme（以 `theme_light()` 为底）
- `theme_popgen_minimal(...)` —— 白底无网格（LD 衰减、XP-CLR 用）
- `theme_popgen_classic(...)` —— 白底带坐标轴线（XP-EHH 用）
- `save_popgen(plot, file, width, height, dpi)` —— 统一导出
- `save_popgen_both(plot, prefix, ...)` —— 同时导出 PDF + PNG

---

## 数据格式约定

### 分组表（PCA / ADMIXTURE 共用）

```csv
SampleID,Breed
PopA_01,PopA
PopA_02,PopA
```

- 列名统一用 `SampleID` + `Breed`（原始脚本里出现过 `Group` / `breed` / `sampleID` 三种写法，
  现在通过 `group_col` / `sample_col` / `breed_col` 参数兼容，新代码请统一用 `SampleID` + `Breed`）
- PCA 通过 `SampleID` 做 `match()` 对齐，**不依赖行顺序**

### PLINK PCA

- `*.eigenvec`：`FID IID PC1 PC2 ...`，空格分隔，无表头
- `*.eigenval`：一列数值，每行一个 PC 的特征值

### ADMIXTURE

- `*.Q`：每行一个样本，K 列，行和为 1，制表符分隔
- 文件命名：`<prefix><K>.Q`，如 `bestqc.3.Q`

### Fst（CMplot 输入）

```tsv
SNP	CHR	BP	FST
chr1_00001	1	15234	0.0213
```

至少 4 列：SNP 名、染色体、物理位置、值。第 4 列及以后每列是一条性状。

---

## 加一张新图的标准流程

1. 在 `R/` 下**对应的子文件夹**里新建 `plot_xxx.R`（群体结构 / 遗传多样性 / 选择信号），
   写一个函数，**第一个参数是路径，最后几个参数是 width/height/dpi**
2. 颜色用 `popgen_colors()` / `popgen_colors_ramp()`，尺寸用 `POPGEN_*` 常量
3. 导出用 `save_popgen()`，不要直接 `ggsave` / `pdf()` / `png()`
4. **不需要改 `load_all.R`** —— 它会递归扫描 `R/` 下所有 `.R` 文件自动加载
   （`theme_popgen.R` 会被强制排在最前）
5. 在 `tools/run_demo.R` 里加一段调用，输出到 `demo_out/`
6. 跑 `Rscript tools/run_demo.R` 确认图正常
7. 在 `README.md` 的图库里加一行

> 文件名建议全小写 + 下划线（`LD_Decay_Plot.R`），不要用中文和空格。
> R 本身能处理中文文件名，但 URL 引用和 git 历史里会很别扭。

---

## 已知坑

### 1. Windows 上 R 的 locale 必须是 UTF-8

Windows R 默认 `LC_CTYPE` 可能是 `C` 或 GBK，此时 `source(encoding = "UTF-8")`
遇到中文注释会报 `invalid input found on input connection`，**并静默截断文件**——
后面的函数定义会全部丢失，报错却是莫名其妙的「object not found」。

`R/load_all.R` 和 `tools/run_demo.R` 开头已经处理：

```r
for (.loc in c(".UTF-8", "en_US.UTF-8", "Chinese.UTF-8")) {
  if (!is.na(suppressWarnings(Sys.setlocale("LC_CTYPE", .loc)))) break
}
```

自己写新脚本时也要复制这三行。

### 2. CMplot 版本差异

`CMplot` 4.0 起把输出文件标记参数从 `memo` 改成了 `file.name`。
本仓库按 4.5.1 写。用旧版会报 `unused argument (memo = ...)`。

### 3. ADMIXTURE 的 mfrow 与面板数

`mfrow` 默认 `6`，配合 `height = 18` 和 `mar = c(10, 4, 3, 10)`——这是原脚本调好的比例。

如果 `k_range` 超过 6 个（如 `2:9` 共 8 个），PDF 会分成两页。
想把 8 个一页画完，必须**同时**加大 `height`，否则边距会把面板挤爆，
报 `figure margins too large`。建议：

```r
Admixture_Plot(..., k_range = 2:9, mfrow = 8, height = 26)
```

### 4. `family = "GB1"` 不可移植

`GB1` 是 Windows 专有中文字体名，换 Linux/macOS 直接报错。
默认不指定；确实需要中文字体时用 `font_family` 参数显式传。

### 5. 3D 图的图例位置

`scatterplot3d` 会重置 `par(mar)`。图例要放到 3D 框外，必须通过
`scatterplot3d(..., mar = c(5, 4, 4, 10))` 留右边距，不能只在外层 `par(mar=...)`。

### 6. 绘图设备要收干净

两个坑连着出：

- `on.exit` 里必须**先 `par(op)` 再 `dev.off()`**。反过来的话 `par()` 会在设备
  已关闭后执行，R 当场另开一个默认设备，结果是凭空多一个设备 + 工作目录多一个
  `Rplots.pdf`；且文件被设备占用时 Windows 下 `unlink()` 是静默失败的。
- `Rscript` 下顶层表达式的可见返回值会**自动打印**。返回 ggplot 的函数不赋值，
  就会 print，同样触发默认设备。用 `x <- 函数()` 或 `invisible(函数())` 接住。

自检：`names(grDevices::dev.list())` 正常应为 `NULL`。

### 7. 装 ggtree 系列的坑

`BiocManager::install()` 会先访问 `bioconductor.org/config.yaml` 校验版本。
网络不稳时报 `Bioconductor version cannot be validated; no internet connection?`，
但此时 **packages 目录本身可能是通的**，只是那个 yaml 超时。
国内镜像（清华/阿里/USTC）的 `packages/3.20/bioc` 路径目前多数 404。

可用的做法是绕过版本校验，直接从官方 packages 目录装：

```r
install.packages(
  c("treeio", "ggtree", "ggtreeExtra"),
  repos = c("https://bioconductor.org/packages/3.20/bioc",
            "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"),
  dependencies = c("Depends", "Imports")
)
```

`3.20` 对应 R 4.4。其它 R 版本要换对应编号（3.19↔R4.3，3.21↔R4.5）。

### 8. ggtreeExtra + scatterpie 在当前版本下不能一起用

想画「树上挂饼图」时，网上教程几乎都是：

```r
geom_fruit(geom = geom_scatterpie, ...)   # ❌ 本机跑不通
```

实测（ggtreeExtra 1.16.0 + scatterpie 0.2.6 + ggplot2 3.5.2）会在两处失败：

1. `ggtreeExtra:::.convert_to_name()` 里 `is.layer(x)` 的实现是
   `inherits(do.call(x, list()), "Layer")` —— **用空参数试调用 geom**。
   `geom_scatterpie` 的 `cols`（0.2.6）或 `data`（0.2.2）是必需参数，空调用直接报
   `argument "cols" is missing`。而且 scatterpie 返回对象的 class 是
   `"layer_scatterpie"`，不是 `"Layer"`，`is.layer()` 恒为 FALSE。
2. 即便绕过第 1 条，scatterpie 的 geom 在 ggplot2 3.5.x 下还会报
   `could not find function "from_theme"`。

**本仓库的处理**：`Tree_Pie_Plot()` 不用这两个包，改成 ggtree 画树 + 手工算扇形
多边形（`geom_polygon`）+ patchwork 横向拼接，只依赖 `ggtree` / `treeio` / `patchwork`。

自检代码：

```r
ggtreeExtra:::is.layer("geom_col")          # TRUE
ggtreeExtra:::is.layer("geom_scatterpie")   # 报错或 FALSE
```

### 9. 用 VennDiagram 画韦恩图的两个坑

第一，**必须显式关掉欧拉模式**：

```r
draw.triple.venn(..., euler.d = TRUE, scaled = TRUE)   # 默认：按集合大小画欧拉图，形状变形
draw.triple.venn(..., euler.d = FALSE, scaled = FALSE) # Venny 那种标准三角/四椭圆布局
```

第二，**`draw.pairwise.venn` 在 `cat.dist` 传长度为 2 的向量时会出事**：

```r
draw.pairwise.venn(area1, area2, cross.area, fill = c("red","yellow"),
                   cat.dist = c(0.05, 0.05))   # area1 被画到右边、fill 跟着错位、
                                               # 两个集合名还叠在正上方
```

实测同一份数据，`cat.dist` 传标量时位置是对的，传向量就错位。
这个行为比较隐蔽，所以 `Venn_Plot()` 里 **2 集合不走 VennDiagram，直接用自己的
`.venn2Grob()`（grid 原语画两个圆 + 三个数字）**；3 / 4 集合仍用 VennDiagram。

### 6. 绘图设备要按正确顺序收拾

**`on.exit` 里必须先 `par(op)` 再 `dev.off()`**：

```r
on.exit({
  graphics::par(op)                 # 先恢复参数
  if (dev_opened) grDevices::dev.off()   # 再关设备
}, add = TRUE)
```

顺序反了，`par(op)` 会在设备已关闭之后执行，R 当场另开一个默认 pdf 设备，
结果是调用方凭空多一个设备、工作目录多一个 `Rplots.pdf`。
而且此时文件被设备占着，Windows 下 `unlink()` 会静默失败删不掉。

**`Rscript` 下顶层表达式的返回值会自动打印。**
返回 ggplot 对象的函数（`PCA_2D_Plot` / `PCA_Scree_Plot`）如果不赋值，
返回值就被 print → 触发默认设备 → 又是 `Rplots.pdf`：

```r
PCA_2D_Plot(...)              # 错：返回值被打印
fig <- PCA_2D_Plot(...)       # 对
invisible(PCA_2D_Plot(...))   # 也可以
```

Debug 时用这个看有没有漏关：

```r
names(grDevices::dev.list())   # 正常应该是 NULL
```
