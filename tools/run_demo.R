#### ============================================================
#### run_demo.R —— 跑一遍全部作图函数，生成 demo_out/ 下的示例图
####
#### 用途有两个：
####   1. 人看：README 里贴的就是这些图
####   2. AI 看：改完代码后跑这个脚本，图不对就是改坏了
####
#### 运行（在仓库根目录）：Rscript tools/run_demo.R
#### ============================================================

## 先切 UTF-8 locale，否则 source(encoding="UTF-8") 会被中文注释截断
## （原因详见 R/load_all.R 的注释）
for (.loc in c(".UTF-8", "en_US.UTF-8", "Chinese.UTF-8")) {
  if (!is.na(suppressWarnings(Sys.setlocale("LC_CTYPE", .loc)))) break
}
rm(.loc)

source("R/load_all.R")

dir.create("demo_out", showWarnings = FALSE)

## 注意：Rscript 下顶层表达式的可见返回值会被自动打印。返回 ggplot 对象的函数
## 一旦被打印，R 会开默认 pdf 设备并生成 Rplots.pdf。
## 所以凡是返回非 NULL 的调用，一律用变量接住；不接就包 invisible()。
## （返回 invisible(NULL) 的函数无此问题）

#### ---------- 1. 二维 PCA ----------
fig_pca2d <- PCA_2D_Plot(
  eigenvec_path = "demo_data/demo.eigenvec",
  eigenval_path = "demo_data/demo.eigenval",
  group_path    = "demo_data/demo.group.csv",
  out_file      = "demo_out/pca_2d.png"
)

#### ---------- 2. PCA 碎石图 ----------
fig_scree <- PCA_Scree_Plot(
  eigenval_path = "demo_data/demo.eigenval",
  out_file      = "demo_out/pca_scree.png"
)

#### ---------- 3. 三维 PCA ----------
invisible(PCA_3D_Plot(
  eigenvec_path = "demo_data/demo.eigenvec",
  eigenval_path = "demo_data/demo.eigenval",
  group_path    = "demo_data/demo.group.csv",
  out_file      = "demo_out/pca_3d.png"
))

#### ---------- 4. ADMIXTURE ----------
## 正式版 PDF（投稿用矢量图）
Admixture_Plot(
  sname       = "demo_data/demo.admixture.samples.csv",
  file_prefix = "demo_data/demo.K.",
  out         = "demo_out/admixture",
  out_type    = "pdf"
)

## PNG 预览版（仅供 README 在网页里展示，图与 PDF 一致）
Admixture_Plot(
  sname       = "demo_data/demo.admixture.samples.csv",
  file_prefix = "demo_data/demo.K.",
  out         = "demo_out/admixture_preview",
  out_type    = "png",
  mfrow       = 8,
  png_width   = 1300, png_height = 2600, dpi = 120
)

#### ---------- 5. Fst 曼哈顿图 ----------
## 正式版 PDF
FST_Manhattan_Plot(
  fst_file = "demo_data/demo.fst.txt",
  out_dir  = "demo_out",
  memo     = "FST"
)

## PNG 预览版（仅供 README 展示）
FST_Manhattan_Plot(
  fst_file = "demo_data/demo.fst.txt",
  out_dir  = "demo_out",
  memo     = "FST",
  file     = "png"
)

#### ---------- 6. LD 衰减 ----------
fig_ld <- LD_Decay_Plot(
  files = c("25R"      = "demo_data/demo.ld.25R.stat.gz",
            "49W"      = "demo_data/demo.ld.49W.stat.gz",
            "80aixiao" = "demo_data/demo.ld.80aixiao.stat.gz"),
  out_file = "demo_out/ld_decay.png",
  title    = "LD Decay"
)

## 同一份数据按 20 bp 分箱取均值，曲线更顺（原脚本的第二种画法）
fig_ld_bin <- LD_Decay_Plot(
  files = c("25R"      = "demo_data/demo.ld.25R.stat.gz",
            "49W"      = "demo_data/demo.ld.49W.stat.gz",
            "80aixiao" = "demo_data/demo.ld.80aixiao.stat.gz"),
  bin_size = 20,
  out_file = "demo_out/ld_decay_smooth.png",
  title    = "LD Decay (20 bp bins)"
)

#### ---------- 7. SNP 密度 ----------
SNP_Density_Plot("demo_data/demo.snp_pos.txt",
               out_dir = "demo_out", memo = "SNP_density", file = "pdf")
SNP_Density_Plot("demo_data/demo.snp_pos.txt",
               out_dir = "demo_out", memo = "SNP_density", file = "png")

#### ---------- 8. 进化树饼图 ----------
## 依赖 ggtree / treeio / patchwork；没装就跳过，不影响其它图
if (all(vapply(c("ggtree", "treeio", "patchwork"),
               requireNamespace, logical(1), quietly = TRUE))) {
  fig_tree <- Tree_Pie_Plot("demo_data/demo.tree.nwk",
                          "demo_data/demo.tree.pie.csv",
                          out_file = "demo_out/tree_pie.png",
                          width = 7, height = 12, dpi = 200)
} else {
  message("跳过进化树饼图：未安装 ggtree / treeio / patchwork")
}

#### ---------- 9. XP-CLR 曼哈顿图 ----------
fig_xpclr <- XPCLR_Manhattan_Plot(
  xpclr_file = "demo_data/demo.xpclr.txt",
  out_file   = "demo_out/xpclr_manhattan.png",
  top_pct    = 0.05,
  top_out    = "demo_out/xpclr_top5.txt",
  dpi        = 300
)

#### ---------- 10. XP-EHH 曼哈顿图 ----------
fig_xpehh <- XPEHH_Manhattan_Plot(
  xpehh_file = "demo_data/demo.xpehh.txt",
  out_file   = "demo_out/xpehh_manhattan.png",
  dpi        = 300
)

#### ---------- 11. 韦恩图 ----------
## 读进来后把每列的空位（NA / 空串）清掉，再交给 Venn_Plot
vdat  <- utils::read.csv("demo_data/demo.venn.csv", stringsAsFactors = FALSE)
vsets <- lapply(vdat, function(x) {
  x <- as.character(x)
  x[!is.na(x) & nzchar(x)]
})
names(vsets) <- names(vdat)

## 2 集合
Venn_Plot(vsets[1:2],
         out_file = "demo_out/venn_2sets.png",
         title    = "Shared genes (2 lists)",
         width = 7, height = 7, dpi = 200)

## 3 集合
Venn_Plot(vsets[1:3],
         out_file = "demo_out/venn_3sets.png",
         title    = "Shared genes (3 lists)",
         width = 8, height = 8, dpi = 200)

## 4 集合，顺带导出每个区域的元素明细
fig_venn <- Venn_Plot(vsets,
                     out_file   = "demo_out/venn_4sets.png",
                     title      = "Shared genes (4 lists)",
                     detail_out = "demo_out/venn_regions.tsv",
                     width = 8, height = 8, dpi = 200)

rm(fig_pca2d, fig_scree, fig_ld, fig_ld_bin, fig_xpclr, fig_xpehh, fig_venn)
if (exists("fig_tree")) rm(fig_tree)

## 收尾：先关掉所有遗留设备，再清理 Rscript 下可能产生的默认设备文件。
## 顺序不能反 —— 设备开着时 Windows 会锁住 Rplots.pdf，unlink 静默失败。
open_dev <- names(grDevices::dev.list())
if (length(open_dev)) cat("遗留设备:", paste(open_dev, collapse = ", "), "\n")
grDevices::graphics.off()
if (file.exists("Rplots.pdf")) unlink("Rplots.pdf")

cat("\n全部示例图已生成到 demo_out/\n")
print(list.files("demo_out"))
