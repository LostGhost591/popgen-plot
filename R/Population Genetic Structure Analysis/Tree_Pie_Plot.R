#### ============================================================
#### Tree_Pie_Plot.R —— 系统发育树 + 群体祖先成分饼图
####
#### 实现方式（零兼容风险）：
####   · ggtree 画树
####   · 手工把每个样本的成分比例算成扇形多边形，用 geom_polygon 在独立面板里画
####   · patchwork 横向拼两张图，y 轴范围一致所以节点自动对齐
####
#### 为什么不用常见的 ggtreeExtra::geom_fruit + scatterpie：
####   ggtreeExtra 1.16 的 .convert_to_name() 内部用 `is.layer(x)` 空参数试调用
####   geom，而 geom_scatterpie 的 cols/data 是必需参数，直接报错；且 scatterpie
####   返回对象 class 是 "layer_scatterpie" 而非 "Layer"，is.layer() 恒为 FALSE。
####   即便绕过，scatterpie 的 geom 在 ggplot2 3.5.2 下还会因 from_theme() 报错。
####   所以这里干脆手算扇形，只依赖 ggplot2 + ggtree + patchwork。
####
#### 依赖：ggtree、treeio（Bioconductor）、patchwork（CRAN）
####
#### 输入：
####   1) newick 树文件，如 iqtree / raxml 的 .treefile
####   2) 饼图数据（csv，第一列是 tip 名，与树叶节点名完全一致）：
####        Tip,PopA,PopB,PopC
####        Sample01,0.80,0.15,0.05
####  每一行是一个样本在各祖先成分上的比例，行和应约为 1。
#### ============================================================

## 内部工具：把成分比例转成扇形多边形坐标
.makePiePolygons <- function(pie, cols, x0 = 1, r = 0.45, nseg = 36) {
  out <- vector("list", nrow(pie) * length(cols))
  idx <- 0L
  for (i in seq_len(nrow(pie))) {
    vals <- as.numeric(pie[i, cols])
    vals[is.na(vals)] <- 0
    if (sum(vals) <= 0) next
    props <- vals / sum(vals)
    ang   <- c(0, cumsum(props) * 2 * pi)   # 从 12 点方向顺时针起算
    for (k in seq_along(cols)) {
      if (props[k] <= 0) next
      t <- seq(ang[k], ang[k + 1], length.out = nseg)
      idx <- idx + 1L
      out[[idx]] <- data.frame(
        x     = x0 + r * sin(t),
        y     = pie$y[i] + r * cos(t),
        group = paste0(i, "_", k),
        fill  = cols[k],
        stringsAsFactors = FALSE
      )
    }
  }
  if (!idx) stop("所有样本的成分比例都是 0，无法画饼")
  do.call(rbind, out[seq_len(idx)])
}

Tree_Pie_Plot <- function(tree_file,
                        pie_file,
                        out_file    = NULL,
                        tip_col     = "Tip",      # 饼图数据里对应叶节点名的列
                        pie_cols    = NULL,       # 画哪几列做饼；NULL = 除 tip_col 外全部
                        layout      = "rectangular",
                        tiplab      = TRUE,
                        tiplab_size = 2.5,
                        pie_radius  = 0.45,       # 饼半径，0.5 = 正好填满样本间距
                        pie_nseg    = 36,         # 每个扇形分几段，越大越圆滑
                        colors      = NULL,       # NULL = 统一色板
                        legend      = TRUE,
                        width  = 7,
                        height = 12,
                        dpi    = 300) {

  need <- c("ggtree", "treeio", "patchwork")
  miss <- need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss)) {
    stop("缺少依赖包: ", paste(miss, collapse = ", "), "\n",
         "请先运行：\n",
         '  install.packages("patchwork")\n',
         '  install.packages(c("treeio", "ggtree"),\n',
         '                   repos = c("https://bioconductor.org/packages/3.20/bioc",\n',
         '                             "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))\n')
  }

  if (!file.exists(tree_file)) stop("找不到树文件: ", tree_file)
  if (!file.exists(pie_file))  stop("找不到饼图数据: ", pie_file)

  #### ---- 读树 ----
  tree <- treeio::read.tree(tree_file)
  n    <- length(tree$tip.label)

  #### ---- 读饼图数据 ----
  pie <- utils::read.csv(pie_file, stringsAsFactors = FALSE, check.names = FALSE)
  if (!tip_col %in% names(pie)) {
    stop("饼图数据里找不到列 '", tip_col, "'，现有列：", paste(names(pie), collapse = ", "))
  }
  if (is.null(pie_cols)) pie_cols <- setdiff(names(pie), tip_col)
  if (!length(pie_cols)) stop("没有可画饼的列")

  #### ---- 检查 tip 名是否对得上 ----
  missing_in_tree <- setdiff(pie[[tip_col]], tree$tip.label)
  if (length(missing_in_tree)) {
    stop("以下样本在树里找不到同名的叶节点（前 10 个）: ",
         paste(utils::head(missing_in_tree, 10), collapse = ", "))
  }

  #### ---- 叶节点编号：有根树的 tip node 固定为 1..Ntip，同时也是 y 坐标 ----
  pie$y <- match(pie[[tip_col]], tree$tip.label)

  #### ---- 配色 ----
  if (is.null(colors)) colors <- popgen_colors(length(pie_cols))
  colors <- stats::setNames(rep(colors, length.out = length(pie_cols)), pie_cols)

  #### ---- 左：树 ----
  p_tree <- ggtree::ggtree(tree, layout = layout, size = 0.6)
  if (isTRUE(tiplab)) {
    p_tree <- p_tree + ggtree::geom_tiplab(size = tiplab_size)
  }

  #### ---- 右：一列饼（独立面板，用 coord_fixed 保证是正圆）----
  poly <- .makePiePolygons(pie, pie_cols, x0 = 1, r = pie_radius, nseg = pie_nseg)

  p_pie <- ggplot(poly, aes(x = x, y = y, group = group, fill = fill)) +
    geom_polygon(colour = NA) +
    coord_fixed() +
    scale_y_reverse(limits = c(n + 0.5, 0.5)) +   # 与树同向：tip1 在最上面
    scale_x_continuous(limits = c(0.5, 1.5)) +
    scale_fill_manual(
      values = colors,
      name   = if (isTRUE(legend)) "Ancestry" else NULL,
      guide  = if (isTRUE(legend)) "legend" else "none"
    ) +
    theme_void() +
    theme(legend.position = if (isTRUE(legend)) "right" else "none")

  #### ---- 拼接 ----
  p <- p_tree + p_pie + patchwork::plot_layout(widths = c(3, 1))

  if (!is.null(out_file)) save_popgen(p, out_file, width, height, dpi)
  p
}
