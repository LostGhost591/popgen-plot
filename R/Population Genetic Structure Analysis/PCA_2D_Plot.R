#### ============================================================
#### PCA_2D_Plot.R —— 二维 PCA 散点图（带 95% 置信椭圆）
####
#### 画法完全沿用原始脚本：
####   geom_point(size = 2) + coord_equal() + stat_ellipse(type = "norm", geom = "polygon")
#### 唯一改动是把路径、配色、尺寸抽成参数，样式统一从 theme_popgen.R 读。
####
#### 用法：
####   source("R/load_all.R")
####   PCA_2D_Plot("demo_data/demo.eigenvec", "demo_data/demo.eigenval",
####               "demo_data/demo.group.csv", out_file = "demo_out/pca_2d.png")
#### ============================================================

PCA_2D_Plot <- function(eigenvec_path,
                        eigenval_path,
                        group_path,
                        group_col   = "Breed",     # group.csv 里装分组的列
                        sample_col  = "SampleID",  # group.csv 里装样本名的列
                        out_file    = NULL,        # NULL = 不出图，只返回对象
                        colors      = NULL,        # NULL = 用统一色板
                        point_size  = POPGEN_POINT_SIZE,
                        ellipse     = TRUE,
                        ellipse_level = 0.95,
                        ellipse_alpha = 0.2,
                        n_pc_axis   = 10,
                        width       = POPGEN_WIDTH,
                        height      = POPGEN_HEIGHT,
                        dpi         = POPGEN_DPI) {

  #### ---- 读取 PCA 结果 ----
  pca <- readr::read_table(eigenvec_path, col_names = TRUE, show_col_types = FALSE)
  eigenval <- scan(eigenval_path, quiet = TRUE)

  #### ---- 整理列名：去掉 FID 列，第一列叫 ind，其余叫 PC1..PCn ----
  pca <- pca[, -1]
  names(pca)[1] <- "ind"
  names(pca)[2:ncol(pca)] <- paste0("PC", 1:(ncol(pca) - 1))

  #### ---- 方差解释百分比 ----
  n_pc <- min(n_pc_axis, length(eigenval))
  pve <- data.frame(PC  = seq_len(n_pc),
                    pve = eigenval[seq_len(n_pc)] / sum(eigenval) * 100)

  #### ---- 按 SampleID 对齐分组（比按行号绑定安全）----
  grp <- as.data.frame(data.table::fread(group_path))
  if (!sample_col %in% names(grp)) stop("group 文件里找不到列: ", sample_col)
  if (!group_col  %in% names(grp)) stop("group 文件里找不到列: ", group_col)

  pca$group_var <- grp[[group_col]][match(pca$ind, grp[[sample_col]])]
  if (anyNA(pca$group_var)) {
    miss <- pca$ind[is.na(pca$group_var)]
    stop("以下样本在 group 文件里没有匹配（前 10 个）: ",
         paste(utils::head(miss, 10), collapse = ", "))
  }
  pca$group_var <- factor(pca$group_var, levels = unique(pca$group_var))

  #### ---- 配色：默认走统一色板 ----
  lv <- levels(pca$group_var)
  if (is.null(colors)) colors <- popgen_colors(length(lv))
  colors <- stats::setNames(rep(colors, length.out = length(lv)), lv)

  #### ---- 作图（画法与原脚本一致）----
  p <- ggplot(pca, aes(x = PC1, y = PC2,
                       colour = group_var, shape = group_var)) +
    geom_point(size = point_size) +
    scale_colour_manual(values = colors) +
    coord_equal() +
    theme_popgen() +
    labs(colour = group_col, shape = group_col) +
    xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) +
    ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)"))

  if (isTRUE(ellipse)) {
    p <- p +
      stat_ellipse(aes(fill = group_var), type = "norm", geom = "polygon",
                   alpha = ellipse_alpha, colour = NA, level = ellipse_level) +
      scale_fill_manual(values = colors, guide = "none")
  }

  #### ---- 导出 + 附带方差解释数据，方便调用方取用 ----
  attr(p, "pve") <- pve
  if (!is.null(out_file)) save_popgen(p, out_file, width, height, dpi)
  p
}

#### 碎石图：每种 PC 的方差解释占比（原脚本里那张柱状图）
PCA_Scree_Plot <- function(eigenval_path, out_file = NULL,
                           n_pc = 10, width = POPGEN_WIDTH,
                           height = POPGEN_HEIGHT, dpi = POPGEN_DPI) {
  eigenval <- scan(eigenval_path, quiet = TRUE)
  n_pc <- min(n_pc, length(eigenval))
  pve <- data.frame(PC  = seq_len(n_pc),
                    pve = eigenval[seq_len(n_pc)] / sum(eigenval) * 100)

  p <- ggplot(pve, aes(x = PC, y = pve)) +
    geom_bar(stat = "identity", fill = PALETTE_POPGEN[2]) +
    scale_x_continuous(breaks = pve$PC) +
    theme_popgen() +
    xlab("Principal component") +
    ylab("Percentage variance explained")

  if (!is.null(out_file)) save_popgen(p, out_file, width, height, dpi)
  p
}
