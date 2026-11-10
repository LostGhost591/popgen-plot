#### ============================================================
#### PCA_3D_Plot.R —— 三维 PCA 散点图
####
#### 画法沿用原始脚本：scatterplot3d + 按群体分配颜色和 pch 形状。
#### 相对原始脚本做了两处修正（都可用参数关掉）：
####   1. drop_first_row 默认 FALSE —— 原脚本 `eigenvec[-1, ]` 会静默丢掉第一个样本
####   2. 只画一次 legend —— 原脚本写了三遍 legend，图上会叠两层
#### ============================================================

PCA_3D_Plot <- function(eigenvec_path,
                        eigenval_path,
                        group_path,
                        group_col = "Breed",
                        sample_col = "SampleID",
                        out_file = NULL,
                        colors = NULL,
                        drop_first_row = FALSE,
                        legend_pos = "topright",
                        legend_cex = 1,
                        width = 8,
                        height = 6,
                        dpi = 300) {

  if (!requireNamespace("scatterplot3d", quietly = TRUE)) {
    stop("缺少包 scatterplot3d，请先 install.packages('scatterplot3d')")
  }

  #### ---- 读取 ----
  eigenvec <- utils::read.table(eigenvec_path, header = FALSE)
  eigenval <- utils::read.table(eigenval_path, header = FALSE)
  group_info <- as.data.frame(data.table::fread(group_path))

  #### ---- 整理列名 ----
  eigenvec <- eigenvec[, -1]
  names(eigenvec)[1] <- "SampleID"
  names(eigenvec)[2:ncol(eigenvec)] <- paste0("PC", 1:(ncol(eigenvec) - 1))

  #### 若 eigenvec 第一行是 PLINK 输出的表头，设 drop_first_row = TRUE
  if (isTRUE(drop_first_row)) eigenvec <- eigenvec[-1, ]

  #### 统一转成数值，避免 read.table 把整列读成字符
  for (k in 2:ncol(eigenvec)) eigenvec[[k]] <- as.numeric(eigenvec[[k]])

  #### ---- 方差解释百分比 ----
  variance_percent <- eigenval$V1 / sum(eigenval$V1) * 100

  #### ---- 合并分组 ----
  if (!sample_col %in% names(group_info)) stop("group 文件里找不到列: ", sample_col)
  if (!group_col  %in% names(group_info)) stop("group 文件里找不到列: ", group_col)

  pca_data <- merge(eigenvec,
                    group_info[, c(sample_col, group_col)],
                    by.x = "SampleID", by.y = sample_col,
                    all.x = TRUE, sort = FALSE)
  if (anyNA(pca_data[[group_col]])) {
    stop("有样本在 group 文件里匹配不到分组，请检查 SampleID")
  }

  #### ---- 颜色 / 形状：颜色走统一色板，形状沿用 16, 17, 18... ----
  groups <- unique(pca_data[[group_col]])
  if (is.null(colors)) colors <- popgen_colors(length(groups))
  color_map <- stats::setNames(rep(colors, length.out = length(groups)), groups)

  shapes    <- 16:(16 + length(groups) - 1)
  shape_map <- stats::setNames(shapes, groups)

  point_colors <- unname(color_map[as.character(pca_data[[group_col]])])
  point_shapes <- unname(shape_map[as.character(pca_data[[group_col]])])

  #### ---- 画图函数 ----
  draw <- function() {
    op <- graphics::par(xpd = TRUE)
    on.exit(graphics::par(op), add = TRUE)

    scatterplot3d::scatterplot3d(
      pca_data$PC1, pca_data$PC2, pca_data$PC3,
      color = point_colors,
      pch   = point_shapes,
      mar   = c(5, 4, 4, 10),   # 右边留宽，给框外的图例让位
      xlab  = paste0("PC1 (", round(variance_percent[1], 2), "%)"),
      ylab  = paste0("PC2 (", round(variance_percent[2], 2), "%)"),
      zlab  = paste0("PC3 (", round(variance_percent[3], 2), "%)")
    )

    ## 图例放到 3D 框外右侧，避免压住框线
    graphics::legend(legend_pos,
                     inset  = c(-0.32, 0.02),
                     legend = groups,
                     col    = unname(color_map),
                     pch    = unname(shape_map),
                     bty    = "n",
                     cex    = legend_cex,
                     title  = group_col,
                     xpd    = NA)
  }

  #### ---- 导出 ----
  if (!is.null(out_file)) {
    dir.create(dirname(out_file), showWarnings = FALSE, recursive = TRUE)
    if (grepl("\\.pdf$", out_file, ignore.case = TRUE)) {
      grDevices::pdf(out_file, width = width, height = height)
    } else {
      grDevices::png(out_file, width = width, height = height,
                     units = "in", res = dpi)
    }
    on.exit(grDevices::dev.off(), add = TRUE)
  }

  draw()
  invisible(pca_data)
}
