#### ============================================================
#### LD_Decay_Plot.R —— LD 衰减曲线（多群体 r² 随物理距离衰减）
####
#### 画法沿用原始脚本：geom_line + 固定 x/y 范围 + theme_minimal 观感
####
#### 输入是 PLINK 的 LD 统计文件（--r2 输出，常为 *.stat.gz）：
####   第一行是注释/表头，之后每行两列：Distance  R2
####
#### 用法：
####   LD_Decay_Plot(
####     files = c("25R" = "25R.stat.gz", "49W" = "49W.stat.gz",
####               "80aixiao" = "80aixiao.stat.gz"),
####     out_file = "demo_out/ld_decay.png"
####   )
#### ============================================================

LD_Decay_Plot <- function(files,                    # 命名向量 c(群体名 = 文件路径)，或普通向量
                          pop_names  = NULL,        # files 没有 names 时用这个指定群体名
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
                          add_smooth = FALSE,       # 是否再叠一层 geom_smooth 顺滑
                          smooth_span = 0.3,
                          skip       = 1,           # 跳过文件开头几行（PLINK 的表头行数）
                          width      = POPGEN_WIDTH,
                          height     = POPGEN_HEIGHT,
                          dpi        = POPGEN_DPI) {

  if (!length(files)) stop("files 不能为空")

  #### ---- 群体名：命名向量 > pop_names > 文件名 ----
  nm <- names(files)
  if (is.null(nm) || all(!nzchar(nm))) nm <- pop_names
  if (is.null(nm) || length(nm) != length(files)) {
    nm <- tools::file_path_sans_ext(basename(files))
  }

  #### ---- 逐个读入并合并 ----
  ld_list <- lapply(seq_along(files), function(i) {
    if (!file.exists(files[i])) stop("找不到文件: ", files[i])
    d <- readr::read_tsv(
      files[i],
      skip       = skip,
      col_names  = c("Distance", "R2"),
      col_types  = "dd",
      show_col_types = FALSE
    )
    d$Pop <- nm[i]
    d[, c("Pop", "Distance", "R2")]
  })
  ld_all <- dplyr::bind_rows(ld_list)

  if (!nrow(ld_all)) stop("读到的 LD 数据是空的，检查 skip 和列格式是否正确")

  #### ---- 可选：按距离分箱取均值（对应原脚本的 bin_size = 20 那段）----
  if (!is.null(bin_size) && is.finite(bin_size) && bin_size > 0) {
    ld_all <- ld_all |>
      dplyr::mutate(bin = floor(Distance / bin_size) * bin_size + bin_size / 2) |>
      dplyr::group_by(Pop, bin) |>
      dplyr::summarise(R2 = mean(R2, na.rm = TRUE), .groups = "drop") |>
      dplyr::rename(Distance = bin)
    ## 最后一个 bin 的中心可能略超出 xlim 上界，先裁掉，免得 scale 报警告
    ld_all <- ld_all[ld_all$Distance <= max(xlim), ]
  }

  #### ---- 配色：默认走统一色板 ----
  pops <- unique(ld_all$Pop)
  if (is.null(colors)) colors <- popgen_colors(length(pops))
  colors <- stats::setNames(rep(colors, length.out = length(pops)), pops)

  #### ---- 作图 ----
  p <- ggplot(ld_all, aes(x = Distance, y = R2, colour = Pop)) +
    geom_line(linewidth = linewidth)

  if (isTRUE(add_smooth)) {
    p <- p + geom_smooth(se = FALSE, span = smooth_span, linewidth = linewidth * 1.2)
  }

  p <- p +
    scale_colour_manual(values = colors) +
    scale_x_continuous(limits = xlim, breaks = breaks_x) +
    scale_y_continuous(limits = ylim) +
    labs(x = xlab, y = ylab, title = title) +
    theme_popgen_minimal() +
    theme(legend.title = element_blank())

  if (!is.null(out_file)) save_popgen(p, out_file, width, height, dpi)
  p
}
