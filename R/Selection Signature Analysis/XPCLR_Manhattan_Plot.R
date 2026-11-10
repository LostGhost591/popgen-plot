#### ============================================================
#### XPCLR_Manhattan_Plot.R —— XP-CLR 曼哈顿图（多染色体累加坐标）
####
#### 画法沿用手写脚本：
####   · 各染色体首尾相接拼成累加坐标（bp_cum）
####   · 按染色体循环 8 色，highlight_chr（默认 Z）强制用第一色
####   · 画 99% 分位虚线（threshold_q）
####   · x 轴只标 keep_axis 里列出的染色体
####
#### 输入（制表符分隔，带表头，列名可调）：
####   chrom  start  stop  xpclr_norm
####   1      10000  15000 2.31
#### ============================================================

XPCLR_Manhattan_Plot <- function(xpclr_file,
                               out_file = NULL,
                               chr_col   = "chrom",
                               start_col = "start",
                               end_col   = "stop",
                               value_col = "xpclr_norm",
                               chr_levels = NULL,      # 默认为 1..33 + "Z"
                               keep_axis = c("1", "2", "3", "4", "5", "6", "7", "8",
                                             "9", "10", "12", "14", "17", "21", "33", "Z"),
                               colors     = NULL,      # NULL = PALETTE_CHR8
                               highlight_chr = "Z",    # 强制用第一色的染色体
                               threshold_q   = 0.99,
                               point_size    = 1,
                               point_alpha   = 0.8,
                               title = "XP-CLR Manhattan Plot",
                               xlab  = "Chromosome",
                               ylab  = "Normalized XP-CLR Score",
                               top_pct = NULL,         # 如 0.05 = 同时导出前 5% 位点
                               top_out = NULL,
                               width = 12, height = 6, dpi = 300) {

  if (!file.exists(xpclr_file)) stop("找不到文件: ", xpclr_file)

  df <- utils::read.table(xpclr_file, header = TRUE, stringsAsFactors = FALSE,
                          fill = TRUE, quote = "", comment.char = "")

  need <- c(chr_col, start_col, value_col)
  miss <- setdiff(need, names(df))
  if (length(miss)) {
    stop("输入文件缺少列: ", paste(miss, collapse = ", "),
         "\n现有列: ", paste(names(df), collapse = ", "))
  }

  d <- data.frame(
    chrom = as.character(df[[chr_col]]),
    start = as.numeric(df[[start_col]]),
    value = suppressWarnings(as.numeric(df[[value_col]])),
    stringsAsFactors = FALSE
  )
  ## 染色体长度用 stop 列算（与原始脚本一致），没有 stop 列就退回 start
  d$end <- if (end_col %in% names(df)) as.numeric(df[[end_col]]) else d$start

  #### ---- 染色体排序：先按 chr_levels，未列出的排在后面 ----
  if (is.null(chr_levels)) chr_levels <- c(as.character(1:33), "Z")
  lv <- c(intersect(chr_levels, unique(d$chrom)),
          sort(setdiff(unique(d$chrom), chr_levels)))
  d$chrom <- factor(d$chrom, levels = lv)
  d <- d[order(d$chrom, d$start), ]

  #### ---- 累加坐标 ----
  chr_len <- d |>
    dplyr::group_by(chrom) |>
    dplyr::summarise(chr_len = max(end, na.rm = TRUE), .groups = "drop") |>
    dplyr::mutate(tot = cumsum(as.numeric(chr_len)) - chr_len)

  d <- d |>
    dplyr::left_join(chr_len[, c("chrom", "tot")], by = "chrom") |>
    dplyr::arrange(chrom, start) |>
    dplyr::mutate(bp_cum = start + tot)

  #### ---- x 轴：刻度全画，标签只标 keep_axis ----
  axis_df <- d |>
    dplyr::group_by(chrom) |>
    dplyr::summarise(center = (max(bp_cum, na.rm = TRUE) + min(bp_cum, na.rm = TRUE)) / 2,
                     .groups = "drop")
  axis_df$label <- ifelse(as.character(axis_df$chrom) %in% keep_axis,
                          as.character(axis_df$chrom), "")

  #### ---- 配色：按染色体循环取色 ----
  if (is.null(colors)) colors <- PALETTE_CHR8
  col_map <- stats::setNames(rep(colors, length.out = length(lv)), lv)
  if (!is.null(highlight_chr) && highlight_chr %in% lv) {
    col_map[highlight_chr] <- colors[1]
  }

  #### ---- 作图 ----
  p <- ggplot(d, aes(x = bp_cum, y = value, colour = chrom)) +
    geom_point(alpha = point_alpha, size = point_size) +
    scale_colour_manual(values = col_map) +
    scale_x_continuous(labels = axis_df$label, breaks = axis_df$center) +
    geom_hline(yintercept = stats::quantile(d$value, threshold_q, na.rm = TRUE),
               linetype = "dashed", colour = "grey40") +
    labs(title = title, x = xlab, y = ylab) +
    theme_popgen_minimal() +
    theme(
      legend.position     = "none",
      plot.title          = element_text(hjust = 0.5, face = "bold", size = 14),
      panel.grid.major.x  = element_blank(),
      panel.grid.minor.x  = element_blank(),
      axis.text.x         = element_text(angle = 0, size = 9, vjust = 0.5)
    )

  #### ---- 可选：导出超过分位阈值的位点 ----
  if (!is.null(top_pct)) {
    thr <- stats::quantile(d$value, 1 - top_pct, na.rm = TRUE)
    top <- df[d$value >= thr, , drop = FALSE]
    if (!is.null(top_out)) {
      dir.create(dirname(top_out), showWarnings = FALSE, recursive = TRUE)
      utils::write.table(top, top_out, sep = "\t", quote = FALSE, row.names = FALSE)
    }
  }

  if (!is.null(out_file)) save_popgen(p, out_file, width, height, dpi)
  p
}
