#### ============================================================
#### XPEHH_Manhattan_Plot.R —— XP-EHH 曼哈顿图（双尾阈值）
####
#### 画法沿用手写脚本：
####   · 输入 8 列无表头，第 1 列是 "染色体:坐标" 字符串，第 8 列是标准化 XP-EHH
####   · 按染色体序号做累加坐标
####   · 8 色循环着色
####   · 上、下各 5% 分位画两条黑色虚线
####   · x 轴所有染色体都有刻度，只显示 show_chr_labels 里列出的标签
####
#### 输入（空白分隔，无表头）：
####   1:12345  12345  0.1  0.9  0.1  0.8  -1.23  -1.45
#### ============================================================

XPEHH_Manhattan_Plot <- function(xpehh_file,
                               out_file = NULL,
                               col_names = c("pos", "gpos", "p1", "ihh1",
                                             "p2", "ihh2", "xpehh_raw", "xpehh_norm"),
                               value_col = "xpehh_norm",
                               show_chr_labels = c("1", "2", "3", "4", "5", "6", "7",
                                                   "8", "9", "11", "13", "16", "19",
                                                   "22", "27", "33", "Z"),
                               colors     = NULL,       # NULL = PALETTE_CHR8
                               top_pct    = 0.05,       # 上/下各取这个比例作阈值
                               add_threshold = TRUE,
                               point_size  = 0.5,
                               point_alpha = 0.7,
                               xlab = "Chr",
                               ylab = "XP-EHH",
                               width = 14, height = 7, dpi = 300) {

  if (!file.exists(xpehh_file)) stop("找不到文件: ", xpehh_file)

  d <- utils::read.table(xpehh_file, header = FALSE, stringsAsFactors = FALSE,
                         fill = TRUE)
  if (ncol(d) < length(col_names)) {
    stop("输入文件至少需要 ", length(col_names), " 列，实际只有 ", ncol(d), " 列")
  }
  names(d)[seq_along(col_names)] <- col_names

  if (!value_col %in% names(d)) stop("找不到数值列: ", value_col)

  #### ---- 从 "chr:pos" 里拆染色体和物理位置 ----
  parts <- strsplit(as.character(d$pos), ":", fixed = TRUE)
  d$chr <- vapply(parts, function(x) if (length(x) >= 1) x[1] else NA_character_, character(1))
  d$physical_pos <- vapply(parts,
                           function(x) if (length(x) >= 2) suppressWarnings(as.numeric(x[2])) else NA_real_,
                           numeric(1))
  d$value <- suppressWarnings(as.numeric(d[[value_col]]))
  d <- d[!is.na(d$chr) & !is.na(d$physical_pos), ]

  #### ---- 染色体排序：数字升序在前，其它（Z / X / 带字母的）在后 ----
  all_chr  <- unique(d$chr)
  num_chr  <- all_chr[grepl("^[0-9]+$", all_chr)]
  char_chr <- all_chr[!grepl("^[0-9]+$", all_chr)]
  ordered_chr <- c(sort(as.numeric(num_chr)), sort(char_chr))
  d$chr_num <- match(d$chr, ordered_chr)

  #### ---- 累加坐标 ----
  chr_len <- d |>
    dplyr::group_by(chr) |>
    dplyr::summarise(len = max(physical_pos, na.rm = TRUE), .groups = "drop") |>
    dplyr::mutate(ord = match(chr, ordered_chr)) |>
    dplyr::arrange(ord) |>
    dplyr::mutate(bp_add = cumsum(len) - len)

  d <- merge(d, chr_len[, c("chr", "bp_add")], by = "chr", sort = FALSE)
  d <- d[order(d$chr_num, d$physical_pos), ]
  d$pos_abs <- d$physical_pos + d$bp_add

  #### ---- x 轴：刻度全画，标签只标指定的染色体 ----
  axis_centers <- d |>
    dplyr::group_by(chr) |>
    dplyr::summarise(center = mean(pos_abs, na.rm = TRUE), .groups = "drop") |>
    dplyr::mutate(ord = match(chr, ordered_chr)) |>
    dplyr::arrange(ord)
  all_labels <- ifelse(axis_centers$chr %in% show_chr_labels, axis_centers$chr, "")

  #### ---- 双尾阈值：上 top_pct 与下 top_pct ----
  sv <- sort(d$value, na.last = NA)
  k  <- max(1L, ceiling(length(sv) * top_pct))
  threshold_lower <- sv[k]
  threshold_upper <- sv[length(sv) - k + 1L]

  #### ---- 8 色循环 ----
  if (is.null(colors)) colors <- PALETTE_CHR8
  d$color_group <- (d$chr_num - 1) %% length(colors) + 1

  #### ---- 作图 ----
  p <- ggplot(d, aes(x = pos_abs, y = value)) +
    geom_point(aes(colour = factor(color_group)),
               alpha = point_alpha, size = point_size) +
    scale_colour_manual(values = colors) +
    scale_x_continuous(breaks = axis_centers$center,
                       labels = all_labels,
                       expand = expansion(mult = 0.02)) +
    scale_y_continuous(expand = expansion(mult = 0.05)) +
    labs(x = xlab, y = ylab) +
    theme_popgen_classic() +
    theme(
      legend.position   = "none",
      axis.text.x       = element_text(angle = 0, hjust = 0.5, size = 10, colour = "black"),
      axis.text.y       = element_text(size = 10, colour = "black"),
      axis.title.x      = element_text(size = 12, face = "bold", margin = margin(t = 10)),
      axis.title.y      = element_text(size = 12, face = "bold", margin = margin(r = 10)),
      axis.line         = element_line(colour = "black", linewidth = 0.5),
      axis.ticks        = element_line(colour = "black", linewidth = 0.5),
      axis.ticks.length = unit(0.2, "cm")
    )

  if (isTRUE(add_threshold)) {
    p <- p + geom_hline(yintercept = c(threshold_upper, threshold_lower),
                        linetype = "dashed", colour = "black",
                        linewidth = 0.8, alpha = 0.8)
  }

  if (!is.null(out_file)) save_popgen(p, out_file, width, height, dpi)
  attr(p, "thresholds") <- c(lower = threshold_lower, upper = threshold_upper)
  p
}
