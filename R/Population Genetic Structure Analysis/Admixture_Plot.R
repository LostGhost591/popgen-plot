#### ============================================================
#### Admixture_Plot.R —— ADMIXTURE 堆叠柱状图（K = 2..K 多面板）
####
#### 画法沿用原始脚本：mfrow 多面板 + barplot(space = 0) + 右侧 K 标签
####                + 底部两行 axis 画群体分隔线和群体名
####
#### 相对原脚本的修正（都可用参数退回原行为）：
####   1. mfrow 默认 6，与原脚本完全一致。
####      注意：k_range 超过 6 个（如 2:9 共 8 个）时 PDF 会分成两页 —— 这是原脚本
####      就有的行为。想一页画完，传 mfrow = length(k_range) 并把 height 同步加大，
####      否则 mar = c(10,4,3,10) 会把面板挤到报 "figure margins too large"
####   2. dev.off() 只在真的开过设备时才调用 —— 原代码无条件关闭，
####      out.type 不设时会把调用方的绘图设备一起关掉
####   3. family = "GB1" 改成参数 font_family，默认不指定 —— "GB1"
####      是 Windows 专有中文字体名，换机器/换系统会直接报错
#### ============================================================

Admixture_Plot <- function(sname,                    # 含 breed + sampleID 的 csv
                           file_prefix,              # Q 文件前缀，如 "admixture/bestqc."
                           out,                      # 输出文件名前缀（不带后缀）
                           k_range   = 2:9,          # 画哪几个 K
                           height    = 18,
                           width     = 14,
                           cexlab    = 3,
                           out_type  = c("pdf", "png", "none"),
                           liner     = 3,
                           adj       = -0.25,
                           mfrow     = 6,            # 与原脚本一致；传 NULL 则自动 = K 的个数
                           font_family = "",         # PDF 中文字体；Windows 可试 "GB1"
                           breed_col  = "breed",
                           sample_col = "sampleID",
                           dpi        = 900,
                           png_width  = 2000 * 10,   # 原脚本写死的像素尺寸
                           png_height = 1500 * 10) {

  out_type <- match.arg(out_type)

  #### ---- 读样本与分组 ----
  sn <- utils::read.csv(file = sname, header = TRUE, as.is = TRUE)
  if (!breed_col  %in% names(sn)) stop("csv 里找不到列: ", breed_col)
  if (!sample_col %in% names(sn)) stop("csv 里找不到列: ", sample_col)

  group <- unique(sn[[breed_col]])

  #### ---- 算每个群体的柱宽与标签位置（原逻辑，未改）----
  len <- vector(); div <- vector(); lab <- vector()
  for (i in seq_along(group)) {
    len <- c(len, nrow(sn[sn[[breed_col]] == group[i], ]))
    div <- c(div, sum(len))
    lab <- c(lab, sum(len) - len[i] / 2 + 0.5)
  }

  #### ---- 打开设备 ----
  dev_opened <- FALSE
  if (out_type == "pdf") {
    dir.create(dirname(paste0(out, ".pdf")), showWarnings = FALSE, recursive = TRUE)
    if (nzchar(font_family)) {
      grDevices::pdf(paste0(out, ".pdf"), height = height, width = width,
                     family = font_family)
    } else {
      grDevices::pdf(paste0(out, ".pdf"), height = height, width = width)
    }
    dev_opened <- TRUE
  }
  if (out_type == "png") {
    dir.create(dirname(paste0(out, ".png")), showWarnings = FALSE, recursive = TRUE)
    grDevices::png(paste0(out, ".png"),
                   height = png_height, width = png_width,
                   type = "cairo", res = dpi)
    dev_opened <- TRUE
  }

  #### ---- 面板布局 ----
  if (is.null(mfrow)) mfrow <- length(k_range)
  op <- graphics::par(mar = c(10, 4, 3, 10), las = 2, mfrow = c(mfrow, 1))

  ## 收尾顺序必须是「先恢复 par，再关设备」。
  ## 反过来的话 par(op) 会在设备已关闭时执行，R 会当场另开一个默认 pdf 设备，
  ## 结果是调用方凭空多一个设备、工作目录多一个 Rplots.pdf。
  on.exit({
    graphics::par(op)
    if (dev_opened) grDevices::dev.off()
  }, add = TRUE)

  #### ---- 逐个 K 画堆叠柱 ----
  for (k in k_range) {
    qfile <- paste0(file_prefix, k, ".Q")
    if (!file.exists(qfile)) {
      warning("跳过 K=", k, "：找不到文件 ", qfile)
      next
    }
    q <- utils::read.table(file = qfile, as.is = TRUE)

    graphics::barplot(t(as.matrix(q)),
                      col        = popgen_colors_ramp(k),
                      border     = FALSE,
                      space      = 0,
                      names.arg  = sn[[sample_col]],
                      cex.names  = 0.946,
                      axes       = FALSE)

    graphics::legend("right", paste0("K=", k),
                     cex = cexlab + 1,
                     x.intersp = 100000, y.intersp = 100000,
                     bty = "n", adj = adj)

    graphics::axis(2, line = -3)
  }

  #### ---- 底部：群体分隔线 + 群体名（原逻辑，未改）----
  graphics::axis(1, line = liner, at = c(0, div), labels = NA,
                 lwd.tick = 2, lwd = 0, tck = -0.03)
  graphics::axis(1, line = liner, at = c(0, lab, div[length(group)]),
                 labels = c("", group, ""),
                 cex.axis = cexlab + 0.5, las = 0, padj = 0.5, tck = 0)

  invisible(NULL)
}

#### 兼容原函数名，老脚本的 admixture(...) 仍然能跑
admixture <- Admixture_Plot
