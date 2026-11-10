#### ============================================================
#### FST_Manhattan_Plot.R —— Fst 曼哈顿图（CMplot）
####
#### 画法沿用原始脚本：CMplot(plot.type = "m") + 6 色染色体交替 + 阈值线
#### 染色体交替色从 theme_popgen.R 的 PALETTE_CHR 取，值与原脚本完全一致。
####
#### 输入文件格式（制表符分隔，带表头，至少 4 列）：
####   SNP        CHR   BP        FST
####   rs123      1     12345     0.021
####   ...（第 4 列及以后每列都是一条要画的性状）
#### ============================================================

FST_Manhattan_Plot <- function(fst_file,
                               out_dir = NULL,      # 输出目录；NULL = 当前目录
                               memo    = "FST",     # 输出文件名标记
                               plot_type = "m",
                               ylab      = "FST",
                               LOG10     = FALSE,
                               cex       = 0.5,
                               signal_cex = 0.5,
                               threshold     = 0.3,
                               threshold_col = "red",
                               threshold_lty = 1,
                               threshold_lwd = 1,
                               width  = 14,
                               height = 6,
                               file   = "pdf",
                               col    = PALETTE_CHR) {

  if (!requireNamespace("CMplot", quietly = TRUE)) {
    stop("缺少包 CMplot，请先 install.packages('CMplot')")
  }

  dat <- utils::read.table(fst_file, header = TRUE, check.names = FALSE)

  if (ncol(dat) < 4) {
    stop("Fst 文件至少需要 4 列：SNP, CHR, BP, 值。当前只有 ", ncol(dat), " 列")
  }

  #### 切到输出目录（CMplot 只能往工作目录写）
  if (!is.null(out_dir)) {
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
    old_wd <- setwd(out_dir)
    on.exit(setwd(old_wd), add = TRUE)
  }

  CMplot::CMplot(dat,
                 plot.type  = plot_type,
                 LOG10      = LOG10,
                 cex        = cex,
                 col        = col,
                 ylab       = ylab,
                 width      = width,
                 height     = height,
                 file       = file,
                 file.name  = memo,   # CMplot >= 4.0 用 file.name（旧版叫 memo）
                 signal.cex = signal_cex,
                 threshold     = threshold,
                 threshold.col = threshold_col,
                 threshold.lty = threshold_lty,
                 threshold.lwd = threshold_lwd)

  invisible(NULL)
}
