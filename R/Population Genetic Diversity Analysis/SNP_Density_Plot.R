#### ============================================================
#### SNP_Density_Plot.R —— 全基因组 SNP 密度图（CMplot，plot.type = "d"）
####
#### 原理：按染色体长度（该染色体上最大 SNP 坐标）把染色体切成 bin.size 大小的窗口，
####       统计每个窗口里的 SNP 数，用颜色深浅表示密度。
####
#### 输入（制表符分隔，至少 2 列，不含表头）：
####   CHR   BP
####   1     15234
####   1     18902
####   ...
#### 第 1 列染色体（绘制顺序 = 首次出现的顺序），第 2 列 SNP 的基因组坐标。
####
#### 用法：
####   SNP_Density_Plot("demo_data/snp_pos.txt", out_dir = "demo_out")
#### ============================================================

SNP_Density_Plot <- function(snp_file,
                             out_dir  = NULL,      # 输出目录；NULL = 当前目录
                             memo     = "SNP_density",
                             bin_size = 1e6,       # 密度窗口大小（bp），1e6 = 1Mb
                             col      = c("darkgreen", "yellow", "red"),  # 低 -> 中 -> 高
                             file     = "pdf",     # "pdf" / "png" / "jpg" / "tiff"
                             width    = 14,
                             height   = 6,
                             dpi      = 300,
                             chr_labels = NULL,    # 自定义染色体标签，NULL = 用原值
                             verbose  = TRUE) {

  if (!requireNamespace("CMplot", quietly = TRUE)) {
    stop("缺少包 CMplot，请先 install.packages('CMplot')")
  }
  if (!file.exists(snp_file)) stop("找不到文件: ", snp_file)

  dat <- utils::read.table(snp_file, header = FALSE,
                           colClasses = c("character", "numeric"),
                           sep = "", check.names = FALSE)

  if (ncol(dat) < 2) {
    stop("SNP 位置文件至少需要 2 列：CHR, BP。当前只有 ", ncol(dat), " 列")
  }
  names(dat) <- c("CHR", "BP")[seq_len(ncol(dat))]

  #### CMplot 画密度图只需要 2 列，但内部仍按 Pmap 处理，
  #### 保险起见补齐 4 列（SNP, CHR, BP, 值）
  if (ncol(dat) == 2) {
    dat <- data.frame(
      SNP = paste0("snp", seq_len(nrow(dat))),
      CHR = dat$CHR,
      BP  = dat$BP,
      stringsAsFactors = FALSE
    )
  }

  #### 切到输出目录（CMplot 只能往工作目录写）
  if (!is.null(out_dir)) {
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
    old_wd <- setwd(out_dir)
    on.exit(setwd(old_wd), add = TRUE)
  }

  CMplot::CMplot(dat,
                 plot.type  = "d",
                 bin.size   = bin_size,
                 col        = col,
                 file       = file,
                 file.name  = memo,
                 width      = width,
                 height     = height,
                 dpi        = dpi,
                 chr.labels = chr_labels,
                 verbose    = verbose)

  invisible(NULL)
}
