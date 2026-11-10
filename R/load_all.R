#### ============================================================
#### load_all.R —— 一次载入 R/ 下全部作图函数（含子文件夹）
####
#### 用法（在仓库根目录）：
####   source("R/load_all.R")
####
#### 之后可直接调用：
####   群体结构：PCA_2D_Plot()  PCA_Scree_Plot()  PCA_3D_Plot()
####             Admixture_Plot()  Tree_Pie_Plot()
####   遗传多样性：LD_Decay_Plot()  SNP_Density_Plot()
####   选择信号：FST_Manhattan_Plot()
#### ============================================================

local({
  ## ---- 先切到 UTF-8 locale ----
  ## Windows 上 R 默认 LC_CTYPE 常是 "C" 或 GBK，此时 source(encoding="UTF-8")
  ## 遇到中文注释会报 "invalid input found" 并静默截断文件（后面的函数定义丢失）。
  ## 切不动也不影响执行，只是中文注释可能显示为乱码。
  for (loc in c(".UTF-8", "en_US.UTF-8", "Chinese.UTF-8")) {
    if (!is.na(suppressWarnings(Sys.setlocale("LC_CTYPE", loc)))) break
  }

  ## ---- 定位本文件所在目录（= R/），不依赖工作目录 ----
  ofile <- tryCatch(normalizePath(sys.frame(1)$ofile), error = function(e) NULL)
  if (is.null(ofile) || !file.exists(ofile)) {
    ofile <- tryCatch(normalizePath("R/load_all.R"), error = function(e) NULL)
  }
  base <- if (!is.null(ofile) && file.exists(ofile)) dirname(ofile) else "R"

  ## ---- 递归收集所有 .R 文件（支持子文件夹）----
  all_r <- list.files(base, pattern = "\\.R$", recursive = TRUE,
                      full.names = TRUE, ignore.case = TRUE)

  ## 排除本文件自己，否则会无限递归
  all_r <- all_r[basename(all_r) != "load_all.R"]

  ## 排除 legacy/ 下的旧版手写脚本（不是函数，且写着绝对路径）
  all_r <- all_r[!grepl("(^|/)legacy/", gsub("\\\\", "/", all_r))]

  if (!length(all_r)) stop("在 ", base, " 下没找到任何 .R 文件")

  ## theme_popgen.R 必须最先加载：其它函数都依赖其中的样式常量
  theme_file <- all_r[basename(all_r) == "theme_popgen.R"]
  rest <- sort(setdiff(all_r, theme_file))

  for (f in c(theme_file, rest)) {
    source(f, encoding = "UTF-8")
  }

  assign(".popgen_dir", base, envir = globalenv())
  rel <- sub(paste0("^", gsub("\\\\", "/", base), "/"), "", gsub("\\\\", "/", c(theme_file, rest)))
  message("popgen-plot 已载入 ", length(rest) + length(theme_file), " 个文件（样式版本 ", popgen_version, "）：")
  for (x in rel) message("  - ", x)
})
