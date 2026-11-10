#### ============================================================
#### make_demo_data.R —— 生成脱敏的示例数据
####
#### 全部为随机生成的假样本，样本名、群体名、位点名均为虚构，
#### 不含任何真实个体信息。格式与 PLINK / ADMIXTURE / CMplot 输出一致。
####
#### 运行：Rscript tools/make_demo_data.R
#### ============================================================

set.seed(20261110)

out_dir <- "demo_data"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

#### ---------- 基本设定：3 个群体 × 20 个样本 = 60 ----------
pops    <- c("PopA", "PopB", "PopC")
n_per   <- 20
samples <- paste0(rep(pops, each = n_per), "_", sprintf("%02d", rep(seq_len(n_per), 3)))
pop_of  <- rep(pops, each = n_per)

#### ---------- 1. 分组文件（PCA 二维/三维/ADMIXTURE 共用）----------
group_df <- data.frame(SampleID = samples, Breed = pop_of, stringsAsFactors = FALSE)
write.csv(group_df, file.path(out_dir, "demo.group.csv"), row.names = FALSE, quote = FALSE)

#### ---------- 2. PCA eigenvec / eigenval ----------
n_pc    <- 20
centers <- list(PopA = c(-0.055,  0.020, -0.010),
                PopB = c( 0.062, -0.028,  0.008),
                PopC = c( 0.008,  0.072,  0.032))

pcs <- matrix(0, nrow = length(samples), ncol = n_pc)
for (i in seq_along(pops)) {
  idx <- which(pop_of == pops[i])
  pcs[idx, 1:3] <- matrix(rep(centers[[pops[i]]], length(idx)),
                          nrow = length(idx), byrow = TRUE)
}
## 前 3 个 PC 加群体内噪声，第 4 个以后纯随机
pcs[, 1:3] <- pcs[, 1:3] + matrix(rnorm(length(samples) * 3, 0, 0.012),
                                  nrow = length(samples))
pcs[, 4:n_pc] <- matrix(rnorm(length(samples) * (n_pc - 3), 0, 0.01),
                        nrow = length(samples))

eigenvec <- data.frame(FID = pop_of, IID = samples, pcs, check.names = FALSE)
names(eigenvec)[3:ncol(eigenvec)] <- paste0("PC", seq_len(n_pc))
utils::write.table(eigenvec, file.path(out_dir, "demo.eigenvec"),
                   row.names = FALSE, col.names = FALSE, quote = FALSE, sep = " ")

## 特征值：单调递减，PC1 最大
eigenval <- round(seq(18, 1.2, length.out = n_pc) * exp(-seq(0, 1.8, length.out = n_pc)), 6)
utils::write.table(eigenval, file.path(out_dir, "demo.eigenval"),
                   row.names = FALSE, col.names = FALSE, quote = FALSE)

#### ---------- 3. Fst 文件（SNP CHR BP FST）----------
n_chr  <- 10
n_snp  <- 200
fst_df <- do.call(rbind, lapply(seq_len(n_chr), function(ch) {
  data.frame(
    SNP = sprintf("chr%d_%05d", ch, seq_len(n_snp)),
    CHR = ch,
    BP  = sort(sample(1:30000000, n_snp)),
    FST = round(abs(rnorm(n_snp, 0.05, 0.04)), 4),
    stringsAsFactors = FALSE
  )
}))
## 造几个超过阈值 0.3 的信号点，让阈值线看得出来
sig <- sample(seq_len(nrow(fst_df)), 12)
fst_df$FST[sig] <- round(runif(12, 0.31, 0.48), 4)
utils::write.table(fst_df, file.path(out_dir, "demo.fst.txt"),
                   row.names = FALSE, quote = FALSE, sep = "\t")

#### ---------- 4. ADMIXTURE 的 Q 文件（K = 2..9）----------
## 每个 K 重新给三个群体分配主成分（互不相同），
## 且 K 越大主成分占比越低，模拟真实数据里 K 增大后成分逐渐混合的现象
for (K in 2:9) {
  q <- matrix(0, nrow = length(samples), ncol = K)
  mains <- sample(rep(seq_len(K), length.out = length(pops)))  # 该 K 下三个群体的主成分
  for (i in seq_along(pops)) {
    idx   <- which(pop_of == pops[i])
    main  <- mains[i]
    share <- runif(length(idx), 0.72 - 0.05 * (K - 2), 0.92 - 0.05 * (K - 2))
    q[cbind(idx, main)] <- share
    others <- setdiff(seq_len(K), main)
    for (j in idx) {
      w <- if (length(others)) {
        as.numeric(stats::rmultinom(1, 1000, rep(1, length(others)))) / 1000
      } else numeric(0)
      if (length(others)) q[j, others] <- (1 - q[j, main]) * w
    }
  }
  q <- q / rowSums(q)                           # 归一化，每行和为 1
  utils::write.table(round(q, 4), file.path(out_dir, paste0("demo.K.", K, ".Q")),
                     row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
}

#### ---------- 5. ADMIXTURE 需要的样本表（列名与函数默认一致）----------
adm_df <- data.frame(sampleID = samples, breed = pop_of, stringsAsFactors = FALSE)
write.csv(adm_df, file.path(out_dir, "demo.admixture.samples.csv"),
          row.names = FALSE, quote = FALSE)

#### ---------- 6. SNP 位置文件（SNP 密度图用，两列 CHR / BP）----------
## 每条染色体的 SNP 数不等，并在随机某 1Mb 区段加倍，形成密度热点
snp_list <- lapply(seq_len(n_chr), function(ch) {
  n  <- sample(2000:6000, 1)
  bp <- sample.int(30e6, n)
  hot <- sample.int(29, 1) * 1e6
  bp <- c(bp, sample(hot:(hot + 1e6), n %/% 5, replace = TRUE))
  data.frame(CHR = ch, BP = sort(as.integer(bp)))
})
snp_df <- do.call(rbind, snp_list)
utils::write.table(snp_df, file.path(out_dir, "demo.snp_pos.txt"),
                   row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")

#### ---------- 7. LD 衰减统计文件（PLINK --r2 输出格式，*.stat.gz）----------
## 第一行是注释，之后每行两列：Distance, R2
## 三个群体给不同的衰减速度，曲线能看出区别
ld_pops <- c("25R" = 0.30, "49W" = 0.22, "80aixiao" = 0.15)
for (p in names(ld_pops)) {
  d  <- seq(0, 500, by = 1)
  r2 <- ld_pops[[p]] * exp(-d / (350 + 600 * ld_pops[[p]])) +
        stats::rnorm(length(d), 0, 0.012)
  r2 <- pmin(pmax(r2, 0), 1)
  txt <- c("#Distance(bp)\tR2",
           paste(format(d, scientific = FALSE), sprintf("%.5f", r2), sep = "\t"))
  con <- gzfile(file.path(out_dir, paste0("demo.ld.", p, ".stat.gz")), "w")
  writeLines(txt, con)
  close(con)
}

#### ---------- 8. 进化树 + 饼图数据（树饼图用）----------
## 随机树只是为了跑通代码，不代表真实系统发育关系
if (requireNamespace("ape", quietly = TRUE)) {
  tree <- ape::rtree(length(samples), tip.label = samples)
  ape::write.tree(tree, file = file.path(out_dir, "demo.tree.nwk"))

  ## 饼图数据：保留群体结构但不极端（主成分 0.45~0.80），
  ## 这样每个饼都能看出明显扇形；直接拿 K=3 的 Q 会几乎全是一个颜色
  w <- matrix(stats::runif(length(samples) * 3, 0.05, 0.35), ncol = 3)
  for (i in seq_along(pops)) {
    idx <- which(pop_of == pops[i])
    w[idx, i] <- stats::runif(length(idx), 0.45, 0.80)
  }
  w <- w / rowSums(w)
  pie_df <- data.frame(Tip  = samples,
                       PopA = round(w[, 1], 3),
                       PopB = round(w[, 2], 3),
                       PopC = round(w[, 3], 3),
                       stringsAsFactors = FALSE)
  write.csv(pie_df, file.path(out_dir, "demo.tree.pie.csv"),
            row.names = FALSE, quote = FALSE)
} else {
  message("未安装 ape，跳过进化树示例数据")
}

#### ---------- 9. XP-CLR 输入（chrom / start / stop / xpclr_norm）----------
xp_chr <- c(as.character(1:33), "Z")
xpclr_list <- lapply(xp_chr, function(ch) {
  n_win <- sample(150:300, 1)
  st <- sort(sample.int(25e6, n_win))
  data.frame(chrom      = ch,
             start      = st,
             stop       = st + 10000L,
             xpclr_norm = round(abs(stats::rnorm(n_win, 0.8, 0.7)), 3),
             stringsAsFactors = FALSE)
})
xpclr_df <- do.call(rbind, xpclr_list)
## 造一批高分窗口，让 99% 阈值线以上有信号点
sig_xp <- sample(seq_len(nrow(xpclr_df)), 40)
xpclr_df$xpclr_norm[sig_xp] <- round(stats::runif(40, 3.5, 6), 3)
utils::write.table(xpclr_df, file.path(out_dir, "demo.xpclr.txt"),
                   row.names = FALSE, quote = FALSE, sep = "\t")

#### ---------- 10. XP-EHH 输入（8 列无表头，第 1 列是 "chr:pos"）----------
xpehh_list <- lapply(xp_chr, function(ch) {
  n_win <- sample(150:300, 1)
  p     <- sort(sample.int(25e6, n_win))
  data.frame(pos        = paste0(ch, ":", p),
             gpos       = p,
             p1         = round(stats::runif(n_win), 3),
             ihh1       = round(stats::runif(n_win, 0, 2), 3),
             p2         = round(stats::runif(n_win), 3),
             ihh2       = round(stats::runif(n_win, 0, 2), 3),
             xpehh_raw  = round(stats::rnorm(n_win, 0, 1), 3),
             xpehh_norm = round(stats::rnorm(n_win, 0, 1), 3),
             stringsAsFactors = FALSE)
})
xpehh_df <- do.call(rbind, xpehh_list)
## 上、下各造一批极值，让双尾阈值线外有信号
sig_up <- sample(seq_len(nrow(xpehh_df)), 30)
xpehh_df$xpehh_norm[sig_up] <- round(stats::runif(30, 3, 5), 3)
sig_dn <- sample(setdiff(seq_len(nrow(xpehh_df)), sig_up), 30)
xpehh_df$xpehh_norm[sig_dn] <- round(-stats::runif(30, 3, 5), 3)
utils::write.table(xpehh_df, file.path(out_dir, "demo.xpehh.txt"),
                   row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")

#### ---------- 11. 韦恩图基因集（最多 4 列，每列长度不等，空位留空）----------
venn_sets <- list(
  List1 = sprintf("gene%03d",   1:150),
  List2 = sprintf("gene%03d", 100:250),
  List3 = sprintf("gene%03d", 120:280),
  List4 = sprintf("gene%03d", 130:320)
)
max_n   <- max(lengths(venn_sets))
venn_df <- as.data.frame(
  lapply(venn_sets, function(x) c(x, rep(NA_character_, max_n - length(x)))),
  stringsAsFactors = FALSE
)
names(venn_df) <- names(venn_sets)
utils::write.csv(venn_df, file.path(out_dir, "demo.venn.csv"),
                 row.names = FALSE, na = "")

cat("示例数据已生成到", normalizePath(out_dir), "\n")
cat("文件清单一览：\n")
print(list.files(out_dir))
