# WGS 群体遗传学分析流程：从原始数据到成图

本文档回答三个问题，每一项分析都按这个顺序走：

1. **这一步在做什么**（步骤）
2. **用什么代码算出来**（分析命令 → 中间文件）
3. **用什么代码画图**（`popgen-plot` 里的函数 → 成图）

分析命令按通用写法整理，路径用变量代替；作图函数全部在本仓库 `demo_data/` 上实测跑通，`Rscript tools/run_demo.R` 可复现。

文中每个绘图小节给出的 R 代码，都是本仓库 `R/` 下函数的调用，不是外部脚本；每个作图小节末尾附了对应源码文件的链接。

> **状态说明**：第 1 部分的 WGS 命令来自本机实际跑过的脚本（`重测序数据分析.sh`），已通用化为 [`scripts/00_wgs_align_gatk.sh`](../scripts/00_wgs_align_gatk.sh)，软件参数一字未改；本机没有 bwa / gatk，因此只做了语法检查与依赖检查，未重跑全流程。第 2–4 部分的作图代码全部实测通过（见 `demo_out/`）。

---

## 目录

- [总览](#overview)
- [第一部分：WGS 原始数据 → 基因型](#part1)
  - [0.1 配置区：唯一需要改的地方](#s0-1)
  - [0.2 阶段 1：比对与去重](#s0-2)
  - [0.3 阶段 2：深度与覆盖度报告](#s0-3)
  - [0.4 阶段 3：变异检测与过滤](#s0-4)
  - [0.5 转 PLINK（下游分析的起点）](#s0-5)
  - [0.6 相对原脚本改了什么](#s0-6)
- [第二部分：群体结构](#part2)
  - [1.1 主成分分析 PCA](#s1-1)
  - [1.2 祖先成分 ADMIXTURE](#s1-2)
  - [1.3 系统发育树 + 祖先成分饼图](#s1-3)
- [第三部分：遗传多样性](#part3)
  - [2.1 LD 衰减](#s2-1)
  - [2.2 SNP 密度](#s2-2)
- [第四部分：选择信号](#part4)
  - [3.1 群体分化 FST](#s3-1)
  - [3.2 XP-CLR](#s3-2)
  - [3.3 XP-EHH](#s3-3)
- [第五部分：通用](#part5)
  - [4.1 候选基因交集（韦恩图）](#s4-1)
- [附：一次跑完的顺序](#run)
- [待补](#todo)

---

## <a id="overview"></a>总览

| 环节 | 分析 | 软件 | 中间文件 | 作图函数 | 示例图 |
|---|---|---|---|---|---|
| 0 | 比对与变异检测 | bwa / samtools / qualimap / GATK | `snp_pass.vcf.gz` → `QC.bed` | — | — |
| 1.1 | 群体结构 PCA | plink `--pca` / GCTA | `QC.pca.eigenvec` `.eigenval` | `PCA_2D_Plot` `PCA_3D_Plot` `PCA_Scree_Plot` | `pca_2d.png` `pca_3d.png` `pca_scree.png` |
| 1.2 | 祖先成分 | ADMIXTURE | `QC.<K>.Q` | `Admixture_Plot` | `admixture_preview.png` |
| 1.3 | 系统发育树 + 成分 | plink IBS / FastTree | `.nwk` + `pie.csv` | `Tree_Pie_Plot` | `tree_pie.png` |
| 2.1 | LD 衰减 | PopLDdecay | `<pop>.stat.gz` | `LD_Decay_Plot` | `ld_decay.png` |
| 2.2 | SNP 密度 | plink bim | `snp_pos.txt` | `SNP_Density_Plot` | `Marker_Density.SNP_density.png` |
| 3.1 | 群体分化 FST | vcftools | `fst.txt`（4 列） | `FST_Manhattan_Plot` | `Rect_Manhtn.FST.png` |
| 3.2 | 选择信号 XP-CLR | xpclr | `xpclr_norm.txt` | `XPCLR_Manhattan_Plot` | `xpclr_manhattan.png` |
| 3.3 | 选择信号 XP-EHH | selscan | `xpehh.txt`（8 列） | `XPEHH_Manhattan_Plot` | `xpehh_manhattan.png` |
| 4 | 候选基因交集 | — | 基因列表 | `Venn_Plot` | `venn_3sets.png` |

全文约定：

```
data/    原始与外部输入（参考基因组、样本表、分组表）
bam/     比对结果
vcf/     变异检测结果
out/     各步骤中间文件 + 成图
```

**所有分析共用一个前缀 `QC`**，即 `QC.bed/.bim/.fam` 与 `QC.vcf`。后面每一步都从这个前缀出发，不再重复过滤。

---

# <a id="part1"></a>第一部分：WGS 原始数据 → 基因型

这一部分是纯命令行，来自本机实际跑过的脚本（`重测序数据分析.sh`）改造而成。

起点是 clean data（`${s}.Clean.R1.fq.gz`）。从下机数据到 clean data 的质控（fastp / Trimmomatic）
不在这个脚本里——按你们实验室的流程跑完，把 clean data 放进 `$DATADIR` 即可。

> **完整可执行脚本**：[`scripts/00_wgs_align_gatk.sh`](../scripts/00_wgs_align_gatk.sh)
>
> 只需要改脚本顶部的「配置区」，其余路径全部自动派生。
> 软件参数（`--max-alternate-alleles 4`、`--sample-ploidy 2`、`-ERC GVCF`、两条过滤表达式、
> `-M`、`-nr 10000`）与原脚本完全一致，没有调整。

运行方式：

```bash
bash scripts/00_wgs_align_gatk.sh              # 按配置区跑全流程
bash scripts/00_wgs_align_gatk.sh A1 A2 A3     # 临时指定样本（覆盖配置区）
RUN_STAGE=1 bash scripts/00_wgs_align_gatk.sh  # 只跑阶段 1（1=比对 2=质控 3=GATK all=全跑）
```

---

## <a id="s0-1"></a>0.1 配置区：唯一需要改的地方

原来的「设置工作目录和参考文件路径」写死在脚本里，换个人、换个物种就跑不了。
现在集中成一段，改完自动派生所有子路径：

```bash
# 项目与参考
WORKDIR=/data/sda/zhou/aixiao            # 项目根目录
DATADIR=$WORKDIR/A                       # clean data 所在目录
REFDIR=/data/sda/zhou/chickengenomeindex6a/refbwaindex

REF_NAME=Gallus_gallus.GRCg6a.dna.toplevel.fa          # 参考基因组文件名
GFF_NAME=Gallus_gallus.GRCg6a.104.chr.gff3             # GFF（下游注释）
GTF_NAME=Gallus_gallus.GRCg6a.104.chr.gtf              # GTF（下游注释）
DICT_NAME=Gallus_gallus.GRCg6a.dna.toplevel.dict       # Picard 字典（GATK 必需）

# 样本
SAMPLES=(A1 A2 A3 A4 A5 A6 A7 A8 A9)     # 样本名，数量任意
CLEAN_R1=.Clean.R1.fq.gz                 # clean data 命名：$DATADIR/${s}${CLEAN_R1}
CLEAN_R2=.Clean.R2.fq.gz

# 资源与阶段
THREADS_ALN=1        # bwa mem -t
THREADS_SORT=4       # samtools sort -@
MEM_SORT=4G          # samtools sort -m
JAVA_MEM=4g          # gatk --java-options -Xmx
THREADS_QC=2         # qualimap -nt
QC_DEPTH=10000       # qualimap -nr
RUN_STAGE=all        # 1 / 2 / 3 / all
KEEP_BAM=1           # 0 = 排序后删中间 bam，只留去重版
```

派生规则：`REF=$REFDIR/$REF_NAME`、`FAI=$REF.fai`、`BWAOUT=$DATADIR/3.map/bwa`、
`QCOUT=$DATADIR/3.map/result`、`VARDIR=$DATADIR/4.snp_indel/GATK`。

**换物种**只改四个 `*_NAME`；**换项目**只改前三行；**样本数变了**改 `SAMPLES`
（或取消 `SAMPLES_FILE` 那行的注释，改成一行一个样本的文件）。

---

## <a id="s0-2"></a>0.2 阶段 1：比对与去重

**步骤**：bwa mem 比对（带 read group）→ 转 BAM → 排序 → 去 PCR 重复 → 建索引。

```bash
bwa mem "$REF" "$R1" "$R2" -t "$THREADS_ALN" -M \
    -R "@RG\tID:${s}\tLB:${s}\tPL:ILLUMINA\tSM:${s}" \
  | samtools view -bS -h - > "${s}.bam"

samtools sort -@ "$THREADS_SORT" -m "$MEM_SORT" "${s}.bam" -O bam -o "${s}.sorted.bam"

samtools rmdup "${s}.sorted.bam" "${s}.rmdup.bam"
samtools index "${s}.rmdup.bam"
```

产出 `$DATADIR/3.map/bwa/${s}.rmdup.bam`。read group 里的 `SM` 必须与后续样本表一致，
否则 GATK 合并时会按样本名报错。

> `samtools rmdup` 在 samtools ≥ 1.21 已移除。脚本注释里给了等价写法：
> `sort -n` → `fixmate -m` → `sort` → `markdup -r`。

---

## <a id="s0-3"></a>0.3 阶段 2：深度与覆盖度报告

**步骤**：每个样本出一份 PDF，看平均深度、覆盖度、插入片段长度，判断测序质量是否达标。

```bash
qualimap bamqc --java-mem-size="$JAVA_MEM" -nt "$THREADS_QC" -nr "$QC_DEPTH" \
    -bam "$BAM" -outfile report.pdf -outformat PDF \
    -outdir "$QCOUT/${s}.depth_and_coverage_report"
```

一个样本一个目录，输出在 `$DATADIR/3.map/result/`。

---

## <a id="s0-4"></a>0.4 阶段 3：变异检测与过滤

**步骤**：单样本 gVCF → 合并 → 统一基因分型 → SNP / INDEL 分别硬过滤。

```bash
## 单样本 gVCF
gatk --java-options "-Xmx${JAVA_MEM}" HaplotypeCaller -R "$REF" \
    -I "$BAM" -O "${s}.g.vcf.gz" \
    --max-alternate-alleles 4 --sample-ploidy 2 \
    -ERC GVCF --tmp-dir "$TMPDIR"

## 合并（--variant 由脚本按样本列表自动拼，样本数任意）→ 基因分型
gatk --java-options "-Xmx${JAVA_MEM}" CombineGVCFs -R "$REF" \
    "${VARGVCFS[@]}" -O all.g.vcf.gz --tmp-dir "$TMPDIR"
gatk --java-options "-Xmx${JAVA_MEM}" GenotypeGVCFs -R "$REF" \
    -V all.g.vcf.gz -O all.raw.vcf.gz --tmp-dir "$TMPDIR"

## SNP 硬过滤
SNP_FILTER="QD < 2.0 || MQ < 40.0 || FS > 60.0 || SOR > 3.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0"
gatk SelectVariants -select-type SNP -V all.raw.vcf.gz -O snp.vcf.gz
gatk VariantFiltration -V snp.vcf.gz \
    --filter-expression "$SNP_FILTER" --filter-name "Filter" \
    -O snp.filter.vcf.gz
gatk SelectVariants -R "$REF" -V snp.filter.vcf.gz --exclude-filtered -O snp_pass.vcf.gz
tabix -p vcf snp_pass.vcf.gz

## INDEL 硬过滤（阈值比 SNP 宽）
INDEL_FILTER="QD < 2.0 || FS > 200.0 || SOR > 10.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0"
gatk SelectVariants -select-type INDEL -V all.raw.vcf.gz -O indel.vcf.gz
gatk VariantFiltration -V indel.vcf.gz \
    --filter-expression "$INDEL_FILTER" --filter-name "Filter" \
    -O indel.filter.vcf.gz
gatk SelectVariants -R "$REF" -V indel.filter.vcf.gz --exclude-filtered -O indel_pass.vcf.gz
tabix -p vcf indel_pass.vcf.gz
```

产出：`$DATADIR/4.snp_indel/GATK/snp_pass.vcf.gz`、`indel_pass.vcf.gz`。

---

## <a id="s0-5"></a>0.5 转 PLINK（下游分析的起点）

原脚本到 VCF 为止。后面所有群体遗传分析都从 PLINK 前缀出发，补这一步：

```bash
plink --vcf snp_pass.vcf.gz --allow-extra-chr --chr-set 29 \
      --double-id --make-bed --out out/QC.raw

## 群体层面过滤：缺失率 5%、MAF 3%、HWE 1e-6
plink --bfile out/QC.raw --geno 0.05 --maf 0.03 --hwe 1e-6 \
      --make-bed --out out/QC
```

需要 VCF 的分析（PopLDdecay / xpclr / selscan）再转回去：

```bash
plink --bfile out/QC --recode vcf-iid --out out/QC
bgzip out/QC.vcf && tabix -p vcf out/QC.vcf.gz
```

---

## <a id="s0-6"></a>0.6 相对原脚本改了什么

| 原写法 | 问题 | 现在 |
|---|---|---|
| `\  # 注释` | 续行符后面跟注释，shell 把整行当注释，命令被从中截断 | 注释挪到命令行上方 |
| `CombineGVCFs --variant ${samples[0]} … ${samples[8]}` 写死 9 行 | 相对路径找不到 gVCF；样本数一变就得手改 | 数组自动拼 `--variant`，数量任意 |
| `qualimap -bam .../${i}.sorted.dedup.bam` | 实际生成的是 `${i}.rmdup.bam`，路径和文件名都不对 | 指向真实文件 |
| `mkdir $workdir/data/4.snp_indel/GATK` + `B=$datadir/4.snp_indel/GATK` | 两个目录不是同一个（`$workdir/data` vs `$workdir/A`） | 统一到 `$DATADIR/4.snp_indel/GATK` |
| 路径、样本、线程、内存全部硬编码 | 换人 / 换物种 / 换机器都要全文搜改 | 顶部配置区集中管理 |

> 软件参数一个都没动。改的只有路径组织方式和 shell 语法。

# <a id="part2"></a>第二部分：群体结构

## <a id="s1-1"></a>1.1 主成分分析 PCA

**步骤**

1. LD 剪枝（PCA 前必须做，否则近缘位点会主导主成分）
2. 计算前 20 个主成分
3. 用特征值算各 PC 的解释方差
4. 出图：二维散点（带 95% 置信椭圆）、三维散点、碎石图

```bash
plink --bfile out/QC --indep-pairwise 50 5 0.2 --out out/QC.prune
plink --bfile out/QC --extract out/QC.prune.prune.in \
      --pca 20 --out out/QC.pca
```

产出 `out/QC.pca.eigenvec`（`FID IID PC1..PC20`，空格分隔无表头）、`out/QC.pca.eigenval`（每行一个特征值）。

分组表 `data/group.csv`：

```csv
SampleID,Breed
PopA_01,PopA
PopA_02,PopA
PopB_01,PopB
```

**作图**

```r
source("R/load_all.R")

## 二维散点 + 95% 置信椭圆
PCA_2D_Plot(
  eigenvec_path = "out/QC.pca.eigenvec",
  eigenval_path = "out/QC.pca.eigenval",
  group_path    = "data/group.csv",
  group_col     = "Breed",
  out_file      = "out/PCA_2D.png"
)

## 三维散点
PCA_3D_Plot(
  eigenvec_path = "out/QC.pca.eigenvec",
  eigenval_path = "out/QC.pca.eigenval",
  group_path    = "data/group.csv",
  out_file      = "out/PCA_3D.png"
)

## 碎石图（前 10 个 PC 的解释方差）
PCA_Scree_Plot("out/QC.pca.eigenval", out_file = "out/PCA_Scree.png")
```

**源码**：[`PCA_2D_Plot.R`](../R/Population%20Genetic%20Structure%20Analysis/PCA_2D_Plot.R) · [`PCA_3D_Plot.R`](../R/Population%20Genetic%20Structure%20Analysis/PCA_3D_Plot.R)（`PCA_Scree_Plot` 与 `PCA_2D_Plot` 在同一个文件）

常用开关：`ellipse = FALSE` 去椭圆，`colors = c(...)` 换配色，`n_pc_axis` 控制取前几个 PC 算解释率。

---

## <a id="s1-2"></a>1.2 祖先成分 ADMIXTURE

**步骤**

1. 对 K = 2…9 各跑一次，记 CV error
2. 取 CV error 最小的 K 作为最优祖先数
3. 每个 K 输出一个 `.Q`（每行一个样本，K 列，行和为 1）
4. 出图：多层堆叠柱，一行一个 K

```bash
for K in $(seq 2 9); do
  admixture --cv out/QC.bed $K > out/admixture/log${K}.out 2>&1
done

## 看哪个 K 的 CV error 最小
grep -h "CV error" out/admixture/log*.out
```

产出 `out/admixture/QC.2.Q … QC.9.Q`。

样本表 `data/admixture.samples.csv`：

```csv
sampleID,breed
PopA_01,PopA
PopB_01,PopB
```

**作图**

```r
Admixture_Plot(
  sname       = "data/admixture.samples.csv",
  file_prefix = "out/admixture/QC.",   # 对应 QC.2.Q … QC.9.Q
  out         = "out/ADMIXTURE",
  k_range     = 2:9,
  out_type    = "pdf"                  # 投稿用矢量图
)
```

配色默认走统一色板（`popgen_colors_ramp()`）。**要跟以前出的图完全一致**，配色和旧尺寸都要显式传：

```r
Admixture_Plot(
  sname       = "data/admixture.samples.csv",
  file_prefix = "out/admixture/QC.",
  out         = "out/ADMIXTURE_rainbow",
  height = 15, width = 22,          # 旧图尺寸
  cexlab = 2, liner = 5,            # 旧字号与分隔线位置
  adj = -0.25, mfrow = 4,           # 旧图例位置与每页面板数
  out_type    = "pdf",              # 注意：旧脚本里叫 out.type
  colors      = function(k) rainbow(k),
  font_family = "GB1"               # 旧脚本写死的 Windows 中文字体，非 Windows 去掉
)
```

`colors` 三种写法：

| 传什么 | 效果 |
|---|---|
| 不传 / `NULL` | 统一色板，全仓库配色一致（默认） |
| `function(k) rainbow(k)` | 复刻旧图的 `rainbow()` |
| `c("#E31A1C", "#1F78B4", ...)` | 固定色板，按 K 循环取色 |

**源码**：[`Admixture_Plot.R`](../R/Population%20Genetic%20Structure%20Analysis/Admixture_Plot.R)

要 PNG 预览版时加 `out_type = "png"`，并用 `mfrow` 控制每页行数：

```r
Admixture_Plot(sname = "data/admixture.samples.csv",
               file_prefix = "out/admixture/QC.",
               out = "out/ADMIXTURE_preview", out_type = "png",
               mfrow = 8, png_width = 2600, png_height = 3600, dpi = 120)
```

---

## <a id="s1-3"></a>1.3 系统发育树 + 祖先成分饼图

**步骤**

1. 由 `QC.vcf` 生成多序列比对（phy）
2. 建树（ML 用 FastTree/RAxML；NJ 可用 PLINK IBS 距离）
3. 取 ADMIXTURE 最优 K 的 Q 值转成宽表 `pie.csv`
4. 出图：树 + 叶节点挂成分饼图

```bash
## 建树
python tools/vcf2phylip.py -i out/QC.vcf -o out/QC
FastTree -nt out/QC.min4.phy > out/QC.tree.nwk

## NJ 路线（替代方案）：PLINK 算 IBS 距离
plink --bfile out/QC --genome --out out/QC.ibs
```

**Q 值转饼图表**（R 里一行）：

```r
q   <- read.table("out/admixture/QC.4.Q", header = FALSE)
ids <- read.csv("data/admixture.samples.csv")
pie <- data.frame(Tip = ids$sampleID, q)
names(pie)[-1] <- paste0("K", seq_len(ncol(q)))
write.csv(pie, "out/tree_pie.csv", row.names = FALSE)
```

`pie.csv` 格式：第一列叶节点名，后面每列一个祖先成分比例。

```csv
Tip,K1,K2,K3
PopA_01,0.80,0.15,0.05
PopB_01,0.10,0.85,0.05
```

**作图**

```r
Tree_Pie_Plot(
  tree_file = "out/QC.tree.nwk",
  pie_file  = "out/tree_pie.csv",
  tip_col   = "Tip",
  out_file  = "out/Tree_Pie.png",
  width = 7, height = 12
)
```

**源码**：[`Tree_Pie_Plot.R`](../R/Population%20Genetic%20Structure%20Analysis/Tree_Pie_Plot.R)

依赖 `ggtree` / `treeio` / `patchwork`，缺任何一个会直接报错并提示安装，不影响其它图。

---

# <a id="part3"></a>第三部分：遗传多样性

## <a id="s2-1"></a>2.1 LD 衰减

**步骤**

1. 按群体拆样本
2. 每个群体单独算 r² 随距离的衰减
3. 多群体叠一张图，比较衰减速度（快 = 重组历史长或有效群体小）

```bash
for pop in $(cat data/pops.list); do
  PopLDdecay -InVCF out/QC.vcf.gz \
             -SubPop data/${pop}.txt \
             -MaxDist 500 -OutStat out/ld/${pop}.stat
done
## 输出默认 gzip 压缩：out/ld/<pop>.stat.gz
```

**作图**

```r
LD_Decay_Plot(
  files = c(PopA = "out/ld/PopA.stat.gz",
            PopB = "out/ld/PopB.stat.gz",
            PopC = "out/ld/PopC.stat.gz"),
  out_file = "out/LD_Decay.png",
  title    = "LD Decay"
)
```

**源码**：[`LD_Decay_Plot.R`](../R/Population%20Genetic%20Diversity%20Analysis/LD_Decay_Plot.R)

曲线毛刺多时按距离分箱取均值：

```r
LD_Decay_Plot(files, bin_size = 20, out_file = "out/LD_Decay_smooth.png")
```

`xlim` 默认 `c(0, 500)`（单位 bp），算全基因组时按 `-MaxDist` 同步调大。

---

## <a id="s2-2"></a>2.2 SNP 密度

**步骤**

1. 从 `QC.bim` 取染色体与物理位置
2. 按窗口（默认 1 Mb）统计位点数
3. 出图：染色体着色的密度图，找低密度区（可能是组装缺口或低复杂度区）

```bash
awk '{print $1"\t"$4}' out/QC.bim > out/snp_pos.txt
```

`out/snp_pos.txt` 两列，无表头：

```
1	15234
1	15987
2	10233
```

**作图**

```r
SNP_Density_Plot("out/snp_pos.txt",
                 out_dir  = "out",
                 memo     = "SNP_density",
                 bin_size = 1e6,        # 1 Mb 窗口
                 file     = "pdf")
```

**源码**：[`SNP_Density_Plot.R`](../R/Population%20Genetic%20Diversity%20Analysis/SNP_Density_Plot.R)

---

# <a id="part4"></a>第四部分：选择信号

## <a id="s3-1"></a>3.1 群体分化 FST

**步骤**

1. 选两个对比群体（各写一个只有样本 ID 的单列表）
2. 滑窗计算 FST（常用 50 kb 窗口 / 50 kb 步长）
3. 转成 4 列格式（位点名 / 染色体 / 位置 / FST）
4. 出图：曼哈顿图 + 阈值线

```bash
vcftools --vcf out/QC.vcf \
         --weir-fst-pop data/popA.txt \
         --weir-fst-pop data/popB.txt \
         --fst-window-size 50000 --fst-window-step 50000 \
         --out out/FST_AB

## 窗口表 → CMplot 需要的 4 列（跳过表头行）
awk 'NR>1 {print "chr"$1"_"$2"\t"$1"\t"$2"\t"$3}' \
    out/FST_AB.windowed.weir.fst > out/fst.txt
```

`out/fst.txt`：

```
SNP	CHR	BP	FST
chr1_50000	1	50000	0.0213
```

阈值取全基因组 99% 分位更客观：

```bash
awk 'NR>1{print $4}' out/FST_AB.windowed.weir.fst | sort -g \
  | awk '{a[NR]=$1} END{print a[int(NR*0.99)]}'
```

**作图**

```r
FST_Manhattan_Plot(
  fst_file  = "out/fst.txt",
  out_dir   = "out",
  memo      = "FST",
  threshold = 0.3          # 换成上面算出的 99% 分位值
)
```

**源码**：[`FST_Manhattan_Plot.R`](../R/Selection%20Signature%20Analysis/FST_Manhattan_Plot.R)

`file = "png"` 出位图版；`LOG10 = TRUE` 转 -log10 纵轴。

---

## <a id="s3-2"></a>3.2 XP-CLR

**步骤**

1. 拆分染色体，准备每条染色体的 vcf 与遗传距离 map
2. 逐染色体跑 xpclr（50 kb 窗口 / 50 kb 步长），A 为参考群、B 为目标群
3. 合并各染色体，取 `norm_xpclr` 列（空值按 0 处理）
4. 出图：多染色体拼接的曼哈顿图 + 99% 分位线，顺带导出前 5% 位点

```bash
for k in $(seq 1 29) Z; do
  xpclr --out out/xpclr/chr${k} --format vcf \
        --input out/chr/chr${k}.vcf \
        --samplesA data/popA.txt --samplesB data/popB.txt \
        --map out/map/chr${k}.map.distance --chr ${k} \
        --phased --size 50000 --step 50000
done

## 合并：chrom / start / stop / xpclr_norm
{
  echo -e "chrom\tstart\tstop\txpclr_norm"
  for k in $(seq 1 29) Z; do
    awk -v c=${k} 'BEGIN{OFS="\t"} {v=($12==""?0:$12); print c,$2,$3,v}' \
        out/xpclr/chr${k}.xpclr.txt
  done
} > out/xpclr_norm.txt
```

**作图**

```r
XPCLR_Manhattan_Plot(
  xpclr_file = "out/xpclr_norm.txt",   # 需带表头：chrom/start/stop/xpclr_norm
  out_file   = "out/XPCLR.png",
  top_pct    = 0.05,                   # 同时导出前 5% 位点
  top_out    = "out/xpclr_top5.txt",
  highlight_chr = "Z"                  # Z 染色体强制用第一色
)
```

**源码**：[`XPCLR_Manhattan_Plot.R`](../R/Selection%20Signature%20Analysis/XPCLR_Manhattan_Plot.R)

染色体顺序默认 `1..33 + Z`，非哺乳动物改 `chr_levels`；轴标签密度改 `keep_axis`。

---

## <a id="s3-3"></a>3.3 XP-EHH

**步骤**

1. beagle 定相（phased vcf）
2. 按群体提取样本、拆染色体
3. 每条染色体算 xpehh（B 为目标群，A 为参考群）
4. `norm` 加 50 kb 窗口并标准化
5. 合并成 8 列，第一列 `chr:pos`
6. 出图：上下各 5% 作双尾阈值

```bash
## 清理 + 定相
java -jar beagle.jar gt=out/QC.vcf.gz out=out/QC.phased ne=200

## 拆染色体（示例 1..29 + Z）
for k in $(seq 1 29) Z; do
  bcftools view -S data/popB.txt out/QC.phased.vcf.gz -Ou \
    | bcftools view -r ${k} -Ov -o out/chr/B.chr${k}.vcf
  bcftools view -S data/popA.txt out/QC.phased.vcf.gz -Ou \
    | bcftools view -r ${k} -Ov -o out/chr/A.chr${k}.vcf
done

## 遗传距离 map（cM = bp / 1e6）
for k in $(seq 1 29) Z; do
  plink --vcf out/chr/A.chr${k}.vcf --out out/map/chr${k} --recode --allow-extra-chr
  awk 'BEGIN{OFS=" "}{print 1,".",$4/1000000,$4}' out/map/chr${k}.map \
      > out/map/chr${k}.map.distance
done

## 计算 + 窗口标准化
for k in $(seq 1 29) Z; do
  selscan --xpehh --vcf out/chr/B.chr${k}.vcf --vcf-ref out/chr/A.chr${k}.vcf \
          --map out/map/chr${k}.map.distance --threads 10 \
          --out out/xpehh/chr${k}
  norm --xpehh --files out/xpehh/chr${k}.xpehh.out --bp-win --winsize 50000
done

## 合并成 8 列：第 1 列 chr:pos，第 8 列是标准化值
for k in $(seq 1 29) Z; do
  awk -v c=${k} 'BEGIN{OFS="\t"}{$1=c":"$1; print}' \
      out/xpehh/chr${k}.xpehh.out.100bins.norm \
    > out/xpehh/chr${k}.8col.txt
done
cat out/xpehh/chr*.8col.txt > out/xpehh_8col.txt

## 校验：应为 8 列，且第 8 列有正负号
awk 'NR==1{print "列数:", NF}' out/xpehh_8col.txt
head -2 out/xpehh_8col.txt
```

> `.100bins.norm` 的列数随 selscan 版本略有差异。若 `NF` 不是 8，按实际列号改上面的 `awk`，
> 保证「输出 8 列、第 1 列 `chr:pos`、第 8 列标准化值」即可——函数只依赖这三点。

**作图**

```r
XPEHH_Manhattan_Plot(
  xpehh_file = "out/xpehh_8col.txt",   # 8 列无表头，第 1 列 "chr:pos"
  out_file   = "out/XP_EHH.png",
  top_pct    = 0.05                    # 上下各 5% 作阈值
)
```

**源码**：[`XPEHH_Manhattan_Plot.R`](../R/Selection%20Signature%20Analysis/XPEHH_Manhattan_Plot.R)

值的正负带方向：`> 2` 表示目标群受选择，`< -2` 表示参考群受选择。

---

# <a id="part5"></a>第五部分：通用

## <a id="s4-1"></a>4.1 候选基因交集（韦恩图）

**步骤**

1. 每个分析各出一份候选基因列表（如 FST 前 5%、XP-CLR 前 5%、XP-EHH 前 5% 的区间注释结果）
2. 取交集，看共同受选择的基因
3. 出图：2/3/4 集合韦恩图，顺带导出每个区域的元素明细

```r
source("R/load_all.R")

genesA <- readLines("out/FST_top5pct.genes.txt")
genesB <- readLines("out/XPCLR_top5pct.genes.txt")
genesC <- readLines("out/XPEHH_top5pct.genes.txt")

Venn_Plot(list(FST = genesA, `XP-CLR` = genesB, `XP-EHH` = genesC),
          out_file   = "out/Venn_3sets.png",
          title      = "Candidate genes",
          detail_out = "out/venn_regions.tsv")
```

**源码**：[`Venn_Plot.R`](../R/General%20Plots/Venn_Plot.R)

`detail_out` 会写出每个区域具体包含哪些基因，相当于 Venny 里点数字看内容。返回值的 `regions` 元素也能直接在 R 里取用。

---

# <a id="run"></a>附：一次跑完的顺序

```r
source("R/load_all.R")

## 1. 群体结构
PCA_2D_Plot("out/QC.pca.eigenvec", "out/QC.pca.eigenval",
            "data/group.csv", out_file = "out/PCA_2D.png")
PCA_3D_Plot("out/QC.pca.eigenvec", "out/QC.pca.eigenval",
            "data/group.csv", out_file = "out/PCA_3D.png")
PCA_Scree_Plot("out/QC.pca.eigenval", out_file = "out/PCA_Scree.png")
Admixture_Plot("data/admixture.samples.csv", "out/admixture/QC.",
               "out/ADMIXTURE", k_range = 2:9)
Tree_Pie_Plot("out/QC.tree.nwk", "out/tree_pie.csv", out_file = "out/Tree_Pie.png")

## 2. 遗传多样性
LD_Decay_Plot(c(PopA = "out/ld/PopA.stat.gz",
                PopB = "out/ld/PopB.stat.gz"),
              out_file = "out/LD_Decay.png")
SNP_Density_Plot("out/snp_pos.txt", out_dir = "out", bin_size = 1e6)

## 3. 选择信号
FST_Manhattan_Plot("out/fst.txt", out_dir = "out", memo = "FST")
XPCLR_Manhattan_Plot("out/xpclr_norm.txt", out_file = "out/XPCLR.png", top_pct = 0.05)
XPEHH_Manhattan_Plot("out/xpehh_8col.txt", out_file = "out/XP_EHH.png")
```

想先看效果不碰自己的数据，直接跑样例：

```bash
Rscript tools/run_demo.R
```

图在 `demo_out/`，输入假的但格式与真实数据完全一致（`demo_data/`）。

---

# <a id="todo"></a>待补

以下分析已有成熟软件流程，但本仓库还没有对应作图函数，暂列在这里，后续按需补：

- 系统发育之外的群体历史：Treemix、SFS（`fastsimcoal2` / `∂a∂i`）
- 基因流检验：f3 / f4 / D-statistics（AdmixTools、Dsuite）
- 渗入片段：rIBD（beagle）、Loter、RFMix
- 多样性统计：ROH、近交系数、核苷酸多样性 π、亲缘系数（IBD）
- 其他选择信号：Tajima's D、iHS、CLR（SweeD）

每补一个，按本文档的三段式（步骤 → 分析代码 → 作图函数）追加一节。
