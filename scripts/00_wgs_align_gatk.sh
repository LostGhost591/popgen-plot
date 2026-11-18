#!/usr/bin/env bash
# =============================================================================
# 00_wgs_align_gatk.sh —— WGS 上游：比对 → 去重 → 质控报告 → GATK 变异检测
#
# 来源：由「重测序数据分析.sh」通用化改造，软件参数与原脚本完全一致，
#       只把写死的路径 / 样本列表 / 输出目录换成了变量。
#
# 用法：
#   bash 00_wgs_align_gatk.sh                # 用配置区的 SAMPLES 或 SAMPLES_FILE
#   bash 00_wgs_align_gatk.sh A1 A2 A3       # 临时指定样本（覆盖配置区）
#   RUN_STAGE=1 bash 00_wgs_align_gatk.sh    # 只跑第 1 阶段（1=比对 2=质控 3=GATK all=全跑）
#
# 设计约定：
#   1. 只有「配置区」需要改；下面所有路径都从这几个变量派生，不出现物种名 / 项目名
#   2. 每一步失败立即退出，不留半成品（要跳过就用 RUN_STAGE 分阶段跑）
#   3. 样本数任意，CombineGVCFs 的 --variant 列表自动生成
# =============================================================================

set -euo pipefail

# ------------------------------- 配置区（只改这里）------------------------------
WORKDIR=/data/sda/zhou/aixiao            # 项目根目录
DATADIR=$WORKDIR/A                       # 样本 clean data 所在目录
REFDIR=/data/sda/zhou/chickengenomeindex6a/refbwaindex

REF_NAME=Gallus_gallus.GRCg6a.dna.toplevel.fa          # 参考基因组文件名
GFF_NAME=Gallus_gallus.GRCg6a.104.chr.gff3             # GFF（下游注释用）
GTF_NAME=Gallus_gallus.GRCg6a.104.chr.gtf              # GTF（下游注释用）
DICT_NAME=Gallus_gallus.GRCg6a.dna.toplevel.dict       # Picard 字典（GATK 用）

SAMPLES=(A1 A2 A3 A4 A5 A6 A7 A8 A9)     # 样本名列表
# SAMPLES_FILE=$DATADIR/samples.list     # 也可以改成从文件读（一行一个，取消注释即可）

CLEAN_R1=.Clean.R1.fq.gz                 # clean data 命名：$DATADIR/${s}${CLEAN_R1}
CLEAN_R2=.Clean.R2.fq.gz

THREADS_ALN=1        # bwa mem -t（原脚本写死 1，机器够就调大）
THREADS_SORT=4       # samtools sort -@
MEM_SORT=4G          # samtools sort -m
JAVA_MEM=4g          # gatk --java-options -Xmx
THREADS_QC=2         # qualimap -nt
QC_DEPTH=10000       # qualimap -nr

RUN_STAGE=all        # 1 / 2 / 3 / all
KEEP_BAM=1           # 1=保留中间 bam；0=排序后删掉未排序与排序版，只留 rmdup
# ------------------------------------------------------------------------------

# ------------------------------- 派生变量（不用改）------------------------------
REF=$REFDIR/$REF_NAME
FAI=$REF.fai
DICT=$REFDIR/$DICT_NAME
GFF=$REFDIR/$GFF_NAME
GTF=$REFDIR/$GTF_NAME

TMPDIR=$WORKDIR/tmp
MAPDIR=$DATADIR/3.map
BWAOUT=$MAPDIR/bwa
QCOUT=$MAPDIR/result
VARDIR=$DATADIR/4.snp_indel/GATK

# ------------------------------------------------------------------------------

die() { echo "[ERROR] $*" >&2; exit 1; }
say() { echo ">>> $*"; }

# 样本：命令行参数优先，其次 SAMPLES_FILE，最后用配置区数组
if [ $# -gt 0 ]; then
  SAMPLES=("$@")
elif [ -n "${SAMPLES_FILE:-}" ]; then
  [ -f "$SAMPLES_FILE" ] || die "样本文件不存在：$SAMPLES_FILE"
  mapfile -t SAMPLES < <(grep -v '^[[:space:]]*$' "$SAMPLES_FILE")
fi
[ ${#SAMPLES[@]} -gt 0 ] || die "样本列表为空"

# 依赖检查：缺软件直接报错，不往后跑
for t in bwa samtools gatk qualimap tabix bgzip; do
  command -v "$t" >/dev/null 2>&1 || die "缺少软件：$t"
done

# 参考文件检查
for f in "$REF" "$FAI" "$DICT"; do
  [ -s "$f" ] || die "参考文件缺失：$f"
done

mkdir -p "$TMPDIR"

echo "============================================================"
echo " 样本数 : ${#SAMPLES[@]}  (${SAMPLES[*]})"
echo " 参考   : $REF"
echo " 输出   : $BWAOUT / $QCOUT / $VARDIR"
echo " 阶段   : $RUN_STAGE"
echo "============================================================"

# ============================== 阶段 1：比对 ==================================
run_align() {
  mkdir -p "$BWAOUT"
  cd "$BWAOUT" || die "进不去 $BWAOUT"

  for s in "${SAMPLES[@]}"; do
    R1=$DATADIR/${s}${CLEAN_R1}
    R2=$DATADIR/${s}${CLEAN_R2}
    [ -s "$R1" ] || { echo "[跳过] $s 缺 R1：$R1" >&2; continue; }
    [ -s "$R2" ] || { echo "[跳过] $s 缺 R2：$R2" >&2; continue; }

    say "[$s] bwa mem"
    bwa mem "$REF" "$R1" "$R2" -t "$THREADS_ALN" -M \
        -R "@RG\tID:${s}\tLB:${s}\tPL:ILLUMINA\tSM:${s}" \
      | samtools view -bS -h - > "${s}.bam"

    say "[$s] sort"
    samtools sort -@ "$THREADS_SORT" -m "$MEM_SORT" "${s}.bam" -O bam -o "${s}.sorted.bam"

    say "[$s] dedup"
    # 注：samtools >= 1.21 已移除 rmdup。若你的版本没有 rmdup，改用：
    #   samtools sort -n -@ $THREADS_SORT "${s}.bam" -o "${s}.qsort.bam"
    #   samtools fixmate -m -@ $THREADS_SORT "${s}.qsort.bam" "${s}.fixmate.bam"
    #   samtools sort -@ $THREADS_SORT "${s}.fixmate.bam" -o "${s}.pos.bam"
    #   samtools markdup -r -@ $THREADS_SORT "${s}.pos.bam" "${s}.rmdup.bam"
    samtools rmdup "${s}.sorted.bam" "${s}.rmdup.bam"
    samtools index "${s}.rmdup.bam"

    if [ "$KEEP_BAM" -eq 0 ]; then
      rm -f "${s}.bam" "${s}.sorted.bam"
    fi
  done
}

# ============================== 阶段 2：质控报告 ==============================
run_qualimap() {
  mkdir -p "$QCOUT"

  for s in "${SAMPLES[@]}"; do
    BAM=$BWAOUT/${s}.rmdup.bam
    [ -s "$BAM" ] || { echo "[跳过] $s 缺去重 bam：$BAM" >&2; continue; }

    say "[$s] qualimap bamqc"
    qualimap bamqc --java-mem-size="$JAVA_MEM" -nt "$THREADS_QC" -nr "$QC_DEPTH" \
        -bam "$BAM" -outfile report.pdf -outformat PDF \
        -outdir "$QCOUT/${s}.depth_and_coverage_report"
  done
}

# ============================== 阶段 3：GATK ==================================
run_gatk() {
  mkdir -p "$VARDIR"
  cd "$VARDIR" || die "进不去 $VARDIR"

  # 3.1 单样本 gVCF
  for s in "${SAMPLES[@]}"; do
    BAM=$BWAOUT/${s}.rmdup.bam
    [ -s "$BAM" ] || die "$s 缺去重 bam：$BAM（先跑阶段 1）"

    say "[$s] HaplotypeCaller"
    gatk --java-options "-Xmx${JAVA_MEM}" HaplotypeCaller -R "$REF" \
        -I "$BAM" -O "${s}.g.vcf.gz" \
        --max-alternate-alleles 4 --sample-ploidy 2 \
        -ERC GVCF --tmp-dir "$TMPDIR"
  done

  # 3.2 合并 gVCF —— --variant 列表按样本自动生成，样本数任意
  VARGVCFS=()
  for s in "${SAMPLES[@]}"; do
    VARGVCFS+=(--variant "${s}.g.vcf.gz")
  done

  say "CombineGVCFs（${#SAMPLES[@]} 个样本）"
  gatk --java-options "-Xmx${JAVA_MEM}" CombineGVCFs -R "$REF" \
      "${VARGVCFS[@]}" -O all.g.vcf.gz --tmp-dir "$TMPDIR"

  say "GenotypeGVCFs"
  gatk --java-options "-Xmx${JAVA_MEM}" GenotypeGVCFs -R "$REF" \
      -V all.g.vcf.gz -O all.raw.vcf.gz --tmp-dir "$TMPDIR"

  # 3.3 SNP 硬过滤（表达式与原脚本一致）
  SNP_FILTER="QD < 2.0 || MQ < 40.0 || FS > 60.0 || SOR > 3.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0"
  INDEL_FILTER="QD < 2.0 || FS > 200.0 || SOR > 10.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0"

  say "SNP 过滤"
  gatk SelectVariants -select-type SNP -V all.raw.vcf.gz -O snp.vcf.gz
  gatk VariantFiltration -V snp.vcf.gz \
      --filter-expression "$SNP_FILTER" \
      --filter-name "Filter" \
      -O snp.filter.vcf.gz
  gatk SelectVariants -R "$REF" -V snp.filter.vcf.gz --exclude-filtered -O snp_pass.vcf.gz
  tabix -p vcf snp_pass.vcf.gz

  say "INDEL 过滤"
  gatk SelectVariants -select-type INDEL -V all.raw.vcf.gz -O indel.vcf.gz
  gatk VariantFiltration -V indel.vcf.gz \
      --filter-expression "$INDEL_FILTER" \
      --filter-name "Filter" \
      -O indel.filter.vcf.gz
  gatk SelectVariants -R "$REF" -V indel.filter.vcf.gz --exclude-filtered -O indel_pass.vcf.gz
  tabix -p vcf indel_pass.vcf.gz

  echo "------------------------------------------------------------"
  echo " 完成：$VARDIR"
  echo "   SNP   : snp_pass.vcf.gz"
  echo "   INDEL : indel_pass.vcf.gz"
  echo "------------------------------------------------------------"
}

# ============================== 入口 ==========================================
case "$RUN_STAGE" in
  1)   run_align ;;
  2)   run_qualimap ;;
  3)   run_gatk ;;
  all) run_align; run_qualimap; run_gatk ;;
  *)   die "RUN_STAGE 只能是 1 / 2 / 3 / all，当前：$RUN_STAGE" ;;
esac
