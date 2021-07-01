#!/bin/bash
chrom=$1
cd /home/wuzhongzi/hip_str/03DNA_output/project525/all556/eSTR_1MB_F6/01.fdr/02.anova
python ./RunESTRAnova.py \
--sigsnps ZSNP.tsv \
--sigstrs ZSTR.tsv \
--samples /home/wuzhongzi/hip_str/RNAdata/14.rpkm/RNA556/f6_260.id \
--strgt /home/wuzhongzi/hip_str/03DNA_output/project525/all556/fa6_norm_str.vcf.gz \
--snpgt /home/wuzhongzi/hip_str/03DNA_output/project525/all556/SNP_data/fa6_c${chrom}_snpOnly_norm.vcf.gz \
--chrom chr${chrom} \
--mingt 3 \
--out f6_chr${chrom}_anova.tab \
--expr /home/wuzhongzi/hip_str/RNAdata/14.rpkm/RNA556/PEER/F6/Corr_Expr.csv_f20
