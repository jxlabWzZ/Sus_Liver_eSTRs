#!/bin/bash
chrom=$1
batches=$(ls /home/wuzhongzi/hip_str/03DNA_output/project525/all556/eSTR_1MB_F6/01.fdr/03.caviar/f6/chr${chrom}/batch* | awk -F"/" '{print $NF}')
for batch in $batches
do
mkdir -p /home/wuzhongzi/hip_str/03DNA_output/project525/all556/eSTR_1MB_F6/01.fdr/03.caviar/tmp/f6/${chrom}-${batch}/
#cd /home/wuzhongzi/hip_str/03DNA_output/project525/all556/eSTR_1MB_F6/01.fdr/03.caviar
python ./RunCaviarGTEx.py \
--zsnp ZSNP.tsv \
--zstr ZSTR.tsv \
--tissue f6 \
--samples /home/wuzhongzi/hip_str/RNAdata/14.rpkm/RNA556/f6_260.id \
--strgt /home/wuzhongzi/hip_str/03DNA_output/project525/all556/fa6_norm_str.vcf.gz \
--snpgt /home/wuzhongzi/hip_str/03DNA_output/project525/all556/SNP_data/fa6_c${chrom}_snpOnly_norm.vcf.gz \
--mingt 3 \
--recompute-z \
--tmpdir /home/wuzhongzi/hip_str/03DNA_output/project525/all556/eSTR_1MB_F6/01.fdr/03.caviar/tmp/f6/${chrom}-${batch}/ \
--num-causal 1 \
--use-topn-strs 1 \
--genes-file /home/wuzhongzi/hip_str/03DNA_output/project525/all556/eSTR_1MB_F6/01.fdr/03.caviar/f6/chr${chrom}/${batch} \
--out f6_chr${chrom}_${batch}_caviar.tab \
--expr /home/wuzhongzi/hip_str/RNAdata/14.rpkm/RNA556/PEER/F6/Corr_Expr.csv_f20
done
