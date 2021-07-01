# 提取指定个体 --keep
# 删除样本 --remove 
# get f6 samples STR
/home/yangbin/bin/vcftools_0.1.13/bin/vcftools --vcf myliver_merged_filter.vcf --recode --recode-INFO-all --stdout --keep f6.id > f6_str.vcf
# get f7 samples STR 
/home/yangbin/bin/vcftools_0.1.13/bin/vcftools --vcf myliver_merged_filter.vcf --recode --recode-INFO-all --stdout --keep f7.id > f7_str.vcf

## to norm STR GB (f6)
/home/wuzhongzi/miniconda3/envs/py27/bin/python 1.2_STR_normalization.py \
--vcf f6_str.vcf \
--out f6_norm_str.vcf \
--gtfield "GB" \
--minsamples 50 \
--mincount 3 \
--mingenotypes 2 

## to norm STR GB (f7)
/home/wuzhongzi/miniconda3/envs/py27/bin/python 1.2_STR_normalization.py \
--vcf f7_str.vcf \
--out f7_norm_str.vcf \
--gtfield "GB" \
--minsamples 50 \
--mincount 3 \
--mingenotypes 2 


