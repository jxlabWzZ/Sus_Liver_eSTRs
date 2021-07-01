###
### step one: genotype STR
###

#!/bin/bash
#PBS -N STR_hip2
#PBS -l nodes=1:ppn=10
#PBS -e /home/wuzhongzi/qsublog
#PBS -o /home/wuzhongzi/qsublog
#PBS -q cu
#PBS -t 0-200

wd=/home/wuzhongzi/hip_str/DNA_output/
cd $wd

myregion=(`ls /home/wuzhongzi/hip_str/pip_reference_bed/*.bed`)
## ref_STR region were divided into 200 subfiles

chrom=`echo ${myregion[$PBS_ARRAYID]} | awk '{split($0,arra,"tmp2/"); print arra[2]}'`

/home/wuzhongzi/biosoft/HipSTR/HipSTR --bams ${bam_list} \
--fasta /home/wuzhongzi/hip_str/pip_reference_bed/all_chroms.fa \
--regions ${myregion[$PBS_ARRAYID]} \
--str-vcf myliver_${chrom}.vcf.gz \
--max-str-len 150 \
--log out_${chrom}.log

## bam_list=BAM1,BAM2,BAM3,...,BAM548

###
### step two: filter STR
###

/home/wuzhongzi/miniconda3/envs/py27/bin/python ~/biosoft/HipSTR/scripts/filter_vcf.py --vcf myliver_merged.vcf \
                             --min-call-depth 5 \
                             --min-call-qual 0.9 \
                             --max-call-flank-indel 0.15 \
                             --max-call-stutter 0.15 \
                             --min-call-allele-bias -2 \
                             --min-call-strand-bias -2 \
                             --min-loc-calls 100 > myliver_merged_filter.vcf
							 
## 548 samples and approximately 1.15 million STR loci
## Illumina ~7x, 150bp pair end reads.
