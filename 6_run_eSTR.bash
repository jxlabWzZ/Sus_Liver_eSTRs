cd /home/wuzhongzi/hip_str/03DNA_output/project525/all556/eSTR_1MB_F6

for i in {1..18}
do
/home/wuzhongzi/miniconda3/envs/py27/bin/python ~/myproject/scripts/LinRegAssociationTest_v2_wzz.py \
--expr /home/wuzhongzi/hip_str/RNAdata/14.rpkm/RNA556/PEER/F6/Corr_Expr.csv_f20 \
--exprannot ~/hip_str/RNAdata/15.eSTR/GTex_eQTL/Liver548_Expr_annot.txt2 \
--strgt /home/wuzhongzi/hip_str/03DNA_output/project525/all556/fa6_norm_str.vcf.gz \
--norm --out ./estr1Mb_chr${i}_f0 \
--distfromgene 1000000 --chrom chr${i}
done

