chrom=$1
cd /home/wuzhongzi/hip_str/03DNA_output/project525/all556/eSTR_1MB_F6/01.fdr/04.heritability

python2 /home/wuzhongzi/biosoft/gtex-estrs/Scripts/Heritability/strsnp_h2_gcta2.py \
--expr /home/wuzhongzi/hip_str/RNAdata/14.rpkm/RNA556/PEER/F6/Corr_Expr.csv_f20 \
--exprannot ~/hip_str/RNAdata/15.eSTR/GTex_eQTL/Liver548_Expr_annot.txt2 \
--distfromgene 1000000 \
--strgt /home/wuzhongzi/hip_str/03DNA_output/project525/all556/fa6_norm_str.vcf.gz \
--snpgt /home/wuzhongzi/hip_str/03DNA_output/project525/all556/SNP_data/fa6_c${chrom}_snpOnly_norm.vcf.gz \
--chrom ${chrom} \
--out f6ok_all5_${chrom}.out.txt \
--tmpdir ./tmp \
--method GCTA \
--include-str FE \
--estr-results fdr1.txt \
--esnps esnp2.txt \
--samples id6.txt
