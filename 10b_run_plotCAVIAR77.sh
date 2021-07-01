
for a in `cat $1`
do
i=`echo $a |sed 's/_/\t/g' | cut -f 1`

mkdir ${i}
grep "$i" ~/hip_str/RNAdata/15.eSTR/GTex_eQTL/Liver548_Expr_annot.txt2 | sed 's/,/\t/g' |awk '{if($2>=1000000) print $1"\t"$2-1000000"\t"$3+1000000; else print $1"\t0\t"$3+1000000}' > ${i}/g.bed

m=`cat ${i}/g.bed | cut -f 1 | sed 's/chr//g'`
n=`cat ${i}/g.bed | cut -f 2`
p=`cat ${i}/g.bed | cut -f 3`

q=`echo $a |sed 's/_/\t/g' | cut -f 2`
d=`echo $a |sed 's/_/\t/g' | cut -f 3`

awk -v d=$d -v q=$q -v i=$i '{if($1==i && $2==q && $4==d) print $1"\t"$2"\t"$4"\t"$9}' ~/hip_str/03DNA_output/project525/all556/eSTR_1MB_F7/01.fdr/tmp > ${i}/str.bed7
#awk -v q=$q -v n=$n -v p=$p '{if($2==q && $4>=n && $4<=p) print $1"\t"$2"\t"$4"\t"$9}' ~/hip_str/03DNA_output/project525/all556/eSTR_1MB_F7/01.fdr/03.caviar/aasnp.txt | awk -v i=$i '{if($1==i) print $0}' > ${i}/snp.bed7
cat ~/hip_str/03DNA_output/project525/all556/eSTR_1MB_F6/01.fdr/08.plot/SNPeqtl/f7/out_chr* | awk -v q=$q -v n=$n -v p=$p '{if($2==q && $4>=n && $4<=p) print $1"\t"$2"\t"$4"\t"$9}' | awk -v i=$i '{if($1==i) print $0}' > ${i}/snp.bed7
cat ${i}/str.bed7 ${i}/snp.bed7 > ${i}/ttt_7
less -S ~/hip_str/03DNA_output/project525/all556/eSTR_1MB_F7/01.fdr/03.caviar_rep/tmp/f7/*/${i}/CAVIAR_post | grep -v "SNP_ID" | sed 's/:/\t/g' | cut -f 2,4 > ${i}/cs7.bed
#awk '{print $3"\t"$0}' ${i}/ttt | awk -F\\t 'FNR==NR {a[$1]=$0; next}; $1 in a {print a[$1]"\t"$0}' - ${i}/cs6.bed | cut -f 2-5,7 > ${i}/aaaa
bedtools intersect -a gene_list_gff3.bed -b ${i}/g.bed -wa | grep -v "protein" | awk '{print $4","$2","$3","$5}' > ${i}/g.csv
Rscript rr_plotCAVIAR77.R /home/wuzhongzi/hip_str/03DNA_output/project525/all556/eSTR_1MB_F6/01.fdr/01.topSTR/tmp $i 
done



