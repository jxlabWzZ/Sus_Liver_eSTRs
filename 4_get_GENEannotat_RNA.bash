### get gene annotation from gtf and expr data

# 1- Remove genes Ids from the Expression data ## gene names (symbol name)
less -S /home/wuzhongzi/hip_str/RNAdata/15.eSTR/rpkm_single_phe.bed.gz | awk '{print $1}' > ids
 
# 2- from symbol to search ensemble gene name

awk -F\\t 'FNR==NR {a[$2]=$0; next}; $1 in a {print a[$1]}' ~/hip_str/RNAdata/09.stringtie_merge/gene_1.95_all.bed ids | awk '{print $1}' > ids2

# 3- all gene from the gencode annotation version 1.95

less -S /home/wuzhongzi/sy_20191008/susref/STAR/Sus_scrofa.Sscrofa11.1.95.gtf | awk '$3 == "gene" {print $1","$4","$5","$10","$7","$1","$4","$5","$12","$7}' > Sus1.95.annot.expr.txt

# Then remove unwanted char ; " and "
sed 's/"//g' Sus1.95.annot.expr.txt > tmp
sed 's/;//g' tmp > Sus1.95.annot.expr.txt

# 4- Then we can get a well formatted Expression annotation file with only genes that are present in the expression file and in that same order

awk -F, 'FNR==NR {a[$4]=$0; next}; $1 in a {print a[$1]}' Sus1.95.annot.expr.txt ids2 > Liver548_Expr_annot.txt

## 行转列 awk -F "+" '{for(i=1;i<=NF;i++) a[i,NR]=$i}END{for(i=1;i<=NF;i++) {for(j=1;j<=NR;j++) printf a[i,j] ",";print ""}}' ids2 > ids3

## awk -F\\t 'FNR==NR {a[$2]=$0; next}; $1 in a {print a[$1]}' ../09.stringtie_merge/gene_1.95_all.bed rpkm548.gene | less -S