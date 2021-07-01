#%pylab inline

# Params
DATADIR = "/home/wuzhongzi/hip_str/03DNA_output/project525/all556/eSTR_1MB_F6/01.fdr/08.plot/"
RESULTSDIR = "/home/wuzhongzi/hip_str/03DNA_output/project525/all556/eSTR_1MB_F6/01.fdr/08.plot/"

# Allow us to edit fonts in Illustrator
import matplotlib
matplotlib.rcParams['ps.useafm'] = True
matplotlib.rcParams['pdf.use14corefonts'] = True
matplotlib.rcParams['text.usetex'] = True

# Import libraries
import os
import pandas as pd
import scipy.stats
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

def ExtractData(gene, chrom, start, tissue, newdata=True):
    if "(" in tissue: tissue = tissue.split("(")[0]
    vcf = "/home/wuzhongzi/hip_str/03DNA_output/project525/all556/a200_fa7.vcf.gz"
    # Pull out STR genotypes
    cmd = """bcftools query -r %s:%s-%s -f"[%%SAMPLE\\t%%GB\\n]" %s | \
    grep -v "\." | sed 's/|/\\t/' | awk '{print $1 "\\t" $2 "\\t" $3}' > str_genotypes.tab"""%(chrom, start, start, vcf)
    os.system(cmd)
    # Pull out gene expression
    if newdata:
        expr="/home/wuzhongzi/hip_str/RNAdata/14.rpkm/RNA556/PEER/%s/Corr_Expr.csv_f20"%tissue
    else: expr="/home/wuzhongzi/hip_str/RNAdata/14.rpkm/RNA556/PEER/%s/Corr_Expr.csv_f20"%tissue
    allgenes = open(expr,"r").readline().split(",")
    colnum = allgenes.index('"' + gene + '"')+2
    cmd = """cat %s | cut -d',' -f 1,%s | grep -v ENSG > expr.tab"""%(expr, colnum)
    os.system(cmd) 
 
 ####### Great example ########
gene = "ENSSSCG00000039506"
chrom = "chr5"
start = 21981597
tissue = "F7"
reflen = 0 # TODO
period = 1 # TODO

# linear p: 0.15069674861460525     quadratic p: 3.235125855915583e-11
###############




# Extract data
ExtractData(gene, chrom, start, tissue)
strgt = pd.read_csv("str_genotypes.tab", sep="\t", names=["sample","str1","str2"])
strgt["sample"] = strgt["sample"].apply(lambda x: "-".join(x.split("-")[0:2]))
expr = pd.read_csv("expr.tab", names=["sample","expr"])
data = pd.merge(strgt, expr)
data["str"] = data["str1"]+data["str2"]
data['expr'] = data['expr'].apply(pd.to_numeric)


########## Linear #############
boxcol = "gray"
fig = plt.figure()
fig.set_size_inches((8,5))
ax = fig.add_subplot(111)

sns.swarmplot(x="str", y="expr", ax=ax, data=data, color="black", zorder=0)
sns.boxplot(x="str", y="expr", ax=ax, data=data, color="white", linewidth=0.5, 
        boxprops={'facecolor':'None', 'edgecolor': boxcol}, showfliers=False)
# Set box properties
for i,artist in enumerate(ax.artists):
    artist.set_edgecolor(boxcol)
    artist.set_facecolor('None')

    # Each box has 6 associated Line2D objects (to make the whiskers, fliers, etc.)
    # Loop over them here, and use the same colour as above
    x = 5
    for j in range(i*x,i*x+x):
        line = ax.lines[j]
        line.set_color(boxcol)
        line.set_mfc(boxcol)
        line.set_mec(boxcol)
totals = data.groupby("str", as_index = False).agg({"expr": len})
print(totals.sort_values("str"))
means = data.groupby("str", as_index = False).agg({"expr": ['mean']})
means = means.sort_values("str")
means["num"] = range(means.shape[0])
    
ax.plot(means["num"], means["expr"], color="red", marker="o", linewidth=2)
ax.set_xlabel("Mean STR Dosage.", size=15)
ax.set_ylabel("Normalized Expression - %s"%tissue, size=15)
ax.set_xticklabels(["%.1f"%((item*0.5+reflen)/period) for item in sorted(list(set(data["str"])))], size=12)
ax.set_yticklabels(["%.1f"%(item) for item in ax.get_yticks()], size=12)
ax.set_title("")
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.yaxis.set_ticks_position('left')
ax.xaxis.set_ticks_position('bottom')
ax.set_title("%s %s %s"%(chrom, start,gene))
plt.suptitle("");
plt.show()
plt.savefig("/home/wuzhongzi/hip_str/03DNA_output/project525/all556/eSTR_1MB_F6/01.fdr/08.plot/chipseq77/f7_%s.pdf"%gene)

########## Non-linear (hold one constant) #############
#alleles = set(list(data["str1"])+list(data["str2"]))
#fig = plt.figure()
#ax = fig.add_subplot(111)
#for al in alleles:
#    adata = data[(data["str1"]==al) | (data["str2"]==al)].copy()
#    adata["x"] = adata.apply(lambda x: [x["str1"],x["str2"]][int(x["str1"]==al)], 1) # if str1=al, use str2, else use str2
#    xx = adata.groupby("x", as_index=False).agg({"expr": np.mean,"str1": len})
#    xx = xx[xx["str1"]>=3]
#    ax.plot(xx["x"],xx["expr"], label=al)
#fig.legend(title="Allele 1")
#ax.set_xlabel("Allele 2")
#ax.set_ylabel("Mean Expression")
#plt.show()
plt.savefig("/home/wuzhongzi/hip_str/03DNA_output/project525/all556/eSTR_1MB_F6/01.fdr/08.plot/chipseq77/f7_%s.png"%gene)

