import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

tissues = ['F6','F7']
colors = ['green','grey']
genes = ['ENSSSCG00000024827','ENSSSCG00000030395', 'ENSSSCG00000003653','ENSSSCG00000007982']

for G in genes:
    fig = plt.figure(figsize=(6, 4))
    ax = fig.add_subplot(111)
    i=0

    for tissu in tissues[:2]:
        #print(tissu)
        File="/home/wuzhongzi/hip_str/03DNA_output/project525/all556/eSTR_1MB_"+tissu+"/01.fdr/tmp"   #'PQValues.txt'
        File1=pd.read_csv(File, sep='\t')
        SET = File1.loc[File1['gene']==G]
        ax.scatter(SET['str.start']/1e6, -np.log10(SET['p.wald']),color=colors[i], label=tissues[i])
        #print(SET['str.start','p.wald'])
        i=i+1
        #print(tissu)
    ax.set_xlabel('Genomic position on '+G+' (Mb)')
    ax.set_ylabel('-log10(pvalue)')
    plt.legend(colors, labels = tissues,loc='upper right',fontsize=8)
    plt.show()
    plt.savefig("/home/wuzhongzi/hip_str/03DNA_output/project525/all556/eSTR_1MB_F6/01.fdr/tmp_"+G+"_a.png")
    
