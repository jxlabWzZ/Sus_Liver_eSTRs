import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

tissues = ['F6']
tissu = tissues[0]

#OUT3=pd.read_csv("/storage/szfeupe/Runs/GTEx_estr/Analysis_by_Tissue/SUMMARY/LR_SummaryTest_Table.tsv", sep='\t')

for tissu in tissues:
    OUT = pd.read_csv("/home/wuzhongzi/hip_str/03DNA_output/project525/all556/eSTR_1MB_"+tissu+"/01.fdr/tmp", '\t')
    #CHH = OUT.loc[OUT['chrom']=='chr'+str(i)]['str.start']
    n=0; p=1
    colors=['black','grey']
    #colors = ['green','blue','cyan','red','pink','black','grey','yellow']
    for i in range(1,19):
        CHH = OUT.loc[OUT['chrom']=='chr'+str(i)]
        #X=[x for x in range(n,CHH.shape[0]+n)]
        XP=np.array(CHH['str.start'])/1e6+n
        plt.plot(XP, -np.log10(CHH['p.wald']), ls='', marker='.', linewidth=8)    #many colors
        #plt.scatter(X, -np.log10(CHH['p.wald']),color=colors[p], marker='.') #alternate 2 colors
        plt.xticks([])
        #n=X[-1]
        n=XP[-1]
        if p==0: 
            p=1
        else:
            p=0
    plt.xlabel('Chromosomes')
    plt.ylabel('-log10(pvalue)')
    plt.title('Manhattan plot for '+tissu)
    plt.show()
    plt.savefig("/home/wuzhongzi/hip_str/03DNA_output/project525/all556/eSTR_1MB_"+tissu+"/01.fdr/tmp6.png")
    plt.show()
    plt.savefig("/home/wuzhongzi/hip_str/03DNA_output/project525/all556/eSTR_1MB_"+tissu+"/01.fdr/tmp6.pdf")
