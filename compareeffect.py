tissu = ['F7']

import matplotlib.pyplot as plt
import pandas as pd

##Old study
oldestr = pd.read_csv('/home/wuzhongzi/hip_str/03DNA_output/project525/all556/eSTR_1MB_F6/01.fdr/fdr1.out',sep=',')
c=['','grey','b','g','black','r','purple','blue' ]; i=0

for Tissue in tissu:
    i=i+1
    ##Get current study
    newestr = pd.read_csv("/home/wuzhongzi/hip_str/03DNA_output/project525/all556/eSTR_1MB_"+Tissue+"/01.fdr/fdr1.out",sep='\t')
        ##Reformat gene name without version to match old
    G2 = [x.split('.')[0] for x in list(newestr['gene'])]           
    newestr['gene']=G2

    ##only select gene present in both dataset
    NEW = newestr.loc[newestr['gene'].isin(oldestr['gene'])]     #New data in old
    OLD = oldestr.loc[oldestr['gene'].isin(NEW['gene'])]         #Old data in new

    ##adjusting tart position in old data
    #OLD['str.start'] = OLD['str.start']-1 
    #OLD['str.start'] = list(OLD['start'])
    
    #Select 3 col of interest 
    Int=['gene', 'str.start', 'beta','p.wald']
    G0= OLD.loc[:,Int]      #X_o.loc[:,Int]
    GN= NEW.loc[:,Int]      #X_n.loc[:,Int]

    #G0=  G0[G0['signif.estr']==True]
    newrows={x['gene']+'-'+str(int(x['str.start'])) : x['beta'] for y,x in GN.iterrows()}
    oldrows={x['gene']+'-'+str(int(x['str.start'])): x['beta'] for y,x in G0.iterrows()}
    #print(len(newrows.keys()), '\t',len(oldrows.keys()))
    Key = [x for x in newrows.keys() if x in oldrows.keys()]
    #print (len(Key))
    Xn = [newrows[x] for x in Key]
    Yo = [oldrows[x] for x in Key]

    #Tissue vs nature paper 
    plt.scatter(Xn,Yo,color=c[i],marker='.')#,edgecolors='b')
    #plt.xlim(-1, 1)
    #plt.ylim(-1, 1)
    plt.ylabel('eSTR effect size from F6')
    plt.xlabel('eSTR effect size from'+Tissue)
    plt.title('Effect sizes agreement between '+Tissue+' tissue \n and F6. Restricted at eSTRs')
    plt.axis('equal')
    plt.axhline(y=0)
    plt.axvline(x=0)
    plt.grid()
    plt.show()
    plt.savefig("/home/wuzhongzi/hip_str/03DNA_output/project525/all556/eSTR_1MB_F6/01.fdr/tmpf6f7.png")
    CC=[Xn[i]*Yo[i] for i in range(len(Xn))]
    print(len([x for x in CC if x>0]),"=====",len(Xn), '+++', len(Yo))
    print()
