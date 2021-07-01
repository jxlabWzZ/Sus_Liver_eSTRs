#Libraries
#/home/wuzhongzi/miniconda3/envs/py27/bin/python
import argparse
import math
import numpy as np
import os
from collections import Counter
import pandas as pd
import random
import shutil
import statsmodels.api as sm
import statsmodels.formula.api as smf
import sys
#Variables
EXPRFILE = None
EXPRANNOTFILE = None
DISTFROMGENE = None
STRGTFILE = ""
OUTFILE = ""
SCRIPTDIR = "/home/wuzhongzi/hip_str/RNAdata/15.eSTR/GTex_eQTL/"
TMPDIR = "/home/wuzhongzi/hip_str/RNAdata/15.eSTR/GTex_eQTL/"
CHECKCHROM = False
PERMUTE_EXPR = False
DEBUG = False
MINSAMPLES = 50
MINGENOTYPE = 3
LINEAR_ONLY = False
LINEAR=False
QUAD=False
ALLELE=False
No_ANOVA=False
NORM=True

#Functions
def PROGRESS(msg):
    sys.stderr.write("%s\n"%msg.strip())

def ZNorm(vals,m,sd):
    if m is None:
        m = np.mean(vals)
        sd = math.sqrt(np.var(vals))
    if sd == 0: return None
    return [(item-m)/sd for item in vals]

def GetOutlierLimit(X):  
    # Remove outlier beyond 2xSDs
    M = np.mean(X)
    sd= math.sqrt(np.var(X))
    return M+(5*sd), M-(5*sd)

def LinearRegression(data, Y, norm=False, minsamples=0, alleles=False):
    """
    Perform linear regression, return beta, beta_se, p
    """
    print alleles, norm
    if norm:
        if alleles:
            data['x1']=data['x1'].astype(float) ; data['x2']=data['x2'].astype(float)
            data['x1+x2'] = data[['x1', 'x2']].sum(axis=1)
            X = ZNorm(data['x1+x2'], None, None)
        else:
            X = ZNorm(data, None, None) 
            Y = ZNorm(Y, None, None)
        if X is None or Y is None: return None, None, None
        if np.var(X)==0: return None, None, None
        if len(X) <= minsamples: return None, None, None
    else:
        X = data
        print 'Here!!!!!!'
    X = sm.add_constant(X)
    mod_ols = sm.OLS(Y, X, missing='drop')
    res_ols = mod_ols.fit()
    pval = res_ols.pvalues[1]
    print 'P-value: ', pval 
    slope = res_ols.params[1]
    err = res_ols.bse[1]
    return res_ols, slope, err, pval

def QuadraticRegression(data,norm=True, minsamples=120):
    """
    Perform non linear regression, return beta, p_beta, alpha, p_alpha, se_alpha, se_beta
    data.columns has to have [allele1='x1', allele2='x2', gene_expression='expr'] not normalized
    """
    if data.shape[0] <= minsamples: 
        return None, None, None, None, None, None, None
    m = np.mean(list(locus_str["x1"].astype(int))+list(locus_str["x2"].astype(int))) 
    sd= math.sqrt(np.var(locus_str["x1"].astype(int) + locus_str["x2"].astype(int) ))
    if norm:
        #data['x1'] = ZNorm(data['x1'].astype(int), sd, m)  
        #data['x2'] = ZNorm(data['x2'].astype(int), sd, m)
        #data['expr'] = ZNorm(data['expr'].astype(float), None, None)
        data.loc[:,'x1'] =ZNorm(data['x1'].astype(int), sd, m)
        data.loc[:,'x2'] =ZNorm(data['x2'].astype(int), sd, m)       
        data.loc[:,'expr'] =ZNorm(data['expr'].astype(int), sd, m)
        
        if data['x1'].isnull().all()  or data['x2'].isnull().all() or data['expr'].isnull().all(): 
            return None, None, None, None, None, None, None  
        
    data['X']=data['x1']+data['x2']
    data['X2']=data['x1']**2 + data['x2']**2
    #mod_ols = sm.OLS('expr ,X ,X2', data = data).fit()
    mod_ols = smf.ols(formula = 'expr ~ X + X2', data = data).fit()

    alpha = list(mod_ols.params)[1]
    beta = list(mod_ols.params)[2]
    alpha_p = list(mod_ols.pvalues)[1]
    beta_p  = list(mod_ols.pvalues)[2]
    alpha_se=list(mod_ols.bse)[1]
    beta_se=list(mod_ols.bse)[2]
    return mod_ols, alpha, alpha_se, alpha_p, beta, beta_se, beta_p

def Modelcompare(model1, model2):
    """
    Performs ANOVA test to compares 2 models
    Here, model1 = Linear regression ,  model2=Quad regression
    Outputs rsquares , delta AIC (aic_model1- aic_model2) and pvalue
    """
    rsq1 = model1.rsquared
    rsq2 = model2.rsquared
    anova_output = sm.stats.anova_lm(model1 , model2)
    anova_pval = anova_output["Pr(>F)"].values[1]
    delta_aic = model1.aic - model2.aic 
    delta_bic = model1.bic - model2.bic
    return rsq1, rsq2, delta_aic,delta_bic, anova_pval

#arguments help
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="get datasets for LMM assocation tests with STRs and SNPs")
    parser.add_argument("--expr", help="Normalized expression residuals", type=str, required=True)
    parser.add_argument("--exprannot", help="Expression annotation file", type=str, required=True)
    parser.add_argument("--chrom", help="Restrict analysis to this chromosome", type=str, required=False)
    parser.add_argument("--strgt", help="File with noramlized STR genotypes", type=str, required=True)
    parser.add_argument("--scriptdir", help="Directory containing scripts", type=str, required=False)
    parser.add_argument("--checkchrom", help="Only load exons for relevant chromosome", action="store_true")
    parser.add_argument("--tmpdir", help="Tmp directory", type=str)
    parser.add_argument("--permute", help="Permute expression values", action="store_true")
    parser.add_argument("--norm", help="Normalize genotypes before doing association", action="store_true")
    parser.add_argument("--min-samples", help="Require data for this many samples", type=int, default=0)
    parser.add_argument("--min-genotypes", help="Require this many genotypes per test", type=int, default=3)
    parser.add_argument("--out", help="Write data files to this file", type=str, required=True)
    parser.add_argument("--debug", help="Print debug info", action="store_true")
    parser.add_argument("--alleles", help="Input genotypes entered are in the allele1,allele2 form", action="store_true")
    parser.add_argument("--linear_only", help="Linear regression will be performed", action="store_true")
    parser.add_argument("--linear", help="Linear regression will be performed", action="store_true")
    parser.add_argument("--quadratic", help="Non linear regression will also be performed", action="store_true")
    parser.add_argument("--distfromgene", help="Look at variants within this distance of gene boundaries", type=int, required=True)
    parser.add_argument("--noanova", help="In case of quadratic option ANOVA test will not be performed to determine the best model ", action="store_true")
    
#arguments 
    args = parser.parse_args()
    EXPRFILE = args.expr
    EXPRANNOTFILE = args.exprannot
    CHROM = args.chrom
    DISTFROMGENE = args.distfromgene
    STRGTFILE = args.strgt
    OUTFILE = args.out
    MINSAMPLES = args.min_samples
    MINGENOTYPE = args.min_genotypes
    if "chr" not in str(CHROM): CHROM="chr%s"%CHROM
    if args.norm: NORM = True
    if args.scriptdir is not None: SCRIPTDIR = args.scriptdir
    if args.checkchrom: CHECKCHROM = True
    if args.tmpdir is not None: TMPDIR = args.tmpdir
    if args.permute: PERMUTE_EXPR = True
    if args.linear_only: LINEAR_ONLY = True
    if args.linear: LINEAR = True
    if args.debug: DEBUG = True
    if args.alleles: 
            ALLELE = True
    if args.quadratic: 
        QUAD = True
        if not ALLELE:
            sys.stderr.write("ERROR: Input genotypes entered should be in the form (allele1,allele2) for quadratic regression\n OR --alleles option missing \n")  ;   sys.exit(1)            
        if args.noanova: No_ANOVA = True
###
    if not QUAD:
        LINEAR_ONLY = True
    # Load expression values
    PROGRESS("Load expression : ALLELE %s"%ALLELE)
    if CHECKCHROM:
        x = list(pd.read_csv(EXPRFILE, nrows=1).columns.values)
        #x = [item for item in x if item == "Unnamed: 0" or item == CHROM]
        expr = pd.read_csv(EXPRFILE, usecols=x) # reading all the columns takes way too much memory, pull out only ones needed
    else:
        expr = pd.read_csv(EXPRFILE)
    if "Unnamed: 0" in expr.columns:            ##Often seen in csv created with Pandas
        expr.index = expr["Unnamed: 0"].values
        expr = expr.drop("Unnamed: 0", 1)
        
    # Load expression annotation
    PROGRESS("Load annotation")
    expr_annot = pd.read_csv(EXPRANNOTFILE)
    expr_annot.index = expr_annot["gene"].values
    expr_annot = expr_annot.loc[[item for item in expr.columns if item in expr_annot.index],:]
    expr_annot = expr_annot[expr_annot["chr"] == CHROM]

    # Load STR genotypes
    PROGRESS("Load STRs")
    strgt = pd.read_csv(STRGTFILE, sep="\t", low_memory=False)
    strgt = strgt[strgt["chrom"] == CHROM]

    # Restrict to STR samples
    str_samples = list(set(strgt.columns[2:].values))
    samples_to_remove = []
    for item in str_samples:
        if item not in expr.index: samples_to_remove.append(item) #str_samples.remove(item)
    for item in samples_to_remove: str_samples.remove(item)
    expr = expr.loc[str_samples,:]
    PROGRESS("There are %s samples"%str(strgt.shape))

    #Set up output file
    f = open(OUTFILE, "w")
    if LINEAR_ONLY:
        f.write("\t".join([ "chrom", "gene","str.id", "str.start", "n.miss", "beta", "beta.se", "lambda.remel", "p.wald"])+"\n")
    if QUAD:
        if LINEAR and ALLELE:
            if No_ANOVA:
                f.write("\t".join(["chrom","gene",  "str.id", "str.start","alpha", "alpha.se", "alpha.pval", "beta", "beta.se", "beta.pval", "linear.beta", "linear.beta.se", "linear.pval"])+"\n")
            else:
                f.write("\t".join(["chrom","gene",  "str.id", "str.start","alpha", "alpha.se", "alpha.pval", "beta", "beta.se", "beta.pval", "linear.beta", "linear.beta.se", "linear.pval","quad_rsq", "lin_rsq", "delta_aic[Lin-Quad]", "delta_bic", "anova_pva"])+"\n")
        else:
            f.write("\t".join(["chrom","gene",  "str.id", "str.start","alpha", "alpha.se", "alpha.pval", "beta", "beta.se", "beta.pval"])+"\n")
    
    # For each gene:
    # Pull out STRs within distance of gene ends
    # For each STR - gene pair, get expr, str for samples with data and do linreg
    PROGRESS("Expression annotation size %s "%str(expr_annot.shape))
    for i in range(expr_annot.shape[0]):
        gene = expr_annot.index.values[i]
        PROGRESS(" Getting data for %s"%gene)
        start = expr_annot["start"].values[i]
        end = expr_annot["stop"].values[i]
        cis_strs = strgt[(strgt["start"] >= (start-DISTFROMGENE)) & (strgt["start"] <= (end+DISTFROMGENE))]
        PROGRESS("%s STRs tested \n"%str(cis_strs.shape[0]))
       
        if PERMUTE_EXPR:
            expr[gene] = random.sample(list(expr[gene].values), expr.shape[0])
        y = pd.DataFrame({"expr":list(expr.loc[:, gene])})
        y.index = str_samples

        for j in range(cis_strs.shape[0]):
            # cis STR data
            locus_str = cis_strs.iloc[[j],:][str_samples].transpose()
            locus_str.index = str_samples
            locus_str.columns = ["STR_%s"%(cis_strs["start"].values[j])]
            test_str=locus_str.columns[0]
            str_start = cis_strs["start"].values[j]
            if ALLELE:
                    locus_str['x1'] = locus_str[test_str].apply(lambda x: x.split(',')[0] )
                    locus_str['x2'] = locus_str[test_str].apply(lambda x: x.split(',')[1] )
                    samples_to_keep = [str_samples[k] for k in range(len(str_samples)) if (str(locus_str['x1'].values[k]) != "NA")and(str(locus_str['x2'].values[k]) != "NA")]
            else:
                samples_to_keep = [str_samples[k] for k in range(len(str_samples)) if str(locus_str.iloc[:,0].values[k]) != "None" ]
                if QUAD:
                    sys.exit(1)      #This should not happen
                    
            locus_str_t = locus_str.loc[samples_to_keep,:]
            locus_str = locus_str_t[test_str].astype(float)
            #Expression
            locus_y = y.loc[samples_to_keep,:]
            Locus_data = locus_str_t.join(locus_y)
            
## Run regression
            if LINEAR_ONLY:
                res_ols,beta, beta_se, p = LinearRegression(locus_str, locus_y["expr"].values, norm=NORM, minsamples=MINSAMPLES, alleles=ALLELE)
                if beta is not None:
                    f.write("\t".join(map(str, [gene, CHROM, "STR_%s"%str_start, str_start, len(str_samples)-locus_str.shape[0], beta, beta_se, -1, p]))+"\n")
            
            if QUAD:
                data = Locus_data[['expr','x1','x2']]
                quad_model, alpha, alpha_se, alpha_pval, beta, beta_se, beta_pval=QuadraticRegression(data,norm=NORM, minsamples=120)
                if LINEAR:
                    data = Locus_data[['x1','x2']]
                    lin_model, slope, err, pval = LinearRegression(data, locus_y["expr"].values, norm=NORM, minsamples=MINSAMPLES, alleles=ALLELE)
                    if No_ANOVA:
                        f.write("\t".join(map(str, [gene, CHROM, test_str, str_start,alpha, alpha_se, alpha_pval, beta, beta_se, beta_pval, slope, err, pval] ) )+"\n")
                        pass
                    else:
                        if quad_model==None or lin_model==None:
                            lin_rsq, quad_rsq, delta_aic,delta_bic, anova_pval= None, None, None, None, None
                        else:
                            lin_rsq, quad_rsq, delta_aic,delta_bic, anova_pval = Modelcompare(lin_model , quad_model)
                        
                        f.write("\t".join(map(str, [gene, CHROM, test_str, str_start,alpha, alpha_se, alpha_pval, beta, beta_se, beta_pval, slope, err, pval,quad_rsq, lin_rsq, delta_aic,delta_bic, anova_pval] ) )+"\n")
            #   #   #
                else:
                    f.write("\t".join(map(str, [gene, CHROM, test_str, str_start, alpha, alpha_se, alpha_pval, beta, (beta_se), beta_pval] ) )+"\n")
        #   #   #          
    f.close()
