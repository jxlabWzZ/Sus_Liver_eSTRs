import argparse
import math
import numpy as np
import os
import pandas as pd
import random
import shutil
import sys
import gzip


"""
Analyze heritability of gene expression due to SNPs/STRs
Notes:
** The reported SE for STRs when treating the STR as fixed effect are the SE on *beta*, not on *beta^2* **
"""

GCTA_MIN_VE = 0.000001 # Lowest VE that GCTA reports when constrained

def PROGRESS(msg, printit=True):
    if printit: # false for some messages when not in debug mode
        sys.stderr.write("%s\n"%msg.strip())

def GetMAF(x):
    """
    Get SNP MAF
    """
    x=x.convert_objects(convert_numeric=True)
    vals = x.values
#    print(len(vals))
    maf = sum(vals)*1.0/(2*len(vals))
    return min([maf, 1-maf])

def GRM(snpdata):
    """
    Return n by n GRM, scaled to have mean diagonal element 1
    """
#    snpdata.apply(pd.to_numeric, errors='coerce')
    p, n = snpdata.shape
    print(p,'  ',n)
    K = np.zeros((n, n))
    for i in range(p):
#        print(snpdata.iloc[i,:])
        gt = snpdata.iloc[i,:].apply(lambda x: float(x))
        var = np.var(gt)
        #print(var)
        if var == 0: continue
        m = np.mean(gt)
        x = np.reshape((np.array(gt)-m).transpose(), (n, 1))
        gt_m = 1/(var*p)*x.dot(x.transpose())
        K = K + gt_m
#    print(K )
    # Make sure diagonal had mean 1 (should anyway)
    diag_mean = np.mean(np.diagonal(K))
#    print(np.diagonal(K))
    K = K/diag_mean
    return K

def WriteGCTAGRM(K, grmfile, p):
    """
    Calculate GRM and output file in GCTA format (.gz)
    Need to write:
    $grmfile.grm.gz: ind1, ind2, num nonmissing SNPs, relatedness (space-separated) (indices start at 1, rows into $grmfile.ind)
      only includes lower triangle of GRM
    $grmfile.ind: family ID, individual ID (space-separated)
    """
    n = K.shape[0]
    # Calculate GRM
    f = gzip.open("%s.grm.gz"%grmfile, "wb")
    for i in range(n):
        for j in range(i, n):
            val = K[i, j]
            f.write(" ".join(map(str, [j+1, i+1, p, val]))+"\n")
#            f.write(bytes(" ".join([str(j+1), str(i+1), str(p), str(val)]))+"\n" , 'UTF-8'))
    f.close()
    # Write ind file
    f = open("%s.grm.id"%grmfile, "w")
    for i in range(n):
        f.write(" ".join(map(str, [i, i]))+"\n")
    f.close()

def ParseGCTAResults(gctafile, include_str):
    """
    return cis_snp_h2, cis_snp_h2_se, cis_str_h2, cis_str_h2_se, logL
    """
    f = open(gctafile, "r")
    lines = f.readlines()
    cis_str_h2 = None
    cis_str_h2_se = None
    for line in lines:
        items = line.strip().split()
        if len(items) < 1: continue
        if include_str == "NO" or include_str == "FE" or include_str == "SAMPLES":
            if items[0] == "V(G)/Vp":
                cis_snp_h2 = items[1]
                cis_snp_h2_se = items[2]
        else:
            if items[0] == "V(G1)/Vp":
                cis_snp_h2 = items[1]
                cis_snp_h2_se = items[2]
            if items[0] == "V(G2)/Vp":
                cis_str_h2 = items[1]
                cis_str_h2_se = items[2]
        if items[0] == "logL":
            logL = items[1]
        if include_str == "FE":
            if items[0] == "Fix_eff":
                ind = lines.index(line)+2
                items = lines[ind].strip().split()
                cis_str_h2 = float(items[0])**2
                cis_str_h2_se = items[1]
        if items[0] == "Pval":
            Pval = items[1]
    if include_str == "FE":
        return cis_snp_h2, cis_snp_h2_se, cis_str_h2, cis_str_h2_se, logL, Pval
    else:
        return cis_snp_h2, cis_snp_h2_se, cis_str_h2, cis_str_h2_se, logL, 'N/A'

def GetPermutedLocusSTRs(locus_str):
    """
    Return locus_str dataframe with STR genotypes permuted
    """
    gts = list(locus_str.iloc[:,0])
    random.shuffle(gts)
    locus_str_perm = pd.DataFrame({locus_str.columns[0]: gts})
    locus_str_perm.index = locus_str.index
    return locus_str_perm

def z(vals):
    vals = list(map(float, list(vals)))
    m = np.mean(vals)
    s = math.sqrt(np.var(vals))
    return [(item-m)*1.0/s for item in vals]

def ZNorm(locus_vars):
    """
    Znormalize variants
    """
    columns = locus_vars.columns
    for c in columns:
        locus_vars[c] = z(locus_vars[c])
    return locus_vars


def WriteGCTAPhenotypeFile(locus, exprfile):
    """
    Write the phenotype file 
    """
    n=len(locus)
    f = open(exprfile, "w")
    for i in range(n):
        f.write(' '.join([str(i),str(i),str(locus[i]),'\n']))
    f.close()
    
def WriteGCTACovarFile(locus, strcovarfile):
    """ Write GCTA covariable file using normalized genotype
    """
    f = open(strcovarfile, "w")
    n=locus.shape[0]
    for i in range(n):
        N_geno = list(locus.iloc[i,:].values)
        f.write(" ".join([str(i), str(i)]+[str(m) for m in N_geno])+"\n")
    f.close()
    
DISTFROMGENE = 100000    
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="heritability of gene expression due to SNPs vs STRs")
    parser.add_argument("--expr", help="Normalized expression residuals", type=str, required=True)
    parser.add_argument("--samples", help="Restric to samples in this file", type=str, required=False)
    parser.add_argument("--exprannot", help="Expression annotation file", type=str, required=True)
    parser.add_argument("--chrom", help="Restrict analysis to this chromosome", type=str, required=True)
    parser.add_argument("--distfromgene", help="Look at STRs/SNPs within this distance of gene boundaries", type=int, required=True)
    parser.add_argument("--strgt", help="File with noramlized STR genotypes", type=str, required=True)
    parser.add_argument("--snpgt", help="File with normalized SNP genotypes", type=str, required=True)
    parser.add_argument("--esnps", help="File with regression of SNPs with LR pvalues", type=str, required=False)
    parser.add_argument("--out", help="Write results to this file", type=str, required=True)
    parser.add_argument("--tmpdir", help="Use this directory instead of /tmp for temporary files", type=str, default="/tmp")
    parser.add_argument("--snp-maf", help="Only use SNPs with this MAF in our samples", type=float, default=0.01)
    parser.add_argument("--method", help="Which program to run. GEMMA or GCTA", type=str, required=True)
    parser.add_argument("--include-str", help="Include STR in analysis? Should be either NO, FE (fixed effect), or RE (random effect), or SAMPLES (use the samples that have a call at the STR, but don't use the STR)", type=str, required=True)
    parser.add_argument("--estr-results", help="eSTR LMM results (used to choose which STR to include", type=str, required=False)
    parser.add_argument("--reml-no-constrain", help="Don't constrain var comp. ests to be between 0 and 1 (--reml-no-constrain for GCTA)", action="store_true")
    parser.add_argument("--switch-order", help="Switch order of variance components (STR first, SNP second)", action="store_true")
    parser.add_argument("--estr-genes-only", help="Only run on eSTR genes (below specified FDR)", type=float, default=1, required=False)
    parser.add_argument("--permute-strs", help="Generate empirical p-value by permuting STR genotypes this many times", type=int, default=0, required=False)
    parser.add_argument("--ctrl", help="This is an unlinked control analysis. Use random STR instead of best", action="store_true")
    parser.add_argument("--debug", help="Print debug status messages", action="store_true")
    args = parser.parse_args()
    EXPRFILE = args.expr
    EXPRANNOTFILE = args.exprannot
    CHROM = args.chrom
    if "chr" not in str(CHROM): CHROM="chr%s"%CHROM
    DISTFROMGENE = args.distfromgene
    STRGTFILE = args.strgt
    SNPGTFILE = args.snpgt
    OUTFILE = args.out
    TMPDIR = args.tmpdir
    SNPMAF = args.snp_maf
    LMM_METHOD = args.method
    INCLUDE_STR = args.include_str
    ESTR_RESULTS_FILE = args.estr_results
    REML_NO_CONSTRAIN = args.reml_no_constrain
    ESTR_GENES_ONLY = args.estr_genes_only
    SAMPLES = args.samples
    PERMUTE_STRS = args.permute_strs
    UNLINKED_CTRL = args.ctrl
    DEBUG = args.debug
    if INCLUDE_STR == "RE" and LMM_METHOD == "GEMMA":
        PROGRESS("ERROR: Cannot run GEMMA with STR as random effect")
        sys.exit(1)
    if args.esnps:
        ESNPS = pd.read_csv(args.esnps, sep='\t')
    if LMM_METHOD not in ["GCTA", "GEMMA"]:
        PROGRESS("ERROR: Command must be one of GCTA, GEMMA")
        sys.exit(1)
    if LMM_METHOD == "GEMMA":
        PROGRESS("GEMMA not implemented yet")
        sys.exit(1)
    if INCLUDE_STR not in ["NO", "FE", "RE","SAMPLES"]:
        PROGRESS("ERROR: Command must be one of NO, FE, RE","SAMPLES")
        sys.exit(1)
    if INCLUDE_STR != "NO" and ESTR_RESULTS_FILE is None:
        PROGRESS("ERROR: Must supply --estr-results file if --include-str is not NO")
        sys.exit(1)
    if INCLUDE_STR == "RE" and REML_NO_CONSTRAIN:
        PROGRESS("ERROR: can't use --include-str RE and --reml-no-constrain. GCTA won't converge")
        sys.exit(1)
    if PERMUTE_STRS > 0 and ESTR_RESULTS_FILE is None:
        PROGRESS("ERROR: Must supply --estr-results file if --permute-strs is > 0")
        sys.exit(1)
    if ESTR_GENES_ONLY < 1 and ESTR_RESULTS_FILE is None:
        PROGRESS("ERROR: Must supply --estr-results-file if --estr-genes-only is < 1")
        sys.exit(1)
    if PERMUTE_STRS > 0 and INCLUDE_STR != "RE":
        PROGRESS("ERROR: Permutation analysis only for STR as random effect")
        sys.exit(1)
    if SAMPLES is not None:
        samples = pd.read_csv(SAMPLES, sep=' ')
        sample = [x.strip() for x in list(samples['x'])]
        samples = ['-'.join(s.split('.')).strip('\n') for s in sample]   #Using Ids.txt file created from expression cleaning
        
    # Load expression and annotation
    PROGRESS("\n\nLoad expression", printit=DEBUG)
    expr = pd.read_csv(EXPRFILE)
    PROGRESS("Load annotation", printit=DEBUG)
    expr_annot = pd.read_csv(EXPRANNOTFILE)
    expr_annot.index = expr_annot["gene"].values
    #expr_annot = expr_annot.loc[expr.columns[map(lambda x: x in expr_annot.index, expr.columns)],:]
    expr_annot = expr_annot.loc[[item for item in expr.columns if item in expr_annot.index],:]
    expr_annot = expr_annot[expr_annot["chr"] == CHROM]
 
    # Load SNP genotypes
    PROGRESS("Load SNPs", printit=DEBUG)
    snpgt = pd.read_csv(SNPGTFILE, sep="\t")
    snpgt = snpgt.loc[snpgt['#chrom']==CHROM]
    esnps = ESNPS.loc[ESNPS['chrom']==CHROM]
# Load STR genotypes
    PROGRESS("Load STRs", printit=DEBUG)
    strgt = pd.read_csv(STRGTFILE, sep="\t",low_memory=False)
    strgt = strgt.loc[strgt['#chrom']==CHROM].reindex()
    print "**",strgt.shape
# Restrict to samples
    if SAMPLES is None:
        str_samples = list(set(strgt.columns[2:].values).intersection(set(snpgt.columns[2:].values)))    
    str_samples = [ s for s in samples if s in expr.index and s in list(strgt.columns)] #keep only genotypes+expression data
    #str_samples = [ s for s in str_samples if s in in list(strgt.columns)]
####
    print len(str_samples)
    expr = expr.loc[str_samples,:]
    snpgt = snpgt[["#chrom","start"] + str_samples]
    #print(strgt['GTEX-ZZPU'])
    strgt = strgt[["#chrom","start"] + str_samples]
    print "**Samples*********************** ",len(str_samples), strgt.shape
# Load STR results
    if ESTR_RESULTS_FILE is not None:
        estr_results = pd.read_csv(ESTR_RESULTS_FILE, sep="\t")

# Set up output file
    outf = open(OUTFILE, "w")
    if INCLUDE_STR == "NO" or INCLUDE_STR == "SAMPLES":
        outf.write("\t".join(["chrom","gene","num_snps","cis_snp_h2","cis_snp_h2_se","logL","nsamp", "pval"])+"\n")
    else:
        outf.write("\t".join(["chrom","gene","str_start","num_snps","cis_snp_h2","cis_snp_h2_se",\
                                  "cis_str_h2","cis_str_h2_se","logL","nsamp","cis_str_h2_pval"])+"\n") #nperm for cis strs

# Restrict to eSTR genes
    if ESTR_GENES_ONLY < 1:
        genelist = set(estr_results[estr_results["qvalue"]<=ESTR_GENES_ONLY].gene)
        genelist = [item for item in genelist if item in expr_annot.index]
        expr_annot = expr_annot.loc[genelist]
        if ESTR_RESULTS_FILE is not None:
            estr_results = estr_results.loc[estr_results['gene'].isin(genelist)]
            expr = expr[genelist]
            
# For each gene, pull out data and perform specified method
    for i in range(expr_annot.shape[0]):
        PROGRESS("\t\t Starting the loop the %s th gene"% str(i))
        #gene = 'ENSSSCG00000004008'#expr_annot.index.values[i] 
        #ensgene = 'ENSSSCG00000004008' #expr_annot["gene.id"].values[i]
        gene = expr_annot.index.values[i]
        ensgene = expr_annot["gene"].values[i]
        PROGRESS("Getting data for %s"%gene, printit=DEBUG)
        genedir = os.path.join(TMPDIR,"%s/"%gene)
        if not os.path.exists(genedir):
            os.mkdir(genedir)
        start = expr_annot["start"].values[i]
        end = expr_annot["stop"].values[i]
        y = expr.loc[:,gene]       ; print gene                      ##*** on print
        
# Pull out STRs
        samples_to_keep = list(y.dropna(axis=0, how='any').index) ## str_samples
        best_str_start = None
        if INCLUDE_STR != "NO":
            try:
                if UNLINKED_CTRL:
                    print gene, ' .Unlinked' #don't want this now ***
                    possible_starts = list(strgt[(strgt["start"] >= (start-DISTFROMGENE)) & (strgt["start"] <= (end+DISTFROMGENE))].start)
                    best_str_start = random.sample(possible_starts, 1)[0]
                else: 
#make sure to match on Ensembl gene (gene is ILMN if using array)
                    best_str_start = estr_results[estr_results["gene"]==ensgene].sort_values("p.wald")["str.start"].values[0]
            except:
                print sys.exc_info()
                PROGRESS("[%s]\tERROR: couldn't find STR LMM results"%gene)
                continue #break
            try:
                locus_str = strgt[(strgt["start"] == best_str_start)].iloc[[0],:][str_samples].transpose()
                PROGRESS("Transpose locus str %s"% str(locus_str.shape))
            except:
                PROGRESS("[%s]\tERROR: couldn't find STR genotypes for position %s"%(gene, best_str_start))
                continue
            locus_str.index = str_samples
            locus_str.columns = ["STR_%s"%best_str_start]
            samples_to_keep = [str_samples[k] for k in range(len(str_samples)) if str(locus_str.iloc[:,0].values[k]) != "None"]
            locus_str = locus_str.loc[samples_to_keep,:]
            print locus_str.columns  #****
            # Make sure STRs are normalized
            try:
                locus_str = ZNorm(locus_str)
            except:
                PROGRESS("[%s]\tERROR: couldn't Z normalize STR genotypes"%(gene))
                continue
        # Pull out SNPs
        PROGRESS("Getting SNPs data")
        
#        s = list(esnps.loc[esnps['gene']==gene]['str.start'])[0] ######################
#        cis_snps = snpgt.loc[snpgt["start"] == s] #####################
        cis_snps = snpgt[(snpgt["start"] >= (start-DISTFROMGENE)) & (snpgt["start"] <= (end+DISTFROMGENE))]
        
        locus_snp = cis_snps[samples_to_keep].transpose()
        locus_snp.index = samples_to_keep
        locus_snp.columns = cis_snps["start"].apply(lambda x: "SNP_%s"%x)
        locus_snp_maf = locus_snp.apply(lambda x: GetMAF(x), 0)
        print locus_snp.columns, ' cis SNPs', locus_snp.shape
        if locus_snp.shape[1] == 0:
            PROGRESS("[%s]\tERROR: no common SNPs in region\t %s"%(gene, locus_snp.shape))
            continue
        # Get expression
        y = pd.DataFrame({"expr":list(expr.loc[:,gene])})
        y.index = str_samples
        locus_y = y.loc[samples_to_keep,["expr"]]
        # Z normalize
        locus_y = (locus_y - np.mean(locus_y))/math.sqrt(np.var(locus_y))
        # Make SNP GRM
        locus_snp=locus_snp.apply(pd.to_numeric, errors='coerce')
        locus_snp=locus_snp.dropna(axis=1, how='any')
        K = GRM(locus_snp.transpose())
        if str(np.mean(K)) == "nan":
            PROGRESS("[%s]\tERROR: nans in GRM"%gene)
            sys.exit(1) # This shouldn't happen
        # Write GRM
        if LMM_METHOD == "GCTA":
            exprfile = os.path.join(genedir, "expr.pheno")
            mgrmfile = os.path.join(genedir, "mgrm.txt")
            snpgrmfile = os.path.join(genedir, "snp.grm.txt")
            if REML_NO_CONSTRAIN:
                reml_command = "--reml-no-constrain"
            else: reml_command = "--reml"
            gcta_cmd = "/home/yangbin/bin/gcta_1.92.1beta6/gcta64 %s --mgrm-gz %s --pheno %s --out %s/gcta --thread-num 2 "%(reml_command, mgrmfile, exprfile, genedir)
            g = open(mgrmfile, "w")
            g.write(snpgrmfile+"\n")
            WriteGCTAGRM(K, snpgrmfile, p=locus_snp.shape[1])
            if INCLUDE_STR == "FE": # --qcovar
                strcovarfile = os.path.join(genedir, "str.qcovar")
                WriteGCTACovarFile(locus_str, strcovarfile)
                gcta_cmd += " --qcovar %s --reml-est-fix"%strcovarfile
            if INCLUDE_STR == "RE": 
                K_str = GRM(locus_str.transpose())
                strgrmfile = os.path.join(genedir, "str.grm.txt")
                g.write(strgrmfile+"\n")
                WriteGCTAGRM(K_str, strgrmfile, p=locus_str.shape[1])
            g.close()
            locus_y["expr"].fillna("NA", inplace=True) #locus_y.dropna(axis=0, how='any')#
            WriteGCTAPhenotypeFile(locus_y["expr"].values, exprfile)
            print gcta_cmd
#           gcta_cmd += " > /dev/null 2>&1"              ###get output  to identify error
            os.system(gcta_cmd)
            # Parse results
            gctafile = os.path.join(genedir, "gcta.hsq")
            if not os.path.exists(gctafile):
                PROGRESS("[%s]\tERROR: GCTA could not analyze this gene"%gene)
                continue
            cis_snp_h2, cis_snp_h2_se, cis_str_h2, cis_str_h2_se, logL, pval = ParseGCTAResults(gctafile, INCLUDE_STR)

        # Output results
        if INCLUDE_STR == "NO" or INCLUDE_STR == "SAMPLES":
            outf.write("\t".join(map(str, [CHROM, gene, locus_snp.shape[1], cis_snp_h2, cis_snp_h2_se, logL, len(samples_to_keep)]),pval)+"\n")
        else:
            outf.write("\t".join(map(str, [CHROM, gene, best_str_start, locus_snp.shape[1], cis_snp_h2, cis_snp_h2_se,\
                                               cis_str_h2, cis_str_h2_se, logL, len(samples_to_keep), pval]))+"\n")
            ##used if we consider all cis_strs to genes as FE --> len(cis_str_h2_null)
       
        outf.flush()
    outf.close()
                                                                                  
