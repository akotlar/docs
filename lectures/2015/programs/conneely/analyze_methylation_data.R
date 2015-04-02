################################################################################################
###Call CpGassoc library
################################################################################################
library(CpGassoc)

################################################################################################
###Read in and examine methylation data from Illumina 27k array
###This data was downloaded from GEO (accession number GSE23789)
################################################################################################
load('titration_GSE23789.rda')
###What objects are available in the environment?
ls()

###Lets inspect the objects to see what we have
dim(beta)
beta[1:5,]
###beta is a dataframe of beta values for 27578 CpG sites (rows) and 8 samples (columns)
###The 8 samples represent different mixtures of two cell types from the same individual
###The label indicates cell type A (colon cancer) vs. B (B cell)

dim(siga)
siga[1:5,]
dim(sigb)
sigb[1:5,]
###siga and sigb are dataframes containing the unmethylated (A) and methylated (B) signals
###for the same 27578 CpG sites and the same 8 samples.
###They have the same dimension as the beta (in fact, Illumina computes beta = sigb/(siga+sigb+100))

dim(pval)
pval[1:5,]
###pval is a dataframe containing detection p-values for the sample 27578 sites and the same 8 samples
###These p-values assess whether there is significant detection of signal (rather than just noise)
###and are useful for quality control. 

dim(annot)
annot[1:5,]
###annot is a dataframe of basic annotation for all 27578 CpG sites

###Next, verify that samples and CpGs are sorted in the same order in all dataframes
###This is an important step in any analysis involving multiple objects

###Want to make sure the row names match for all 27578 rows
table(row.names(siga)==row.names(beta))
table(row.names(sigb)==row.names(beta))
table(row.names(pval)==row.names(beta))
table(annot[,1]==row.names(beta))
###and the column names match for all 8 samples
table(names(siga)==names(beta))
table(names(sigb)==names(beta))
table(names(pval)==names(beta))


################################################################################################
###Now we are ready to perform some basic QC of the data
################################################################################################

###We will use the cpg.qc function from the package CpGassoc
###Type help(cpg.qc) to explore the function

meth_norm=cpg.qc(beta,siga,sigb,pval,p.cutoff=0.001,cpg.miss=0.15,sample.miss=0.05,constant100=FALSE,sig.return=TRUE)
###This command will perform several steps at once
###   p.cutoff=0.001 -> set to NA all data points with detection p-values > .001 
###   Remove samples with low total signal (less than half of experiment-wide median or 2000 units)  
###   cpg.miss=0.15 -> remove CpG sites with > 15% of data points missing (NA)
###        Note: normally we would set this to 10% or 5%, but that is too restrictive with only 8 samples
###   sample.miss=0.05 -> remove samples with > 5% of data points missing (NA)
###   constant100=FALSE -> recompute beta values as sigb/(siga+sigb)
###   sig.return=TRUE -> returns cleaned signal A & B data as well as cleaned beta values (will be useful for next step)

###If we had used sig.return=FALSE, meth_norm would be a matrix of beta values 
###Since we used sig.return=TRUE, it is a list containing 3 matrices (beta values, signal A, and signal B)

###We can retrieve and inspect these matrices
beta_qc = meth_norm[[1]]
siga_qc = meth_norm[[2]]
sigb_qc = meth_norm[[3]]
dim(beta_qc)
beta_qc[1:5,]
dim(siga_qc)
siga_qc[1:5,]
dim(sigb_qc)
sigb_qc[1:5,]
###The QCed data matrices have fewer rows than the originals to reflect the 302 CpG sites that were removed

################################################################################################
###Using QCed data, we can quantile normalize the signal data
################################################################################################

###We will use the normalizeQuantiles function from the package limma
library(limma)

#Stack the A and B signals into a single matrix
stacked_signal = rbind(siga_qc,sigb_qc)
n=dim(siga_qc)[1]

#quantile normalize the data
after_qn=normalizeQuantiles(stacked_signal,ties=TRUE)

#Retrieve normalized signals
siga_qn=after_qn[1:n,]
sigb_qn=after_qn[(n+1):(2*n),]

#beta=M/U+M or B/A+B
beta_qn=sigb_qn/(siga_qn+sigb_qn)

#We can compare beta_qc and beta_qn
beta_qc[1:5,]
beta_qn[1:5,]

#Plot all of the beta values for a single sample for the unnormalized vs normalized data
plot(beta_qc[,1],beta_qn[,1],cex=.2)

#We can do this for the other 7 samples as well
par(ask=TRUE)
for (i in 2:8) {
   plot(beta_qc[,i],beta_qn[,i],cex=.2)
}

#Inspect the pairwise relationship between beta values 
#(use just the first 1000 CpG sites to keep this manageable)
pairs(beta_qn[1:1000,],cex=.2)


################################################################################################
### We next use the function cpg.assoc to run a differential methylation analysis
################################################################################################

###Type help(cpg.assoc) to explore the function

###Create cell type proportion variable for analysis, based on the known titration proportions
cellprop <- c(1, 0, .9, 1, 0, .9, .75, .5)
print(cellprop)

### Use CpGassoc to perform regressions of beta values on cellprop, for each CpG site
ewas <- cpg.assoc(beta_qn,cellprop)

### View summary of ewas
print(ewas)

### Useful to know: Holm = step-down Bonferroni
### Useful to know: BH = Benjamini-Hochberg FDR adjustment

### How many CpG sites were FDR-significant?  Holm-significant?

################################################################################################
### We can visualize our results with plots
################################################################################################

###Manhattan plot
###Type help(manhattan) to see syntax and options

###Check out annot object to see what information is available
names(annot)
annot[1:5,]

manhattan(ewas,annot$Name,annot$Chr,annot$MapInfo)
### Here we use the annotation information to get chromosome and position of each site
### Note that the manhattan function will merge the annotation information for you
### so there is no need to check that the sort order is the same
### (In fact, the annot file includes some CpG sites that did not pass QC and thus were not analyzed)

### Plot shows that CpG sites across the genome associate with cell type 
### (solid line = Holm-signifiance, dotted line = FDR-significance)

###QQ plot 
plot(ewas)
### This shows the same excess of significant results we saw in the Manhattan plot

################################################################################################
### Now we'll try extracting information from the ewas for use in some downstream analyses
################################################################################################

### What information is available?  Lets investigate:
ewas
names(ewas)

dim(ewas$results)

ewas$results[1:10,]

table(ewas$results[,1]==row.names(beta_qn))
### ewas$results provides a matrix of results with CpGs sorted in the original order

### Lets use this information to annotate our earlier scatterplot:
### Color all Holm-significant CpG sites red
plot(beta_qn[,c("A1","B1")],col=1+1*ewas$results$Holm.sig)
pairs(beta_qn[1:5000,],cex=.2,col=1+1*ewas$results$Holm.sig)

### Color all FDR-significant CpG sites red
plot(beta_qn[,c("A1","B1")],col=1+1*(ewas$results$FDR<.05))
pairs(beta_qn[1:5000,],cex=.2,col=1+1*(ewas$results$FDR<.05))

################################################################################################
### Enrichment analysis for biological features
################################################################################################

### We talked briefly about CpG islands in class.  Are the FDR-significant sites more likely to be on CpG islands?

annot[1:5,]
### annot has a field (logical TRUE/FALSE) indicating CpG island status 
### However, we will need to subset it to include just the CpG sites that were in the ewas (those that survived QC)
keep = (annot$Name%in%row.names(beta_qn))
annot_qn=annot[keep,]
table(annot_qn$Name==row.names(beta_qn))

table(ewas$results$FDR<.05)
table(annot_qn$CPG_ISLAND)

table(ewas$results$FDR<.05,annot_qn$CPG_ISLAND)

### It appears that our significant sites are less likely to be on CpG islands
### Is this more extreme than we would expect under random chance?  
### We can check with a simple enrichment analysis:

chisq.test(table(ewas$results$FDR<.05,annot_qn$CPG_ISLAND))

### or better yet, use Fishers exact test:

fisher.test(table(ewas$results$FDR<.05,annot_qn$CPG_ISLAND))


