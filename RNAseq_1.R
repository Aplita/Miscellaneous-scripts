################################################################################
# RNA-Seq 1             #
#########################
# This script performs basic Differential Expression analysis of test MTB sequences.
# Uses edgeR and includes Volcano Plots.
#
# Written by Ana Paula Vargas.
# Bioinformatics lab - UPCH. Lima, Peru.
# Last updated: 14-ago-18.
###############################################################################


#source("https://bioconductor.org/biocLite.R")
#biocLite("edgeR")
#biocLite("locfit")
#biocLite("Rsubread")
library("Rsubread")
library("edgeR")



bams = list.files(pattern="*.bam", full.names=TRUE)
gtf = ".gtf"

fc     = featureCounts(bams, annot.ext=gtf, isGTFAnnotationFile=TRUE, isPaired=FALSE)
fc.dge = DGEList (counts=fc$counts, genes=fc$annotation)

condiciones = c("pH2PZA-","pH2PZA-","pH2PZA100","pH2PZA100","pH7PZA100","pH7PZA100")
fc.dge$samples$group = as.factor(condiciones)

#cpm(fc.dge) counts per million
#cpm(fc.dge) > 1 For each gene and sample, true or false if cpm is more than 1

#re-analyze, keeping genes that have more than 1 cpm in at least 2 samples
keep = rowSums(cpm(fc.dge) > 1) >= 2
fc.dge = fc.dge[keep, , keep.lib.sizes=FALSE]

fc.dge.norm = calcNormFactors(fc.dge)
fc.dge.disp = estimateCommonDisp(fc.dge.norm)
fc.dge.disp2 = estimateTagwiseDisp(fc.dge.disp)

group = as.factor(condiciones)
design = model.matrix(~group)
fc.dge.disp3 = estimateDisp(fc.dge.norm, design)
fit = glmFit(fc.dge.disp3, design)

#Differencial expression
#All conditions
lrt = glmLRT(fit,coef = c(2,3))

#pH2PZA- vs pH2PZA100
lrt = glmLRT(fit,coef = c(2))

#pH2PZA- vs pH7PZA100
lrt = glmLRT(fit,coef = c(3))

#pH2PZA100 vs pH7PZA100
lrt = glmLRT(fit,contrast = c(0,-1,1))


tt  = topTags(lrt,n=1000,sort.by="PValue",p.value=0.05)
#tt$comparison
#write.csv(tt$table,file="all.csv")
#write.csv(tt$table,file="pH2PZA-vspH7PZA100.csv")
#write.csv(tt$table,file="pH2PZA100vspH7PZA100.csv")




# Volcano Plot
res <- read.table("pH2PZA100vspH7PZA100__purged.csv", header=TRUE)
head(res)
with(res, plot(logFC, -log10(PValue), pch=20, main="Volcano plot", xlim=c(-2.5,2)))
with(subset(res, FDR<.05 ), points(logFC, -log10(PValue), pch=20, col="red"))
with(subset(res, abs(logFC)>1), points(logFC, -log10(PValue), pch=20, col="orange"))
with(subset(res, FDR<.05 & abs(logFC)>1), points(logFC, -log10(PValue), pch=20, col="green"))
library(calibrate)
with(subset(res, FDR<.05 & abs(logFC)>1), textxy(logFC, -log10(PValue), labs=GeneID, cex=.8))



with(res, plot(logFC.grouppH2PZA100, -log10(PValue), pch=20, main="Dispersión de los niveles de expresión por efecto de pirazinamida", xlim=c(-5.5,5)))
with(subset(res, logFC.grouppH2PZA100>2.0 ), points(logFC.grouppH2PZA100, -log10(PValue), pch=20, col="red"))
with(subset(res, logFC.grouppH2PZA100< -1.5 ), points(logFC.grouppH2PZA100, -log10(PValue), pch=20, col="green"))
with(subset(res, logFC.grouppH2PZA100< -1.5 ), points(logFC.grouppH2PZA100, -log10(PValue), pch=20, col="brown"))
with(subset(res, logFC.grouppH2PZA100>2.0), textxy(logFC.grouppH2PZA100, -log10(PValue), labs=GeneID, cex=.5))
with(subset(res, logFC.grouppH2PZA100< -1.5), textxy(logFC.grouppH2PZA100, -log10(PValue), labs=GeneID, cex=.5))
