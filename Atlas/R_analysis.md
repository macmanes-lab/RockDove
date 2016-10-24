#### Rstats

```R
source("http://bioconductor.org/biocLite.R")
biocLite(c("Rgraphviz", "topGO", "edgeR", "tximport", "readr", "locfit", "statmod", "gplots", "AnnotationDbi", "impute", "GO.db", "preprocessCore", "KEGG.db", "topGO", "packageNames", "hgu133a.db", "hgu95av2.db", "annotate", "hgu133plus2.db", "SNPlocs.Hsapiens.dbSNP.20100427", "minet", "OrderedList", "mygene"))
install.packages(c("matrixStats", "Hmisc", "foreach", "doParallel", "fastcluster", "dynamicTreeCut", "survival", "WGCNA", "mygene"))
orgCodes = c("Hs", "Mm", "Rn", "Pf", "Sc", "Dm", "Bt", "Ce", "Cf", "Dr", "Gg");
orgExtensions = c(rep(".eg", 4), ".sgd", rep(".eg", 6));
packageNames = paste("org.", orgCodes, orgExtensions, ".db", sep="");
library(edgeR)
library(tximport)
library(readr)
library(locfit)
library(statmod)
library(gplots)
library(topGO)
library(Rgraphviz)

dir <- setwd('/Users/macmanes/Dropbox/RockDove_working')
tx2gene <- read.table(file.path(dir, "rockdove.chicken.tximport"), header=T)

all_samples <- read.table(file.path(dir, "mapping/atlas.list"), header=F)
all_files <- file.path(dir, "mapping/files/atlas", all_samples$V1, "quant.sf")
names(all_files) <- paste0(all_samples$V1)
txi.all <- tximport(all_files, type = "salmon", tx2gene = tx2gene, reader = read_tsv)

cts <- txi.all$counts
normMat <- txi.all$length
normMat <- normMat/exp(rowMeans(log(normMat)))
o <- log(calcNormFactors(cts/normMat)) + log(colSums(cts/normMat))

atlasphenotype <- factor(c(
'mh',
'mp',
'mt',
'mh',
'mp',
'mt',
'fh',
'fo',
'fp',
'mh',
'mp',
'mt',
'mh',
'mp',
'mt',
'mh',
'mp',
'mt',
'mh',
'mp',
'mt',
'mh',
'mp',
'mt',
'fh',
'fo',
'fp',
'fh',
'fo',
'fp',
'fh',
'fo',
'fp',
'fh',
'fo',
'fp',
'mh',
'mp',
'mt',
'fh',
'fo',
'fp',
'fo',
'fh',
'fo',
'fp',
'fh',
'fo',
'fp',
'fh',
'fp',
'mh',
'mp',
'mt',
'mp',
'fo',
'mp',
'mt',
'mh',
'mp',
'mt',
'mh',
'mp',
'mt',
'mp',
'mt'
))

atlasdesign <- model.matrix( ~0+atlasphenotype )

for_analysis <- read.csv("~/Dropbox/RockDove_working/for_analysis_chicken_atlas_entrez.csv")

atlasobject <- DGEList(counts=for_analysis[,5:70], genes=for_analysis[,1:4], group=atlasphenotype)

atlasobject$offset <- t(t(log(normMat)) + o)
atlasobject <- calcNormFactors(atlasobject)
atlasobject <- estimateCommonDisp(atlasobject)
atlasobject <- estimateTagwiseDisp(atlasobject)
atlasobject <- estimateDisp(atlasobject, atlasdesign)
atlasobject <- estimateGLMCommonDisp(atlasobject, atlasdesign, verbose=TRUE)
atlasobject <- estimateGLMTrendedDisp(atlasobject, atlasdesign)
atlasobject <- estimateGLMTagwiseDisp(atlasobject, atlasdesign)

fit <- glmQLFit( atlasobject, atlasdesign, robust=T, lfc=1 )
```

#### contrast_fp_v_mp
```R
contrast_fp_v_mp <- glmQLFTest( fit,
                      contrast=makeContrasts(
                      atlasphenotypefp - atlasphenotypemp,
                      levels=atlasdesign ) )

summary(decideTestsDGE(contrast_fp_v_mp, adjust.method="fdr", p.value=0.01, lfc=1))

detags_fp_v_mp <- rownames(contrast_fp_v_mp)[as.logical(decideTestsDGE(contrast_fp_v_mp, adjust.method="fdr", p.value=0.01, lfc=1))]
size = summary(decideTestsDGE(contrast_fp_v_mp, adjust.method="fdr", p.value=0.01, lfc=1))[3,1] + summary(decideTestsDGE(contrast_fp_v_mp, adjust.method="fdr", p.value=0.01, lfc=1))[1,1]
output <- topTags(contrast_fp_v_mp, n=size)
write.table(output, file="diffexp_detags_fp_v_mp.csv", sep = "," , row.names = TRUE)

png("fp_v_mp.png", width = 8, height = 6, units = 'in', res = 300)
par(bty = 'n')
plotSmear(atlasobject, de.tags=detags_fp_v_mp, , xlim=c(0,20), smooth.scatter=T, pair = c("fp","mp"))
abline( h=c(-2,2) , col="dodgerblue" )
dev.off()
plotSmear(atlasobject, de.tags=detags_fp_v_mp, , xlim=c(0,20), smooth.scatter=T, pair = c("fp","mp"))
abline( h=c(-2,2) , col="dodgerblue" )
```
#### contrast_fh_v_mh
```R
contrast_fh_v_mh <- glmQLFTest( fit,
                      contrast=makeContrasts(
                      atlasphenotypefh - atlasphenotypemh,
                      levels=atlasdesign ) )

summary(decideTestsDGE(contrast_fh_v_mh, adjust.method="fdr", p.value=0.01, lfc=1))

detags_fh_v_mh <- rownames(contrast_fh_v_mh)[as.logical(decideTestsDGE(contrast_fh_v_mh, adjust.method="fdr", p.value=0.01, lfc=1))]
size = summary(decideTestsDGE(contrast_fh_v_mh, adjust.method="fdr", p.value=0.01, lfc=1))[3,1] + summary(decideTestsDGE(contrast_fh_v_mh, adjust.method="fdr", p.value=0.01, lfc=1))[1,1]
output <- topTags(contrast_fh_v_mh, n=size)
write.table(output, file="diffexp_detags_fh_v_mh.csv", sep = "," , row.names = TRUE)

png("fh_v_mh.png", width = 8, height = 6, units = 'in', res = 300)
par(bty = 'n')
plotSmear(atlasobject, de.tags=detags_fh_v_mh, , xlim=c(0,20), smooth.scatter=T, pair = c("fp","mp"))
abline( h=c(-2,2) , col="dodgerblue" )
dev.off()
plotSmear(atlasobject, de.tags=detags_fh_v_mh, , xlim=c(0,20), smooth.scatter=T, pair = c("fp","mp"))
abline( h=c(-2,2) , col="dodgerblue" )
```
