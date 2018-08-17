#Cell Quality workflow

library(ggplot2)
library(dplyr)
library(parallel)
library(reshape2)
library(dropestr)
library(Matrix)

args <- commandArgs(TRUE)
inputFile = args[1]
mtChrom = args[2]
output = args[3]

outputDir = unlist(strsplit(output,"matrix"))[1]

dir.create(outputDir)

holder <- readRDS(inputFile)

umis_per_cell <- sort(Matrix::colSums(holder$cm_raw), decreasing=T)
est_cell_num <- EstimateCellsNumber(umis_per_cell)

#Quality scores
scores <- ScorePipelineCells(holder, mit.chromosome.name=mtChrom, predict.all=T, verbose=T)[names(umis_per_cell)]

#dropEst threshold
intersect_cbs <- names(scores[1:est_cell_num$expected])
intersect_cbs <- intersect_cbs[scores[intersect_cbs] > 0.9]
#intersect_cbs <- intersect_cbs[scores[intersect_cbs] > 0.5]

unknown_cell_scores <- scores[(est_cell_num$expected + 1):length(scores)]
rescued_cbs <- names(unknown_cell_scores)[unknown_cell_scores > 0.5]
#rescued_cbs <- names(unknown_cell_scores)[unknown_cell_scores > 0.9]

unknown_cell_scores <- scores[1:est_cell_num$expected]
filtered_cbs <- names(unknown_cell_scores)[unknown_cell_scores < 0.1]

c(Unchanged=length(intersect_cbs), Rescued=length(rescued_cbs), 
  Filtered=length(filtered_cbs))

r_cm_rescued <- holder$cm_raw[, c(names(umis_per_cell)[1:est_cell_num$expected],rescued_cbs)]
r_cm_rescued <- r_cm_rescued[, -which(colnames(r_cm_rescued) %in% filtered_cbs)]
r_cm_rescued <- r_cm_rescued[grep("^[^;]+$", rownames(r_cm_rescued)),]

holder$cm_rescued = r_cm_rescued

writeMM(r_cm_rescued,file = paste0(outputDir,'matrix.mtx'))
write.table(colnames(r_cm_rescued),file = paste0(outputDir,'barcodes.tsv'),sep = "\t",quote=F,row.names = F,col.names = F)
write.table(cbind(rownames(r_cm_rescued),rownames(r_cm_rescued)),file = paste0(outputDir,'genes.tsv'),sep = "\t",quote=F,row.names = F,col.names = F)