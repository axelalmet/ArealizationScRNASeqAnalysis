library(ggplot2)
library(DESeq2)
library(dplyr)
library(edgeR)
library(SingleCellExperiment)
library(zellkonverter)
library(plyr)
library(EnhancedVolcano)
library(scuttle)

data_directory <- "../scrnaseqdata/" # Change to whever the data is

# sce <- readH5AD(paste0(data_directory, "WatanabeLab_ArealizationOrganoids/cortical_organoid_merged.h5ad"))

### Pseudobulking step
# pbulk_sce <- aggregateAcrossCells(sce, id=colData(sce)[,c("orig_sample", "condition", "samples", "day", "celltype", "condition_celltypes")])
# colData(pbulk_sce) <- colData(pbulk_sce)[, c(1:24, 31)]
# writeH5AD(pbulk_sce, paste0(data_directory, "WatanabeLab_ArealizationOrganoids/cortical_organoid_pbulk.h5ad"), compression="gzip")

# Load the counts
pbulk_sce <- readH5AD(paste0(data_directory, "WatanabeLab_ArealizationOrganoids/cortical_organoid_pbulk.h5ad"))

pbulk_counts <- assay(pbulk_sce, "counts")
meta <- data.frame(colData(pbulk_sce))

# Relevel so that CTR and D35 are the control
meta$day <- relevel(meta$day, ref="D35")
meta$condition <- relevel(meta$condition, ref="ctr")

### edgeR analysis for a single cell type
celltype_oi <- "aRG/oRG" # Specify cell type label
celltype_oi_label <- gsub("/", "-", gsub(" ", "-", celltype_oi))


counts_celltype <- pbulk_counts[, meta$celltype %in% celltype_oi]
meta_celltype <- meta[meta$celltype %in% celltype_oi,]

meta_celltype$samples <- droplevels(meta_celltype$samples)
meta_celltype$day <- droplevels(meta_celltype$day)
meta_celltype$celltype <- droplevels(meta_celltype$celltype)
meta_celltype$condition_celltypes <- droplevels(meta_celltype$condition_celltypes)

design_celltype <- model.matrix(~0 + samples, meta_celltype)
celltype_edger <- DGEList(counts=counts_celltype, genes=rownames(counts_celltype), design=design_celltype)

keep_celltype <- filterByExpr(celltype_edger, design_celltype)
celltype_edger <- celltype_edger[keep_celltype, , keep.lib.sizes=FALSE]

# Normalise
celltype_edger <- calcNormFactors(celltype_edger)

# Estimate dispersion
celltype_edger <- estimateDisp(celltype_edger, design_celltype, robust=TRUE)

# Perform GLM fit
celltype_fit <- glmQLFit(celltype_edger, design_celltype)


#### DLN analysis (only present at D56 and D97)
# D56 analysis
celltype_qlf_d56_ctr_vs_bmp <- glmQLFTest(celltype_fit, contrast=c(1, 0, -1, 0, 0, 0))
celltype_qlf_d56_ctr_vs_fgf <- glmQLFTest(celltype_fit, contrast=c(0, 0, -1, 0, 1, 0))
celltype_qlf_d56_bmp_vs_fgf <- glmQLFTest(celltype_fit, contrast=c(-1, 0, 0, 0, 1, 0))

# D97 analysis
celltype_qlf_d97_ctr_vs_bmp <- glmQLFTest(celltype_fit, contrast=c(0, 1, 0, -1, 0, 0))
celltype_qlf_d97_ctr_vs_fgf <- glmQLFTest(celltype_fit, contrast=c(0, 0, 0, -1, 0, 1))
celltype_qlf_d97_bmp_vs_fgf <- glmQLFTest(celltype_fit, contrast=c(0, -1, 0, 0, 0, 1))

degs_d56_ctr_vs_bmp <- data.frame(topTags(celltype_qlf_d56_ctr_vs_bmp, n=dim(celltype_qlf_d56_ctr_vs_bmp)[1]))
degs_d56_ctr_vs_bmp <- degs_d56_ctr_vs_bmp[order(-degs_d56_ctr_vs_bmp$logFC),]
degs_d56_ctr_vs_fgf <- data.frame(topTags(celltype_qlf_d56_ctr_vs_fgf, n=dim(celltype_qlf_d56_ctr_vs_fgf)[1]))
degs_d56_ctr_vs_fgf <- degs_d56_ctr_vs_fgf[order(-degs_d56_ctr_vs_fgf$logFC),]
degs_d56_bmp_vs_fgf <- data.frame(topTags(celltype_qlf_d56_bmp_vs_fgf, n=dim(celltype_qlf_d56_bmp_vs_fgf)[1]))
degs_d56_bmp_vs_fgf <- degs_d56_bmp_vs_fgf[order(-degs_d56_bmp_vs_fgf$logFC),]

write.csv(degs_d56_ctr_vs_bmp,
          file=paste0(data_directory, "WatanabeLab_ArealizationOrganoids/organoids_pbulk_edgeR_qlf_degs_D56_BMPvsCTR_", celltype_oi ,".csv"))
write.csv(degs_d56_ctr_vs_fgf,
          file=paste0(data_directory, "WatanabeLab_ArealizationOrganoids/organoids_pbulk_edgeR_qlf_degs_D56_FGFvsCTR_", celltype_oi,".csv"))
write.csv(degs_d56_bmp_vs_fgf,
          file=paste0(data_directory, "WatanabeLab_ArealizationOrganoids/organoids_pbulk_edgeR_qlf_degs_D56_FGFvsBMP_", celltype_oi,".csv"))

degs_d97_ctr_vs_bmp <- data.frame(topTags(celltype_qlf_d97_ctr_vs_bmp, n=dim(celltype_qlf_d97_ctr_vs_bmp)[1]))
degs_d97_ctr_vs_bmp <- degs_d97_ctr_vs_bmp[order(-degs_d97_ctr_vs_bmp$logFC),]
degs_d97_ctr_vs_fgf <- data.frame(topTags(celltype_qlf_d97_ctr_vs_fgf, n=dim(celltype_qlf_d97_ctr_vs_fgf)[1]))
degs_d97_ctr_vs_fgf <- degs_d97_ctr_vs_fgf[order(-degs_d97_ctr_vs_fgf$logFC),]
degs_d97_bmp_vs_fgf <- data.frame(topTags(celltype_qlf_d97_bmp_vs_fgf, n=dim(celltype_qlf_d97_bmp_vs_fgf)[1]))
degs_d97_bmp_vs_fgf <- degs_d97_bmp_vs_fgf[order(-degs_d97_bmp_vs_fgf$logFC),]

write.csv(degs_d97_ctr_vs_bmp,
          file=paste0(data_directory, "WatanabeLab_ArealizationOrganoids/organoids_pbulk_edgeR_qlf_degs_D97_BMPvsCTR_", celltype_oi ,".csv"))
write.csv(degs_d97_ctr_vs_fgf,
          file=paste0(data_directory, "WatanabeLab_ArealizationOrganoids/organoids_pbulk_edgeR_qlf_degs_D97_FGFvsCTR_", celltype_oi,".csv"))
write.csv(degs_d97_bmp_vs_fgf,
          file=paste0(data_directory, "WatanabeLab_ArealizationOrganoids/organoids_pbulk_edgeR_qlf_degs_D97_FGFvsBMP_", celltype_oi,".csv"))

keyvals <- ifelse(
  degs_d56_ctr_vs_bmp$logFC < -0.5 & degs_d56_ctr_vs_bmp$FDR < 0.05, '#33A02C',
  ifelse(degs_d56_ctr_vs_bmp$FDR < 0.05 & degs_d56_ctr_vs_bmp$logFC > 0.5, '#1F78B4',
         'black'))
names(keyvals)[keyvals == '#33A02C'] <- 'CTR'
names(keyvals)[keyvals == '#1F78B4'] <- 'BMP'
volcano_ctr_vs_bmp <- EnhancedVolcano(degs_d56_ctr_vs_bmp,
                                      lab = rownames(degs_d56_ctr_vs_bmp),
                                      selectLab = rownames(degs_d56_ctr_vs_bmp)[which(names(keyvals) %in% c('BMP', 'CTR'))],
                                      x = 'logFC',
                                      y = 'PValue',
                                      title = 'BMP vs. CTR (DLN, D56)',
                                      pCutoff = 0.05,
                                      FCcutoff = 0.5,
                                      cutoffLineType="blank",
                                      colCustom = keyvals,
                                      gridlines.major = FALSE,
                                      gridlines.minor = FALSE)
volcano_ctr_vs_bmp
ggsave(file=paste0(data_directory, "WatanabeLab_ArealizationOrganoids/", celltype_oi_label, "/D56/organoids_pbulk_allcells_edgeR_qlf_degs_D56_BMPvsCTR_", celltype_oi_label, "_volcano.pdf"))

keyvals <- ifelse(
  degs_d56_ctr_vs_fgf$logFC < -0.5 & degs_d56_ctr_vs_fgf$FDR < 0.05, '#33A02C',
  ifelse(degs_d56_ctr_vs_fgf$FDR < 0.05 & degs_d56_ctr_vs_fgf$logFC > 0.5, '#E31A1C',
         'black'))
names(keyvals)[keyvals == '#33A02C'] <- 'CTR'
names(keyvals)[keyvals == '#E31A1C'] <- 'FGF'
volcano_ctr_vs_fgf <- EnhancedVolcano(degs_d56_ctr_vs_fgf,
                                      lab = rownames(degs_d56_ctr_vs_fgf),
                                      selectLab = rownames(degs_d56_ctr_vs_fgf)[which(names(keyvals) %in% c('FGF', 'CTR'))],
                                      x = 'logFC',
                                      y = 'PValue',
                                      title = 'FGF vs. CTR (DLN, D56)',
                                      pCutoff = 0.05,
                                      FCcutoff = 0.5,
                                      cutoffLineType="blank",
                                      colCustom = keyvals,
                                      gridlines.major = FALSE,
                                      gridlines.minor = FALSE)
volcano_ctr_vs_fgf
ggsave(file=paste0(data_directory, "WatanabeLab_ArealizationOrganoids/", celltype_oi_label, "/D56/organoids_pbulk_allcells_edgeR_qlf_degs_D56_FGFvsCTR_", celltype_oi_label, "_volcano.pdf"))

keyvals <- ifelse(
  degs_d56_bmp_vs_fgf$logFC < -0.5 & degs_d56_bmp_vs_fgf$FDR < 0.05, '#1F78B4',
  ifelse(degs_d56_bmp_vs_fgf$FDR < 0.05 & degs_d56_bmp_vs_fgf$logFC > 0.5, '#E31A1C',
         'black'))
names(keyvals)[keyvals == '#1F78B4'] <- 'BMP'
names(keyvals)[keyvals == '#E31A1C'] <- 'FGF'
volcano_fgf_vs_bmp <- EnhancedVolcano(degs_d56_bmp_vs_fgf,
                                      lab = rownames(degs_d56_bmp_vs_fgf),
                                      selectLab = rownames(degs_d56_bmp_vs_fgf)[which(names(keyvals) %in% c('FGF', 'BMP'))],
                                      x = 'logFC',
                                      y = 'PValue',
                                      title = 'FGF vs. BMP (DLN, D56)',
                                      pCutoff = 0.05,
                                      FCcutoff = 0.5,
                                      cutoffLineType="blank",
                                      colCustom = keyvals,
                                      gridlines.major = FALSE,
                                      gridlines.minor = FALSE)
volcano_fgf_vs_bmp
ggsave(file=paste0(data_directory, "WatanabeLab_ArealizationOrganoids/", celltype_oi_label, "/D56/organoids_pbulk_allcells_edgeR_qlf_degs_D56_FGFvsBMP_", celltype_oi_label, "_volcano.pdf"))

degs_d97_ctr_vs_bmp <- data.frame(topTags(celltype_qlf_d97_ctr_vs_bmp, n=dim(celltype_qlf_d97_ctr_vs_bmp)[1]))
degs_d97_ctr_vs_bmp <- degs_d97_ctr_vs_bmp[order(-degs_d97_ctr_vs_bmp$logFC),]
degs_d97_ctr_vs_fgf <- data.frame(topTags(celltype_qlf_d97_ctr_vs_fgf, n=dim(celltype_qlf_d97_ctr_vs_fgf)[1]))
degs_d97_ctr_vs_fgf <- degs_d97_ctr_vs_fgf[order(-degs_d97_ctr_vs_fgf$logFC),]
degs_d97_bmp_vs_fgf <- data.frame(topTags(celltype_qlf_d97_bmp_vs_fgf, n=dim(celltype_qlf_d97_bmp_vs_fgf)[1]))
degs_d97_bmp_vs_fgf <- degs_d97_bmp_vs_fgf[order(-degs_d97_bmp_vs_fgf$logFC),]

write.csv(degs_d97_ctr_vs_bmp,
          file=paste0(data_directory, "WatanabeLab_ArealizationOrganoids/organoids_pbulk_edgeR_qlf_degs_D97_BMPvsCTR_", celltype_oi ,".csv"))
write.csv(degs_d97_ctr_vs_fgf,
          file=paste0(data_directory, "WatanabeLab_ArealizationOrganoids/organoids_pbulk_edgeR_qlf_degs_D97_FGFvsCTR_", celltype_oi,".csv"))
write.csv(degs_d97_bmp_vs_fgf,
          file=paste0(data_directory, "WatanabeLab_ArealizationOrganoids/organoids_pbulk_edgeR_qlf_degs_D97_FGFvsBMP_", celltype_oi,".csv"))

keyvals <- ifelse(
  degs_d97_ctr_vs_bmp$logFC < -0.5 & degs_d97_ctr_vs_bmp$FDR < 0.05, '#33A02C',
  ifelse(degs_d97_ctr_vs_bmp$FDR < 0.05 & degs_d97_ctr_vs_bmp$logFC > 0.5, '#1F78B4',
         'black'))
names(keyvals)[keyvals == '#33A02C'] <- 'CTR'
names(keyvals)[keyvals == '#1F78B4'] <- 'BMP'
volcano_ctr_vs_bmp <- EnhancedVolcano(degs_d97_ctr_vs_bmp,
                                      lab = rownames(degs_d97_ctr_vs_bmp),
                                      selectLab = rownames(degs_d97_ctr_vs_bmp)[which(names(keyvals) %in% c('BMP', 'CTR'))],
                                      x = 'logFC',
                                      y = 'PValue',
                                      title = 'BMP vs. CTR (DLN, D97)',
                                      pCutoff = 0.05,
                                      FCcutoff = 0.5,
                                      cutoffLineType="blank",
                                      colCustom = keyvals,
                                      gridlines.major = FALSE,
                                      gridlines.minor = FALSE)
volcano_ctr_vs_bmp
ggsave(file=paste0(data_directory, "WatanabeLab_ArealizationOrganoids/", celltype_oi_label, "/D97/organoids_pbulk_allcells_edgeR_qlf_degs_D97_BMPvsCTR_", celltype_oi_label, "_volcano.pdf"))

keyvals <- ifelse(
  degs_d97_ctr_vs_fgf$logFC < -0.5 & degs_d97_ctr_vs_fgf$FDR < 0.05, '#33A02C',
  ifelse(degs_d97_ctr_vs_fgf$FDR < 0.05 & degs_d97_ctr_vs_fgf$logFC > 0.5, '#E31A1C',
         'black'))
names(keyvals)[keyvals == '#33A02C'] <- 'CTR'
names(keyvals)[keyvals == '#E31A1C'] <- 'FGF'
volcano_ctr_vs_fgf <- EnhancedVolcano(degs_d97_ctr_vs_fgf,
                                      lab = rownames(degs_d97_ctr_vs_fgf),
                                      selectLab = rownames(degs_d97_ctr_vs_fgf)[which(names(keyvals) %in% c('FGF', 'CTR'))],
                                      x = 'logFC',
                                      y = 'PValue',
                                      title = 'FGF vs. CTR (DLN, D97)',
                                      pCutoff = 0.05,
                                      FCcutoff = 0.5,
                                      cutoffLineType="blank",
                                      colCustom = keyvals,
                                      gridlines.major = FALSE,
                                      gridlines.minor = FALSE)
volcano_ctr_vs_fgf
ggsave(file=paste0(data_directory, "WatanabeLab_ArealizationOrganoids/", celltype_oi_label, "/D97/organoids_pbulk_allcells_edgeR_qlf_degs_D97_FGFvsCTR_", celltype_oi_label, "_volcano.pdf"))

keyvals <- ifelse(
  degs_d97_bmp_vs_fgf$logFC < -0.5 & degs_d97_bmp_vs_fgf$FDR < 0.05, '#1F78B4',
  ifelse(degs_d97_bmp_vs_fgf$FDR < 0.05 & degs_d97_bmp_vs_fgf$logFC > 0.5, '#E31A1C',
         'black'))
names(keyvals)[keyvals == '#1F78B4'] <- 'BMP'
names(keyvals)[keyvals == '#E31A1C'] <- 'FGF'
volcano_fgf_vs_bmp <- EnhancedVolcano(degs_d97_bmp_vs_fgf,
                                      lab = rownames(degs_d97_bmp_vs_fgf),
                                      selectLab = rownames(degs_d97_bmp_vs_fgf)[which(names(keyvals) %in% c('FGF', 'BMP'))],
                                      x = 'logFC',
                                      y = 'PValue',
                                      title = 'FGF vs. BMP (DLN, D97)',
                                      pCutoff = 0.05,
                                      FCcutoff = 0.5,
                                      cutoffLineType="blank",
                                      colCustom = keyvals,
                                      gridlines.major = FALSE,
                                      gridlines.minor = FALSE)
volcano_fgf_vs_bmp
ggsave(file=paste0(data_directory, "WatanabeLab_ArealizationOrganoids/", celltype_oi_label, "/D97/organoids_pbulk_allcells_edgeR_qlf_degs_D97_FGFvsBMP_", celltype_oi_label, "_volcano.pdf"))

# ULN + DLN analysis ( present at D56 and D97)
# D56 analysis
celltype_qlf_d56_ctr_vs_bmp <- glmQLFTest(celltype_fit, contrast=c(1, 0, -1, 0, 0, 0))
celltype_qlf_d56_ctr_vs_fgf <- glmQLFTest(celltype_fit, contrast=c(0, 0, -1, 0, 1, 0))
celltype_qlf_d56_bmp_vs_fgf <- glmQLFTest(celltype_fit, contrast=c(-1, 0, 0, 0, 1, 0))

# D97 analysis
celltype_qlf_d97_ctr_vs_bmp <- glmQLFTest(celltype_fit, contrast=c(0, 1, 0, -1, 0, 0))
celltype_qlf_d97_ctr_vs_fgf <- glmQLFTest(celltype_fit, contrast=c(0, 0, 0, -1, 0, 1))
celltype_qlf_d97_bmp_vs_fgf <- glmQLFTest(celltype_fit, contrast=c(0, -1, 0, 0, 0, 1))

degs_d56_ctr_vs_bmp <- data.frame(topTags(celltype_qlf_d56_ctr_vs_bmp, n=dim(celltype_qlf_d56_ctr_vs_bmp)[1]))
degs_d56_ctr_vs_bmp <- degs_d56_ctr_vs_bmp[order(-degs_d56_ctr_vs_bmp$logFC),]
degs_d56_ctr_vs_fgf <- data.frame(topTags(celltype_qlf_d56_ctr_vs_fgf, n=dim(celltype_qlf_d56_ctr_vs_fgf)[1]))
degs_d56_ctr_vs_fgf <- degs_d56_ctr_vs_fgf[order(-degs_d56_ctr_vs_fgf$logFC),]
degs_d56_bmp_vs_fgf <- data.frame(topTags(celltype_qlf_d56_bmp_vs_fgf, n=dim(celltype_qlf_d56_bmp_vs_fgf)[1]))
degs_d56_bmp_vs_fgf <- degs_d56_bmp_vs_fgf[order(-degs_d56_bmp_vs_fgf$logFC),]

write.csv(degs_d56_ctr_vs_bmp,
          file=paste0(data_directory, "WatanabeLab_ArealizationOrganoids/organoids_pbulk_edgeR_qlf_degs_D56_BMPvsCTR_ULN+DLN.csv"))
write.csv(degs_d56_ctr_vs_fgf,
          file=paste0(data_directory, "WatanabeLab_ArealizationOrganoids/organoids_pbulk_edgeR_qlf_degs_D56_FGFvsCTR_ULN+DLN.csv"))
write.csv(degs_d56_bmp_vs_fgf,
          file=paste0(data_directory, "WatanabeLab_ArealizationOrganoids/organoids_pbulk_edgeR_qlf_degs_D56_FGFvsBMP_ULN+DLN.csv"))

degs_d97_ctr_vs_bmp <- data.frame(topTags(celltype_qlf_d97_ctr_vs_bmp, n=dim(celltype_qlf_d97_ctr_vs_bmp)[1]))
degs_d97_ctr_vs_bmp <- degs_d97_ctr_vs_bmp[order(-degs_d97_ctr_vs_bmp$logFC),]
degs_d97_ctr_vs_fgf <- data.frame(topTags(celltype_qlf_d97_ctr_vs_fgf, n=dim(celltype_qlf_d97_ctr_vs_fgf)[1]))
degs_d97_ctr_vs_fgf <- degs_d97_ctr_vs_fgf[order(-degs_d97_ctr_vs_fgf$logFC),]
degs_d97_bmp_vs_fgf <- data.frame(topTags(celltype_qlf_d97_bmp_vs_fgf, n=dim(celltype_qlf_d97_bmp_vs_fgf)[1]))
degs_d97_bmp_vs_fgf <- degs_d97_bmp_vs_fgf[order(-degs_d97_bmp_vs_fgf$logFC),]

write.csv(degs_d97_ctr_vs_bmp,
          file=paste0(data_directory, "WatanabeLab_ArealizationOrganoids/organoids_pbulk_edgeR_qlf_degs_D97_BMPvsCTR_ULN+DLN.csv"))
write.csv(degs_d97_ctr_vs_fgf,
          file=paste0(data_directory, "WatanabeLab_ArealizationOrganoids/organoids_pbulk_edgeR_qlf_degs_D97_FGFvsCTR_ULN+DLN.csv"))
write.csv(degs_d97_bmp_vs_fgf,
          file=paste0(data_directory, "WatanabeLab_ArealizationOrganoids/organoids_pbulk_edgeR_qlf_degs_D97_FGFvsBMP_ULN+DLN.csv"))
#

keyvals <- ifelse(
  degs_d56_ctr_vs_bmp$logFC < -0.5 & degs_d56_ctr_vs_bmp$FDR < 0.05, '#33A02C',
  ifelse(degs_d56_ctr_vs_bmp$FDR < 0.05 & degs_d56_ctr_vs_bmp$logFC > 0.5, '#1F78B4',
         'black'))
names(keyvals)[keyvals == '#33A02C'] <- 'CTR'
names(keyvals)[keyvals == '#1F78B4'] <- 'BMP'
volcano_ctr_vs_bmp <- EnhancedVolcano(degs_d56_ctr_vs_bmp,
                                      lab = rownames(degs_d56_ctr_vs_bmp),
                                      selectLab = rownames(degs_d56_ctr_vs_bmp)[which(names(keyvals) %in% c('BMP', 'CTR'))],
                                      x = 'logFC',
                                      y = 'PValue',
                                      title = 'BMP vs. CTR (ULN+DLN, D56)',
                                      pCutoff = 0.05,
                                      FCcutoff = 0.5,
                                      cutoffLineType="blank",
                                      colCustom = keyvals,
                                      gridlines.major = FALSE,
                                      gridlines.minor = FALSE)
volcano_ctr_vs_bmp
ggsave(file=paste0(data_directory, "WatanabeLab_ArealizationOrganoids/", celltype_oi_label, "/D56/organoids_pbulk_allcells_edgeR_qlf_degs_D56_BMPvsCTR_", celltype_oi_label, "_volcano.pdf"))

keyvals <- ifelse(
  degs_d56_ctr_vs_fgf$logFC < -0.5 & degs_d56_ctr_vs_fgf$FDR < 0.05, '#33A02C',
  ifelse(degs_d56_ctr_vs_fgf$FDR < 0.05 & degs_d56_ctr_vs_fgf$logFC > 0.5, '#E31A1C',
         'black'))
names(keyvals)[keyvals == '#33A02C'] <- 'CTR'
names(keyvals)[keyvals == '#E31A1C'] <- 'FGF'
volcano_ctr_vs_fgf <- EnhancedVolcano(degs_d56_ctr_vs_fgf,
                                      lab = rownames(degs_d56_ctr_vs_fgf),
                                      selectLab = rownames(degs_d56_ctr_vs_fgf)[which(names(keyvals) %in% c('FGF', 'CTR'))],
                                      x = 'logFC',
                                      y = 'PValue',
                                      title = 'FGF vs. CTR (ULN+DLN, D56)',
                                      pCutoff = 0.05,
                                      FCcutoff = 0.5,
                                      cutoffLineType="blank",
                                      colCustom = keyvals,
                                      gridlines.major = FALSE,
                                      gridlines.minor = FALSE)
volcano_ctr_vs_fgf
ggsave(file=paste0(data_directory, "WatanabeLab_ArealizationOrganoids/", celltype_oi_label, "/D56/organoids_pbulk_allcells_edgeR_qlf_degs_D56_FGFvsCTR_", celltype_oi_label, "_volcano.pdf"))

keyvals <- ifelse(
  degs_d56_bmp_vs_fgf$logFC < -0.5 & degs_d56_bmp_vs_fgf$FDR < 0.05, '#1F78B4',
  ifelse(degs_d56_bmp_vs_fgf$FDR < 0.05 & degs_d56_bmp_vs_fgf$logFC > 0.5, '#E31A1C',
         'black'))
names(keyvals)[keyvals == '#1F78B4'] <- 'BMP'
names(keyvals)[keyvals == '#E31A1C'] <- 'FGF'
volcano_fgf_vs_bmp <- EnhancedVolcano(degs_d56_bmp_vs_fgf,
                                      lab = rownames(degs_d56_bmp_vs_fgf),
                                      selectLab = rownames(degs_d56_bmp_vs_fgf)[which(names(keyvals) %in% c('FGF', 'BMP'))],
                                      x = 'logFC',
                                      y = 'PValue',
                                      title = 'FGF vs. BMP (ULN+DLN, D56)',
                                      pCutoff = 0.05,
                                      FCcutoff = 0.5,
                                      cutoffLineType="blank",
                                      colCustom = keyvals,
                                      gridlines.major = FALSE,
                                      gridlines.minor = FALSE)
volcano_fgf_vs_bmp
ggsave(file=paste0(data_directory, "WatanabeLab_ArealizationOrganoids/", celltype_oi_label, "/D56/organoids_pbulk_allcells_edgeR_qlf_degs_D56_FGFvsBMP_", celltype_oi_label, "_volcano.pdf"))

degs_d97_ctr_vs_bmp <- data.frame(topTags(celltype_qlf_d97_ctr_vs_bmp, n=dim(celltype_qlf_d97_ctr_vs_bmp)[1]))
degs_d97_ctr_vs_bmp <- degs_d97_ctr_vs_bmp[order(-degs_d97_ctr_vs_bmp$logFC),]
degs_d97_ctr_vs_fgf <- data.frame(topTags(celltype_qlf_d97_ctr_vs_fgf, n=dim(celltype_qlf_d97_ctr_vs_fgf)[1]))
degs_d97_ctr_vs_fgf <- degs_d97_ctr_vs_fgf[order(-degs_d97_ctr_vs_fgf$logFC),]
degs_d97_bmp_vs_fgf <- data.frame(topTags(celltype_qlf_d97_bmp_vs_fgf, n=dim(celltype_qlf_d97_bmp_vs_fgf)[1]))
degs_d97_bmp_vs_fgf <- degs_d97_bmp_vs_fgf[order(-degs_d97_bmp_vs_fgf$logFC),]

write.csv(degs_d97_ctr_vs_bmp,
          file=paste0(data_directory, "WatanabeLab_ArealizationOrganoids/organoids_pbulk_edgeR_qlf_degs_D97_BMPvsCTR_", celltype_oi ,".csv"))
write.csv(degs_d97_ctr_vs_fgf,
          file=paste0(data_directory, "WatanabeLab_ArealizationOrganoids/organoids_pbulk_edgeR_qlf_degs_D97_FGFvsCTR_", celltype_oi,".csv"))
write.csv(degs_d97_bmp_vs_fgf,
          file=paste0(data_directory, "WatanabeLab_ArealizationOrganoids/organoids_pbulk_edgeR_qlf_degs_D97_FGFvsBMP_", celltype_oi,".csv"))

keyvals <- ifelse(
  degs_d97_ctr_vs_bmp$logFC < -0.5 & degs_d97_ctr_vs_bmp$FDR < 0.05, '#33A02C',
  ifelse(degs_d97_ctr_vs_bmp$FDR < 0.05 & degs_d97_ctr_vs_bmp$logFC > 0.5, '#1F78B4',
         'black'))
names(keyvals)[keyvals == '#33A02C'] <- 'CTR'
names(keyvals)[keyvals == '#1F78B4'] <- 'BMP'
volcano_ctr_vs_bmp <- EnhancedVolcano(degs_d97_ctr_vs_bmp,
                                      lab = rownames(degs_d97_ctr_vs_bmp),
                                      selectLab = rownames(degs_d97_ctr_vs_bmp)[which(names(keyvals) %in% c('BMP', 'CTR'))],
                                      x = 'logFC',
                                      y = 'PValue',
                                      title = 'BMP vs. CTR (ULN+DLN, D97)',
                                      pCutoff = 0.05,
                                      FCcutoff = 0.5,
                                      cutoffLineType="blank",
                                      colCustom = keyvals,
                                      gridlines.major = FALSE,
                                      gridlines.minor = FALSE)
volcano_ctr_vs_bmp
ggsave(file=paste0(data_directory, "WatanabeLab_ArealizationOrganoids/", celltype_oi_label, "/D97/organoids_pbulk_allcells_edgeR_qlf_degs_D97_BMPvsCTR_", celltype_oi_label, "_volcano.pdf"))

keyvals <- ifelse(
  degs_d97_ctr_vs_fgf$logFC < -0.5 & degs_d97_ctr_vs_fgf$FDR < 0.05, '#33A02C',
  ifelse(degs_d97_ctr_vs_fgf$FDR < 0.05 & degs_d97_ctr_vs_fgf$logFC > 0.5, '#E31A1C',
         'black'))
names(keyvals)[keyvals == '#33A02C'] <- 'CTR'
names(keyvals)[keyvals == '#E31A1C'] <- 'FGF'
volcano_ctr_vs_fgf <- EnhancedVolcano(degs_d97_ctr_vs_fgf,
                                      lab = rownames(degs_d97_ctr_vs_fgf),
                                      selectLab = rownames(degs_d97_ctr_vs_fgf)[which(names(keyvals) %in% c('FGF', 'CTR'))],
                                      x = 'logFC',
                                      y = 'PValue',
                                      title = 'FGF vs. CTR (ULN+DLN, D97)',
                                      pCutoff = 0.05,
                                      FCcutoff = 0.5,
                                      cutoffLineType="blank",
                                      colCustom = keyvals,
                                      gridlines.major = FALSE,
                                      gridlines.minor = FALSE)
volcano_ctr_vs_fgf
ggsave(file=paste0(data_directory, "WatanabeLab_ArealizationOrganoids/", celltype_oi_label, "/d97/organoids_pbulk_allcells_edgeR_qlf_degs_D97_FGFvsCTR_", celltype_oi_label, "_volcano.pdf"))

keyvals <- ifelse(
  degs_d97_bmp_vs_fgf$logFC < -0.5 & degs_d97_bmp_vs_fgf$FDR < 0.05, '#1F78B4',
  ifelse(degs_d97_bmp_vs_fgf$FDR < 0.05 & degs_d97_bmp_vs_fgf$logFC > 0.5, '#E31A1C',
         'black'))
names(keyvals)[keyvals == '#1F78B4'] <- 'BMP'
names(keyvals)[keyvals == '#E31A1C'] <- 'FGF'
volcano_fgf_vs_bmp <- EnhancedVolcano(degs_d97_bmp_vs_fgf,
                                      lab = rownames(degs_d97_bmp_vs_fgf),
                                      selectLab = rownames(degs_d97_bmp_vs_fgf)[which(names(keyvals) %in% c('FGF', 'BMP'))],
                                      x = 'logFC',
                                      y = 'PValue',
                                      title = 'FGF vs. BMP (ULN+DLN, D97)',
                                      pCutoff = 0.05,
                                      FCcutoff = 0.5,
                                      cutoffLineType="blank",
                                      colCustom = keyvals,
                                      gridlines.major = FALSE,
                                      gridlines.minor = FALSE)
volcano_fgf_vs_bmp
ggsave(file=paste0(data_directory, "WatanabeLab_ArealizationOrganoids/", celltype_oi_label, "/D97/organoids_pbulk_allcells_edgeR_qlf_degs_D97_FGFvsBMP_", celltype_oi_label, "_volcano.pdf"))

### aRG analysis (only present at D35 and D56)
# D35 analysis
celltype_qlf_d35_ctr_vs_bmp <- glmQLFTest(celltype_fit, contrast=c(1, 0, -1, 0, 0, 0))
celltype_qlf_d35_ctr_vs_fgf <- glmQLFTest(celltype_fit, contrast=c(0, 0, -1, 0, 1, 0))
celltype_qlf_d35_bmp_vs_fgf <- glmQLFTest(celltype_fit, contrast=c(-1, 0, 0, 0, 1, 0))

degs_d35_ctr_vs_bmp <- data.frame(topTags(celltype_qlf_d35_ctr_vs_bmp, n=dim(celltype_qlf_d35_ctr_vs_bmp)[1]))
degs_d35_ctr_vs_bmp <- degs_d35_ctr_vs_bmp[order(-degs_d35_ctr_vs_bmp$logFC),]
degs_d35_ctr_vs_fgf <- data.frame(topTags(celltype_qlf_d35_ctr_vs_fgf, n=dim(celltype_qlf_d35_ctr_vs_fgf)[1]))
degs_d35_ctr_vs_fgf <- degs_d35_ctr_vs_fgf[order(-degs_d35_ctr_vs_fgf$logFC),]
degs_d35_bmp_vs_fgf <- data.frame(topTags(celltype_qlf_d35_bmp_vs_fgf, n=dim(celltype_qlf_d35_bmp_vs_fgf)[1]))
degs_d35_bmp_vs_fgf <- degs_d35_bmp_vs_fgf[order(-degs_d35_bmp_vs_fgf$logFC),]

write.csv(degs_d35_ctr_vs_bmp,
          file=paste0(data_directory, "WatanabeLab_ArealizationOrganoids/organoids_pbulk_edgeR_qlf_degs_D35_BMPvsCTR_", celltype_oi ,".csv"))
write.csv(degs_d35_ctr_vs_fgf,
          file=paste0(data_directory, "WatanabeLab_ArealizationOrganoids/organoids_pbulk_edgeR_qlf_degs_D35_FGFvsCTR_", celltype_oi,".csv"))
write.csv(degs_d35_bmp_vs_fgf,
          file=paste0(data_directory, "WatanabeLab_ArealizationOrganoids/organoids_pbulk_edgeR_qlf_degs_D35_FGFvsBMP_", celltype_oi,".csv"))

# D56 analysis
celltype_qlf_d56_ctr_vs_bmp <- glmQLFTest(celltype_fit, contrast=c(0, 1, 0, -1, 0, 0))
celltype_qlf_d56_ctr_vs_fgf <- glmQLFTest(celltype_fit, contrast=c(0, 0, 0, -1, 0, 1))
celltype_qlf_d56_bmp_vs_fgf <- glmQLFTest(celltype_fit, contrast=c(0, -1, 0, 0, 0, 1))

degs_d56_ctr_vs_bmp <- data.frame(topTags(celltype_qlf_d56_ctr_vs_bmp, n=dim(celltype_qlf_d56_ctr_vs_bmp)[1]))
degs_d56_ctr_vs_bmp <- degs_d56_ctr_vs_bmp[order(-degs_d56_ctr_vs_bmp$logFC),]
degs_d56_ctr_vs_fgf <- data.frame(topTags(celltype_qlf_d56_ctr_vs_fgf, n=dim(celltype_qlf_d56_ctr_vs_fgf)[1]))
degs_d56_ctr_vs_fgf <- degs_d56_ctr_vs_fgf[order(-degs_d56_ctr_vs_fgf$logFC),]
degs_d56_bmp_vs_fgf <- data.frame(topTags(celltype_qlf_d56_bmp_vs_fgf, n=dim(celltype_qlf_d56_bmp_vs_fgf)[1]))
degs_d56_bmp_vs_fgf <- degs_d56_bmp_vs_fgf[order(-degs_d56_bmp_vs_fgf$logFC),]


write.csv(degs_d56_ctr_vs_bmp,
          file=paste0(data_directory, "WatanabeLab_ArealizationOrganoids/organoids_pbulk_edgeR_qlf_degs_D56_BMPvsCTR_", celltype_oi ,".csv"))
write.csv(degs_d56_ctr_vs_fgf,
          file=paste0(data_directory, "WatanabeLab_ArealizationOrganoids/organoids_pbulk_edgeR_qlf_degs_D56_FGFvsCTR_", celltype_oi,".csv"))
write.csv(degs_d56_bmp_vs_fgf,
          file=paste0(data_directory, "WatanabeLab_ArealizationOrganoids/organoids_pbulk_edgeR_qlf_degs_D56_FGFvsBMP_", celltype_oi,".csv"))

keyvals <- ifelse(
  degs_d35_ctr_vs_bmp$logFC < -0.5 & degs_d35_ctr_vs_bmp$FDR < 0.05, '#33A02C',
  ifelse(degs_d35_ctr_vs_bmp$FDR < 0.05 & degs_d35_ctr_vs_bmp$logFC > 0.5, '#1F78B4',
         'black'))
names(keyvals)[keyvals == '#33A02C'] <- 'CTR'
names(keyvals)[keyvals == '#1F78B4'] <- 'BMP'
volcano_ctr_vs_bmp <- EnhancedVolcano(degs_d35_ctr_vs_bmp,
                                      lab = rownames(degs_d35_ctr_vs_bmp),
                                      selectLab = rownames(degs_d35_ctr_vs_bmp)[which(names(keyvals) %in% c('BMP', 'CTR'))],
                                      x = 'logFC',
                                      y = 'PValue',
                                      title = 'BMP vs. CTR (aRG, D35)',
                                      pCutoff = 0.05,
                                      FCcutoff = 0.5,
                                      cutoffLineType="blank",
                                      colCustom = keyvals,
                                      gridlines.major = FALSE,
                                      gridlines.minor = FALSE)
volcano_ctr_vs_bmp
ggsave(file=paste0(data_directory, "WatanabeLab_ArealizationOrganoids/", celltype_oi_label, "/D35/organoids_pbulk_allcells_edgeR_qlf_degs_D35_BMPvsCTR_", celltype_oi_label, "_volcano.pdf"))

keyvals <- ifelse(
  degs_d35_ctr_vs_fgf$logFC < -0.5 & degs_d35_ctr_vs_fgf$FDR < 0.05, '#33A02C',
  ifelse(degs_d35_ctr_vs_fgf$FDR < 0.05 & degs_d35_ctr_vs_fgf$logFC > 0.5, '#E31A1C',
         'black'))
names(keyvals)[keyvals == '#33A02C'] <- 'CTR'
names(keyvals)[keyvals == '#E31A1C'] <- 'FGF'
volcano_ctr_vs_fgf <- EnhancedVolcano(degs_d35_ctr_vs_fgf,
                                      lab = rownames(degs_d35_ctr_vs_fgf),
                                      selectLab = rownames(degs_d35_ctr_vs_fgf)[which(names(keyvals) %in% c('FGF', 'CTR'))],
                                      x = 'logFC',
                                      y = 'PValue',
                                      title = 'FGF vs. CTR (aRG, D35)',
                                      pCutoff = 0.05,
                                      FCcutoff = 0.5,
                                      cutoffLineType="blank",
                                      colCustom = keyvals,
                                      gridlines.major = FALSE,
                                      gridlines.minor = FALSE)
volcano_ctr_vs_fgf
ggsave(file=paste0(data_directory, "WatanabeLab_ArealizationOrganoids/", celltype_oi_label, "/D35/organoids_pbulk_allcells_edgeR_qlf_degs_D35_FGFvsCTR_", celltype_oi_label, "_volcano.pdf"))

keyvals <- ifelse(
  degs_d35_bmp_vs_fgf$logFC < -0.5 & degs_d35_bmp_vs_fgf$FDR < 0.05, '#1F78B4',
  ifelse(degs_d35_bmp_vs_fgf$FDR < 0.05 & degs_d35_bmp_vs_fgf$logFC > 0.5, '#E31A1C',
         'black'))
names(keyvals)[keyvals == '#1F78B4'] <- 'BMP'
names(keyvals)[keyvals == '#E31A1C'] <- 'FGF'
volcano_fgf_vs_bmp <- EnhancedVolcano(degs_d35_bmp_vs_fgf,
                                      lab = rownames(degs_d35_bmp_vs_fgf),
                                      selectLab = rownames(degs_d35_bmp_vs_fgf)[which(names(keyvals) %in% c('FGF', 'BMP'))],
                                      x = 'logFC',
                                      y = 'PValue',
                                      title = 'FGF vs. BMP (aRG, D35)',
                                      pCutoff = 0.05,
                                      FCcutoff = 0.5,
                                      cutoffLineType="blank",
                                      colCustom = keyvals,
                                      gridlines.major = FALSE,
                                      gridlines.minor = FALSE)
volcano_fgf_vs_bmp
ggsave(file=paste0(data_directory, "WatanabeLab_ArealizationOrganoids/", celltype_oi_label, "/D35/organoids_pbulk_allcells_edgeR_qlf_degs_D35_FGFvsBMP_", celltype_oi_label, "_volcano.pdf"))

degs_d56_ctr_vs_bmp <- data.frame(topTags(celltype_qlf_d56_ctr_vs_bmp, n=dim(celltype_qlf_d56_ctr_vs_bmp)[1]))
degs_d56_ctr_vs_bmp <- degs_d56_ctr_vs_bmp[order(-degs_d56_ctr_vs_bmp$logFC),]
degs_d56_ctr_vs_fgf <- data.frame(topTags(celltype_qlf_d56_ctr_vs_fgf, n=dim(celltype_qlf_d56_ctr_vs_fgf)[1]))
degs_d56_ctr_vs_fgf <- degs_d56_ctr_vs_fgf[order(-degs_d56_ctr_vs_fgf$logFC),]
degs_d56_bmp_vs_fgf <- data.frame(topTags(celltype_qlf_d56_bmp_vs_fgf, n=dim(celltype_qlf_d56_bmp_vs_fgf)[1]))
degs_d56_bmp_vs_fgf <- degs_d56_bmp_vs_fgf[order(-degs_d56_bmp_vs_fgf$logFC),]

write.csv(degs_d56_ctr_vs_bmp,
          file=paste0(data_directory, "WatanabeLab_ArealizationOrganoids/organoids_pbulk_edgeR_qlf_degs_D56_BMPvsCTR_", celltype_oi ,".csv"))
write.csv(degs_d56_ctr_vs_fgf,
          file=paste0(data_directory, "WatanabeLab_ArealizationOrganoids/organoids_pbulk_edgeR_qlf_degs_D56_FGFvsCTR_", celltype_oi,".csv"))
write.csv(degs_d56_bmp_vs_fgf,
          file=paste0(data_directory, "WatanabeLab_ArealizationOrganoids/organoids_pbulk_edgeR_qlf_degs_D56_FGFvsBMP_", celltype_oi,".csv"))

keyvals <- ifelse(
  degs_d56_ctr_vs_bmp$logFC < -0.5 & degs_d56_ctr_vs_bmp$FDR < 0.05, '#33A02C',
  ifelse(degs_d56_ctr_vs_bmp$FDR < 0.05 & degs_d56_ctr_vs_bmp$logFC > 0.5, '#1F78B4',
         'black'))
names(keyvals)[keyvals == '#33A02C'] <- 'CTR'
names(keyvals)[keyvals == '#1F78B4'] <- 'BMP'
volcano_ctr_vs_bmp <- EnhancedVolcano(degs_d56_ctr_vs_bmp,
                                      lab = rownames(degs_d56_ctr_vs_bmp),
                                      selectLab = rownames(degs_d56_ctr_vs_bmp)[which(names(keyvals) %in% c('BMP', 'CTR'))],
                                      x = 'logFC',
                                      y = 'PValue',
                                      title = 'BMP vs. CTR (aRG, D56)',
                                      pCutoff = 0.05,
                                      FCcutoff = 0.5,
                                      cutoffLineType="blank",
                                      colCustom = keyvals,
                                      gridlines.major = FALSE,
                                      gridlines.minor = FALSE)
volcano_ctr_vs_bmp
ggsave(file=paste0(data_directory, "WatanabeLab_ArealizationOrganoids/", celltype_oi_label, "/D56/organoids_pbulk_allcells_edgeR_qlf_degs_D56_BMPvsCTR_", celltype_oi_label, "_volcano.pdf"))

keyvals <- ifelse(
  degs_d56_ctr_vs_fgf$logFC < -0.5 & degs_d56_ctr_vs_fgf$FDR < 0.05, '#33A02C',
  ifelse(degs_d56_ctr_vs_fgf$FDR < 0.05 & degs_d56_ctr_vs_fgf$logFC > 0.5, '#E31A1C',
         'black'))
names(keyvals)[keyvals == '#33A02C'] <- 'CTR'
names(keyvals)[keyvals == '#E31A1C'] <- 'FGF'
volcano_ctr_vs_fgf <- EnhancedVolcano(degs_d56_ctr_vs_fgf,
                                      lab = rownames(degs_d56_ctr_vs_fgf),
                                      selectLab = rownames(degs_d56_ctr_vs_fgf)[which(names(keyvals) %in% c('FGF', 'CTR'))],
                                      x = 'logFC',
                                      y = 'PValue',
                                      title = 'FGF vs. CTR (aRG, D56)',
                                      pCutoff = 0.05,
                                      FCcutoff = 0.5,
                                      cutoffLineType="blank",
                                      colCustom = keyvals,
                                      gridlines.major = FALSE,
                                      gridlines.minor = FALSE)
volcano_ctr_vs_fgf
ggsave(file=paste0(data_directory, "WatanabeLab_ArealizationOrganoids/", celltype_oi_label, "/D56/organoids_pbulk_allcells_edgeR_qlf_degs_D56_FGFvsCTR_", celltype_oi_label, "_volcano.pdf"))

keyvals <- ifelse(
  degs_d56_bmp_vs_fgf$logFC < -0.5 & degs_d56_bmp_vs_fgf$FDR < 0.05, '#1F78B4',
  ifelse(degs_d56_bmp_vs_fgf$FDR < 0.05 & degs_d56_bmp_vs_fgf$logFC > 0.5, '#E31A1C',
         'black'))
names(keyvals)[keyvals == '#1F78B4'] <- 'BMP'
names(keyvals)[keyvals == '#E31A1C'] <- 'FGF'
volcano_fgf_vs_bmp <- EnhancedVolcano(degs_d56_bmp_vs_fgf,
                                      lab = rownames(degs_d56_bmp_vs_fgf),
                                      selectLab = rownames(degs_d56_bmp_vs_fgf)[which(names(keyvals) %in% c('FGF', 'BMP'))],
                                      x = 'logFC',
                                      y = 'PValue',
                                      title = 'FGF vs. BMP (aRG, D56)',
                                      pCutoff = 0.05,
                                      FCcutoff = 0.5,
                                      cutoffLineType="blank",
                                      colCustom = keyvals,
                                      gridlines.major = FALSE,
                                      gridlines.minor = FALSE)
volcano_fgf_vs_bmp
ggsave(file=paste0(data_directory, "WatanabeLab_ArealizationOrganoids/", celltype_oi_label, "/D56/organoids_pbulk_allcells_edgeR_qlf_degs_D56_FGFvsBMP_", celltype_oi_label, "_volcano.pdf"))

# aRG/oRG analysis (only available at D97)
celltype_qlf_ctr_vs_bmp <- glmQLFTest(celltype_fit, coef=2)
celltype_qlf_ctr_vs_fgf <- glmQLFTest(celltype_fit, coef=3)
celltype_qlf_bmp_vs_fgf <- glmQLFTest(celltype_fit, contrast=c(0, -1, 1))

degs_ctr_vs_bmp <- data.frame(topTags(celltype_qlf_ctr_vs_bmp, n=dim(celltype_qlf_ctr_vs_bmp)[1]))
degs_ctr_vs_bmp <- degs_ctr_vs_bmp[order(-degs_ctr_vs_bmp$logFC),]
degs_ctr_vs_fgf <- data.frame(topTags(celltype_qlf_ctr_vs_fgf, n=dim(celltype_qlf_ctr_vs_fgf)[1]))
degs_ctr_vs_fgf <- degs_ctr_vs_fgf[order(-degs_ctr_vs_fgf$logFC),]
degs_bmp_vs_fgf <- data.frame(topTags(celltype_qlf_bmp_vs_fgf, n=dim(celltype_qlf_bmp_vs_fgf)[1]))
degs_bmp_vs_fgf <- degs_bmp_vs_fgf[order(-degs_bmp_vs_fgf$logFC),]
#
# celltype_oi_label <- gsub("/", "-", celltype_oi)
write.csv(degs_ctr_vs_bmp,
          file=paste0(data_directory, "WatanabeLab_ArealizationOrganoids/organoids_pbulk_edgeR_qlf_degs_D97_BMPvsCTR_", celltype_oi_label ,".csv"))
write.csv(degs_ctr_vs_fgf,
          file=paste0(data_directory, "WatanabeLab_ArealizationOrganoids/organoids_pbulk_edgeR_qlf_degs_D97_FGFvsCTR_", celltype_oi_label,".csv"))
write.csv(degs_bmp_vs_fgf,
          file=paste0(data_directory, "WatanabeLab_ArealizationOrganoids/organoids_pbulk_edgeR_qlf_degs_D97_FGFvsBMP_", celltype_oi_label,".csv"))
#
# #
keyvals <- ifelse(
  degs_ctr_vs_bmp$logFC < -0.5 & degs_ctr_vs_bmp$FDR < 0.05, '#33A02C',
  ifelse(degs_ctr_vs_bmp$FDR < 0.05 & degs_ctr_vs_bmp$logFC > 0.5, '#1F78B4',
         'black'))
names(keyvals)[keyvals == '#33A02C'] <- 'CTR'
names(keyvals)[keyvals == '#1F78B4'] <- 'BMP'
volcano_ctr_vs_bmp <- EnhancedVolcano(degs_ctr_vs_bmp,
                                      lab = rownames(degs_ctr_vs_bmp),
                                      selectLab = rownames(degs_ctr_vs_bmp)[which(names(keyvals) %in% c('BMP', 'CTR'))],
                                      x = 'logFC',
                                      y = 'PValue',
                                      title = paste0('BMP vs. CTR (', celltype_oi, ', D97)'),
                                      # title = paste0('BMP vs. CTR (', paste0(celltype_oi, collapse='+'), ', any timepoint)'),
                                      pCutoff = 0.05,
                                      FCcutoff = 0.5,
                                      cutoffLineType="blank",
                                      colCustom = keyvals,
                                      gridlines.major = FALSE,
                                      gridlines.minor = FALSE)
volcano_ctr_vs_bmp

ggsave(file=paste0(data_directory, "WatanabeLab_ArealizationOrganoids/aRG+oRG/D97/organoids_pbulk_edgeR_qlf_degs_BMPvsCTR_", celltype_oi_label, "_volcano.pdf"))

# FGF vs CTR at any timepoint
keyvals <- ifelse(
  degs_ctr_vs_fgf$logFC < -0.5 & degs_ctr_vs_fgf$FDR < 0.05, '#33A02C',
  ifelse(degs_ctr_vs_fgf$FDR < 0.05 & degs_ctr_vs_fgf$logFC > 0.5, '#E31A1C',
         'black'))
names(keyvals)[keyvals == '#33A02C'] <- 'CTR'
names(keyvals)[keyvals == '#E31A1C'] <- 'FGF'
volcano_ctr_vs_fgf <- EnhancedVolcano(degs_ctr_vs_fgf,
                                      lab = rownames(degs_ctr_vs_fgf),
                                      selectLab = rownames(degs_ctr_vs_fgf)[which(names(keyvals) %in% c('FGF', 'CTR'))],
                                      x = 'logFC',
                                      y = 'PValue',
                                      title = paste0('FGF vs. CTR (', celltype_oi, ', D97)'),
                                      # title = paste0('FGF vs. CTR (', paste0(celltype_oi, collapse='+'), ', any timepoint)'),
                                      pCutoff = 0.05,
                                      FCcutoff = 0.5,
                                      cutoffLineType="blank",
                                      colCustom = keyvals,
                                      gridlines.major = FALSE,
                                      gridlines.minor = FALSE)
volcano_ctr_vs_fgf
ggsave(file=paste0(data_directory, "WatanabeLab_ArealizationOrganoids/aRG+oRG/D97/organoids_pbulk_allcells_edgeR_qlf_degs_FGFvsCTR_", celltype_oi_label, "_volcano.pdf"))


# FGF vs BMP at any timepoint
keyvals <- ifelse(
  degs_bmp_vs_fgf$logFC < -0.5 & degs_bmp_vs_fgf$FDR < 0.05, '#1F78B4',
  ifelse(degs_bmp_vs_fgf$FDR < 0.05 & degs_bmp_vs_fgf$logFC > 0.5, '#E31A1C',
         'black'))
names(keyvals)[keyvals == '#1F78B4'] <- 'BMP'
names(keyvals)[keyvals == '#E31A1C'] <- 'FGF'
volcano_fgf_vs_bmp <- EnhancedVolcano(degs_bmp_vs_fgf,
                                      lab = rownames(degs_bmp_vs_fgf),
                                      selectLab = rownames(degs_bmp_vs_fgf)[which(names(keyvals) %in% c('FGF', 'BMP'))],
                                      x = 'logFC',
                                      y = 'PValue',
                                      title = paste0('FGF vs. BMP (', celltype_oi, ', D97)'),
#                                       # title = paste0('FGF vs. BMP (', paste0(celltype_oi, collapse='+'), ', any timepoint)'),
                                      pCutoff = 0.05,
                                      FCcutoff = 0.5,
                                      cutoffLineType="blank",
                                      colCustom = keyvals,
                                      gridlines.major = FALSE,
                                      gridlines.minor = FALSE)
volcano_fgf_vs_bmp
ggsave(file=paste0(data_directory, "WatanabeLab_ArealizationOrganoids/aRG+oRG/D97/organoids_pbulk_allcells_edgeR_qlf_degs_FGFvsBMP_", celltype_oi_label, "_volcano.pdf"))


# PP/CR analysis (only available at D35)
celltype_qlf_ctr_vs_bmp <- glmQLFTest(celltype_fit, contrast=c(1, 0,-1))
celltype_qlf_ctr_vs_fgf <- glmQLFTest(celltype_fit, contrast=c(0, 1,-1))
celltype_qlf_bmp_vs_fgf <- glmQLFTest(celltype_fit, contrast=c(-1, 1,0))

degs_ctr_vs_bmp <- data.frame(topTags(celltype_qlf_ctr_vs_bmp, n=dim(celltype_qlf_ctr_vs_bmp)[1]))
degs_ctr_vs_bmp <- degs_ctr_vs_bmp[order(-degs_ctr_vs_bmp$logFC),]
degs_ctr_vs_fgf <- data.frame(topTags(celltype_qlf_ctr_vs_fgf, n=dim(celltype_qlf_ctr_vs_fgf)[1]))
degs_ctr_vs_fgf <- degs_ctr_vs_fgf[order(-degs_ctr_vs_fgf$logFC),]
degs_bmp_vs_fgf <- data.frame(topTags(celltype_qlf_bmp_vs_fgf, n=dim(celltype_qlf_bmp_vs_fgf)[1]))
degs_bmp_vs_fgf <- degs_bmp_vs_fgf[order(-degs_bmp_vs_fgf$logFC),]

celltype_oi_label <- gsub("/", "-", celltype_oi)
write.csv(degs_ctr_vs_bmp,
          file=paste0(data_directory, "WatanabeLab_ArealizationOrganoids/organoids_pbulk_edgeR_qlf_degs_D35_BMPvsCTR_", celltype_oi_label ,".csv"))
write.csv(degs_ctr_vs_fgf,
          file=paste0(data_directory, "WatanabeLab_ArealizationOrganoids/organoids_pbulk_edgeR_qlf_degs_D35_FGFvsCTR_", celltype_oi_label,".csv"))
write.csv(degs_bmp_vs_fgf,
          file=paste0(data_directory, "WatanabeLab_ArealizationOrganoids/organoids_pbulk_edgeR_qlf_degs_D35_FGFvsBMP_", celltype_oi_label,".csv"))


## IP analysis (present at all timepoints)
# D35 analysis
celltype_qlf_d35_ctr_vs_bmp <- glmQLFTest(celltype_fit, contrast=c(1, 0, 0, -1, 0, 0, 0, 0, 0))
celltype_qlf_d35_ctr_vs_fgf <- glmQLFTest(celltype_fit, contrast=c(0, 0, 0, -1, 0, 0, 1, 0, 0))
celltype_qlf_d35_bmp_vs_fgf <- glmQLFTest(celltype_fit, contrast=c(-1, 0, 0, 0, 0, 0, 1, 0, 0))

# D56 analysis
celltype_qlf_d56_ctr_vs_bmp <- glmQLFTest(celltype_fit, contrast=c(0, 1, 0, 0, -1, 0, 0, 0, 0))
celltype_qlf_d56_ctr_vs_fgf <- glmQLFTest(celltype_fit, contrast=c(0, 0, 0, 0, -1, 0, 0, 1, 0))
celltype_qlf_d56_bmp_vs_fgf <- glmQLFTest(celltype_fit, contrast=c(0, -1, 0, 0, 0, 0, 0, 1, 0))

# D97 analysis
celltype_qlf_d97_ctr_vs_bmp <- glmQLFTest(celltype_fit, contrast=c(0, 0, 1, 0, 0, -1, 0, 0, 0))
celltype_qlf_d97_ctr_vs_fgf <- glmQLFTest(celltype_fit, contrast=c(0, 0, 0, 0, 0, -1, 0, 0, 1))
celltype_qlf_d97_bmp_vs_fgf <- glmQLFTest(celltype_fit, contrast=c(0, 0, -1, 0, 0, 0, 0, 0, 1))

# Extract the DEGs
degs_d35_ctr_vs_bmp <- data.frame(topTags(celltype_qlf_d35_ctr_vs_bmp, n=dim(celltype_qlf_d35_ctr_vs_bmp)[1]))
degs_d35_ctr_vs_bmp <- degs_d35_ctr_vs_bmp[order(-degs_d35_ctr_vs_bmp$logFC),]
degs_d35_ctr_vs_fgf <- data.frame(topTags(celltype_qlf_d35_ctr_vs_fgf, n=dim(celltype_qlf_d35_ctr_vs_fgf)[1]))
degs_d35_ctr_vs_fgf <- degs_d35_ctr_vs_fgf[order(-degs_d35_ctr_vs_fgf$logFC),]
degs_d35_bmp_vs_fgf <- data.frame(topTags(celltype_qlf_d35_bmp_vs_fgf, n=dim(celltype_qlf_d35_bmp_vs_fgf)[1]))
degs_d35_bmp_vs_fgf <- degs_d35_bmp_vs_fgf[order(-degs_d35_bmp_vs_fgf$logFC),]
celltype_oi_label <- gsub("/", "-", celltype_oi)

write.csv(degs_d35_ctr_vs_bmp,
          file=paste0(data_directory, "WatanabeLab_ArealizationOrganoids/", celltype_oi_label, "/D35/organoids_pbulk_edgeR_qlf_degs_D35_BMPvsCTR_", celltype_oi ,".csv"))
write.csv(degs_d35_ctr_vs_fgf,
          file=paste0(data_directory, "WatanabeLab_ArealizationOrganoids/organoids_pbulk_edgeR_qlf_degs_D35_FGFvsCTR_", celltype_oi,".csv"))
write.csv(degs_d35_bmp_vs_fgf,
          file=paste0(data_directory, "WatanabeLab_ArealizationOrganoids/organoids_pbulk_edgeR_qlf_degs_D35_FGFvsBMP_", celltype_oi,".csv"))

keyvals <- ifelse(
  degs_d35_ctr_vs_bmp$logFC < -0.5 & degs_d35_ctr_vs_bmp$FDR < 0.05, '#33A02C',
  ifelse(degs_d35_ctr_vs_bmp$FDR < 0.05 & degs_d35_ctr_vs_bmp$logFC > 0.5, '#1F78B4',
         'black'))
names(keyvals)[keyvals == '#33A02C'] <- 'CTR'
names(keyvals)[keyvals == '#1F78B4'] <- 'BMP'
volcano_ctr_vs_bmp <- EnhancedVolcano(degs_d35_ctr_vs_bmp,
                                      lab = rownames(degs_d35_ctr_vs_bmp),
                                      selectLab = rownames(degs_d35_ctr_vs_bmp)[which(names(keyvals) %in% c('BMP', 'CTR'))],
                                      x = 'logFC',
                                      y = 'PValue',
                                      title = 'BMP vs. CTR (IP, D35)',
                                      pCutoff = 0.05,
                                      FCcutoff = 0.5,
                                      cutoffLineType="blank",
                                      colCustom = keyvals,
                                      gridlines.major = FALSE,
                                      gridlines.minor = FALSE)
volcano_ctr_vs_bmp
ggsave(file=paste0(data_directory, "WatanabeLab_ArealizationOrganoids/", celltype_oi_label, "/D35/organoids_pbulk_allcells_edgeR_qlf_degs_D35_BMPvsCTR_", celltype_oi_label, "_volcano.pdf"))

keyvals <- ifelse(
  degs_d35_ctr_vs_fgf$logFC < -0.5 & degs_d35_ctr_vs_fgf$FDR < 0.05, '#33A02C',
  ifelse(degs_d35_ctr_vs_fgf$FDR < 0.05 & degs_d35_ctr_vs_fgf$logFC > 0.5, '#E31A1C',
         'black'))
names(keyvals)[keyvals == '#33A02C'] <- 'CTR'
names(keyvals)[keyvals == '#E31A1C'] <- 'FGF'
volcano_ctr_vs_fgf <- EnhancedVolcano(degs_d35_ctr_vs_fgf,
                                      lab = rownames(degs_d35_ctr_vs_fgf),
                                      selectLab = rownames(degs_d35_ctr_vs_fgf)[which(names(keyvals) %in% c('FGF', 'CTR'))],
                                      x = 'logFC',
                                      y = 'PValue',
                                      title = 'FGF vs. CTR (IP, D35)',
                                      pCutoff = 0.05,
                                      FCcutoff = 0.5,
                                      cutoffLineType="blank",
                                      colCustom = keyvals,
                                      gridlines.major = FALSE,
                                      gridlines.minor = FALSE)
volcano_ctr_vs_fgf
ggsave(file=paste0(data_directory, "WatanabeLab_ArealizationOrganoids/", celltype_oi_label, "/D35/organoids_pbulk_allcells_edgeR_qlf_degs_D35_FGFvsCTR_", celltype_oi_label, "_volcano.pdf"))

keyvals <- ifelse(
  degs_d35_bmp_vs_fgf$logFC < -0.5 & degs_d35_bmp_vs_fgf$FDR < 0.05, '#1F78B4',
  ifelse(degs_d35_bmp_vs_fgf$FDR < 0.05 & degs_d35_bmp_vs_fgf$logFC > 0.5, '#E31A1C',
         'black'))
names(keyvals)[keyvals == '#1F78B4'] <- 'BMP'
names(keyvals)[keyvals == '#E31A1C'] <- 'FGF'
volcano_fgf_vs_bmp <- EnhancedVolcano(degs_d35_bmp_vs_fgf,
                                      lab = rownames(degs_d35_bmp_vs_fgf),
                                      selectLab = rownames(degs_d35_bmp_vs_fgf)[which(names(keyvals) %in% c('FGF', 'BMP'))],
                                      x = 'logFC',
                                      y = 'PValue',
                                      title = 'FGF vs. BMP (IP, D35)',
                                      pCutoff = 0.05,
                                      FCcutoff = 0.5,
                                      cutoffLineType="blank",
                                      colCustom = keyvals,
                                      gridlines.major = FALSE,
                                      gridlines.minor = FALSE)
volcano_fgf_vs_bmp
ggsave(file=paste0(data_directory, "WatanabeLab_ArealizationOrganoids/", celltype_oi_label, "/D35/organoids_pbulk_allcells_edgeR_qlf_degs_D35_FGFvsBMP_", celltype_oi_label, "_volcano.pdf"))

degs_d56_ctr_vs_bmp <- data.frame(topTags(celltype_qlf_d56_ctr_vs_bmp, n=dim(celltype_qlf_d56_ctr_vs_bmp)[1]))
degs_d56_ctr_vs_bmp <- degs_d56_ctr_vs_bmp[order(-degs_d56_ctr_vs_bmp$logFC),]
degs_d56_ctr_vs_fgf <- data.frame(topTags(celltype_qlf_d56_ctr_vs_fgf, n=dim(celltype_qlf_d56_ctr_vs_fgf)[1]))
degs_d56_ctr_vs_fgf <- degs_d56_ctr_vs_fgf[order(-degs_d56_ctr_vs_fgf$logFC),]
degs_d56_bmp_vs_fgf <- data.frame(topTags(celltype_qlf_d56_bmp_vs_fgf, n=dim(celltype_qlf_d56_bmp_vs_fgf)[1]))
degs_d56_bmp_vs_fgf <- degs_d56_bmp_vs_fgf[order(-degs_d56_bmp_vs_fgf$logFC),]

write.csv(degs_d56_ctr_vs_bmp,
          file=paste0(data_directory, "WatanabeLab_ArealizationOrganoids/organoids_pbulk_edgeR_qlf_degs_D56_BMPvsCTR_", celltype_oi ,".csv"))
write.csv(degs_d56_ctr_vs_fgf,
          file=paste0(data_directory, "WatanabeLab_ArealizationOrganoids/organoids_pbulk_edgeR_qlf_degs_D56_FGFvsCTR_", celltype_oi,".csv"))
write.csv(degs_d56_bmp_vs_fgf,
          file=paste0(data_directory, "WatanabeLab_ArealizationOrganoids/organoids_pbulk_edgeR_qlf_degs_D56_FGFvsBMP_", celltype_oi,".csv"))

keyvals <- ifelse(
  degs_d56_ctr_vs_bmp$logFC < -0.5 & degs_d56_ctr_vs_bmp$FDR < 0.05, '#33A02C',
  ifelse(degs_d56_ctr_vs_bmp$FDR < 0.05 & degs_d56_ctr_vs_bmp$logFC > 0.5, '#1F78B4',
         'black'))
names(keyvals)[keyvals == '#33A02C'] <- 'CTR'
names(keyvals)[keyvals == '#1F78B4'] <- 'BMP'
volcano_ctr_vs_bmp <- EnhancedVolcano(degs_d56_ctr_vs_bmp,
                                      lab = rownames(degs_d56_ctr_vs_bmp),
                                      selectLab = rownames(degs_d56_ctr_vs_bmp)[which(names(keyvals) %in% c('BMP', 'CTR'))],
                                      x = 'logFC',
                                      y = 'PValue',
                                      title = 'BMP vs. CTR (IP, D56)',
                                      pCutoff = 0.05,
                                      FCcutoff = 0.5,
                                      cutoffLineType="blank",
                                      colCustom = keyvals,
                                      gridlines.major = FALSE,
                                      gridlines.minor = FALSE)
volcano_ctr_vs_bmp
ggsave(file=paste0(data_directory, "WatanabeLab_ArealizationOrganoids/", celltype_oi_label, "/D56/organoids_pbulk_allcells_edgeR_qlf_degs_D56_BMPvsCTR_", celltype_oi_label, "_volcano.pdf"))

keyvals <- ifelse(
  degs_d56_ctr_vs_fgf$logFC < -0.5 & degs_d56_ctr_vs_fgf$FDR < 0.05, '#33A02C',
  ifelse(degs_d56_ctr_vs_fgf$FDR < 0.05 & degs_d56_ctr_vs_fgf$logFC > 0.5, '#E31A1C',
         'black'))
names(keyvals)[keyvals == '#33A02C'] <- 'CTR'
names(keyvals)[keyvals == '#E31A1C'] <- 'FGF'
volcano_ctr_vs_fgf <- EnhancedVolcano(degs_d56_ctr_vs_fgf,
                                      lab = rownames(degs_d56_ctr_vs_fgf),
                                      selectLab = rownames(degs_d56_ctr_vs_fgf)[which(names(keyvals) %in% c('FGF', 'CTR'))],
                                      x = 'logFC',
                                      y = 'PValue',
                                      title = 'FGF vs. CTR (IP, D56)',
                                      pCutoff = 0.05,
                                      FCcutoff = 0.5,
                                      cutoffLineType="blank",
                                      colCustom = keyvals,
                                      gridlines.major = FALSE,
                                      gridlines.minor = FALSE)
volcano_ctr_vs_fgf
ggsave(file=paste0(data_directory, "WatanabeLab_ArealizationOrganoids/", celltype_oi_label, "/D56/organoids_pbulk_allcells_edgeR_qlf_degs_D56_FGFvsCTR_", celltype_oi_label, "_volcano.pdf"))

keyvals <- ifelse(
  degs_d56_bmp_vs_fgf$logFC < -0.5 & degs_d56_bmp_vs_fgf$FDR < 0.05, '#1F78B4',
  ifelse(degs_d56_bmp_vs_fgf$FDR < 0.05 & degs_d56_bmp_vs_fgf$logFC > 0.5, '#E31A1C',
         'black'))
names(keyvals)[keyvals == '#1F78B4'] <- 'BMP'
names(keyvals)[keyvals == '#E31A1C'] <- 'FGF'
volcano_fgf_vs_bmp <- EnhancedVolcano(degs_d56_bmp_vs_fgf,
                                      lab = rownames(degs_d56_bmp_vs_fgf),
                                      selectLab = rownames(degs_d56_bmp_vs_fgf)[which(names(keyvals) %in% c('FGF', 'BMP'))],
                                      x = 'logFC',
                                      y = 'PValue',
                                      title = 'FGF vs. BMP (IP, D56)',
                                      pCutoff = 0.05,
                                      FCcutoff = 0.5,
                                      cutoffLineType="blank",
                                      colCustom = keyvals,
                                      gridlines.major = FALSE,
                                      gridlines.minor = FALSE)
volcano_fgf_vs_bmp
ggsave(file=paste0(data_directory, "WatanabeLab_ArealizationOrganoids/", celltype_oi_label, "/D56/organoids_pbulk_allcells_edgeR_qlf_degs_D56_FGFvsBMP_", celltype_oi_label, "_volcano.pdf"))

degs_d97_ctr_vs_bmp <- data.frame(topTags(celltype_qlf_d97_ctr_vs_bmp, n=dim(celltype_qlf_d97_ctr_vs_bmp)[1]))
degs_d97_ctr_vs_bmp <- degs_d97_ctr_vs_bmp[order(-degs_d97_ctr_vs_bmp$logFC),]
degs_d97_ctr_vs_fgf <- data.frame(topTags(celltype_qlf_d97_ctr_vs_fgf, n=dim(celltype_qlf_d97_ctr_vs_fgf)[1]))
degs_d97_ctr_vs_fgf <- degs_d97_ctr_vs_fgf[order(-degs_d97_ctr_vs_fgf$logFC),]
degs_d97_bmp_vs_fgf <- data.frame(topTags(celltype_qlf_d97_bmp_vs_fgf, n=dim(celltype_qlf_d97_bmp_vs_fgf)[1]))
degs_d97_bmp_vs_fgf <- degs_d97_bmp_vs_fgf[order(-degs_d97_bmp_vs_fgf$logFC),]
#
write.csv(degs_d97_ctr_vs_bmp,
          file=paste0(data_directory, "WatanabeLab_ArealizationOrganoids/organoids_pbulk_edgeR_qlf_degs_D97_BMPvsCTR_", celltype_oi ,".csv"))
write.csv(degs_d97_ctr_vs_fgf,
          file=paste0(data_directory, "WatanabeLab_ArealizationOrganoids/organoids_pbulk_edgeR_qlf_degs_D97_FGFvsCTR_", celltype_oi,".csv"))
write.csv(degs_d97_bmp_vs_fgf,
          file=paste0(data_directory, "WatanabeLab_ArealizationOrganoids/organoids_pbulk_edgeR_qlf_degs_D97_FGFvsBMP_", celltype_oi,".csv"))

keyvals <- ifelse(
  degs_d97_ctr_vs_bmp$logFC < -0.5 & degs_d97_ctr_vs_bmp$FDR < 0.05, '#33A02C',
  ifelse(degs_d97_ctr_vs_bmp$FDR < 0.05 & degs_d97_ctr_vs_bmp$logFC > 0.5, '#1F78B4',
         'black'))
names(keyvals)[keyvals == '#33A02C'] <- 'CTR'
names(keyvals)[keyvals == '#1F78B4'] <- 'BMP'
volcano_ctr_vs_bmp <- EnhancedVolcano(degs_d97_ctr_vs_bmp,
                                      lab = rownames(degs_d97_ctr_vs_bmp),
                                      selectLab = rownames(degs_d97_ctr_vs_bmp)[which(names(keyvals) %in% c('BMP', 'CTR'))],
                                      x = 'logFC',
                                      y = 'PValue',
                                      title = 'BMP vs. CTR (IP, D97)',
                                      pCutoff = 0.05,
                                      FCcutoff = 0.5,
                                      cutoffLineType="blank",
                                      colCustom = keyvals,
                                      gridlines.major = FALSE,
                                      gridlines.minor = FALSE)
volcano_ctr_vs_bmp
ggsave(file=paste0(data_directory, "WatanabeLab_ArealizationOrganoids/", celltype_oi_label, "/D97/organoids_pbulk_allcells_edgeR_qlf_degs_D97_BMPvsCTR_", celltype_oi_label, "_volcano.pdf"))

keyvals <- ifelse(
  degs_d97_ctr_vs_fgf$logFC < -0.5 & degs_d97_ctr_vs_fgf$FDR < 0.05, '#33A02C',
  ifelse(degs_d97_ctr_vs_fgf$FDR < 0.05 & degs_d97_ctr_vs_fgf$logFC > 0.5, '#E31A1C',
         'black'))
names(keyvals)[keyvals == '#33A02C'] <- 'CTR'
names(keyvals)[keyvals == '#E31A1C'] <- 'FGF'
volcano_ctr_vs_fgf <- EnhancedVolcano(degs_d97_ctr_vs_fgf,
                                      lab = rownames(degs_d97_ctr_vs_fgf),
                                      selectLab = rownames(degs_d97_ctr_vs_fgf)[which(names(keyvals) %in% c('FGF', 'CTR'))],
                                      x = 'logFC',
                                      y = 'PValue',
                                      title = 'FGF vs. CTR (IP, D97)',
                                      pCutoff = 0.05,
                                      FCcutoff = 0.5,
                                      cutoffLineType="blank",
                                      colCustom = keyvals,
                                      gridlines.major = FALSE,
                                      gridlines.minor = FALSE)
volcano_ctr_vs_fgf
ggsave(file=paste0(data_directory, "WatanabeLab_ArealizationOrganoids/", celltype_oi_label, "/D97/organoids_pbulk_allcells_edgeR_qlf_degs_D97_FGFvsCTR_", celltype_oi_label, "_volcano.pdf"))

keyvals <- ifelse(
  degs_d97_bmp_vs_fgf$logFC < -0.5 & degs_d97_bmp_vs_fgf$FDR < 0.05, '#1F78B4',
  ifelse(degs_d97_bmp_vs_fgf$FDR < 0.05 & degs_d97_bmp_vs_fgf$logFC > 0.5, '#E31A1C',
         'black'))
names(keyvals)[keyvals == '#1F78B4'] <- 'BMP'
names(keyvals)[keyvals == '#E31A1C'] <- 'FGF'
volcano_fgf_vs_bmp <- EnhancedVolcano(degs_d97_bmp_vs_fgf,
                                      lab = rownames(degs_d97_bmp_vs_fgf),
                                      selectLab = rownames(degs_d97_bmp_vs_fgf)[which(names(keyvals) %in% c('FGF', 'BMP'))],
                                      x = 'logFC',
                                      y = 'PValue',
                                      title = 'FGF vs. BMP (IP, D97)',
                                      pCutoff = 0.05,
                                      FCcutoff = 0.5,
                                      cutoffLineType="blank",
                                      colCustom = keyvals,
                                      gridlines.major = FALSE,
                                      gridlines.minor = FALSE)
volcano_fgf_vs_bmp
ggsave(file=paste0(data_directory, "WatanabeLab_ArealizationOrganoids/", celltype_oi_label, "/D97/organoids_pbulk_allcells_edgeR_qlf_degs_D97_FGFvsBMP_", celltype_oi_label, "_volcano.pdf"))
