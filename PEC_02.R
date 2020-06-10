# PEC_02

set.seed(280520)

library(readr)
counts <- read_delim("C:/Users/rubenuoc/Desktop/UOC/MO157/PEC2/archivos targets y counts/counts.csv", 
                     ";", escape_double = FALSE, trim_ws = TRUE)


targets <- read_csv("C:/Users/rubenuoc/Desktop/UOC/MO157/PEC2/archivos targets y counts/targets.csv")


targets_10_g <- lapply(split(targets, targets$Grupo_analisis),
                       function(subdf) subdf[sample(1:nrow(subdf), 10),])


targets_10_g <- do.call('rbind', targets_10_g)

df <- counts[, targets_10_g$Sample_Name]

# No los añadimos ensembles

df <- cbind(df, X1 = counts$X1)


# Cambiamos orden para tenerla primera

df <- df[, c(31, 1:30)]

countdata <- df[,-1]
rownames(countdata) <- df[,1]

# Usamos el paquete DESeq2

library(DESeq2)

# Read counting steps

library("GenomicAlignments")

library("BiocParallel")

register(SerialParam())

# Pasos buenos: 3.5.2

dds <- DESeqDataSetFromMatrix(countData = countdata, colData = targets_10_g, design = ~Group)


# 3.6 Exploratory analysis

nrow(dds)

dds <- dds[ rowSums(counts(dds)) > 1, ]
nrow(dds)

# 3.6.2

vsd <- vst(dds, blind = FALSE)
head(assay(vsd), 3)

colData(vsd)

rld <- rlog(dds, blind = FALSE)
head(assay(rld), 3)


library("dplyr")
library("ggplot2")

dds <- estimateSizeFactors(dds)

df <- bind_rows(
  as_data_frame(log2(counts(dds, normalized=TRUE)[, 1:2]+1)) %>%
    mutate(transformation = "log2(x + 1)"),
  as_data_frame(assay(vsd)[, 1:2]) %>% mutate(transformation = "vst"),
  as_data_frame(assay(rld)[, 1:2]) %>% mutate(transformation = "rlog"))


colnames(df)[1:2] <- c("x", "y")  

ggplot(df, aes(x = x, y = y)) + geom_hex(bins = 80) +
  coord_fixed() + facet_grid( . ~ transformation) 

#3.6.3

sampleDists <- dist(t(assay(vsd)))
sampleDists

# Heatmap

library("pheatmap")
library("RColorBrewer")

sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( vsd$dex, vsd$cell, sep = " - " )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors)

# 3.6.4 PCA plot

plotPCA(vsd, intgroup = c("Group"))

# 3.6.5 MDS Plot (no sale)

library(ggplot2)

mds <- as.data.frame(colData(vsd)) 
cbind(cmdscale(sampleDistMatrix))

ggplot(mds, aes(x = "1", y = "2", color = Group)) +
  geom_point(size = 3) + coord_fixed()

# 3.7. Differential Expression Analysis

dds <- DESeq(dds, parallel =TRUE)

# 1a Comparación

res <- results(dds)
res

res <- results(dds, contrast = c("Group","ELI","NIT"))

res


mcols(res, use.names = TRUE)


summary(res)


res.05 <- results(dds, alpha = 0.05)
table(res.05$padj < 0.05)


resLFC1 <- results(dds, lfcThreshold=1)
table(resLFC1$padj < 0.1)


sum(res$padj < 0.1, na.rm=TRUE)


resSig <- subset(res, padj < 0.1)
head(resSig[ order(resSig$log2FoldChange), ])


head(resSig[ order(resSig$log2FoldChange, decreasing = TRUE), ])


# 3.8 Plotting Results

# Count plot

topGene <- rownames(res)[which.min(res$padj)]


library("ggbeeswarm")

geneCounts <- plotCounts(dds, gene = topGene, intgroup = c("Group"),
                         returnData = TRUE)
ggplot(geneCounts, aes(x = Group, y = count, color = Group, group = Group)) +
  scale_y_log10() + geom_point(size = 3) + geom_line()


# MA Plot


library("apeglm")
resultsNames(dds) # Mirar que salgan los 3

# Canvi results para cada uno

res <- lfcShrink(dds, coef="Group_NIT_vs_ELI", type="apeglm")

plotMA(res, ylim = c(-5, 5))

plotMA(res, ylim = c(-5,5))
topGene <- rownames(res)[which.min(res$padj)]
with(res[topGene, ], {
  points(baseMean, log2FoldChange, col="dodgerblue", cex=2, lwd=2)
  text(baseMean, log2FoldChange, topGene, pos=2, col="dodgerblue")
})

#3.8.3 Gene Clustering (problema)

library("genefilter")
topVarGenes <- head(order(rowVars(assay(vsd)), decreasing = TRUE), 20)

mat  <- assay(vsd)[topVarGenes, ]
mat  <- mat - rowMeans(mat)
anno <- as.data.frame(colData(vsd)[, c("Group")])
rownames(anno) <- colnames(mat) #Añadimos el nombre de las filas al data frame

library("pheatmap")

pheatmap(mat, annotation_col = anno)

# 3.8.4 (summarizeOverLaps)

resGR <- results(dds, name="Group_NIT_vs_ELI", format="GRangesList") # Cambiado a GRangesList


resGR$log2FoldChange <- res$log2FoldChange
resGR



library("org.Hs.eg.db")
resGR$symbol <- mapIds(org.Hs.eg.db, names(resGR), "SYMBOL", "ENSEMBL")


library("Gviz")

window <- resGR[topGene] + 1e6
strand(window) <- "*"
resGRsub <- resGR[resGR %over% window]
naOrDup <- is.na(resGRsub$symbol) | duplicated(resGRsub$symbol)
resGRsub$group <- ifelse(naOrDup, names(resGRsub), resGRsub$symbol)


status <- factor(ifelse(resGRsub$padj < 0.1 & !is.na(resGRsub$padj), "sig", "notsig"))


options(ucscChromosomeNames = FALSE)
g <- GenomeAxisTrack()
a <- AnnotationTrack(resGRsub, name = "gene ranges", feature = status)
d <- DataTrack(resGRsub, data = "log2FoldChange", baseline = 0,
               type = "h", name = "log2 fold change", strand = "+")
plotTracks(list(g, d, a), groupAnnotation = "group",
           notsig = "grey", sig = "hotpink")

# 3.9 Anotación

library("AnnotationDbi")
columns(org.Hs.eg.db)

# Para obtener las keys

tmp=gsub("\\..*","",row.names(res)) # borrar . de ENSEMBL

res$symbol <- mapIds(org.Hs.eg.db,
                     keys=tmp,
                     column="SYMBOL",
                     keytype="ENSEMBL",
                     multiVals="first")

res$entrez <- mapIds(org.Hs.eg.db,
                     keys=tmp,
                     column="ENTREZID",
                     keytype="ENSEMBL",
                     multiVals="first")

resOrdered <- res[order(res$pvalue),]
head(resOrdered)
