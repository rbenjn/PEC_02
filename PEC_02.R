## ----setup, include=FALSE---------------------
knitr::opts_chunk$set(echo = TRUE)


## ----4.1, include=FALSE-----------------------
setwd(".")
dir.create("data")
dir.create("results")


## ----include=FALSE----------------------------
library(knitr)
library(readr)
library(colorspace)
library(gplots)
library(ggplot2)
library(DESeq2)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(hugene21sttranscriptcluster.db)
library(clusterProfiler)
library(dplyr)
library(limma)
library(hexbin)
library(pheatmap)
library(RColorBrewer)


## ----include=FALSE----------------------------
counts <- read_delim("C:/Users/rubenuoc/Desktop/UOC/MO157/PEC2/PEC_02/data/counts.csv", 
                     ";", escape_double = FALSE, trim_ws = TRUE)


targets <- read_csv("C:/Users/rubenuoc/Desktop/UOC/MO157/PEC2/PEC_02/data/targets.csv")



## ----include=FALSE----------------------------
set.seed(280520)

targets_10_g <- lapply(split(targets, targets$Grupo_analisis),
                       function(subdf) subdf[sample(1:nrow(subdf), 10),])


targets_10_g <- do.call('rbind', targets_10_g)



## ----echo=FALSE-------------------------------
head(targets_10_g)


## ----include=FALSE----------------------------
df <- counts[, targets_10_g$Sample_Name]

# Para que podamos realizar con mayor falicidad el análisis posterior, cambiaremos el nombre de cada filas: 1, 2, 3, etc por el nombre de la columna X1 del archivo 'counts'.

# Añadimos la columna X1 a nuestro dataframe

df <- cbind(df, X1 = counts$X1)


# Cambiamos orden para tenerla primera y poder cambiar el

df <- df[, c(31, 1:30)]

countdata <- df[,-1]
rownames(countdata) <- df[,1]


## ----include=FALSE----------------------------
dds <- DESeqDataSetFromMatrix(countData = countdata, colData = targets_10_g, design = ~Group)



## ----echo=FALSE-------------------------------
nrow(dds)

dds <- dds[ rowSums(counts(dds)) > 1, ]
nrow(dds)


## ----echo=FALSE-------------------------------
vsd <- vst(dds, blind = FALSE) # blind = FALSE para evitar que las variables en el diseño del objeto (Group) interfieran en la varianza media esperada
head(assay(vsd), 3)


## ----echo=FALSE-------------------------------
colData(vsd)


## ----echo=FALSE-------------------------------
rld <- rlog(dds, blind = FALSE)
head(assay(rld), 3)


## ----include=FALSE----------------------------
dds <- estimateSizeFactors(dds)

df <- bind_rows(
  as_data_frame(log2(counts(dds, normalized=TRUE)[, 1:2]+1)) %>%
         mutate(transformation = "log2(x + 1)"),
  as_data_frame(assay(vsd)[, 1:2]) %>% mutate(transformation = "vst"),
  as_data_frame(assay(rld)[, 1:2]) %>% mutate(transformation = "rlog"))


## ----echo=FALSE-------------------------------
colnames(df)[1:2] <- c("x", "y")  

ggplot(df, aes(x = x, y = y)) + geom_hex(bins = 80) +
  coord_fixed() + facet_grid( . ~ transformation) 


## ----include=FALSE----------------------------
sampleDists <- dist(t(assay(vsd)))
sampleDists


## ----echo=FALSE-------------------------------

sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( vsd$dex, vsd$cell, sep = " - " )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors)


## ----echo=FALSE-------------------------------
plotPCA(vsd, intgroup = c("Group"))


## ----include=FALSE----------------------------
dds <- DESeq(dds)


## ----echo=FALSE-------------------------------
# Realizamos las distintas comparaciones usando 'constrast' para cada grupo, con su respectivo factor: 'Group' en primer lugar, y los 2 grupos a comparar.

# NIT vs ELI

res_NIT_vs_ELI <- results(dds, contrast = c("Group", "NIT", "ELI"))

NIT_vs_ELI_rownames <- row.names(res_NIT_vs_ELI)

summary(res_NIT_vs_ELI)



## ----echo=FALSE-------------------------------
# NIT vs SFI

res_NIT_vs_SFI <- results(dds, contrast = c("Group", "NIT", "SFI"))

NIT_vs_SFI_rownames <- row.names(res_NIT_vs_SFI)

summary(res_NIT_vs_SFI)


## ----echo=FALSE-------------------------------
# ELI vs SFI

res_ELI_vs_SFI <- results(dds, contrast = c("Group", "ELI", "SFI"))

ELI_vs_SFI_rownames <- row.names(res_ELI_vs_SFI)

summary(res_ELI_vs_SFI)


## ----include=FALSE----------------------------
sum(res_NIT_vs_ELI$padj < 0.1, na.rm=TRUE)

sum(res_NIT_vs_SFI$padj < 0.1, na.rm=TRUE)

sum(res_ELI_vs_SFI$padj < 0.1, na.rm=TRUE)


## ----include=FALSE----------------------------
res_NIT_vs_ELI_Sig <- subset(res_NIT_vs_ELI, padj < 0.1)

res_NIT_vs_SFI_Sig <- subset(res_NIT_vs_SFI, padj < 0.1)

res_ELI_vs_SFI_Sig <- subset(res_ELI_vs_SFI, padj < 0.1)


## ----echo=FALSE-------------------------------
head(res_NIT_vs_ELI_Sig[ order(res_NIT_vs_ELI_Sig$log2FoldChange), ])

head(res_NIT_vs_SFI_Sig[ order(res_NIT_vs_SFI_Sig$log2FoldChange), ])

head(res_ELI_vs_SFI_Sig[ order(res_ELI_vs_SFI_Sig$log2FoldChange), ])


## ----echo=FALSE-------------------------------
head(res_NIT_vs_ELI_Sig[ order(res_NIT_vs_ELI_Sig$log2FoldChange, decreasing = TRUE), ])

head(res_NIT_vs_SFI_Sig[ order(res_NIT_vs_SFI_Sig$log2FoldChange, decreasing = TRUE), ])

head(res_ELI_vs_SFI_Sig[ order(res_ELI_vs_SFI_Sig$log2FoldChange, decreasing = TRUE), ])


## ----include=FALSE----------------------------
columns(org.Hs.eg.db)


## ----echo=FALSE-------------------------------
# NIT vs ELI Sig

tmp_NIT_vs_ELI_Sig=gsub("\\..*","",row.names(res_NIT_vs_ELI_Sig)) # borrar '.' de ENSEMBL

res_NIT_vs_ELI_Sig$symbol <- mapIds(org.Hs.eg.db,
                     keys=tmp_NIT_vs_ELI_Sig,
                     column="SYMBOL",
                     keytype="ENSEMBL",
                     multiVals="first")

res_NIT_vs_ELI_Sig$entrez <- mapIds(org.Hs.eg.db,
                     keys=tmp_NIT_vs_ELI_Sig,
                     column="ENTREZID",
                     keytype="ENSEMBL",
                     multiVals="first")

res_NIT_vs_ELI_Sig_Ordered <- res_NIT_vs_ELI_Sig[order(res_NIT_vs_ELI_Sig$pvalue),]

head(res_NIT_vs_ELI_Sig_Ordered)

# NIT vs SFI Sig

tmp_NIT_vs_SFI_Sig=gsub("\\..*","",row.names(res_NIT_vs_SFI_Sig)) # borrar '.' de ENSEMBL

res_NIT_vs_SFI_Sig$symbol <- mapIds(org.Hs.eg.db,
                     keys=tmp_NIT_vs_SFI_Sig,
                     column="SYMBOL",
                     keytype="ENSEMBL",
                     multiVals="first")

res_NIT_vs_SFI_Sig$entrez <- mapIds(org.Hs.eg.db,
                     keys=tmp_NIT_vs_SFI_Sig,
                     column="ENTREZID",
                     keytype="ENSEMBL",
                     multiVals="first")

res_NIT_vs_SFI_Sig_Ordered <- res_NIT_vs_SFI_Sig[order(res_NIT_vs_SFI_Sig$pvalue),]

head(res_NIT_vs_SFI_Sig_Ordered)

# ELI vs SFI Sig

tmp_ELI_vs_SFI_Sig=gsub("\\..*","",row.names(res_ELI_vs_SFI_Sig)) # borrar '.' de ENSEMBL

res_ELI_vs_SFI_Sig$symbol <- mapIds(org.Hs.eg.db,
                     keys=tmp_ELI_vs_SFI_Sig,
                     column="SYMBOL",
                     keytype="ENSEMBL",
                     multiVals="first")

res_ELI_vs_SFI_Sig$entrez <- mapIds(org.Hs.eg.db,
                     keys=tmp_ELI_vs_SFI_Sig,
                     column="ENTREZID",
                     keytype="ENSEMBL",
                     multiVals="first")

res_ELI_vs_SFI_Sig_Ordered <- res_ELI_vs_SFI_Sig[order(res_ELI_vs_SFI_Sig$pvalue),]

head(res_ELI_vs_SFI_Sig_Ordered)


## ----echo=FALSE-------------------------------
# Para obtener las keys

# NIT vs ELI

tmp_NIT_vs_ELI=gsub("\\..*","",row.names(res_NIT_vs_ELI)) # borrar '.' de ENSEMBL

res_NIT_vs_ELI$symbol <- mapIds(org.Hs.eg.db,
                     keys=tmp_NIT_vs_ELI,
                     column="SYMBOL",
                     keytype="ENSEMBL",
                     multiVals="first")

res_NIT_vs_ELI$entrez <- mapIds(org.Hs.eg.db,
                     keys=tmp_NIT_vs_ELI,
                     column="ENTREZID",
                     keytype="ENSEMBL",
                     multiVals="first")

res_NIT_vs_ELI_Ordered <- res_NIT_vs_ELI[order(res_NIT_vs_ELI$pvalue),]

head(res_NIT_vs_ELI_Ordered)

# NIT vs SFI

res_NIT_vs_SFI_DF <- as.data.frame(res_NIT_vs_SFI)

tmp_NIT_vs_SFI=gsub("\\..*","",row.names(res_NIT_vs_SFI)) # borrar '.' de ENSEMBL

res_NIT_vs_SFI$symbol <- mapIds(org.Hs.eg.db,
                     keys=tmp_NIT_vs_SFI,
                     column="SYMBOL",
                     keytype="ENSEMBL",
                     multiVals="first")

res_NIT_vs_SFI$entrez <- mapIds(org.Hs.eg.db,
                     keys=tmp_NIT_vs_SFI,
                     column="ENTREZID",
                     keytype="ENSEMBL",
                     multiVals="first")

res_NIT_vs_SFI_Ordered <- res_NIT_vs_SFI[order(res_NIT_vs_SFI$pvalue),]

head(res_NIT_vs_SFI_Ordered)

# ELI vs SFI

res_ELI_vs_SFI_DF <- as.data.frame(res_ELI_vs_SFI)

tmp_ELI_vs_SFI=gsub("\\..*","",row.names(res_ELI_vs_SFI_DF)) # borrar '.' de ENSEMBL

res_ELI_vs_SFI$symbol <- mapIds(org.Hs.eg.db,
                     keys=tmp_ELI_vs_SFI,
                     column="SYMBOL",
                     keytype="ENSEMBL",
                     multiVals="first")

res_ELI_vs_SFI$entrez <- mapIds(org.Hs.eg.db,
                     keys=tmp_ELI_vs_SFI,
                     column="ENTREZID",
                     keytype="ENSEMBL",
                     multiVals="first")

res_ELI_vs_SFI_Ordered <- res_ELI_vs_SFI[order(res_ELI_vs_SFI$pvalue),]

head(res_ELI_vs_SFI_Ordered)



## ----Significativos, include=FALSE------------
# Exportamos los resultados significativos

# NIT_vs_ELI

res_NIT_vs_ELI_Sig_Ordered_DF <- as.data.frame(res_NIT_vs_ELI_Sig_Ordered)

write.csv(res_NIT_vs_ELI_Sig_Ordered, file = "./results/NIT_vs_ELI_Sig.csv")

# NIT_vs_SFI

res_NIT_vs_SFI_Sig_Ordered_DF <- as.data.frame(res_NIT_vs_SFI_Sig_Ordered)

write.csv(res_NIT_vs_SFI_Sig_Ordered, file = "./results/NIT_vs_SFI_Sig.csv")

# ELI_vs_SFI

res_ELI_vs_SFI_Sig_Ordered_DF <- as.data.frame(res_ELI_vs_SFI_Sig_Ordered)

write.csv(res_ELI_vs_SFI_Sig_Ordered, file = "./results/ELI_vs_SFI_Sig.csv")


## ----include=FALSE----------------------------
# Exportamos los resultados

# NIT_vs_ELI

res_NIT_vs_ELI_Ordered_DF <- as.data.frame(res_NIT_vs_ELI_Ordered)

write.csv(res_NIT_vs_ELI_Ordered, file = "./results/NIT_vs_ELI.csv")

# NIT_vs_SFI

res_NIT_vs_SFI_Ordered_DF <- as.data.frame(res_NIT_vs_SFI_Ordered)

write.csv(res_NIT_vs_SFI_Ordered, file = "./results/NIT_vs_SFI.csv")

# ELI_vs_SFI

res_ELI_vs_SFI_Ordered_DF <- as.data.frame(res_ELI_vs_SFI_Ordered)

write.csv(res_ELI_vs_SFI_Ordered, file = "./results/ELI_vs_SFI.csv")



## ----echo=FALSE-------------------------------
# Combinamos las 3 comparaciones para el posterior Diagrama de Venn

comb <- c(rownames(res_NIT_vs_ELI_Sig_Ordered_DF), rownames(res_NIT_vs_ELI_Sig_Ordered_DF), rownames(res_ELI_vs_SFI_Sig_Ordered_DF))

res_NIT_vs_ELI_venn <- comb %in% rownames(res_NIT_vs_ELI_Sig_Ordered_DF)
res_NIT_vs_SFI_venn <- comb %in% rownames(res_NIT_vs_ELI_Sig_Ordered_DF)  
res_ELI_vs_SFI_venn <- comb %in% rownames(res_ELI_vs_SFI_Sig_Ordered_DF)

counts_venn <- cbind(res_NIT_vs_ELI_venn, res_NIT_vs_SFI_venn, res_ELI_vs_SFI_venn)

results_venn <- vennCounts(counts_venn)

vennDiagram (results_venn, cex=1, names = c("NIT vs ELI", "NIT vs SFI", "ELI vs SFI"), circle.col = c("red", "blue", "yellow"))

title("Genes en común en las 3 comparaciones")


## ----echo=FALSE-------------------------------
# Con el paquete 'ClusterProfiler'

tmp_NIT_vs_ELI=gsub("\\..*","",row.names(res_NIT_vs_ELI_Ordered_DF)) # borrar '.' de ENSEMBL

ego_NIT_vs_ELI <- enrichGO(gene         = res_NIT_vs_ELI_Ordered_DF$entrez,
                OrgDb         = org.Hs.eg.db,
                keyType       = 'ENTREZID',
                ont           = "CC",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05)
head(summary(ego_NIT_vs_ELI))

dotplot(ego_NIT_vs_ELI, showCategory = 15)


## ----echo=FALSE-------------------------------
ego_NIT_vs_ELI_Sig <- enrichGO(gene         = res_NIT_vs_ELI_Sig_Ordered_DF$entrez,
                OrgDb         = org.Hs.eg.db,
                keyType       = 'ENTREZID',
                ont           = "CC",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05)
head(summary(ego_NIT_vs_ELI_Sig))

dotplot(ego_NIT_vs_ELI_Sig, showCategory = 15)


## ----echo=FALSE-------------------------------
tmp_NIT_vs_SFI=gsub("\\..*","",row.names(res_NIT_vs_SFI_Ordered_DF)) # borrar '.' de ENSEMBL

ego_NIT_vs_SFI <- enrichGO(gene         = res_NIT_vs_SFI_Ordered_DF$entrez,
                OrgDb         = org.Hs.eg.db,
                keyType       = 'ENTREZID',
                ont           = "CC",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05)
head(summary(ego_NIT_vs_SFI))

dotplot(ego_NIT_vs_SFI, showCategory = 15)


## ----echo=FALSE-------------------------------
ego_NIT_vs_SFI_Sig <- enrichGO(gene         = res_NIT_vs_SFI_Sig_Ordered_DF$entrez,
                OrgDb         = org.Hs.eg.db,
                keyType       = 'ENTREZID',
                ont           = "CC",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05)
head(summary(ego_NIT_vs_SFI_Sig))

dotplot(ego_NIT_vs_SFI_Sig, showCategory = 15)


## ----echo=FALSE-------------------------------
# En este caso utilizaremos como keyType 'ENSEMBL' en lugar de ENTREZID, para ver si hay diferencias.

tmp_ELI_vs_SFI=gsub("\\..*","",row.names(res_ELI_vs_SFI_Ordered_DF)) # borrar '.' de ENSEMBL

ego_ELI_vs_SFI <- enrichGO(gene         = tmp_ELI_vs_SFI,
                OrgDb         = org.Hs.eg.db,
                keyType       = 'ENSEMBL',
                ont           = "CC",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05)
head(summary(ego_ELI_vs_SFI))

dotplot(ego_ELI_vs_SFI, showCategory = 15)


## ----echo=FALSE-------------------------------
tmp_ELI_vs_SFI_Sig=gsub("\\..*","",row.names(res_ELI_vs_SFI_Sig_Ordered_DF)) # borrar '.' de ENSEMBL

ego_ELI_vs_SFI_Sig <- enrichGO(gene         = tmp_ELI_vs_SFI_Sig,
                OrgDb         = org.Hs.eg.db,
                keyType       = 'ENSEMBL',
                ont           = "CC",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05)
head(summary(ego_ELI_vs_SFI_Sig))

dotplot(ego_ELI_vs_SFI_Sig, showCategory = 15)


## ----echo=FALSE-------------------------------

listOfTables <- list(NIT_vs_ELI = res_NIT_vs_ELI_Sig_Ordered_DF,                    NIT_vs_SFI = res_NIT_vs_SFI_Sig_Ordered_DF, ELI_vs_SFI = res_ELI_vs_SFI_Sig_Ordered_DF)

listOfSelected <- list()

for (i in 1:length(listOfTables)){
  
  topTab <- listOfTables[[i]]
  
  whichGenes<-topTab$padj<0.1
  selectedIDs <- gsub("\\..*","",row.names(topTab)[whichGenes])
  EntrezIDs <- topTab$entrez
  listOfSelected[[i]] <- EntrezIDs
  names(listOfSelected)[i] <- names(listOfTables)[i]
}

sapply(listOfSelected, length)



## ----include=FALSE----------------------------

mapped_genes2GO <- mappedkeys(org.Hs.egGO)
mapped_genes2KEGG <- mappedkeys(org.Hs.egPATH)
mapped_genes <- union(mapped_genes2GO , mapped_genes2KEGG)



## ----echo=FALSE-------------------------------

library(clusterProfiler)


listOfData <- listOfSelected[1:3]
comparisonsNames <- names(listOfData)
universe <- mapped_genes

for (i in 1:length(listOfData)){
  genesIn <- listOfData[[i]]
  comparison <- comparisonsNames[i]
  enrich.result <- enrichKEGG(gene = genesIn,
                                 pvalueCutoff = 0.10,
                                 pAdjustMethod = "BH",
                                 organism = "human",
                                 universe = universe)
  cat("##################################")
  cat("\nComparison: ", comparison,"\n")
  print(head(enrich.result))
  
  if (length(rownames(enrich.result@result)) != 0) {
    write.csv(as.data.frame(enrich.result), 
              file =paste0("./results/","cluster.Results.",comparison,".csv"), 
              row.names = FALSE)
    
    pdf(file=paste0("./results/","clusterBarplot.",comparison,".pdf"))
    print(barplot(enrich.result, showCategory = 15, font.size = 4, 
                  title = paste0("Cluster Pathway Analysis for ", comparison,". Barplot")))
    dev.off()
    
    pdf(file = paste0("./results/","clustercnetplot.",comparison,".pdf"))
    print(cnetplot(enrich.result, categorySize = "geneNum", schowCategory = 15, 
                   vertex.label.cex = 0.75))
    dev.off()
  } 
}



## ----echo=FALSE-------------------------------
tabla_NIT_vs_ELI <- read.csv2("./results/cluster.Results.NIT_vs_ELI.csv", header = TRUE, sep = ",") 

knitr::kable(head(tabla_NIT_vs_ELI), booktabs = TRUE, caption = 'Tabla NIT_vs_ELI')


## ----echo=FALSE-------------------------------
tabla_NIT_vs_SFI <- read.csv2("./results/cluster.Results.NIT_vs_SFI.csv", header = TRUE, sep = ",") 

knitr::kable(head(tabla_NIT_vs_SFI), booktabs = TRUE, caption = 'Tabla NIT_vs_SFI')


## ----echo=FALSE-------------------------------
tabla_ELI_vs_SFI <- read.csv2("./results/cluster.Results.ELI_vs_SFI.csv", header = TRUE, sep = ",") 

knitr::kable(head(tabla_ELI_vs_SFI), booktabs = TRUE, caption = 'Tabla ELI_vs_SFI')


## ----echo=FALSE-------------------------------
knitr::include_graphics("./results/barplot_NIT_vs_ELI.png")

knitr::include_graphics("./results/cnetplot_NIT_vs_ELI.png")



## ----echo=FALSE-------------------------------
knitr::include_graphics("./results/barplot_NIT_vs_SFI.png")

knitr::include_graphics("./results/cnetplot_NIT_vs_SFI.png")



## ----echo=FALSE-------------------------------
knitr::include_graphics("./results/barplot_ELI_vs_SFI.png")

knitr::include_graphics("./results/cnetplot_ELI_vs_SFI.png")


