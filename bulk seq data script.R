library(DESeq2)
library(pheatmap)
library(dplyr)
library(RColorBrewer)
library(ggplot2)
library(ggrepel)
#setworkingditrectory
#setwd("Desktop")

getwd()
#load count data

count_data <- read.csv('', header= TRUE,  sep = ",",row.names = 1)
head(count_data)
class(count_data)


#no. of gene ids inside the data
nrow(count_data)


#get number of duplicated geneid
sum(duplicated(count_data))

#get NA position in data
which(is.na(count_data))


#performing data fram samples 

#- to perform a DESEq2 analysis a mapping of the columns of the count data to a corresponding group or condition is needed

#- deseq2 uses a data frame where rownames of data frame needs to be the column names of count data 
#- the colmns of new data frame provide the mapping of samples from specific group 
# creation of data frame for the gene names

genenames <- count_data$
  
  head(genenames)



# teh data fails to satisfy the criteria for a poisson distribution and typical RNA Seq data will do the same
#thats why we will go for negative bionomial distribution


#load sample data
coldata <- read.csv("coldata.csv", header = TRUE, row.names = 1)
colnames(coldata)
head(coldata)

#creation of count data

#count_data <- count_data[,2:9]
#View(count_data)

row.names(count_data) <- genenames

class(count_data)

count_data <- as.matrix(count_data)

# coldata dataframe



#pre filtering data
#while pre filtering is not required step, but additionally benifical for two reasons: by removing reads with very few reads, we reduce the memory of ds dataobject to speed up the transformation and testing function deseq2 additionally since feature without any data for difefrnatial expression are not plotted it can be improve visualization

#genes with total count >= 10 are kept for defferentail expression analysis

count_data <- count_data[rowSums(count_data) >= 10,]

head(count_data)
View(count_data)



dds <- DESeqDataSetFromMatrix(countData = count_data, colData = coldata, design = ~Group)

dds
#dds$Group <- relevel(dds$Group, ref = "")

dds <- DESeq(dds)


normalized_counts <- counts(dds, normalized = TRUE)
normalized_counts <- as.matrix(normalized_counts)
head(normalized_counts)
write.csv(normalized_counts, "normalized_count 22.csv")
head(count_data)
res <- results((dds), alpha = 0.05)
res <- results(dds, contrast = c("Group", "Compare Grp", "Ref Grp"), alpha = 0.05)


write.csv(res3, "")

summary(res)


#plotma

plotMA(res, cex = 0.7, ylim= c(-10,10))
abline(h=c(-1,1), col = "red", lwd=3)



resultsNames(res)


#shrinkage of effect size (LFC estimate) is useful for visualization and ranking of genes, to shrink the LFC we pass the dds objectto function lfcshrink, below we specify to use apeglm method for effect size shrinkage  which improve previsios estimator
library(apeglm)
resLFC <- lfcShrink(dds, coef = "condition_RIVA_vs_DM", type = "apeglm")
plotMA(resLFC, cex = 0.7, ylim =c(-10,10))
abline(h=c(-1,1), col="red", lwd=3)

BiocManager::install("apeglm")

#dispersion plot 

plotDispEsts(dds, main= "Dispersion plot")




#PCA plot

rld <- rlogTransformation(dds, blind = FALSE)

head(assay(rld))

hist(assay(rld))

PCAA <- plotPCA(rld, intgroup = "Group")

PCAA + geom_text(aes(label = name), size = 4) + ggtitle("PCA Plot")



#volcano plot
library(EnhancedVolcano)
EnhancedVolcano(res, 
                lab= rownames(res),
                x= "log2FoldChange",
                y= "padj",
                title = "",
                pCutoff = 0.05, 
                FCcutoff = 0.5,
                pointSize = 2,
                labSize = 4,
                selectLab = c(),)
EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'pvalue',
                selectLab = c(),
                title = "",
                pCutoff = 0.05,
                FCcutoff = 0.5,
                pointSize = 2.0,
                labSize = 4.0,
                labCol = 'black',
                labFace = 'bold',
                boxedLabels = FALSE,
                subtitle = "DEGs",
                colAlpha = 4/5,
                legendPosition = 'right',
                legendLabSize = 13,
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 0.5,
                colConnectors = 'black')
EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'pvalue',
                pCutoff = 0.05,
                FCcutoff = 0.5,
                pointSize = 1.5,
                labSize = 3,
                selectLab = c(),
                title = '',
                subtitle = "DEGs",
                #caption = 'FC cutoff, 1.333; p-value cutoff, 10e-4',
                legendPosition = "right",
                legendLabSize = 13,
                col = c('grey30', 'forestgreen', 'royalblue', 'red2'),
                colAlpha = 0.9,
                drawConnectors = FALSE,
                hline = c(10e-8),
                widthConnectors = 0.5)
# Gather Log-fold change and FDR-corrected pvalues from DESeq2 results
## - Change pvalues to -log10 (1.3 = 0.05)
data <- data.frame(gene = row.names(res1),
                   pval = -log10(res1$padj), 
                   lfc = res1$log2FoldChange)

# Remove any rows that have NA as an entry
data <- na.omit(data)


BiocManager::install("EnhancedVolcano")
lib
#distance matrix  based on normalized count

sampleDists <- dist(t(assay(rld)))

library(RColorBrewer)
sampleDistMatrix <- as.matrix(sampleDists)


colnames(sampleDistMatrix)

colors <- colorRampPalette(rev(brewer.pal(9,"Blues")))(255)
Summary(res1)

#generate hetamap
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDistMatrix,
         clustering_distance_cols = sampleDistMatrix, color = colors)
#selected genes
sel_genes <- c()


#heatmap of log transformed normalized counts

top_genes <- res3[order(res3$padj), ][1:100,]


class(top_genes)

top_genes <- row.names(top_genes)
top_genes
summary(res2)
pheatmap(assay(rld)[top_genes,], cluster_rows = FALSE, show_rownames = TRUE, cluster_cols = FALSE, fontsize = 8,
         annotation_col = coldata)



#heatmap fo z score for top 20
#a z score of zero indicates that the gene expression level is the same as the mean expression level acresso
#all samples while possitive z score indicaes that the gene is expressed at higher level than the mean, and negative z score means that gene is expressed at lower level than the mean


cal_z_score <- function(x){(x-mean(x))/sd(x)}



z_score_all <- t(apply(normalized_counts, 1, cal_z_score))

gene <- "Nnt"

z_score_subset <- z_score_all[top_genes,]

pheatmap(z_score_subset, cluster_rows = TRUE, show_rownames = TRUE,show_colnames = TRUE,
         cluster_cols = FALSE, fontsize = 7.5, annotation_col = coldata)


z_score_subset <- z_score_all[rownames(z_score_all) %in% sel_genes,]


plotCounts(dds, "", intgroup = "condition")
#load sample data
sample_info <- read.csv("", header = TRUE, row.names = 1)
colnames(sample_info)
head(sample_info)
res


resSig <- subset(res3, padj < 0.05)
resSig
View(resSig)

resOrdered <- resSig[order(resSig$log2FoldChange),]

top_genes <- resOrdered[order(resOrdered$log2FoldChange, decreasing = TRUE), ][1:50,]

res <- readRDS("")
library(EnhancedVolcano)
EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'pvalue',
                selectLab = c("STS","CGB","PlGF","FLT-1","TNF","IL-1B","CCL5"),
                xlab = bquote(~Log[2]~ 'fold change'),
                pCutoff = 10e-14,
                FCcutoff = 2.0,
                pointSize = 4.0,
                labSize = 6.0,
                labCol = 'black',
                labFace = 'bold',
                boxedLabels = TRUE,
                colAlpha = 4/5,
                legendPosition = 'right',
                legendLabSize = 14,
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 1.0,
                colConnectors = 'black')




