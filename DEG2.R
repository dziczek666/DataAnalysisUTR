# load libraries
library(DESeq2)
library(ggplot2)
library(dplyr)
# read .csv file with counts
counts <- read.delim("../counts_cancer.txt", header=TRUE, row.names = 1)
counts
counts <- counts[which(rowSums(counts) > 50),]
counts <- select(counts, -Chr, -Start, -End, -Length, -Strand)
counts <- select(counts, SRR19969208.bam, SRR19969210.bam, SRR19969211.bam, SRR19969212.bam, SRR19969216.bam, SRR19969217.bam, SRR19969218.bam)
head(counts)

# set conditions
condition <- factor(c("test", "test", "test", "test", "control", "control", "control"))

# create df with metadata
coldata <- data.frame(row.names = colnames(counts), condition)

dds <- DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = ~condition)

# run DESeq2
dds <- DESeq(dds)
dds

vsdata <- vst(dds, blind=FALSE)
vsdata
# PCA
PCA <- plotPCA(vsdata, intgroup = "condition")

png("PCA_trip.png",res = 100)
print(PCA)
dev.off()

# dispersion
dispersion <- plotDispEsts(dds)

png("disp_trip.png", res=300, width = 1000, height = 2000)
print(dispersion)
dev.off()

# create df with DESeq2 results
res <- results(dds)
res

# omit NA
sigs <- na.omit(res)

# filter padj < 0.05
sigs <- sigs[sigs$padj < 0.05,]


write.csv(res, "all_trip_results_res.csv")

# create df with filtered baseMean and log2FC
df1 <- as.data.frame(sigs)
df1.top <- df1[(df1$baseMean > 50) & (abs(df1$log2FoldChange) > 1),]
df1.top <- read.csv('df1.nowytop.csv', header=TRUE, row.names = 1)
df1.top
write.csv(df1.top, 'df1.top.utr_trip.csv')
# sort from the highest to lowest log2FC
df1.top <- df1.top[order(df1.top$log2FoldChange, decreasing = TRUE),]
df1.top
# Get normalized count data
rlog_out <- rlog(dds, blind = FALSE) 
mat <- assay(rlog_out)[rownames(df1.top), rownames(coldata)]
colnames(mat) <- rownames(coldata)
base_mean <- rowMeans(mat)
mat.scaled <- t(apply(mat, 1, scale))
colnames(mat.scaled) <- colnames(mat)
mat

mat.scaled
# keep 50 most significant genes
num_keep <- 25
rows_keep <- c(seq(1:num_keep), seq((nrow(mat.scaled)-num_keep), nrow(mat.scaled)))

# create matrix with log2FC for kept genes
l2_val <- as.matrix(df1.top[rows_keep,]$log2FoldChange)
colnames(l2_val) <- "logFC"

#create matrix with mean for kept genes
mean <- as.matrix(df1.top[rows_keep,]$baseMean)
colnames(mean) <- "avrExpr"

#load libraries
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)

# Map values between colors for min and max l2 _alues
col_logFC <- colorRamp2(c(min(l2_val),0, max (l2_val)), c("blue", "red", "white"))

# Maps between 0% quantile and 75% quantile of mean values
col_avrExpr <- colorRamp2(c(quantile(mean)[1], quantile(mean)[4]), c("white", "red"))

# Set summary for heatmap plot
ha <- HeatmapAnnotation(summary = anno_summary(gp = gpar(fill = 2),
                                               height = unit(2, "cm")))

h1 <- Heatmap(mat.scaled[rows_keep,], cluster_rows = F,
              column_labels = colnames(mat.scaled), name = "Z-score",
              cluster_columns = T)
df1.top_rownames <- rownames(df1.top)

h2 <- Heatmap(l2_val, row_labels = df1.top_rownames[rows_keep],
              cluster_rows = F, name = "logFC", top_annotation = ha, col = col_logFC,
              cell_fun = function(j, i, x, y, w, h, col) {# add text to each grid
                grid.text(round(l2_val[i, j],2), x, y)
              })

h3 <- Heatmap(mean, row_labels = df1.top_rownames[rows_keep],
              cluster_rows = F, name = "avrExpr", col = col_avrExpr,
              cell_fun = function(j, i, x, y, w, h, col) {# add text to each grid
                grid.text(round(mean[i,j],2), x, y)
              })

h <- h1+h2+h3
h

png("heatmap_utr_trip.png", res=300, width = 3000, height = 5500)
print(h)
dev.off()

# load library
library(EnhancedVolcano)

# Create volcano plot
volcano <- EnhancedVolcano(res, x = "log2FoldChange", y = "padj", lab = row.names(res))

png("volcano_trip_UTR.png", res=300, width = 3000, height = 3500)
print(volcano)
dev.off()
volcano
volcano2 <- EnhancedVolcano(res, x = "log2FoldChange", y = "padj", lab = row.names(res),
                            pCutoff = 1e-4, FCcutoff = 1)
volcano2

png("volcano2_UTR_her.png", res=300, width = 4000, height = 4500)
print(volcano2)
dev.off()
