BiocManager::install("clusterProfiler")
BiocManager::install("pathview")
BiocManager::install("enrichplot")
library(clusterProfiler)
library(enrichplot)
# we use ggplot2 to add x axis labels (ex: ridgeplot)
library(ggplot2)

# SET THE DESIRED ORGANISM HERE
organism = "org.Hs.eg.db"
BiocManager::install(organism, character.only = TRUE)
library(organism, character.only = TRUE)

# reading in data from deseq2
df = read.table("uniqA.txt", header=TRUE)
# we want the log2 fold change 
original_gene_list <- df$log2FoldChange

name <- df$reference_id
# name
library(tidyr)
library(gprofiler2)
converted_names <- gconvert(query = name, organism = "hsapiens",
         target="ENSG", mthreshold = Inf, filter_na = TRUE)
converted_names
ens_name <- converted_names[,'name']

# name the vector
names(original_gene_list) <- ens_name
original_gene_list <- ens_name
# omit any NA values 
gene_list<-na.omit(original_gene_list)
gene_list
# sort the list in decreasing order (required for clusterProfiler)
gene_list = sort(gene_list, decreasing = TRUE)
gene_list

enrichRes <- enrichGO(gene          = gene_list,
                      OrgDb         = organism,
                      keyType = "SYMBOL",
                      ont           = "BP",  # Specify the ontology ("MF" for molecular function, "BP" for biological process, "CC" for cellular component)
                      pAdjustMethod = "BH",  # Specify the p-value adjustment method ("BH" for Benjamini & Hochberg)
                      pvalueCutoff  = 0.05,  # Set the significance threshold
                      qvalueCutoff  = 0.05)

# View the top enriched terms
head(enrichRes)


gse <- gseGO(geneList=gene_list, 
             ont ="BP", 
             keyType = "SYMBOL", 
             nPerm = 10000, 
             minGSSize = 3, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = organism, 
             pAdjustMethod = "none")

require(DOSE)
dotplot(gse, showCategory=10, split=".sign") + facet_grid(.~.sign)

x2 <- pairwise_termsim(gse)

emapplot(x2, showCategory = 10)

cnetplot(x2, categorySize="pvalue", foldChange=gene_list, showCategory = 7)

# Use the `Gene Set` param for the index in the title, and as the value for geneSetId
gseaplot(gse, by = "all", title = gse$Description[1], geneSetID = 1)


# Convert gene IDs for gseKEGG function
# We will lose some genes here because not all IDs will be converted
converted_names <- gconvert(query = name, organism = "hsapiens",
                            target="ENTREZGENE_ACC", mthreshold = Inf, filter_na = TRUE)
converted_names

entrez_name <- converted_names[,'target']

names(original_gene_list) <- entrez_name

# omit any NA values 
gene_list<-na.omit(original_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
kegg_gene_list = sort(gene_list, decreasing = TRUE)

kegg_organism = "hsa"
kk2 <- gseKEGG(geneList     = kegg_gene_list,
               organism     = kegg_organism,
               nPerm        = 10000,
               minGSSize    = 2,
               maxGSSize    = 800,
               pvalueCutoff = 0.05,
               pAdjustMethod = "none",
               keyType       = "ncbi-geneid")

x2 <- pairwise_termsim(kk2)
dotplot(kk2, showCategory = 10, title = "Enriched Pathways" , split=".sign") + facet_grid(.~.sign)
emapplot(x2)

cnetplot(kk2, categorySize="pvalue", foldChange=gene_list)
install.packages("ggridges")
ridgeplot(kk2) + labs(x = "enrichment distribution")

# Use the `Gene Set` param for the index in the title, and as the value for geneSetId
gseaplot(kk2, by = "all", title = kk2$Description[5], geneSetID = 1)

library(pathview)
kk2
# Produce the native KEGG plot (PNG)
dme <- pathview(gene.data=kegg_gene_list, pathway.id="hsa04062", species = kegg_organism)

# Produce a different plot (PDF) (not displayed here)
dme <- pathview(gene.data=kegg_gene_list, pathway.id="hsa05323", species = kegg_organism, kegg.native = F)
