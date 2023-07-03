getwd()

library(clusterProfiler)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(enrichplot)
library(enrichplot)
library(ggplot2)
deg_res <- read.csv('nowy_sigs.csv')
deg_res
deg_res <- deg_res[(abs(deg_res$log2FoldChange) > 1),]
genelist <- deg_res$refeference_id
gene_entrez <- mapIds(org.Hs.eg.db, keys = genelist, keytype = "REFSEQ", column = "SYMBOL")

GO_results <- enrichGO(gene = genelist, keyType = 'SYMBOL', OrgDb = "org.Hs.eg.db", ont = "ALL")

dotplot(GO_results, title = "Wzbogacone terminy GO dla typu TNBC")
plot(barplot(GO_results, showCategory =20, title = "Wzbogacone terminy GO dla typu TNBC"))

ego2 <- simplify(GO_results, cutoff=0.7, by="p.adjust", select_fun=min)

ego <- pairwise_termsim(GO_results)
ego2 <- pairwise_termsim(ego2)

p1 <- emapplot(ego, cex_label_category=.8, cex_line=.5) + coord_cartesian()
p2 <- emapplot(ego2, cex_label_category=.8, cex_line=.5) + coord_cartesian()


p1 <- p1 + scale_fill_continuous(low = "#e06663", high = "#327eba", name = "p.adjust",
                                 guide = guide_colorbar(reverse = TRUE, order=1), trans='log10')
p2 <- p2 + scale_fill_continuous(low = "#e06663", high = "#327eba", name = "p.adjust",
                                 guide = guide_colorbar(reverse = TRUE, order=1), trans='log10')
cowplot::plot_grid(p1, p2, labels=c("A", "B"), rel_widths=c(1, 1.2))



require(DOSE)

ggo <- groupGO(gene = genelist, OrgDb = "org.Hs.eg.db", ont = "BP", level = 3,  keyType = 'SYMBOL', readable = TRUE)
head(summary(ggo))

ego <- enrichGO(gene = genelist, universe = names(geneList), OrgDb = "org.Hs.eg.db", keyType = 'SYMBOL', ont = "CC", pvalueCutoff = 0.01,readable = TRUE)
head(summary(ego))

kk <- enrichKEGG(gene = gene_entrez, organism = "human", pvalueCutoff = 0.01)
head(summary(kk))

barplot(ggo, drop = TRUE, showCategory = 20)
barplot(ego, showCategory = 10, title = "Wzbogacenie terminów GO dla podtypu TNBC")

cnetplot(ego, categorySize = "pvalue", foldchange = deg_res$log2FoldChange)
cnetplot(kk, categorySize = "geneNum", foldChange = deg_res$log2FoldChange, showCategory = 10,
         node_label = "gene")


ck <- compareCluster(geneCluster = genelist, fun = "enrichKEGG")
plot(ck)
