## Packages Loading 

PACKAGES <- c(
  "clusterProfiler",
  "dplyr",
  "enrichplot",
  "ggplot2",
  "org.Hs.eg.db",
  "pathview",
  "tidyr"
)

invisible(lapply(PACKAGES, library, character.only = TRUE))



IMG_PATH <- "C:/Users/aless/Desktop/script/Images/"
OUT_PATH <- "C:/Users/aless/Desktop/script/output/"


# Data Input
# Load ranked gene list

# up-regulated genes
up_genes <- read.table(paste(OUT_PATH, "rankprod/genes_up_wb.txt", sep = ""))
head(up_genes)
up_genes$id <- rownames(up_genes)
colnames(up_genes)[3] <- "FC"
up_genes$log2FC <- log2(up_genes$FC)

# down-regulated genes
down_genes <- read.table(paste(OUT_PATH, "rankprod/genes_down_wb.txt", sep = ""))
head(down_genes)
down_genes$id <- rownames(down_genes)
colnames(down_genes)[3] <- "FC"
down_genes$log2FC <- log2(down_genes$FC)

# merge the two lists for GSEA
de_genes <- rbind(up_genes, down_genes)


# Adjust p-value with fdr method and filter for adj-p-value < 0.05
de_genes$adjpval <- p.adjust(de_genes$P.value, method = "fdr", n = length(de_genes$P.value))
de_genes <- de_genes[which(de_genes$adjpval < 0.05), ]


# Plot the distribution of the log2FC for up- and down-regulated genes
de_genes <- de_genes %>%
  mutate(
    DEG = case_when(log2FC >= 1.5 ~ "Up-regulated",
                    log2FC <= -1.5 ~ "Down-regulated",
                    TRUE ~ "Unchanged")
  )

ggplot(de_genes, aes(DEG, log2FC)) + 
  geom_boxplot()

table(de_genes$DEG)

summary(de_genes)

# Plot the volcano plot of DEGs
# Differentially expressed genes (DEGs) are usually considered as those with an 
# absolute fold change greater or equal to 2 and a FDR value of 0.05 or less.
ggplot(de_genes, aes(log2FC, -log(adjpval,10))) + # -log10 conversion  
  geom_point(size = 2/5) +
  xlab(expression("log"[2]*"FC")) + 
  ylab(expression("-log"[10]*"FDR"))

ggplot(de_genes, aes(log2FC, -log(adjpval, 10))) +
  geom_point(aes(color = DEG), size = 2/5) +
  xlab(expression("log"[2]*"FC")) + 
  ylab(expression("-log"[10]*"FDR")) +
  scale_color_manual(values = c("dodgerblue3", "gray50", "firebrick3")) +
  guides(colour = guide_legend(override.aes = list(size=1.5)))


# Ranks - log2fold-change
organism <- "org.Hs.eg.db"

ranks <- de_genes$log2FC
names(ranks) <- de_genes$id


# Sort the list in decreasing order (required for clusterProfiler)
ranks <- sort(ranks, decreasing = TRUE)

# GO GSEA
# use ensembl
gse <- gseGO(geneList = ranks,
             keyType = "ENSEMBL",
             ont = "BP",
             pvalueCutoff = 0.05,
             OrgDb = organism,
             pAdjustMethod = "BH")


# Dotplot the results
require(DOSE)
dotplot(gse, showCategory = 10) + 
  ggtitle("GSEA platelets - whole blood")

ggsave(
  "gsea/gsea_wb.jpg",
  path = IMG_PATH,
  device = "jpeg"
)

# plot the running Enrichment Score of the first gene set
# gseaplot(gse, by = "all", title = gse$Description[1], geneSetID = 1)

# The GSEA plot shows only one value for the adjusted p-value because the p-values 
# of the first 18 enriched terms share at least 30 core genes, hence most of them are composed by the same genes. 

gene_sets <- strsplit(gse@result[c(1:10),]$core_enrichment, "\t")
list_gene_set <- sapply(gene_sets, function(x) strsplit(x, "/"))
genes_gsea <- unique(c(list_gene_set[[1]],
                       list_gene_set[[2]],
                       list_gene_set[[3]],
                       list_gene_set[[4]],
                       list_gene_set[[5]],
                       list_gene_set[[6]],
                       list_gene_set[[7]],
                       list_gene_set[[8]],
                       list_gene_set[[9]],
                       list_gene_set[[10]]))

print("gsea genes:")
length(genes_gsea)
write.csv(genes_gsea, paste(OUT_PATH, "gsea/gene_list_gsea_wb.csv", sep = ""), row.names = FALSE)

# sapply(list_gene_set, function(x) length(x))

# ORA
# Filter the genes using a log2FC > 1.5
genes <- names(ranks)[ranks >= 1.5]
print("gene with logFC > 1.5")
length(genes)

write.csv(genes, paste(OUT_PATH, "rankprod/gene_list_logFCthresh_wb.csv", sep = ""), row.names = FALSE)

genes_fc <- up_genes$id[which(up_genes$FC > 1.5)]
write.csv(genes_fc, paste(OUT_PATH, "rankprod/gene_list_FCthresh_wb.csv", sep = ""), row.names = FALSE)

# ORA
go_enrich <- enrichGO(gene = genes,
                      ont = "BP",
                      pvalueCutoff = 0.05,
                      OrgDb = organism, 
                      keyType = "ENSEMBL")


# Plot the results
barplot(go_enrich, 
        drop = TRUE, 
        showCategory = 10, 
        title = "GO pathways ORA plateletes - whole blood",
        font.size = 9)

ggsave(
  "gsea/ora_barplot_wb.jpg",
  path = IMG_PATH,
  device = "jpeg"
)

jpeg(file = paste(IMG_PATH, "gsea/ora_upsetplot_wb.png", sep = ""))
upsetplot(go_enrich)
dev.off()