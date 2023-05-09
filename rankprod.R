## Load Packages

PACKAGES <- c(
  "dplyr",
  "edgeR",
  "factoextra",
  "genefilter",
  "ggplot2",
  "glue",
  "RankProd",
  "stringr",
  "tidyr",
  "tidyverse",
  "viridis"
)

invisible(lapply(PACKAGES, library, character.only = TRUE))


IMG_PATH <- "C:/Users/aless/Desktop/script/Images/"
OUT_PATH <- "C:/Users/aless/Desktop/script/output/"

## Load Data

platelets <- read.table(paste(OUT_PATH, "preprocessing/platelets_norm_counts.txt", sep = ""))
colnames(platelets) <- gsub("X.", "", colnames(platelets))

tissues <- read.table(paste(OUT_PATH, "preprocessing/tissues_norm_counts.txt", sep = ""))

platelets <- platelets[which(rownames(platelets) %in% rownames(tissues)), ]
tissues <- tissues[which(rownames(tissues) %in% rownames(platelets)), ]


# Data sets merge
merged <- cbind(platelets, tissues)



# Load the tissues attributes
sample_attributes <- read.table(paste(OUT_PATH, "preprocessing/tissues_sample_attributes.txt", sep = ""))
sample_attributes <- sample_attributes[which(sample_attributes$SAMPID %in% colnames(merged)),]

# create an attribute data frame for the platelets samples
sample_attributes_2 <- data.frame(SAMPID = colnames(platelets), TS = rep("Platelet", 55), SMUBRID = rep(0, 55))
info_samples <- rbind(sample_attributes, sample_attributes_2)
info_samples <- info_samples[order(info_samples$SAMPID),]


info_genes <- read.table(paste(OUT_PATH, "preprocessing/info_genes.txt", sep = ""))

info_genes <- info_genes[which(info_genes$ensembl_id %in% rownames(merged)), ]
rownames(info_genes) <- info_genes$ensembl_id
info_genes <- info_genes[order(info_genes$ensembl_id),]


rm(sample_attributes, sample_attributes_2)
rm(platelets, tissues)


# Plot the distribution of the counts
# Log transform the matrix of counts
log_merged <- log2(merged +1)

# 
# df <- gather(log_merged[, c(1:150)], key = Samples, value = cmp)
# 
# ggplot(df, aes(Samples, cmp)) + 
#   geom_boxplot(fill = "grey", alpha = 0.6) +
#   theme(element_blank())

NUM_SAMPLES <- dim(info_samples)[1]

## Principal component analysis 
# Transpose matrix and perform PCA
info_samples <- info_samples[order(info_samples$SAMPID),]
TISSUES <- info_samples$TS
TISSUES[TISSUES != "Platelet"] <- "Other"

res_pca <- prcomp(t(log_merged))

# Table of dimension contribution
eig_PCA <- get_eig(res_pca)


# Display and save the scree plot 
fviz_eig(res_pca) + 
  ggtitle("Scree plot platelets - tissues") + 
  ylim(0, 30)
ggsave(
  paste(glue("rankprod/screeplot_platelets_tissues.jpg")),
  path = IMG_PATH,
  device='jpg')

# Save order of genes by model contribution 
ord_pca <- facto_summarize(res_pca, "var") %>%
  arrange(desc(contrib))

# Plot PCA and save it
colors <- viridis(2)

fviz_pca_ind(
  res_pca,
  axes = c(1,2),
  geom = c("point"),
  label = "none",
  habillage = TISSUES,
  palette = colors,
  pointshape = 19,
  pointsize = 2
) + 
  ggtitle("PCA platelets - tissues")

ggsave(
  "rankprod/platelets_tissues_pca2.jpg",
  device = "jpeg",
  path = IMG_PATH
)



# Rank Product
# Create a vector that specifies the class of the samples: 0 for platelets and 1 for tissues.
sample_origin <- data.frame(SAMPID = info_samples$SAMPID, CLASS = c(rep(0, 55), rep(1, 588)))
RP.out <- RankProducts(log_merged, sample_origin$CLASS, logged = TRUE, na.rm = FALSE, gene.names = info_genes$gene_name, calculateProduct = FALSE)

head(RP.out$pfp)
RP.genes <- topGene(RP.out,
                    cutoff = 0.05,
                    method = "pval",
                    logged = TRUE,
                    gene.names = info_genes$ensembl_id)

# genes over expressed in other vs platelets
print("genes over expressed in other vs platelets")
head(RP.genes$Table1)
print("total:")
dim(RP.genes$Table1)[1]

# genes over expressed in platelets vs other
print("genes over expressed in platelets vs other")
head(RP.genes$Table2)
print("total:")
dim(RP.genes$Table2)[1] # 12869

# Write the over-expressed genes into a file
write.table(RP.genes$Table2, paste(OUT_PATH, "rankprod/genes_up.txt", sep = ""), sep = "\t")
write.table(RP.genes$Table1, paste(OUT_PATH, "rankprod/genes_down.txt", sep = ""), sep = "\t")
