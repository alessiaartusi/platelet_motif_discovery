## Packages Loading

PACKAGES <- c(
  "dplyr",
  "edgeR",
  "factoextra",
  "ggplot2",
  "glue",
  "RankProd",
  "stringr",
  "tidyr",
  "tidyverse"
)

invisible(lapply(PACKAGES, library, character.only = TRUE))

IMG_PATH <- "C:/Users/aless/Desktop/script/Images/"
OUT_PATH <- "C:/Users/aless/Desktop/script/output/"



## Load Data
platelets <- read.table(paste(OUT_PATH, "preprocessing/platelets_norm_counts.txt", sep = ""))
colnames(platelets) <- gsub("X.", "", colnames(platelets))

wb <- read.table(paste(OUT_PATH, "preprocessing/whole_blood_norm_counts.txt", sep = ""))
colnames(wb) <- gsub("X", "", colnames(wb))

platelets <- platelets[which(rownames(platelets) %in% rownames(wb)), ]
wb <- wb[which(rownames(wb) %in% rownames(platelets)), ]

# Merge data sets
merged <- cbind(platelets, wb)
rm(platelets, wb)


# Plot the distribution of the counts
# Log transform the matrix of counts
log_merged <- log2(merged +1)

df <- gather(log_merged, key = Samples, value = cmp)

ggplot(df, aes(Samples, cmp)) + 
  geom_boxplot(fill = "grey", alpha = 0.6) +
  theme(element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5))


NUM_SAMPLES <- dim(merged)[2]

# create a data frame of attributes per each sample
info_samples <- data.frame(SAMPID = colnames(merged), TS = c(rep("platelets", 55), rep("wb", 20)))

## Principal component analysis 
# Transpose matrix and perform PCA
TISSUES <- info_samples$TS %>%
  as.factor()

res_pca <- prcomp(t(log_merged))

# Table of dimension contribution
eig_PCA <- get_eig(res_pca)


# Display and save the scree plot 
fviz_eig(res_pca) + 
  ggtitle("Scree plot platelets - whole blood")
ggsave(
  paste(glue("rankprod/screeplot_platelets_wb.jpg")),
  path = IMG_PATH,
  device='jpeg')


# Save order of genes by model contribution 
ord_pca <- facto_summarize(res_pca, "var") %>%
  arrange(desc(contrib))


# Plot PCA and save it
colors <- c("#56CBF9", "#FF729F")

fviz_pca_ind(
  res_pca,
  axes = c(1,2),
  geom = c("point"),
  label = "none",
  habillage = TISSUES,
  palette = colors,
  pointshape = 16,
  pointsize = 2,
) + 
  ggtitle("PCA Platelets - Whole Blood")

ggsave(
  "rankprod/platelets_wb_pca.jpg",
  device = "jpeg",
  path = IMG_PATH
)


# Rank Product
# Create a vector that specifies the class of the samples: 0 for platelets and 1 for mk.
sample_origin <- data.frame(SAMPID = info_samples$SAMPID, CLASS = c(rep(0, 55), rep(1, 20)))
RP.out <- RankProducts(log_merged, sample_origin$CLASS, logged = TRUE, na.rm = FALSE, gene.names = rownames(log_merged))

head(RP.out$pfp)
RP.genes <- topGene(RP.out,
                    cutoff = 0.05,
                    method = "pval",
                    logged = TRUE,
                    gene.names = rownames(log_merged))

# genes over expressed in mk vs platelets
head(RP.genes$Table1)
dim(RP.genes$Table1)[1] # 2505

# genes over expressed in platelets vs wb
head(RP.genes$Table2)
dim(RP.genes$Table2)[1] # 2532

# Write the over-expressed genes into a file
write.table(RP.genes$Table2, paste(OUT_PATH, "rankprod/genes_up_wb.txt", sep = ""), sep = "\t")
write.table(RP.genes$Table1, paste(OUT_PATH, "rankprod/genes_down_wb.txt", sep = ""), sep = "\t")
