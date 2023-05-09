# Loading Packages

PACKAGES <-  c(
  "biomaRt",
  "dplyr",
  "edgeR",
  "genefilter",
  "reshape2",
  "tidyr",
  "tidyverse"
)

invisible(lapply(PACKAGES, library, character.only = TRUE))

IMG_PATH <- "C:/Users/aless/Desktop/script/Images/preprocessing/"
OUT_PATH <- "C:/Users/aless/Desktop/script/output/preprocessing/"
IN_PATH <- "C:/Users/aless/Desktop/script/input/"


# Data Retrival

# Download pre-computed expression table
# gunzip(paste(IN_PATH, "GSE113182_raw_count.txt.gz", sep = ""))
# gunzip(paste(IN_PATH, "GSE113182_series_matrix.txt.gz", sep = ""))

data.matrix <- read.table(paste(IN_PATH, "GSE113182_raw_count.txt", sep = ""), header = TRUE, sep = "\t")
data.matrix <- column_to_rownames(data.matrix, var = "Feature.ID")

# df <- gather(log2(data.matrix  + 1), key = "sample", value = "counts")
# ggplot(df, aes(sample, counts)) +
#   geom_boxplot() + 
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5))


series.matrix <- read.csv(paste(IN_PATH, "GSE113182_series_matrix.csv", sep = ""), header = TRUE)
series.matrix <- as.data.frame(sapply(series.matrix, function(x) gsub("\"", "", x)))



# Select only the samples that refere to MEPs (megakaryocyte-erythrocyte progenitors)
data.matrix <- data.matrix %>% select(starts_with("MEP"))
series.matrix <- series.matrix[grep("MEP", series.matrix$Sample_description), ]


# Retrieve gene length using BiomaRt
ensembl <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

genes <- rownames(data.matrix)

genes_length <- getBM(attributes = c("ensembl_gene_id",
                                     "external_gene_name",
                                     "start_position",
                                     "end_position"),
                      filters = c("ensembl_gene_id"),
                      values = genes,
                      mart = ensembl)

gene_info <- dplyr::transmute(genes_length, ensembl_gene_id, external_gene_name, gene_length = end_position - start_position)


# Order the data sets
data.matrix <- data.matrix[order(rownames(data.matrix)), ]
gene_info <- gene_info[order(gene_info$ensembl_gene_id), ]


# Normalization

rpk_counts <- data.matrix*(10^3)/gene_info$gene_length  

NUM_SAMPLES <- dim(data.matrix)[2]

# Store rpk matrix into a DGEList object
rpk_norm <- DGEList(
  counts = rpk_counts, 
  group = c(rep("Dummy_var", NUM_SAMPLES))
)

# Compute and add the normalization factors (TMM default correction) to DGEList
# (inter-sample variation)
rpk_norm <- calcNormFactors(rpk_norm)

# Re-scale normalized counts in counts per million
counts_GeTMM <- as.data.frame(cpm(rpk_norm))

# Remove non-needed objects
rm(rpk_counts, rpk_norm)


## Filtering
# Filtering parameters grouped for ease of access

EXP_MIN <- 1       # Normalized reads to consider non-zero expression for...
# ... a gene/feature in certain sample
EXP_THRESH <- 0.2  # Max fraction of samples that can have minimal expression...
# ... in order for the gene/feature to be kept


# Remove features with very low expression in a certain amount of samples

counts_clean <- filterfun(kOverA(                  # Create filter function
  EXP_THRESH * NUM_SAMPLES, 
  EXP_MIN
)) %>%
  genefilter(counts_GeTMM, .) %>%    # Obtain indexer 
  filter(counts_GeTMM, .)            # Filter data set

# Plot the distribution of counts
df <- gather(counts_clean, key = "sample", value = "counts")

ggplot(df, aes(sample, log(counts + 1))) +
  geom_boxplot(colour ="olivedrab", fill = "olivedrab", alpha = 0.7) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 8))

ggsave(
  "mk_norm_counts.jpg",
  device = "jpeg",
  path = IMG_PATH,
  dpi = 700
)

# Save the data frame of clean counts
write.table(as.data.frame(counts_clean), paste(OUT_PATH, "mk_norm_counts.txt", sep = ""), row.names = T, sep = "\t")