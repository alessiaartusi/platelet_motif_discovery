## Load Packages

PACKAGES <- c(
  "dplyr",
  "edgeR",
  "genefilter",
  "ggplot2",
  "glue",
  "recount3",
  "reshape2",
  "stringr",
  "tidyr",
  "tidyverse"
)

invisible(lapply(PACKAGES, library, character.only = TRUE))

IMG_PATH <- "C:/Users/aless/Desktop/script/Images/preprocessing/"
OUT_PATH <- "C:/Users/aless/Desktop/script/output/preprocessing/"

# PLATELETS Data Pre Processing

## Retrieve Data

### Download data from recount3; the object contains raw counts data and metadata about the samples and the genes.

options(recount3_url = "https://recount-opendata.s3.amazonaws.com/recount3/release")

recount_data <- recount3::create_rse_manual(
  project = "SRP057500",
  project_home = "data_sources/sra",
  organism = "human",
  annotation = "gencode_v26",
  type = "gene"
)


## Extract Useful Data

# r: genes, c: samples
raw_counts <- compute_read_counts(recount_data)

# r: genes, c: gene info
info_genes  <- as.data.frame(rowData(recount_data)) %>%
  dplyr::select(one_of(
    "bp_length",
    "gene_type",
    "gene_name")
  )

info_genes$ensembl_id <- rownames(info_genes)

# the expression matrix contains every possible RNA type, not only protein coding
sort(table(info_genes$gene_type), decreasing = T)

# Keep only protein coding genes
info_genes <- info_genes[which(info_genes$gene_type == "protein_coding"), ]


# r: samples, c: sample type
info_samples  <- as.data.frame(colData(recount_data)) %>%
  dplyr::select(., sra.experiment_title) %>%
  mutate_at("sra.experiment_title", function(x) str_replace(x, ".*:", "")) %>%
  mutate_at("sra.experiment_title", function(x) str_replace(x, ";.*", ""))

# filter only the healthy samples               
healthy_samples  <- info_samples %>%
  dplyr::filter(grepl("Blood_Platelets_HC", sra.experiment_title))

dim(healthy_samples)


#There are 5 samples more than the reported number (55).
healthy_samples %>%
  group_by_all() %>%
  filter(n()>1) %>%
  distinct()


#Filter the healthy samples in the raw counts matrix and for protein coding genes
raw_counts <- raw_counts[which(rownames(raw_counts) %in% info_genes$ensembl_id), 
                         which(colnames(raw_counts) %in% rownames(healthy_samples))]
dim(raw_counts)


## Data Cleaning
df <- raw_counts
#To merge the samples that come from the same patient, compute the expression mean across columns that have the same name

# substitute SRR sample names with the name sample name 'Blood_Platelets_HC-X' to recognize columns with the same name (SRR are different)
colnames(df) <- healthy_samples$sra.experiment_title 

# compute the row-mean of columns with the same name
df_clean <- sapply(split(seq_len(ncol(df)),colnames(df)),function(cis) rowMeans(df[,cis,drop=F]))

dim(df_clean)
raw_counts <- as.data.frame(df_clean)


#When removing the version specifications (.x) from the ensembl ids, we end up with multiple equal ensembl ids.
#Hence, we need to compute the mean across rows for the occurrences that refer to the same gene.

# cut the .x specification
raw_counts$ensembl_id <- gsub("\\..*", "", rownames(raw_counts))

dim(raw_counts)[1] - length(unique(raw_counts$ensembl_id))
# there are 18 duplicated ensembl ids

# since r does not allow to have duplicated row names, we use the transpose of the matrix and compute the mean across columns
# r: samples, c: genes
rc <- t(raw_counts)
colnames(rc) <- raw_counts$ensembl_id
rc <- rc[-56,]

str(rc)
rc <- apply(rc, 2, function(x) as.numeric(x))

raw_counts_clean <- t(apply(rc, 1, function(x) tapply(x, colnames(rc), mean)))

rownames(raw_counts_clean) <- unique(str_sort(healthy_samples$sra.experiment_title))
raw_counts_clean <- t(raw_counts_clean)

dim(raw_counts)[1] - dim(raw_counts_clean)[1]
# 18 rows have been removed


# Remove not needed objects
rm(recount_data, df, df_clean, raw_counts, rc, info_samples)


# Pre processed data visualization
# df <- gather(as.data.frame(raw_counts_clean), key = sample, value = COUNTS)
# 
# ggplot(data = df, aes(sample, COUNTS+1)) +
#   geom_boxplot(colour ="olivedrab", fill = "olivedrab", alpha = 0.7) +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 8)) +
#   scale_y_log10()
#
# ggsave(
#   paste(glue("platelets_raw_counts.jpg")),
#   path = IMG_PATH,
#   device='jpg',
#   dpi=700)


## Normalization
# Normalize counts using GeTMM (for both inter and intra sample variation). 
# Convert raw counts to reads per kb of transcript (intra-sample variation)
info_genes$ensembl_id <- gsub("\\..*", "", info_genes$ensembl_id)
info_genes <- info_genes[!duplicated(info_genes$ensembl_id), ]
info_genes <- info_genes[which(info_genes$ensembl_id %in% rownames(raw_counts_clean)),]
info_genes <- info_genes[order(rownames(info_genes)), ]

rpk_counts <- raw_counts_clean*(10^3)/info_genes$bp_length  

NUM_SAMPLES <- 55

# Store rpk matrix into a DGEList object
rpk_norm <- DGEList(
  counts=rpk_counts, 
  group=c(rep("Dummy_var", NUM_SAMPLES))
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


# Plot the normalized data
data <- gather(as.data.frame(counts_clean), key = Sample, value = Counts)

ggplot(data = data, aes(Sample, Counts+1)) +
  geom_boxplot(colour ="olivedrab", fill = "olivedrab", alpha = 0.7) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=8)) + 
  scale_y_log10()

ggsave(
  paste(glue("preprocessing/platelets_normalized_counts.jpg")),
  path = IMG_PATH,
  device='jpg',
  dpi=700)



# Save the obtained data for further analysis
write.table(as.data.frame(counts_clean), paste(OUT_PATH, "preprocessing/platelets_norm_counts.txt", sep = ""), row.names = T, sep = "\t")
write.table(info_genes, paste(OUT_PATH, "preprocessing/info_genes.txt", sep = ""), sep = "\t")
# info_genes contains only information about protein coding genes

