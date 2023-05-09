## Packages Loading
PACKAGES <- c(
  "dplyr",
  "edgeR",
  "genefilter",
  "glue",
  "recount3",
  "stringr",
  "tidyverse"
)

invisible(lapply(PACKAGES, library, character.only = TRUE))


# Set seed for reproducibility
SEED <- 1234
set.seed(SEED)

# PATH SETTIGN
IMG_PATH <- "C:/Users/aless/Desktop/script/Images/preprocessing/"
IN_PATH <- "C:/Users/aless/Desktop/script/input/"
OUT_PATH <- "C:/Users/aless/Desktop/script/output/preprocessing/"


# TISSUES PRE PROCESSING

## Data Loading
data <- read.delim(paste(IN_PATH, "GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct.gz", sep = ""), skip = 2)

str(data)
dim(data)

## Metadata Loading
# load the information about the genes(bp_length, gene_biotype, gene_name)
info_genes <- read.table(paste(OUT_PATH, "info_genes.txt", sep = "")) 

# load the info about the samples (tissue, code, ...)
sample_attributes <- read.delim(paste(IN_PATH, "GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt", sep = ""))

sample_attributes <- sample_attributes %>%
  select(one_of("SAMPID",
                "SMTSD",
                "SMUBRID")) %>%
  separate(., "SMTSD", " - ", into = c("TS", "SUBTS")) %>% # separate into tissue and subsection
  select(-c("SUBTS"))                                      # as we are only interested in the general tissue information

# manually remove the samples that refer to cell lines, and skin exposed to sun
exclude_tissue <- c("EFO_0002009", "EFO_0000572", "EFO_0002067", "0004264")
sample_attributes <- subset(sample_attributes, !(SMUBRID %in% exclude_tissue))

# the attribute data set contains records of 21166 samples
# exclude the samples that are not present in the expression matrix (which has 17382 records)
sample_attributes$SAMPID <- gsub("-", ".", sample_attributes$SAMPID)
sample_attributes <- sample_attributes[which(sample_attributes$SAMPID %in% colnames(data)), ]
data <- data[, c(1, which(colnames(data) %in% sample_attributes$SAMPID))]


# Reduce the data set by sampling n samples for each tissue
# the lowest number of samples per tissue is 21 (Bladder). 
# Sample 21 samples from each tissue

sample_attributes <- sample_attributes %>%      # aggregate tissues type with very low occurences
  mutate(across("TS", str_replace, "Fallopian Tube|Cervix", "Uterus"))

table(sample_attributes$TS)
N_TISSUES <- length(unique(sample_attributes$TS)) # 28

samples <- sample_attributes %>%
  group_by(TS) %>%
  slice_sample(n = 21)


# Reduce the expression matrix with the selected subset of samples
exprs <- data[, c(1, which(colnames(data) %in% samples$SAMPID))]
rm(data)


# Remove the version information from the ensembl id and compute the mean across rows that (without the version) refer to the same gene
exprs$ensembl_id <- gsub("\\..*", "", exprs$Name)
exprs <- exprs[which(exprs$ensembl_id %in% info_genes$ensembl_id), ]
dim(exprs)[1] - length(unique(exprs$ensembl_id))
# 18 duplicated ensembl ids

# since r does not allow to have duplicated row names, we use the transpose of the matrix and compute the mean across columns
rc <- t(exprs)
colnames(rc) <- exprs$ensembl_id
rc <- rc[- c(1, 590),]

str(rc)
rc <- apply(rc, 2, function(x) as.numeric(x))

exprs_clean <- t(apply(rc, 1, function(x) tapply(x, colnames(rc), mean)))

rownames(exprs_clean) <- str_sort(samples$SAMPID)
exprs_clean <- t(exprs_clean)

dim(exprs)[1] - dim(exprs_clean)[1]
# 18 rows have been removed


# Plotting function
segmented_hist <- function(data, rows){
  
  size = ceiling(dim(data)[1]/rows)
  
  for (seg in (1:rows)){
    
    int_start = (1 + (seg-1)*size)
    int_stop = seg*size
    
    tmp_df <- data[int_start:int_stop,]
    tmp_df <- gather(tmp_df, key = Sample, value = Counts)
    
    print(ggplot(tmp_df, aes(x = Samples, y = log(value))) +
            geom_boxplot(colour ="olivedrab", fill = "olivedrab", alpha = 0.7) +
            theme(axis.text.x = element_blank()) +
            scale_y_log10()
    )

    ggsave(
      paste(glue("tissues_norm_counts{seg}.jpg")),
      path = IMG_PATH,
      device='jpg')
  }
}


# NORMALIZATION
# Filter and order the genes info as gene lengths are needed for normalization
info_genes <- info_genes[which(info_genes$ensembl_id %in% rownames(exprs_clean)), ]
rownames(info_genes) <- info_genes$ensembl_id
info_genes <- info_genes[order(info_genes$ensembl_id), ]

NUM_SAMPLES <- dim(samples)[1] # 588 (21*28)

exprs_clean <- exprs_clean[order(rownames(exprs_clean)), ]

# Convert raw counts to reads per kb of transcript (intra-sample variation)
rpk_counts <- exprs_clean*(10^3)/info_genes$bp_length  

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


# FILTERING

# Filtering parameters
EXP_MIN <- 1       # Normalized reads to consider non-zero expression for a gene/feature in certain sample
EXP_THRESH <- 0.2  # Max fraction of samples that can have minimal expression in order for the gene/feature to be kept


# Remove features with very low expression in a certain amount of samples
counts_clean <- filterfun(kOverA(                  # Create filter function
  EXP_THRESH * NUM_SAMPLES, 
  EXP_MIN
)) %>%
  genefilter(counts_GeTMM, .) %>%    # Obtain indexer 
  filter(counts_GeTMM, .)            # Filter data set


# Plot normalized counts
counts_clean <- as.data.frame(counts_clean)
segmented_hist(counts_clean, 10)

# easier option to visualize only a fraction of the matrix
# df <- gather(as.data.frame(counts_clean[, c(1:500)]), key = Sample, value = Counts)
# 
# ggplot(data = df, aes(Sample, log(Counts+1))) +
#    geom_boxplot(colour ="olivedrab", fill = "olivedrab", alpha = 0.7) +
#    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
#    scale_y_log10()


write.table(counts_clean, paste(OUT_PATH, "tissues_norm_counts.txt", sep = ""), row.names = T, sep = "\t")
write.table(sample_attributes, paste(OUT_PATH, "tissues_sample_attributes.txt", sep = ""), sep = "\t", row.names = T)