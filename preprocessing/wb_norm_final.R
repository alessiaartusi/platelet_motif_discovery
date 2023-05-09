## Load Packages

PACKAGES <- c(
  "dplyr",
  "edgeR",
  "genefilter",
  "ggplot2",
  "glue",
  "recount3",
  "reshape2",
  "R.utils",
  "stringr",
  "tidyr",
  "tidyverse"
)

invisible(lapply(PACKAGES, library, character.only = TRUE))


IMG_PATH <- "C:/Users/aless/Desktop/script/Images/preprocessing/"
IN_PATH <- "C:/Users/aless/Desktop/script/input/"
OUT_PATH <- "C:/Users/aless/Desktop/script/output/preprocessing/"


gunzip(paste(IN_PATH, "GSE120312_Counts_Matrix.txt.gz", sep = ""))
raw_counts <- read.table(paste(IN_PATH, "GSE120312_Counts_Matrix.txt", sep = ""), skip = 1)

col_names <- "geneid	ensembl_gene_id	chr	gene_start_(bp)	gene_end_(bp)	strand	length_gene	length_trans	gene_type	class_types	12003	12004	12005	12006	12008	12009	12015	12016	12018	12044	12083	13103	13104	13105	13115	13116	13117	13140	13141	15002"

col_names <- str_split_1(col_names, "\t")
colnames(raw_counts) <- col_names

raw_counts <- raw_counts[which(raw_counts$gene_type == "protein_coding"), ]


info_genes <- raw_counts[, c(1:10)]
raw_counts <- raw_counts[,c(2, 11:30)]

raw_counts %>%
  group_by_all() %>%
  filter(n()>1) %>%
  distinct()

rownames(raw_counts) <- raw_counts[["ensembl_gene_id"]]
raw_counts <- raw_counts[, -c(1)]


# Normalization 

rpk_counts <- raw_counts*(10^3)/info_genes$length_gene  

NUM_SAMPLES <- dim(raw_counts)[2]

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


# Plot normalized counts

data <- gather(as.data.frame(counts_clean), key = Sample, value = Counts)

ggplot(data = data, aes(Sample, log(Counts+1))) +
  geom_boxplot(colour ="olivedrab", fill = "olivedrab", alpha = 0.7) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 8))

ggsave(
  paste(glue("whole_blood_norm_counts.jpg")),
  path = IMG_PATH,
  device='jpg')



# Save the obtained data for further analysis

write.table(as.data.frame(counts_clean), paste(OUT_PATH, "whole_blood_norm_counts.txt", sep = ""), row.names = T, sep = "\t")