OUT_PATH = "C:/Users/aless/Desktop/script/output/"

# List intersections

## TISSUES and MEP CD34+
ts <- read.csv(file = paste(OUT_PATH, "rankprod/gene_list_logFCthresh.csv", sep = ""))
mep <- read.csv(file = paste(OUT_PATH, "rankprod/gene_list_logFCthresh_mk.csv", sep = ""))

int <- mep$x[mep$x %in% ts$x]
# the intersection is made of 360 genes

## MEP CD34+ and WHOLE BLOOD
wb <- read.csv(file = paste(OUT_PATH, "rankprod/gene_list_logFCthresh_wb.csv", sep = ""))
int_2 <- mep$x[mep$x %in% wb$x]
# the intersection is made of 457 genes

# Intersect the intersections
INT <- int[int %in% int_2]
# 222 genes intersected from tissues, MEP CD34+ and whole blood analyses
write.table(INT, paste(OUT_PATH, "intersected_gene_list.csv", sep = ""), row.names = FALSE)

# Check the similarity of the lists
path <- "C:/Users/aless/Downloads/drive-download-20230427T105031Z-001/"

# TISSUES
ts_old <- read.csv(paste(path, "gene_list_logFCthresh.csv", sep = ""))
same <- mean(ts_old$x %in% ts$x) # 99.74 % (only one gene is not present)
# the new list contains 246 genes more

# MEP CD34+
mep_old <- read.csv(paste(path, "gene_list_logFCthresh_mk.csv", sep = ""))
same <- mean(mep_old$x %in% mep$x) # 100%, the new list contains 313 genes more