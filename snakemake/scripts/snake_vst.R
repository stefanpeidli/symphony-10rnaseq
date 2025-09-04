library('DESeq2')
library('BiocParallel')
register(MulticoreParam(4))

# Snakemake I/O
count_file <- snakemake@input[['counts']]
coldata_file <- snakemake@input[['coldata']]
vst_counts <- snakemake@output[['vst_counts']]

# Load data
df <- read.csv(file = count_file, header = TRUE, sep = "\t", check.names = FALSE)
countData <- df[, c(2, 6:ncol(df))]
rownames(countData) <- countData[[1]]
countData[[1]] <- NULL
colData <- read.csv(file = coldata_file, row.names=1, sep="\t")
colData$perturbation <- factor(colData$perturbation)
colData$perturbation <- droplevels(colData$perturbation)
colData$perturbation <- relevel(colData$perturbation, ref = "NT")
count_cas <- countData[, colData$genotype == "CAS"]
col_cas   <- droplevels(colData[colData$genotype == "CAS", ])

# Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(countData = count_cas,
                              colData = col_cas,
                              design = ~ line + perturbation)
dds <- dds[rowSums(counts(dds)) >= 10, ]  # min counts per gene
# Run DESeq
dds <- DESeq(dds)

vsd <- vst(dds, blind = TRUE)
mat <- assay(vsd)
write.csv(mat, file = vst_counts)
