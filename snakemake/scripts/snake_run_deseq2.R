library('DESeq2')
library('BiocParallel')
register(MulticoreParam(4))

# Snakemake I/O
count_file <- snakemake@input[['counts']]
coldata_file <- snakemake@input[['coldata']]
cell_line <- snakemake@wildcards[['cell_line']]
perturbation <- snakemake@wildcards[['perturbation']]
outpath <- snakemake@output[['deseq2_results']]

# Load data
df <- read.csv(file = count_file, header = TRUE, sep = "\t", check.names = FALSE)
countData <- df[, c(2, 6:ncol(df))]
rownames(countData) <- countData[[1]]
countData[[1]] <- NULL
colData <- read.csv(file = coldata_file, row.names=1, sep="\t")
colData$perturbation <- factor(colData$perturbation)
colData$perturbation <- droplevels(colData$perturbation)
colData$perturbation <- relevel(colData$perturbation, ref = "NT")  # set NT as reference

#TODO: this
count_cas <- countData[, colData$genotype == "CAS"]
col_cas   <- droplevels(colData[colData$genotype == "CAS", ])

# Subset to the specified cell line
count_cas <- count_cas[, col_cas$line == cell_line]
col_cas <- droplevels(col_cas[col_cas$line == cell_line, ])

# Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(countData = count_cas,
                              colData = col_cas,
                              design = ~ perturbation)
dds <- dds[rowSums(counts(dds)) >= 10, ]  # min counts per gene
# Run DESeq
dds <- DESeq(dds)

# Get results
perturbations <- setdiff(levels(colData(dds)$perturbation), "NT")
  
# Construct the coefficient name
coef_name <- paste0("perturbation_", perturbation, "_vs_NT")

# Check if coef exists in resultsNames() to avoid errors
if (!(coef_name %in% resultsNames(dds))) {
message("WARNING: ", perturbation, " — no coefficient named ", coef_name)
}
message("Processing: ", perturbation)
# Get DESeq2 results
res <- results(dds, name = coef_name)
# Shrink LFCs (currently using "normal"; swap to apeglm if available)
# TODO: replace this with apeglm or ashr if possible
res_shr <- lfcShrink(dds, coef = coef_name, type = "normal")
# Order by p-value
resOrdered <- res_shr[order(res_shr$pvalue), ]
# Write results to CSV
write.csv(as.data.frame(resOrdered), file = outpath)
