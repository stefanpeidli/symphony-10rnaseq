library('DESeq2')
library('BiocParallel')
register(MulticoreParam(4))

# Snakemake I/O
count_file <- snakemake@input[['counts']]
coldata_file <- snakemake@input[['coldata']]
cell_line <- snakemake@wildcards[['cell_line']]
perturbation <- snakemake@wildcards[['perturbation']]

# Load data
df <- read.csv(file = count_file, header = TRUE, sep = "\t", check.names = FALSE)
countData <- df[, c(2, 6:ncol(df))]
rownames(countData) <- countData[[1]]
countData[[1]] <- NULL
colData <- read.csv(file = coldata_file, row.names=1, sep="\t")

# Factors
colData$line         <- factor(colData$line)
colData$genotype     <- relevel(factor(colData$genotype), "WT")  # WT = baseline for genotype
colData$perturbation <- relevel(factor(colData$perturbation), "NT")  # NT = baseline control

# Subset to the specified cell line
count <- countData[, colData$line == cell_line]
col <- droplevels(colData[colData$line == cell_line, ])

# Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(countData = count,
                              colData = col,
                              design = ~ genotype + perturbation + genotype:perturbation)
dds <- dds[rowSums(counts(dds)) >= 10, ]  # min counts per gene
dds <- DESeq(dds)

# Also test purely in WT background
keep_wt <- col$genotype == "WT"
dds_wt  <- DESeqDataSetFromMatrix(count[, keep_wt], droplevels(col[keep_wt,]), design = ~ perturbation)
dds_wt  <- DESeq(dds_wt)

# Cas-dependence of results (positive control)
res_genotype_effect <- results(dds, name = paste0("perturbation_", perturbation, ".genotypeCAS"))
write.csv(as.data.frame(res_genotype_effect), file = snakemake@output[['res_genotype_effect']])

# This one should be flat since we are in WT
res_wt_vim_eq <- results(dds_wt, name=paste0("perturbation_", perturbation, "_vs_NT"),
                         lfcThreshold=0.25, altHypothesis="lessAbs")
write.csv(as.data.frame(res_wt_vim_eq), file = snakemake@output[['res_wt_vim_eq']])

# Negative controls, those should technically all look like NT control in BOTH genotypes
res_mruby_vs_nt_cas <- results(dds, contrast=c("perturbation","RUBY","NT"))
write.csv(as.data.frame(res_mruby_vs_nt_cas), file = snakemake@output[['res_mruby_vs_nt_cas']])
res_mock_vs_nt_cas  <- results(dds, contrast=c("perturbation","MOCK","NT"))
write.csv(as.data.frame(res_mock_vs_nt_cas), file = snakemake@output[['res_mock_vs_nt_cas']])
res_mruby_eq_wt <- results(dds_wt, contrast=c("perturbation","RUBY","NT"),
                           lfcThreshold=0.25, altHypothesis="lessAbs")
write.csv(as.data.frame(res_mruby_eq_wt), file = snakemake@output[['res_mruby_eq_wt']])
res_mock_eq_wt  <- results(dds_wt, contrast=c("perturbation","MOCK","NT"),
                           lfcThreshold=0.25, altHypothesis="lessAbs")
write.csv(as.data.frame(res_mock_eq_wt), file = snakemake@output[['res_mock_eq_wt']])

# Global Cas vs WT effect
res_genotype_main <- results(dds, name="genotype_CAS_vs_WT")
write.csv(as.data.frame(res_genotype_main), file = snakemake@output[['res_genotype_main']])