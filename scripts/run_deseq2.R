library('DESeq2')
library('BiocParallel')
register(MulticoreParam(4))

# Load data
df <- read.csv(file = "prelim_combined_counts.tsv", header = TRUE, sep = "\t", check.names = FALSE)
countData <- df[, c(2, 6:ncol(df))]
rownames(countData) <- countData[[1]]
countData[[1]] <- NULL
colData <- read.csv(file = "prelim_combined_coldata.tsv", row.names=1, sep="\t")
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

# Get results
perturbations <- setdiff(levels(colData(dds)$perturbation), "NT")
outdir <- "deseq2_results"
dir.create(outdir, showWarnings = FALSE)
for (pert in perturbations) {
  
  # Construct the coefficient name
  coef_name <- paste0("perturbation_", pert, "_vs_NT")
  
  # Check if coef exists in resultsNames() to avoid errors
  if (!(coef_name %in% resultsNames(dds))) {
    message("Skipping ", pert, " — no coefficient named ", coef_name)
    next
  }
  message("Processing: ", pert)
  # Get DESeq2 results
  res <- results(dds, name = coef_name)
  # Shrink LFCs (currently using "normal"; swap to apeglm if available)
  # TODO: replace this with apeglm or ashr if possible
  res_shr <- lfcShrink(dds, coef = coef_name, type = "normal")
  # Order by p-value
  resOrdered <- res_shr[order(res_shr$pvalue), ]
  # Write results to CSV
  outfile <- file.path(outdir, paste0(pert, "_deseq2_results.csv"))
  write.csv(as.data.frame(resOrdered), file = outfile)
}

# # Plot
# vsd <- vst(dds, blind = TRUE)
# png("PCA_plot.png", width = 1600, height = 1600, res = 200)
# plotPCA(vsd, intgroup = c("line"))  # Lines should cluster together
# dev.off()  # close device and save file

# lines <- levels(droplevels(colData(dds)$line))

# for (ln in lines) {
#   keep <- colData(dds)$line == ln
#   vsd_ln <- vst(dds[, keep], blind = TRUE)
  
#   # open PNG device
#   png(sprintf("PCA_%s_by_perturbation.png", ln), width = 1600, height = 1200, res = 200)
  
#   # explicitly print the ggplot object
#   print(plotPCA(vsd_ln, intgroup = "perturbation"))
  
#   # close device
#   dev.off()
# }