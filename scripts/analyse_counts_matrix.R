library(DESeq2)
library(tidyverse)
library(edgeR)
library(ComplexHeatmap)
library(circlize)
library(biomaRt)
library(ggrepel)

counts_table <- read_tsv("results/8094_Sym101_bulk_counts_matrix.tsv") 

counts_table <- counts_table |> 
dplyr::filter(grepl(pattern = "ENSG", GENE))


# Annotate Ensembl IDs with gene symbols via biomaRt
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
gene_annot <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters = "ensembl_gene_id",
  values = counts_table$GENE,
  mart = ensembl
) |> as_tibble() |> 
dplyr::rename(gene_id = ensembl_gene_id, gene_symbol = hgnc_symbol)

counts_matrix <- column_to_rownames(counts_table,var = "GENE") |> as.matrix()

counts_long <- counts_table |>
                pivot_longer(cols = -GENE,
                             values_to = "counts",
                             names_to = "sample")


controls <- c("Mock","MOCK", "NT1","NT2","NT3","Ruby1","RUBY1", "Ruby2","RUBY2","Ruby3","RUBY3","Ruby4","RUBY4")

sample_metadata <- counts_long |> 
  dplyr::select(sample) |> 
  distinct() |>
  separate(sample, sep = "_", into = c("Line", "Gene", "Guide", "Replicate"), remove = FALSE) |> 
  mutate(
    is_control = Gene %in% controls,
    original_gene = Gene,
    Replicate = if_else(is_control, Guide, Replicate),
    Guide = case_when(
      !is_control ~ Guide,
      original_gene == "Mock|MOCK" ~ NA_character_,
      str_detect(original_gene, "^NT|^Ruby|^RUBY") ~ str_replace(original_gene, "RUBY", "Ruby")
    ),
    Pertubation = case_when(
      !is_control ~ paste0(Gene, "_", Guide),
       TRUE  ~ Guide
    ),
    Gene = if_else(is_control, "Mock", Gene)
  ) |>
  dplyr::select(-original_gene, -is_control)  |> 
  mutate(Cas13 = case_when(str_detect(Line, pattern = "CAS") ~ "Positive",
                          str_detect(Line, pattern = "WT") ~ "Negative", 
                          TRUE ~ NA_character_),
        Cell_line = case_when(str_detect(Line, pattern = "99") ~ "H1299",
                              str_detect(Line, pattern = "75") ~ "H1975",
                              str_detect(Line, pattern = "^A") ~ "A549",
                              TRUE ~ NA_character_
                              
                              )
                ) 

sample_metadata |> 
write_tsv("results/tables/sample_metadata.tsv")


counts_long <- counts_long |> 
left_join(sample_metadata, by = "sample")

png("results/figures/counts_distribution.png")
    ggplot(counts_long, aes(x=log2(counts+1))) +
    geom_density(aes(colour = Cell_line), alpha = 0.1) +
    theme(legend.position = "none")
dev.off()



png("results/figures/counts_distribution_trimmed.png")
    ggplot(counts_long |> filter(counts > 10), aes(x=log2(counts+1))) +
    geom_density(aes(colour = Cell_line), alpha = 0.1) +
    theme(legend.position = "none")
dev.off()


# dds <- DESeqDataSetFromMatrix(
# countData = counts_matrix,
# colData = sample_metadata,
# design = ~1
# )
# min_count <- 10 
# min_samples <- 3
# keep <- rowSums(counts(dds) >= min_count) >= min_samples
# dds <- dds[keep, ]

# vsd <- vst(dds, blind = TRUE)


# pca_data <- plotPCA(vsd, intgroup = "Gene", returnData = FALSE)


# --- 1) PCA on raw counts (log-transformed) ---
# Filter low counts and log-transform for PCA
keep <- rowSums(counts_matrix > 10) >= 3
counts_filtered <- counts_matrix[keep, ]

# Log2 CPM for visualisation (add pseudocount)
log_cpm_raw <- log2(cpm(counts_filtered) + 1)

# PCA
pca_raw <- prcomp(t(log_cpm_raw), scale. = TRUE, center = TRUE)

pca_raw_df <- pca_raw$x |> 
  as.data.frame() |> 
  rownames_to_column("sample") |> 
  left_join(sample_metadata, by = "sample")

# Variance explained
var_explained <- summary(pca_raw)$importance[2, ] * 100

ggplot(pca_raw_df, aes(x = PC1, y = PC2, colour = Cell_line, shape = Cas13)) +
  geom_point(size = 3, alpha = 0.8) +
  labs(
    x = paste0("PC1 (", round(var_explained[1], 1), "%)"),
    y = paste0("PC2 (", round(var_explained[2], 1), "%)"),
    title = "PCA - Raw counts (log2 CPM)"
  ) +
  theme_bw()

ggsave("pca_raw.pdf", width = 8, height = 6)


# --- 2) Correlation heatmap with ComplexHeatmap ---
cor_matrix <- cor(log_cpm_raw, method = "pearson")

# Annotation for heatmap
ha <- HeatmapAnnotation(
  Cell_line = sample_metadata$Cell_line,
  Gene = sample_metadata$Gene,
  Cas13 = sample_metadata$Cas13,
  Guide = sample_metadata$Pertubation,
  col = list(
    Cell_line = c("A549" = "#E41A1C", "H1975" = "#377EB8", "H1299" = "#4DAF4A"),
    Cas13 = c("Positive" = "#984EA3", "Negative" = "#FF7F00")
  )
)

# Colour scale for correlation
col_fun <- colorRamp2(c(0.85, 0.925, 1), c("blue", "white", "red"))

pdf("results/figures/correlation_heatmap.pdf", width = 20, height = 20)
Heatmap(
  cor_matrix,
  name = "Pearson r",
  col = col_fun,
  top_annotation = ha,
  show_row_names = FALSE,
  show_column_names = FALSE,
  clustering_method_rows = "ward.D2",
  clustering_method_columns = "ward.D2"
)
dev.off()


# --- 3) edgeR + TMM normalisation ---
dge <- DGEList(counts = counts_matrix, samples = sample_metadata)

# Filter low expression genes
keep <- filterByExpr(dge, group = dge$samples$Cell_line)
dge <- dge[keep, , keep.lib.sizes = FALSE]

# TMM normalisation
dge <- calcNormFactors(dge, method = "TMM")

# Get normalised log2 CPM
log_cpm_norm <- cpm(dge, log = TRUE, prior.count = 1)


# --- 4) PCA on normalised data ---
pca_norm <- prcomp(t(log_cpm_norm), scale. = TRUE, center = TRUE)

pca_norm_df <- pca_norm$x |> 
  as.data.frame() |> 
  rownames_to_column("sample") |> 
  left_join(sample_metadata, by = "sample")

var_explained_norm <- summary(pca_norm)$importance[2, ] * 100

# Colour by cell line, shape by Cas13
p1 <- ggplot(pca_norm_df, aes(x = PC1, y = PC2, colour = Cell_line, shape = Cas13)) +
  geom_point(size = 3, alpha = 0.8) +
  labs(
    x = paste0("PC1 (", round(var_explained_norm[1], 1), "%)"),
    y = paste0("PC2 (", round(var_explained_norm[2], 1), "%)"),
    title = "PCA - TMM normalised (log2 CPM)"
  ) +
  theme_bw()

# Alternative: colour by Gene target
p2 <- ggplot(pca_norm_df, aes(x = PC1, y = PC2, colour = Gene, shape = Cell_line)) +
  geom_point(size = 3, alpha = 0.7) +
  labs(
    x = paste0("PC1 (", round(var_explained_norm[1], 1), "%)"),
    y = paste0("PC2 (", round(var_explained_norm[2], 1), "%)"),
    title = "PCA - TMM normalised (log2 CPM)"
  ) +
  theme_bw()

ggsave("results/figures/pca_normalised_cellline.pdf", p1, width = 8, height = 6)
ggsave("results/figures/pca_normalised_gene.pdf", p2, width = 10, height = 6)


# --- Quick MDS plot (edgeR's built-in, good sanity check) ---
pdf("mds_plot.pdf", width = 8, height = 6)
plotMDS(dge, col = as.numeric(factor(dge$samples$Cell_line)), 
        pch = ifelse(dge$samples$Cas13 == "Positive", 16, 1))
legend("topright", legend = levels(factor(dge$samples$Cell_line)), 
       col = 1:3, pch = 16)
dev.off()



sample_metadata <- sample_metadata |>
mutate(Pertubation = replace_na(Pertubation, "NoGuide")) |>
mutate(Group = interaction(Cell_line, Pertubation, Cas13, drop = TRUE, sep = "_")) 




# --- Build the model ---
design <- model.matrix(~ 0 + Group, data = sample_metadata)
colnames(design) <- str_remove(colnames(design), "^Group")

# Create DGEList and normalise
dge <- DGEList(counts = counts_matrix[, sample_metadata$sample], 
               samples = sample_metadata)
keep <- filterByExpr(dge, design = design)
dge <- dge[keep, , keep.lib.sizes = FALSE]
dge <- calcNormFactors(dge, method = "TMM")

# Voom transformation
v <- voom(dge, design, plot = TRUE)

# Fit the model
fit <- lmFit(v, design)


# --- Define contrasts ---
# Let's see what coefficients we have
colnames(design)

# Question 1: Mock_NoGuide vs Mock with NT guides (within Cas13+)
# We'll do this per cell line first, then can meta-analyse

# Build contrast strings programmatically: NoGuide vs each NT and NT average per cell line
q1_contrast_strings <- unlist(map(c("A549", "H1299", "H1975"), function(cl) {
  nts <- paste0("NT", 1:3)
  nt_avg <- sprintf("(%s_NT1_Positive + %s_NT2_Positive + %s_NT3_Positive)/3", cl, cl, cl)
  c(
    # Each NT vs NoGuide
    setNames(
      sprintf("%s_%s_Positive - %s_NoGuide_Positive", cl, nts, cl),
      sprintf("%s_NoGuide_vs_%s_Pos", cl, nts)
    ),
    # NT average vs NoGuide
    setNames(
      sprintf("%s - %s_NoGuide_Positive", nt_avg, cl),
      sprintf("%s_NoGuide_vs_NTavg_Pos", cl)
    )
  )
}))

contrasts_q1 <- makeContrasts(contrasts = q1_contrast_strings, levels = design)
colnames(contrasts_q1) <- names(q1_contrast_strings)


# Question 2: NT controls vs Gene targets (Cas13+ only, where knockdown happens)
# The key contrast: Target_Cas13+ vs Control_Cas13+
# Using average of NTs as the control baseline

# Build contrast strings programmatically from design matrix columns
# Automatically picks up all target/version combos per cell line
control_patterns <- c("NT1", "NT2", "NT3", "NoGuide", "Ruby")
design_cols <- colnames(design)

q2_contrast_strings <- unlist(map(set_names(c("A549", "H1299", "H1975")), function(cl) {
  nt_avg <- sprintf("(%s_NT1_Positive + %s_NT2_Positive + %s_NT3_Positive)/3", cl, cl, cl)

  # Find all Cas13+ target columns for this cell line, excluding controls
  cl_pos_cols <- design_cols[grepl(sprintf("^%s_.*_Positive$", cl), design_cols)]
  target_cols <- cl_pos_cols[!grepl(paste(control_patterns, collapse = "|"), cl_pos_cols)]

  # Name each contrast: strip _Positive suffix, add _vs_NT
  contrast_names <- paste0(str_remove(target_cols, "_Positive$"), "_vs_NT")
  setNames(sprintf("%s - %s", target_cols, nt_avg), contrast_names)
}))

contrasts_q2 <- makeContrasts(contrasts = q2_contrast_strings, levels = design)
colnames(contrasts_q2) <- names(q2_contrast_strings)


# Question 3: Gene targets vs Mock (NoGuide) per cell line (Cas13+ only)
q3_contrast_strings <- unlist(map(set_names(c("A549", "H1299", "H1975")), function(cl) {
  mock_col <- sprintf("%s_NoGuide_Positive", cl)

  # Find all Cas13+ target columns for this cell line, excluding controls
  cl_pos_cols <- design_cols[grepl(sprintf("^%s_.*_Positive$", cl), design_cols)]
  target_cols <- cl_pos_cols[!grepl(paste(control_patterns, collapse = "|"), cl_pos_cols)]

  contrast_names <- paste0(str_remove(target_cols, "_Positive$"), "_vs_Mock")
  setNames(sprintf("%s - %s", target_cols, mock_col), contrast_names)
}))

contrasts_q3 <- makeContrasts(contrasts = q3_contrast_strings, levels = design)
colnames(contrasts_q3) <- names(q3_contrast_strings)


# --- Fit contrasts and get results ---
fit_q1 <- contrasts.fit(fit, contrasts_q1)
fit_q1 <- eBayes(fit_q1)

fit_q2 <- contrasts.fit(fit, contrasts_q2)
fit_q2 <- eBayes(fit_q2)

fit_q3 <- contrasts.fit(fit, contrasts_q3)
fit_q3 <- eBayes(fit_q3)

# Summary of significant genes
summary(decideTests(fit_q1, p.value = 0.01, lfc =2))
summary(decideTests(fit_q2, p.value = 0.01, lfc =2))
summary(decideTests(fit_q3, p.value = 0.01, lfc =2))

# --- Volcano plots: one per target, faceted by cell line x guide version ---

# Extract all q2 results into a single data frame, with gene symbol annotations
q2_results <- map_dfr(colnames(fit_q2$coefficients), function(coef) {
  topTable(fit_q2, coef = coef, number = Inf) |>
    rownames_to_column("gene_id") |>
    mutate(contrast = coef)
}) |>
  left_join(gene_annot, by = "gene_id")

# Parse contrast names to extract cell line, target, and version
q2_results <- q2_results |>
  mutate(
    cell_line = str_extract(contrast, "^[^_]+"),
    perturbation = str_remove(contrast, "^[^_]+_") |> str_remove("_vs_NT$"),
    target = str_remove(perturbation, "[._][^._]+$"),
    version = str_extract(perturbation, "[^._]+$")
  ) |>
  mutate(
    sig = case_when(
      adj.P.Val < 0.05 & logFC > 1  ~ "Up",
      adj.P.Val < 0.05 & logFC < -1 ~ "Down",
      TRUE ~ "NS"
    )
  )

# Map short target names in contrasts to full gene symbols for labelling
target_to_symbol <- c(
  MALAT = "MALAT1", NEAT = "NEAT1", LINC = "LINC-PINT",
  FBL = "FBL", FTL = "FTL", VIM = "VIM",
  SNHG15 = "SNHG15", SNHG16 = "SNHG16", SNHG6 = "SNHG6",
  PP = "PPIA"
)

# One volcano per target, with the perturbed gene labelled
targets <- unique(q2_results$target)

volcano_plots <- map(set_names(targets), function(tgt) {
  df <- q2_results |> filter(target == tgt)

  # Find the target gene symbol to label
  tgt_symbol <- target_to_symbol[[tgt]] %||% tgt
  df <- df |> mutate(is_target = str_detect(gene_symbol, fixed(tgt_symbol)))

  ggplot(df, aes(x = logFC, y = -log10(adj.P.Val), colour = sig)) +
    geom_point(size = 0.5, alpha = 0.5) +
    geom_point(data = df |> filter(is_target), size = 2, shape = 21,
               fill = "yellow", colour = "black", stroke = 0.5) +
    geom_text_repel(
      data = df |> filter(is_target),
      aes(label = gene_symbol), colour = "black",
      size = 3, fontface = "bold", max.overlaps = 20,
      box.padding = 0.5, min.segment.length = 0
    ) +
    scale_colour_manual(values = c(Down = "blue", NS = "grey70", Up = "red")) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "grey40") +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", colour = "grey40") +
    facet_grid(cell_line ~ version) +
    labs(
      title = sprintf("Volcano: %s vs NT average", tgt),
      x = "log2 Fold Change",
      y = "-log10(adj. P)"
    ) +
    theme_bw() +
    theme(legend.position = "bottom")
})

# Save each target as a separate PDF
dir.create("results/volcanos", showWarnings = FALSE, recursive = TRUE)
iwalk(volcano_plots, function(p, tgt) {
  df <- q2_results |> filter(target == tgt)
  n_versions <- n_distinct(df$version)
  n_celllines <- n_distinct(df$cell_line)
  ggsave(
    sprintf("results/volcanos/volcano_%s.pdf", tgt),
    p, width = 3 + 3 * n_versions, height = 2 + 3 * n_celllines
  )
})


# --- Q3 Volcano plots: knockdowns vs Mock (NoGuide) ---

q3_results <- map_dfr(colnames(fit_q3$coefficients), function(coef) {
  topTable(fit_q3, coef = coef, number = Inf) |>
    rownames_to_column("gene_id") |>
    mutate(contrast = coef)
}) |>
  left_join(gene_annot, by = "gene_id")

q3_results <- q3_results |>
  mutate(
    cell_line = str_extract(contrast, "^[^_]+"),
    perturbation = str_remove(contrast, "^[^_]+_") |> str_remove("_vs_Mock$"),
    target = str_remove(perturbation, "[._][^._]+$"),
    version = str_extract(perturbation, "[^._]+$")
  ) |>
  mutate(
    sig = case_when(
      adj.P.Val < 0.05 & logFC > 1  ~ "Up",
      adj.P.Val < 0.05 & logFC < -1 ~ "Down",
      TRUE ~ "NS"
    )
  )

targets_q3 <- unique(q3_results$target)

volcano_plots_q3 <- map(set_names(targets_q3), function(tgt) {
  df <- q3_results |> filter(target == tgt)

  tgt_symbol <- target_to_symbol[[tgt]] %||% tgt
  df <- df |> mutate(is_target = str_detect(gene_symbol, fixed(tgt_symbol)))

  ggplot(df, aes(x = logFC, y = -log10(adj.P.Val), colour = sig)) +
    geom_point(size = 0.5, alpha = 0.5) +
    geom_point(data = df |> filter(is_target), size = 2, shape = 21,
               fill = "yellow", colour = "black", stroke = 0.5) +
    geom_text_repel(
      data = df |> filter(is_target),
      aes(label = gene_symbol), colour = "black",
      size = 3, fontface = "bold", max.overlaps = 20,
      box.padding = 0.5, min.segment.length = 0
    ) +
    scale_colour_manual(values = c(Down = "blue", NS = "grey70", Up = "red")) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "grey40") +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", colour = "grey40") +
    facet_grid(cell_line ~ version) +
    labs(
      title = sprintf("Volcano: %s vs Mock", tgt),
      x = "log2 Fold Change",
      y = "-log10(adj. P)"
    ) +
    theme_bw() +
    theme(legend.position = "bottom")
})

dir.create("results/volcanos_vs_mock", showWarnings = FALSE, recursive = TRUE)
iwalk(volcano_plots_q3, function(p, tgt) {
  df <- q3_results |> filter(target == tgt)
  n_versions <- n_distinct(df$version)
  n_celllines <- n_distinct(df$cell_line)
  ggsave(
    sprintf("results/volcanos_vs_mock/volcano_%s_vs_mock.pdf", tgt),
    p, width = 3 + 3 * n_versions, height = 2 + 3 * n_celllines
  )
})