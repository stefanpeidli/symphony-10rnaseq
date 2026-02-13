library(DESeq2)
library(tidyverse)

counts_table <- read_tsv("results/8094_Sym101_bulk_counts_matrix.tsv") 


counts_long <- counts_table |>
                pivot_longer(cols = -GENE,
                             values_to = "counts",
                             names_to = "sample")



# sample_metadata <- counts_long |> 
# select(sample) |> 
# distinct() |>
# separate(sample, sep = "_", into = c("Line", "Gene", "Guide", "Replicate")) |> 
# filter(Gene %in% controls) |>
#   mutate(
#     original_gene = Gene,  # save original for reference
#     Replicate = Guide,     # current Guide becomes Replicate
#     Guide = case_when(
#       original_gene == "Mock" ~ NA_character_,
#       str_detect(original_gene, "^NT") ~ original_gene,  # NT1, NT2, NT3
#       str_detect(original_gene, "Ruby") ~ original_gene,
#       str_detect(original_gene, "RUBY") ~ original_gene,
#       TRUE ~ NA_character_
#     ),
#     Gene = case_when(
#       original_gene == "Mock" ~ "Mock",
#       str_detect(original_gene, "^NT") ~ "Mock",
#       str_detect(original_gene, "^Ruby") ~ "Mock",
#       TRUE ~ original_gene
#     )
#   ) |>
#   select(-original_gene) 




controls <- c("Mock","NT1","NT2","NT3","Ruby1","RUBY1", "Ruby2","RUBY2","Ruby3","RUBY3","Ruby4","RUBY4")

sample_metadata <- counts_long |> 
  select(sample) |> 
  distinct() |>
  separate(sample, sep = "_", into = c("Line", "Gene", "Guide", "Replicate")) |> 
  mutate(
    is_control = Gene %in% controls,
    original_gene = Gene,
    Replicate = if_else(is_control, Guide, Replicate),
    Guide = case_when(
      !is_control ~ Guide,
      original_gene == "Mock" ~ NA_character_,
      str_detect(original_gene, "^NT|^Ruby|^RUBY") ~ str_replace(original_gene, "RUBY", "Ruby")
    ),
    Gene = if_else(is_control, "Mock", Gene)
  ) |>
  select(-original_gene, -is_control)  |> 
  mutate(Cas13 = case_when(str_detect(Line, pattern = "CAS") ~ "Positive",
                          str_detect(Line, pattern = "WT") ~ "Negative", 
                          TRUE ~ NA_character_),
        Cell_line = case_when(str_detect(Line, pattern = "99") ~ "H1299",
                              str_detect(Line, pattern = "75") ~ "H1975",
                              str_detect(Line, pattern = "^A") ~ "A549",
                              TRUE ~ NA_character_
                              
                              )
                ) 

sample_metadata |> write_tsv("sample_metadata.tsv")


png("counts_distribution.png")
    ggplot(counts_long, aes(x=log2(counts+1))) +
    geom_density(aes(colour = sample), alpha = 0.1) +
    theme(legend.position = "none")
dev.off()


    
png("counts_distribution_trimmed.png")
    ggplot(counts_long |> filter(counts > 10), aes(x=log2(counts+1))) +
    geom_density(aes(colour = sample), alpha = 0.1) +
    theme(legend.position = "none")
dev.off()


    
png("counts_boxplots.png")
    ggplot(counts_long, aes(x = sample, y=log2(counts+1))) +
    geom_boxplot(aes(fill = sample), alpha = 0.1) +
    theme(legend.position = "none")
dev.off()



colData <- read_tsv("results/coldata.tsv") |> 
            mutate(pertubation = as.factor(perturbation), 
                    guide = as.factor(guide)) |> 
                    dplyr::rename("SAMPLE" = "...1")


dds <- DESeqDataSetFromMatrix(
countData = counts_mat,
colData = colData,
design = ~1
)
min_count <- 10 
min_samples <- 3
keep <- rowSums(counts(dds) >= min_count) >= min_samples
dds <- dds[keep, ]

vsd <- vst(dds, blind = TRUE)


pca_data <- plotPCA(vsd, intgroup = "sample", returnData = FALSE)
# pca_coords <- pca_data[, c("PC1", "PC2")]
# rownames(pca_coords) <- pca_data[["name"]]
# percentVar <- round(100 * attr(pca_data, "percentVar"))

png("pca_plot.png", width = 3000)
print(pca_data)
dev.off()