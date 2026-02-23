library(ggVennDiagram)
library(patchwork)

get_down_genes <- function(fit, coef, p_thresh = 0.05, lfc_thresh = 0) {
  tt <- topTable(fit, coef = coef, number = Inf)
  tt %>%
    rownames_to_column("gene_id") %>%
    filter(adj.P.Val < p_thresh, logFC < -lfc_thresh) %>%
    pull(gene_id)
}

cell_lines <- c("A549", "H1299", "H1975")
nts <- paste0("NT", 1:3)
venn_colours <- c(A549 = "#E41A1C", H1299 = "#4DAF4A", H1975 = "#377EB8")

# Get down genes for each NT contrast per cell line
down_by_nt <- map(set_names(cell_lines), function(cl) {
  coefs <- set_names(sprintf("%s_NoGuide_vs_%s_Pos", cl, nts), nts)
  map(coefs, ~ get_down_genes(fit_q1, .x, p_thresh = 0.01, lfc_thresh = 2))
})

# Get down genes for NT average contrast per cell line
down_ntavg <- map(set_names(cell_lines), function(cl) {
  get_down_genes(fit_q1, sprintf("%s_NoGuide_vs_NTavg_Pos", cl),  p_thresh = 0.01, lfc_thresh = 2)
})

# --- Venn diagrams: NTs within each cell line ---
venn_within_plots <- imap(down_by_nt, function(gene_lists, cl) {
  ggVennDiagram(gene_lists, label_alpha = 0) +
    scale_fill_gradient(low = "white", high = venn_colours[[cl]]) +
    labs(title = sprintf("%s: Depleted genes (NoGuide vs NTs)", cl)) +
    theme(legend.position = "none")
})

# --- Venn diagram: NT averages across cell lines ---
p_celllines <- ggVennDiagram(down_ntavg, label_alpha = 0) +
  scale_fill_gradient(low = "white", high = "#984EA3") +
  labs(title = "Depleted genes: NT average across cell lines") +
  theme(legend.position = "none")

# --- Combine and save ---
p_within <- wrap_plots(venn_within_plots) +
  plot_annotation(title = "Within cell line: consistency across NT guides")

ggsave("venn_within_cellline_NTs.pdf", p_within, width = 14, height = 5)
ggsave("venn_across_celllines_NTavg.pdf", p_celllines, width = 6, height = 5)

# --- Venn diagrams: each NT across cell lines ---
# Transpose down_by_nt from cell_line -> NT to NT -> cell_line
nt_colours <- c(NT1 = "#FF7F00", NT2 = "#6A3D9A", NT3 = "#A65628")

down_by_nt_across <- map(set_names(nts), function(nt) {
  map(set_names(cell_lines), ~ down_by_nt[[.x]][[nt]])
})

venn_across_plots <- imap(down_by_nt_across, function(gene_lists, nt) {
  ggVennDiagram(gene_lists, label_alpha = 0) +
    scale_fill_gradient(low = "white", high = nt_colours[[nt]]) +
    labs(title = sprintf("%s: DEGs shared across cell lines", nt)) +
    theme(legend.position = "none")
})

p_across_nts <- wrap_plots(venn_across_plots) +
  plot_annotation(title = "Across cell lines: consistency per NT guide")

ggsave("venn_across_celllines_perNT.pdf", p_across_nts, width = 14, height = 5)

