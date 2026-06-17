# =============================================================================
# Plestiodon Selection Analysis
# Tajima's D outlier detection, gene annotation overlap, FST summary, GO enrichment
# Uses: pop_map_ALL.thetas.windowed.pestPG + ragtag.scaffold.agp + GFF3
# =============================================================================

# --- 1. INSTALL / LOAD PACKAGES ----------------------------------------------

packages <- c("ggplot2", "dplyr", "ggrepel", "rtracklayer", "GenomicRanges",
              "data.table", "tidyr", "stringr", "conflicted")

installed <- rownames(installed.packages())
to_install <- packages[!packages %in% installed]
if (length(to_install)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  bioc_pkgs <- c("rtracklayer", "GenomicRanges")
  cran_pkgs <- setdiff(to_install, bioc_pkgs)
  if (length(cran_pkgs)) install.packages(cran_pkgs)
  if (any(bioc_pkgs %in% to_install)) BiocManager::install(bioc_pkgs[bioc_pkgs %in% to_install])
}

lapply(packages, library, character.only = TRUE)

# Resolve namespace conflicts explicitly
conflicts_prefer(dplyr::select)
conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::slice)
conflicts_prefer(dplyr::rename)
conflicts_prefer(dplyr::mutate)
conflicts_prefer(dplyr::first)
conflicts_prefer(dplyr::lag)

# Bioconductor packages
if (!requireNamespace("GO.db", quietly = TRUE)) BiocManager::install("GO.db")
if (!requireNamespace("clusterProfiler", quietly = TRUE)) BiocManager::install("clusterProfiler")
library(GO.db)
library(clusterProfiler)


# --- 2. FILE PATHS -----------------------------------------------------------

base_dir <- file.path(Sys.getenv("USERPROFILE"),
                      "OneDrive - MMU/Documents/PhD Analysis/Bioinformatics/Physalia")

THETAS_FILE <- file.path(base_dir, "selection pressures",
                         "pop_map_ALL.thetas.windowed.pestPG")
AGP_FILE    <- file.path(base_dir, "selection pressures",
                         "ragtag.scaffold.agp")
GFF_FILE    <- file.path(base_dir, "selection pressures",
                         "Plestiodon_longirostris_polished.gff3")
FST_DIR     <- file.path(base_dir, "fst")
OUT_DIR     <- file.path(base_dir, "results", "selection difference")

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)


# --- 3. LOAD TAJIMA'S D DATA -------------------------------------------------

cat("Loading Tajima's D data...\n")

thetas_raw <- fread(THETAS_FILE)

colnames(thetas_raw) <- c("window_info", "Chr", "WinCenter",
                          "tW", "tP", "tF", "tH", "tL",
                          "Tajima", "fuf", "fud", "fayh", "zeng", "nSites")

thetas <- thetas_raw %>%
  dplyr::filter(!is.na(Tajima), !is.nan(Tajima), nSites >= 10) %>%
  dplyr::select(Chr, WinCenter, Tajima, tW, tP, nSites)

cat("Total windows after filtering:", nrow(thetas), "\n")
cat("Unique contigs:", length(unique(thetas$Chr)), "\n")


# --- 4. LOAD AGP FILE AND BUILD CONTIG -> CHROMOSOME MAPPING ----------------

cat("\nLoading AGP mapping file...\n")

agp_raw <- fread(AGP_FILE, skip = 2, header = FALSE,
                 col.names = c("chr", "chr_start", "chr_end", "part_num",
                               "comp_type", "contig", "contig_start",
                               "contig_end", "orientation"))

agp <- as.data.frame(agp_raw) %>%
  dplyr::filter(comp_type == "W") %>%
  dplyr::select(chr, chr_start, chr_end, contig, contig_start, contig_end, orientation) %>%
  dplyr::mutate(
    contig_start = as.numeric(contig_start),
    contig_end   = as.numeric(contig_end)
  )

cat("AGP contigs mapped:", length(unique(agp$contig)), "\n")


# --- 5. TRANSLATE CONTIG COORDINATES TO CHROMOSOME COORDINATES --------------

cat("Translating contig coordinates to chromosome coordinates...\n")

contig_to_chr <- agp %>%
  dplyr::group_by(contig) %>%
  dplyr::summarise(
    chr         = first(chr),
    chr_start   = first(chr_start),
    chr_end     = first(chr_end),
    orientation = first(orientation)
  ) %>%
  dplyr::ungroup()

tajima_chr <- thetas %>%
  dplyr::filter(Chr %in% contig_to_chr$contig) %>%
  dplyr::left_join(contig_to_chr, by = c("Chr" = "contig")) %>%
  dplyr::mutate(
    CHROM = chr,
    POS   = as.integer((chr_start + chr_end) / 2)
  ) %>%
  dplyr::select(CHROM, POS, Tajima, tW, tP, nSites, orig_contig = Chr)

cat("Windows successfully mapped to chromosomes:", nrow(tajima_chr), "\n")
cat("Windows lost (unplaced contigs):", nrow(thetas) - nrow(tajima_chr), "\n")


# --- 6. OUTLIER DETECTION (BOTH TAILS) ---------------------------------------

cat("\nDetecting outliers...\n")

lower_threshold <- quantile(tajima_chr$Tajima, 0.01, na.rm = TRUE)
upper_threshold <- quantile(tajima_chr$Tajima, 0.99, na.rm = TRUE)

cat("Lower threshold (1st percentile):", round(lower_threshold, 4), "\n")
cat("Upper threshold (99th percentile):", round(upper_threshold, 4), "\n")

tajima_chr <- tajima_chr %>%
  dplyr::mutate(
    outlier = Tajima <= lower_threshold | Tajima >= upper_threshold,
    outlier_type = dplyr::case_when(
      Tajima <= lower_threshold ~ "Positive selection",
      Tajima >= upper_threshold ~ "Balancing selection",
      TRUE ~ "Neutral"
    )
  )

cat("Outlier windows (positive selection):", sum(tajima_chr$outlier_type == "Positive selection"), "\n")
cat("Outlier windows (balancing selection):", sum(tajima_chr$outlier_type == "Balancing selection"), "\n")


# --- 7. PREPARE CHROMOSOME COORDINATES FOR MANHATTAN PLOT -------------------

chr_order <- unique(tajima_chr$CHROM)
chr_order <- chr_order[order(as.numeric(gsub("\\D", "", chr_order)))]
tajima_chr$CHROM <- factor(tajima_chr$CHROM, levels = chr_order)

scaffold_lengths <- tajima_chr %>%
  dplyr::group_by(CHROM) %>%
  dplyr::summarise(max_pos = max(POS)) %>%
  dplyr::arrange(CHROM) %>%
  dplyr::mutate(
    offset    = lag(cumsum(as.numeric(max_pos)), default = 0),
    cum_start = offset
  )

tajima_chr <- tajima_chr %>%
  dplyr::left_join(scaffold_lengths %>% dplyr::select(CHROM, cum_start), by = "CHROM") %>%
  dplyr::mutate(cum_pos = POS + cum_start)

scaffold_mids <- tajima_chr %>%
  dplyr::group_by(CHROM) %>%
  dplyr::summarise(mid = (min(cum_pos) + max(cum_pos)) / 2)


# --- 8. LOAD GENE ANNOTATION -------------------------------------------------

cat("\nLoading GFF3 annotation...\n")

gff    <- import(GFF_FILE, format = "gff3")
gff_df <- as.data.frame(gff)

genes <- gff_df %>%
  dplyr::filter(type == "gene") %>%
  dplyr::mutate(
    product_chr = sapply(product, function(x) if (length(x) > 0) x[[1]] else NA_character_),
    gene_id = dplyr::coalesce(ID, Name, product_chr)
  ) %>%
  dplyr::select(seqnames, start, end, gene_id) %>%
  dplyr::rename(CHROM = seqnames) %>%
  dplyr::mutate(CHROM = as.character(CHROM))

cat("Total annotated genes:", nrow(genes), "\n")


# --- 9. FIND GENES OVERLAPPING OUTLIER WINDOWS -------------------------------

cat("Finding gene overlaps...\n")

outlier_windows <- tajima_chr %>% dplyr::filter(outlier)

WINDOW_SIZE <- 50000
outlier_gr <- GRanges(
  seqnames = outlier_windows$CHROM,
  ranges   = IRanges(
    start = pmax(1, outlier_windows$POS - WINDOW_SIZE / 2),
    end   = outlier_windows$POS + WINDOW_SIZE / 2
  ),
  Tajima       = outlier_windows$Tajima,
  outlier_type = outlier_windows$outlier_type
)

gene_gr <- GRanges(
  seqnames = genes$CHROM,
  ranges   = IRanges(start = genes$start, end = genes$end),
  gene_id  = genes$gene_id
)

# Harmonise sequence naming: AGP uses chr1_RagTag, GFF3 uses chr1_RagTag_pilon
if (!any(seqlevels(outlier_gr) %in% seqlevels(gene_gr))) {
  cat("Sequence names differ between ANGSD and GFF3 - attempting to fix...\n")
  seqlevels(outlier_gr) <- paste0(seqlevels(outlier_gr), "_pilon")
}

hits <- findOverlaps(outlier_gr, gene_gr)

overlapping_genes <- data.frame(
  CHROM        = as.character(seqnames(outlier_gr))[queryHits(hits)],
  WIN_CENTER   = outlier_windows$POS[queryHits(hits)],
  Tajima       = outlier_gr$Tajima[queryHits(hits)],
  outlier_type = outlier_gr$outlier_type[queryHits(hits)],
  gene_id      = gene_gr$gene_id[subjectHits(hits)]
) %>% dplyr::distinct()

cat("Outlier windows overlapping annotated genes:", nrow(overlapping_genes), "\n")
cat("Unique genes under selection:", length(unique(overlapping_genes$gene_id)), "\n")

write.csv(overlapping_genes,
          file.path(OUT_DIR, "plestiodon_outlier_genes_TajimaD.csv"),
          row.names = FALSE)


# --- 10. LOAD FST SUMMARY FROM GLOBAL FILES ----------------------------------

cat("\nCompiling population FST summary...\n")

fst_global_files <- list.files(FST_DIR,
                               pattern = "fst_pop_map.*\\.fst\\.global\\.txt",
                               full.names = TRUE,
                               recursive = TRUE)

cat("FST files found:", length(fst_global_files), "\n")

parse_fst_global <- function(f) {
  lines <- readLines(f)
  fst_line <- lines[grepl("FST.Unweight", lines)]
  if (length(fst_line) == 0) return(NULL)
  
  unweighted <- as.numeric(gsub(".*FST.Unweight\\[nObs:(\\d+)\\]:([-0-9.]+).*", "\\2", fst_line))
  weighted   <- as.numeric(gsub(".*Fst.Weight:([-0-9.]+).*", "\\1", fst_line))
  nobs       <- as.numeric(gsub(".*nObs:(\\d+).*", "\\1", fst_line))
  label      <- gsub(".*fst_(pop_map_)?(.+)\\.fst\\.global\\.txt", "\\2", f)
  label      <- gsub("_pop_map_", " vs ", label)
  
  data.frame(comparison = label, FST_unweighted = unweighted,
             FST_weighted = weighted, nSites = nobs)
}

if (length(fst_global_files) > 0) {
  fst_summary <- dplyr::bind_rows(lapply(fst_global_files, parse_fst_global)) %>%
    dplyr::arrange(desc(FST_weighted))
  cat("Population pair FST summary:\n")
  print(fst_summary)
  write.csv(fst_summary,
            file.path(OUT_DIR, "plestiodon_population_FST_summary.csv"),
            row.names = FALSE)
} else {
  cat("No FST global files found - skipping FST summary.\n")
}


# --- 11. MANHATTAN PLOT ------------------------------------------------------

cat("\nGenerating Manhattan plot...\n")

chr_levels <- levels(tajima_chr$CHROM)
tajima_chr <- tajima_chr %>%
  dplyr::mutate(colour_group = factor(match(CHROM, chr_levels) %% 2))

# Build gene labels per chromosome (top 3 genes from most extreme outlier window)
label_df <- overlapping_genes %>%
  dplyr::mutate(CHROM_label = gsub("_pilon$", "", CHROM)) %>%
  dplyr::group_by(CHROM_label) %>%
  dplyr::summarise(
    gene_label = paste(unique(gene_id)[1:min(3, length(unique(gene_id)))], collapse = "/"),
    .groups = "drop"
  )

tajima_chr <- tajima_chr %>%
  dplyr::select(-dplyr::any_of(c("gene_label", "gene_label.x", "gene_label.y", "CHROM_char"))) %>%
  dplyr::mutate(CHROM_char = as.character(CHROM)) %>%
  dplyr::left_join(label_df, by = c("CHROM_char" = "CHROM_label"))

# Only label the single most extreme outlier window per chromosome
tajima_chr <- tajima_chr %>%
  dplyr::group_by(CHROM_char) %>%
  dplyr::mutate(gene_label = ifelse(
    !is.na(gene_label) & outlier & abs(Tajima) == max(abs(Tajima[outlier]), na.rm = TRUE),
    gene_label, NA_character_
  )) %>%
  dplyr::ungroup()

p <- ggplot(tajima_chr, aes(x = cum_pos, y = Tajima, colour = colour_group)) +
  geom_point(size = 0.4, alpha = 0.5, show.legend = FALSE) +
  geom_point(
    data = dplyr::filter(tajima_chr, outlier_type == "Positive selection"),
    aes(x = cum_pos, y = Tajima),
    colour = "#CC0000", size = 1.4, alpha = 0.9
  ) +
  geom_point(
    data = dplyr::filter(tajima_chr, outlier_type == "Balancing selection"),
    aes(x = cum_pos, y = Tajima),
    colour = "#7B2D8B", size = 1.4, alpha = 0.9
  ) +
  geom_hline(yintercept = lower_threshold, linetype = "dashed",
             colour = "#CC0000", linewidth = 0.5) +
  geom_hline(yintercept = upper_threshold, linetype = "dashed",
             colour = "#7B2D8B", linewidth = 0.5) +
  geom_hline(yintercept = 0, linetype = "solid",
             colour = "grey60", linewidth = 0.3) +
  geom_label_repel(
    data = dplyr::filter(tajima_chr, !is.na(gene_label)),
    aes(x = cum_pos, y = Tajima, label = gene_label),
    colour = "black", size = 2.2, max.overlaps = 20,
    box.padding = 0.4, segment.colour = "grey50",
    fill = alpha("white", 0.7), show.legend = FALSE
  ) +
  scale_x_continuous(
    breaks = scaffold_mids$mid,
    labels = gsub("_RagTag.*", "", scaffold_mids$CHROM),
    expand = c(0.01, 0)
  ) +
  scale_colour_manual(values = c("0" = "#4878CF", "1" = "#87AFDA")) +
  labs(
    title    = "Genome-wide Tajima's D — Plestiodon longirostris",
    subtitle = paste0("50kb windows | Outliers: <",
                      round(lower_threshold, 3), " (positive selection, red) | >",
                      round(upper_threshold, 3), " (balancing selection, purple)"),
    x = "Genomic Position (by chromosome)",
    y = "Tajima's D"
  ) +
  theme_classic() +
  theme(
    axis.text.x        = element_text(angle = 45, hjust = 1, size = 7),
    axis.text.y        = element_text(size = 10),
    axis.title         = element_text(size = 12),
    plot.title         = element_text(face = "bold", size = 14),
    plot.subtitle      = element_text(size = 9, colour = "grey40"),
    panel.grid.major.y = element_line(colour = "grey90", linewidth = 0.3)
  )

print(p)

ggsave(file.path(OUT_DIR, "plestiodon_TajimaD_manhattan.pdf"),
       plot = p, width = 16, height = 5, dpi = 300)
ggsave(file.path(OUT_DIR, "plestiodon_TajimaD_manhattan.png"),
       plot = p, width = 16, height = 5, dpi = 300)


# --- 12. GO TERM ANNOTATION --------------------------------------------------

cat("\nExtracting GO terms from GFF3...\n")

# Extract GO terms from mRNA rows, linked to parent gene ID
go_table <- gff_df %>%
  dplyr::filter(type == "mRNA") %>%
  dplyr::mutate(
    gene_id  = sapply(Parent, function(x) x[[1]]),
    go_terms = Ontology_term
  ) %>%
  dplyr::filter(sapply(go_terms, function(x) length(x) > 0)) %>%
  dplyr::select(gene_id, go_terms) %>%
  tidyr::unnest(go_terms) %>%
  dplyr::distinct()

# Look up GO term descriptions - Biological Process only
go_info <- AnnotationDbi::select(GO.db,
                                 keys     = unique(go_table$go_terms),
                                 columns  = c("TERM", "ONTOLOGY"),
                                 keytype  = "GOID") %>%
  dplyr::filter(ONTOLOGY == "BP")

# Join GO terms to outlier genes
outlier_go <- overlapping_genes %>%
  dplyr::left_join(go_table, by = "gene_id", relationship = "many-to-many") %>%
  dplyr::left_join(go_info, by = c("go_terms" = "GOID")) %>%
  dplyr::filter(!is.na(TERM)) %>%
  dplyr::select(CHROM, WIN_CENTER, Tajima, outlier_type, gene_id,
                GO_term = go_terms, GO_description = TERM) %>%
  dplyr::distinct()

cat("Outlier genes with GO biological process terms:", length(unique(outlier_go$gene_id)), "\n")

# Summary table
go_summary <- outlier_go %>%
  dplyr::group_by(outlier_type, GO_description) %>%
  dplyr::summarise(
    n_genes = dplyr::n_distinct(gene_id),
    genes   = paste(unique(gene_id), collapse = "/"),
    .groups = "drop"
  ) %>%
  dplyr::arrange(outlier_type, desc(n_genes))

write.csv(outlier_go, file.path(OUT_DIR, "plestiodon_outlier_genes_GO.csv"),  row.names = FALSE)
write.csv(go_summary, file.path(OUT_DIR, "plestiodon_GO_summary.csv"),        row.names = FALSE)
cat("Saved: plestiodon_outlier_genes_GO.csv\n")
cat("Saved: plestiodon_GO_summary.csv\n")


# --- 13. GO ENRICHMENT TEST --------------------------------------------------

cat("\nRunning GO enrichment analysis...\n")

# Build genome-wide gene-to-GO background
gene2go <- gff_df %>%
  dplyr::filter(type == "mRNA") %>%
  dplyr::mutate(gene_id = sapply(Parent, function(x) x[[1]])) %>%
  dplyr::filter(sapply(Ontology_term, function(x) length(x) > 0)) %>%
  dplyr::select(gene_id, Ontology_term) %>%
  tidyr::unnest(Ontology_term) %>%
  dplyr::distinct() %>%
  dplyr::rename(go_id = Ontology_term)

all_genes <- unique(gene2go$gene_id)
cat("Background genes with GO terms:", length(all_genes), "\n")

for (sel_type in c("Positive selection", "Balancing selection")) {
  
  outlier_genes <- overlapping_genes %>%
    dplyr::filter(outlier_type == sel_type) %>%
    dplyr::pull(gene_id) %>%
    unique()
  
  cat("\nRunning GO enrichment for:", sel_type, "\n")
  cat("Outlier genes:", length(outlier_genes), "| Background:", length(all_genes), "\n")
  
  enrich_result <- enricher(
    gene          = outlier_genes,
    TERM2GENE     = gene2go %>% dplyr::select(go_id, gene_id),
    TERM2NAME     = go_info %>% dplyr::select(GOID, TERM) %>%
      dplyr::rename(go_id = GOID, name = TERM),
    universe      = all_genes,
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.2
  )
  
  label <- gsub(" ", "_", sel_type)
  
  if (!is.null(enrich_result) && nrow(enrich_result@result) > 0) {
    result_df <- as.data.frame(enrich_result@result) %>%
      dplyr::filter(p.adjust < 0.05)
    cat("Significant GO terms (FDR < 0.05):", nrow(result_df), "\n")
    print(result_df %>%
            dplyr::select(Description, GeneRatio, BgRatio, p.adjust) %>%
            head(20))
    write.csv(result_df,
              file.path(OUT_DIR, paste0("plestiodon_GO_enrichment_", label, ".csv")),
              row.names = FALSE)
  } else {
    cat("No significantly enriched GO terms found.\n")
  }
}

#GO Slim
if (!requireNamespace("GSEABase", quietly = TRUE)) BiocManager::install("GSEABase")
library(GSEABase)

# Download the generic GO slim OBO file
download.file("https://current.geneontology.org/ontology/subsets/goslim_generic.obo",
              destfile = file.path(OUT_DIR, "goslim_generic.obo"))

# Read the GO slim OBO file and extract slim term IDs and names
slim_lines <- readLines(file.path(OUT_DIR, "goslim_generic.obo"))

# Parse slim terms - extract ID and name pairs
term_starts <- which(slim_lines == "[Term]")
slim_terms <- data.frame(GOID = character(), slim_term = character(), stringsAsFactors = FALSE)

for (i in term_starts) {
  id_line   <- slim_lines[i + 1]
  name_line <- slim_lines[i + 2]
  goid <- gsub("id: ", "", id_line)
  name <- gsub("name: ", "", name_line)
  slim_terms <- rbind(slim_terms, data.frame(GOID = goid, slim_term = name))
}

cat("GO slim terms loaded:", nrow(slim_terms), "\n")

# Map each outlier gene's GO terms to slim categories using GO.db ancestry
# Get all ancestor terms for each GO term in our outlier set
if (!requireNamespace("GO.db", quietly = TRUE)) BiocManager::install("GO.db")

outlier_go_ids <- unique(outlier_go$GO_term)

# For each outlier GO term, find which slim terms are its ancestors
go_to_slim <- lapply(outlier_go_ids, function(go_id) {
  tryCatch({
    ancestors <- c(go_id, unlist(as.list(GOBPANCESTOR[[go_id]])))
    ancestors <- ancestors[ancestors %in% slim_terms$GOID]
    if (length(ancestors) == 0) return(NULL)
    data.frame(GO_term = go_id,
               slim_id = ancestors,
               stringsAsFactors = FALSE)
  }, error = function(e) NULL)
})

go_to_slim_df <- dplyr::bind_rows(go_to_slim) %>%
  dplyr::left_join(slim_terms, by = c("slim_id" = "GOID")) %>%
  dplyr::filter(!is.na(slim_term), slim_term != "biological_process")  # remove root term

cat("GO terms mapped to slim categories:", length(unique(go_to_slim_df$GO_term)), "\n")

# Join back to outlier genes
outlier_slim <- outlier_go %>%
  dplyr::left_join(go_to_slim_df, by = "GO_term", relationship = "many-to-many") %>%
  dplyr::filter(!is.na(slim_term)) %>%
  dplyr::select(outlier_type, gene_id, slim_term) %>%
  dplyr::distinct()

# Summary: genes per slim category per selection type
slim_summary <- outlier_slim %>%
  dplyr::group_by(outlier_type, slim_term) %>%
  dplyr::summarise(n_genes = dplyr::n_distinct(gene_id), .groups = "drop") %>%
  dplyr::arrange(outlier_type, desc(n_genes))

cat("\nBroad biological categories under selection:\n")
print(slim_summary, n = 50)

write.csv(slim_summary, file.path(OUT_DIR, "plestiodon_GO_slim_summary.csv"), row.names = FALSE)

# --- GO SLIM SUMMARY TABLE ---------------------------------------------------

# Reshape to wide format: one row per slim term, columns for each selection type
slim_table <- slim_summary %>%
  tidyr::pivot_wider(
    names_from  = outlier_type,
    values_from = n_genes,
    values_fill = 0
  ) %>%
  dplyr::mutate(
    Total = `Positive selection` + `Balancing selection`
  ) %>%
  dplyr::arrange(desc(Total)) %>%
  dplyr::rename(
    `Biological Process (GO Slim)` = slim_term,
    `Positive Selection (n genes)` = `Positive selection`,
    `Balancing Selection (n genes)` = `Balancing selection`,
    `Total Genes` = Total
  )

# Save as CSV
write.csv(slim_table,
          file.path(OUT_DIR, "plestiodon_GO_slim_table.csv"),
          row.names = FALSE)
cat("Saved: plestiodon_GO_slim_table.csv\n")

# Save as Word document
if (!requireNamespace("officer", quietly = TRUE)) install.packages("officer")
if (!requireNamespace("flextable", quietly = TRUE)) install.packages("flextable")
library(officer)
library(flextable)

ft <- flextable(slim_table) %>%
  bold(part = "header") %>%
  bg(part = "header", bg = "#2E75B6") %>%
  color(part = "header", color = "white") %>%
  bg(i = seq(2, nrow(slim_table), 2), bg = "#F2F7FC") %>%
  border_outer(part = "all", border = fp_border(color = "#CCCCCC", width = 1)) %>%
  border_inner(part = "all", border = fp_border(color = "#CCCCCC", width = 0.5)) %>%
  set_table_properties(width = 1, layout = "autofit") %>%
  set_caption("Table 1. Broad biological processes under selection in Plestiodon longirostris, identified using GO slim categories. Genes overlapping Tajima's D outlier windows (1st and 99th percentile) were annotated with GO biological process terms and mapped to GO slim categories. n genes = number of unique genes in each category.") %>%
  fontsize(size = 10, part = "all") %>%
  font(fontname = "Arial", part = "all") %>%
  align(j = 2:4, align = "center", part = "all") %>%
  align(j = 1, align = "left", part = "all")

doc <- read_docx() %>%
  body_add_par("Table 1. GO Slim Biological Process Summary", style = "heading 2") %>%
  body_add_flextable(ft)

print(doc, target = file.path(OUT_DIR, "plestiodon_GO_slim_table.docx"))
cat("Saved: plestiodon_GO_slim_table.docx\n")

# --- 14. SUMMARY -------------------------------------------------------------

cat("\n==================================================\n")
cat("Done. Outputs saved to:", OUT_DIR, "\n")
cat("  plestiodon_TajimaD_manhattan.pdf/.png\n")
cat("  plestiodon_outlier_genes_TajimaD.csv\n")
cat("  plestiodon_population_FST_summary.csv\n")
cat("  plestiodon_outlier_genes_GO.csv\n")
cat("  plestiodon_GO_summary.csv\n")
cat("  plestiodon_GO_enrichment_Positive_selection.csv\n")
cat("  plestiodon_GO_enrichment_Balancing_selection.csv\n")
cat("==================================================\n")

cat("\nTop outlier genes (positive selection):\n")
print(dplyr::filter(overlapping_genes, outlier_type == "Positive selection") %>%
        dplyr::arrange(Tajima) %>% head(20))

cat("\nTop outlier genes (balancing selection):\n")
print(dplyr::filter(overlapping_genes, outlier_type == "Balancing selection") %>%
        dplyr::arrange(desc(Tajima)) %>% head(20))

# =============================================================================
# Pairwise FST Outlier Analysis — Site Differences
# Populations: CAI, NS, SB, SI
# =============================================================================

library(tidyr)

# --- FST FILE SETUP ----------------------------------------------------------

POPS <- c("CAI", "NS", "SB", "SI")

# Generate all pairwise combinations
pop_pairs <- combn(POPS, 2, simplify = FALSE)
pair_labels <- sapply(pop_pairs, paste, collapse = "_vs_")

# Map pair labels to file paths
fst_windowed_files <- sapply(pop_pairs, function(pair) {
  f <- file.path(FST_DIR,
                 paste0("fst_pop_map_", pair[1], "_pop_map_", pair[2], ".fst.windowed.tsv"))
  if (!file.exists(f)) {
    # try reversed order
    f <- file.path(FST_DIR,
                   paste0("fst_pop_map_", pair[2], "_pop_map_", pair[1], ".fst.windowed.tsv"))
  }
  f
})

cat("FST windowed files found:\n")
for (i in seq_along(pair_labels)) {
  cat(" ", pair_labels[i], "->", file.exists(fst_windowed_files[i]), "\n")
}


# --- READ AND PARSE FST WINDOWED FILES ---------------------------------------

read_fst_windowed <- function(filepath, pair_label) {
  lines <- readLines(filepath)
  
  # Keep only lines that look like data rows (start with "(" or "region")
  data_lines <- lines[grepl("^\\(|^region", lines)]
  
  if (length(data_lines) <= 1) return(NULL)  # only header, no data
  
  df <- fread(text = paste(data_lines, collapse = "\n"), header = TRUE, fill = TRUE)
  
  colnames(df) <- c("window_info", "Chr", "midPos", "Nsites", "FST")
  
  df %>%
    dplyr::mutate(FST = as.numeric(FST)) %>%
    dplyr::filter(!is.na(FST), !is.nan(FST), FST >= 0, Nsites >= 10) %>%
    dplyr::mutate(comparison = pair_label) %>%
    dplyr::select(Chr, midPos, Nsites, FST, comparison)
}

cat("\nReading FST windowed files...\n")
fst_all <- dplyr::bind_rows(mapply(read_fst_windowed,
                                   fst_windowed_files,
                                   pair_labels,
                                   SIMPLIFY = FALSE))

cat("Total windows loaded:", nrow(fst_all), "\n")
cat("Windows per comparison:\n")
print(table(fst_all$comparison))


# --- TRANSLATE TO CHROMOSOME COORDINATES -------------------------------------

cat("\nTranslating FST windows to chromosome coordinates...\n")

fst_chr <- fst_all %>%
  dplyr::filter(Chr %in% contig_to_chr$contig) %>%
  dplyr::left_join(contig_to_chr, by = c("Chr" = "contig")) %>%
  dplyr::mutate(
    CHROM = chr,
    POS   = as.integer((chr_start + chr_end) / 2)
  ) %>%
  dplyr::select(CHROM, POS, FST, Nsites, comparison, orig_contig = Chr)

cat("Windows mapped to chromosomes:", nrow(fst_chr), "\n")


# --- OUTLIER DETECTION PER PAIR (99th percentile) ----------------------------

cat("\nDetecting FST outliers per population pair...\n")

fst_chr <- fst_chr %>%
  dplyr::group_by(comparison) %>%
  dplyr::mutate(
    threshold = quantile(FST, 0.99, na.rm = TRUE),
    outlier   = FST >= threshold
  ) %>%
  dplyr::ungroup()

outlier_counts <- fst_chr %>%
  dplyr::filter(outlier) %>%
  dplyr::group_by(comparison) %>%
  dplyr::summarise(n_outlier_windows = n(), threshold = first(threshold))

cat("Outlier windows per comparison:\n")
print(outlier_counts)


# --- GENE OVERLAP PER PAIR ---------------------------------------------------

cat("\nFinding genes overlapping FST outlier windows...\n")

WINDOW_SIZE <- 50000

get_overlapping_genes_fst <- function(fst_data, comparison_label) {
  outliers <- fst_data %>% dplyr::filter(outlier)
  if (nrow(outliers) == 0) return(NULL)
  
  # Add _pilon suffix to match GFF3
  chrom_names <- paste0(outliers$CHROM, "_pilon")
  
  fst_gr <- GRanges(
    seqnames = chrom_names,
    ranges   = IRanges(
      start = pmax(1, outliers$POS - WINDOW_SIZE / 2),
      end   = outliers$POS + WINDOW_SIZE / 2
    ),
    FST = outliers$FST
  )
  
  hits <- findOverlaps(fst_gr, gene_gr)
  if (length(hits) == 0) return(NULL)
  
  data.frame(
    comparison   = comparison_label,
    CHROM        = as.character(seqnames(fst_gr))[queryHits(hits)],
    WIN_CENTER   = outliers$POS[queryHits(hits)],
    FST          = fst_gr$FST[queryHits(hits)],
    gene_id      = gene_gr$gene_id[subjectHits(hits)],
    stringsAsFactors = FALSE
  ) %>% dplyr::distinct()
}

fst_genes_list <- lapply(pair_labels, function(lbl) {
  get_overlapping_genes_fst(
    dplyr::filter(fst_chr, comparison == lbl),
    lbl
  )
})

fst_overlapping_genes <- dplyr::bind_rows(fst_genes_list)

cat("Total gene-window overlaps:", nrow(fst_overlapping_genes), "\n")
cat("Unique genes across all comparisons:", length(unique(fst_overlapping_genes$gene_id)), "\n")

# Per pair summary
per_pair_genes <- fst_overlapping_genes %>%
  dplyr::group_by(comparison) %>%
  dplyr::summarise(n_unique_genes = dplyr::n_distinct(gene_id), .groups = "drop")
cat("\nUnique genes per comparison:\n")
print(per_pair_genes)

write.csv(fst_overlapping_genes,
          file.path(OUT_DIR, "plestiodon_pairwise_FST_outlier_genes.csv"),
          row.names = FALSE)


# --- SHARED vs PAIR-SPECIFIC GENES -------------------------------------------

cat("\nIdentifying shared vs pair-specific divergent genes...\n")

gene_pair_counts <- fst_overlapping_genes %>%
  dplyr::group_by(gene_id) %>%
  dplyr::summarise(
    n_comparisons  = dplyr::n_distinct(comparison),
    comparisons    = paste(sort(unique(comparison)), collapse = "; "),
    mean_FST       = mean(FST),
    .groups = "drop"
  ) %>%
  dplyr::arrange(desc(n_comparisons), desc(mean_FST))

shared_genes <- gene_pair_counts %>% dplyr::filter(n_comparisons >= 3)
specific_genes <- gene_pair_counts %>% dplyr::filter(n_comparisons == 1)

cat("Genes diverging in 3+ comparisons (broadly divergent):", nrow(shared_genes), "\n")
cat("Genes diverging in only 1 comparison (pair-specific):", nrow(specific_genes), "\n")

write.csv(gene_pair_counts,
          file.path(OUT_DIR, "plestiodon_FST_gene_sharing.csv"),
          row.names = FALSE)


# --- GO SLIM FOR FST OUTLIER GENES -------------------------------------------

cat("\nRunning GO slim annotation for FST outlier genes...\n")

# Join GO terms to FST outlier genes
fst_go <- fst_overlapping_genes %>%
  dplyr::left_join(go_table, by = "gene_id", relationship = "many-to-many") %>%
  dplyr::left_join(go_info, by = c("go_terms" = "GOID")) %>%
  dplyr::filter(!is.na(TERM)) %>%
  dplyr::select(comparison, gene_id, GO_term = go_terms, GO_description = TERM) %>%
  dplyr::distinct()

# Map to GO slim
fst_go_slim <- fst_go %>%
  dplyr::left_join(go_to_slim_df, by = c("GO_term"), relationship = "many-to-many") %>%
  dplyr::filter(!is.na(slim_term)) %>%
  dplyr::select(comparison, gene_id, slim_term) %>%
  dplyr::distinct()

# Per pair GO slim summary
slim_per_pair <- fst_go_slim %>%
  dplyr::group_by(comparison, slim_term) %>%
  dplyr::summarise(n_genes = dplyr::n_distinct(gene_id), .groups = "drop") %>%
  dplyr::arrange(comparison, desc(n_genes))

# Pooled GO slim summary (across all pairs)
slim_pooled <- fst_go_slim %>%
  dplyr::group_by(slim_term) %>%
  dplyr::summarise(
    n_genes       = dplyr::n_distinct(gene_id),
    n_comparisons = dplyr::n_distinct(comparison),
    comparisons   = paste(sort(unique(comparison)), collapse = "; "),
    .groups = "drop"
  ) %>%
  dplyr::arrange(desc(n_genes))

cat("\nPooled GO slim summary across all pairwise comparisons:\n")
print(slim_pooled, n = 30)

cat("\nTop GO slim categories per comparison:\n")
print(slim_per_pair %>% dplyr::group_by(comparison) %>%
        dplyr::slice_max(n_genes, n = 5), n = 50)

write.csv(slim_per_pair,
          file.path(OUT_DIR, "plestiodon_FST_GO_slim_per_pair.csv"),
          row.names = FALSE)
write.csv(slim_pooled,
          file.path(OUT_DIR, "plestiodon_FST_GO_slim_pooled.csv"),
          row.names = FALSE)


# --- WIDE TABLE FOR THESIS ---------------------------------------------------

# Reshape per-pair slim summary to wide format for easy reading
slim_pair_wide <- slim_per_pair %>%
  tidyr::pivot_wider(
    names_from  = comparison,
    values_from = n_genes,
    values_fill = 0
  ) %>%
  dplyr::mutate(Total = rowSums(dplyr::across(where(is.numeric)))) %>%
  dplyr::arrange(desc(Total)) %>%
  dplyr::rename(`Biological Process (GO Slim)` = slim_term)

write.csv(slim_pair_wide,
          file.path(OUT_DIR, "plestiodon_FST_GO_slim_wide_table.csv"),
          row.names = FALSE)

cat("\nSaved outputs:\n")
cat("  plestiodon_pairwise_FST_outlier_genes.csv\n")
cat("  plestiodon_FST_gene_sharing.csv\n")
cat("  plestiodon_FST_GO_slim_per_pair.csv\n")
cat("  plestiodon_FST_GO_slim_pooled.csv\n")
cat("  plestiodon_FST_GO_slim_wide_table.csv\n")

# =============================================================================
# Pairwise FST Visualisations
# 1. Heatmap - GO slim categories vs population pairs
# 2. UpSet plot - shared vs pair-specific genes
# 3. Per-pair Manhattan plots
# =============================================================================

# Install/load required packages
vis_packages <- c("ggplot2", "dplyr", "tidyr", "stringr", "patchwork")
for (pkg in vis_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
  library(pkg, character.only = TRUE)
}

if (!requireNamespace("UpSetR", quietly = TRUE)) install.packages("UpSetR")
library(UpSetR)


# =============================================================================
# PLOT 1: HEATMAP - GO slim categories vs population pairs
# =============================================================================

cat("Generating heatmap...\n")

# Use per-pair slim summary, filter to terms appearing in at least 2 comparisons
slim_heatmap <- slim_per_pair %>%
  dplyr::group_by(slim_term) %>%
  dplyr::mutate(n_comparisons_present = sum(n_genes > 0)) %>%
  dplyr::ungroup() %>%
  dplyr::filter(n_comparisons_present >= 2) %>%
  # Clean up comparison labels
  dplyr::mutate(
    comparison = gsub("_vs_", " vs ", comparison),
    # Order slim terms by total genes across comparisons
    slim_term  = stringr::str_wrap(slim_term, width = 35)
  )

# Order slim terms by total gene count
term_order <- slim_heatmap %>%
  dplyr::group_by(slim_term) %>%
  dplyr::summarise(total = sum(n_genes)) %>%
  dplyr::arrange(total) %>%
  dplyr::pull(slim_term)

slim_heatmap$slim_term <- factor(slim_heatmap$slim_term, levels = term_order)

# Order comparisons to group islet vs mainland
comp_order <- c("CAI vs NS", "CAI vs SI", "NS vs SI",
                "CAI vs SB", "NS vs SB", "SB vs SI")
slim_heatmap$comparison <- factor(slim_heatmap$comparison, levels = comp_order)

p_heatmap <- ggplot(slim_heatmap, aes(x = comparison, y = slim_term, fill = n_genes)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_text(aes(label = ifelse(n_genes > 0, n_genes, "")),
            size = 3, colour = "white", fontface = "bold") +
  scale_fill_gradientn(
    colours  = c("#EFF3FF", "#6BAED6", "#2171B5", "#084594"),
    name     = "Number\nof genes",
    na.value = "grey95"
  ) +
  geom_vline(xintercept = 3.5, linetype = "dashed",
             colour = "grey40", linewidth = 0.7) +
  scale_x_discrete(
    position = "bottom",
    labels   = function(x) {
      ifelse(x %in% c("CAI vs NS", "CAI vs SI", "NS vs SI"),
             paste0(x, "\n(islet cluster)"),
             paste0(x, "\n(vs mainland)"))
    }
  ) +
  labs(
    title    = "Biological processes diverging between Plestiodon longirostris populations",
    subtitle = "GO slim categories in FST outlier windows (99th percentile) | n = number of divergent genes",
    x        = NULL,
    y        = NULL
  ) +
  theme_minimal() +
  theme(
    axis.text.x        = element_text(angle = 0, hjust = 0.5, size = 9, face = "bold"),
    axis.text.y        = element_text(size = 8, lineheight = 0.9),
    plot.title         = element_text(face = "bold", size = 13),
    plot.subtitle      = element_text(size = 9, colour = "grey40"),
    legend.position    = "right",
    panel.grid         = element_blank(),
    plot.margin        = margin(10, 10, 10, 10)
  )


print(p_heatmap)
ggsave(file.path(OUT_DIR, "plestiodon_FST_heatmap.pdf"),
       plot = p_heatmap, width = 10, height = 8, dpi = 300)
ggsave(file.path(OUT_DIR, "plestiodon_FST_heatmap.png"),
       plot = p_heatmap, width = 10, height = 8, dpi = 300)
cat("Saved: plestiodon_FST_heatmap.pdf/.png\n")


# =============================================================================
# PLOT 2: UPSET PLOT - shared vs pair-specific divergent genes
# =============================================================================

cat("Generating UpSet plot...\n")

# Build binary matrix: genes x comparisons
all_divergent_genes <- unique(fst_overlapping_genes$gene_id)

upset_matrix <- fst_overlapping_genes %>%
  dplyr::select(gene_id, comparison) %>%
  dplyr::distinct() %>%
  dplyr::mutate(
    present    = 1,
    comparison = gsub("_vs_", " vs ", comparison)
  ) %>%
  tidyr::pivot_wider(
    names_from  = comparison,
    values_from = present,
    values_fill = 0
  ) %>%
  as.data.frame()

rownames(upset_matrix) <- upset_matrix$gene_id
upset_matrix <- upset_matrix %>% dplyr::select(-gene_id)

# Reorder columns to group islet vs mainland
col_order <- c("CAI vs NS", "CAI vs SI", "NS vs SI",
               "CAI vs SB", "NS vs SB", "SB vs SI")
col_order <- col_order[col_order %in% colnames(upset_matrix)]
upset_matrix <- upset_matrix[, col_order]

# Colour sets by islet vs mainland involvement
set_colours <- c(
  "CAI vs NS" = "#4878CF",
  "CAI vs SI" = "#6FA8DC",
  "NS vs SI"  = "#9FC5E8",
  "CAI vs SB" = "#CC0000",
  "NS vs SB"  = "#E06666",
  "SB vs SI"  = "#EA9999"
)

pdf(file.path(OUT_DIR, "plestiodon_FST_upset.pdf"), width = 12, height = 6)
upset(
  upset_matrix,
  sets            = col_order,
  sets.bar.color  = set_colours[col_order],
  order.by        = "freq",
  decreasing      = TRUE,
  mb.ratio        = c(0.6, 0.4),
  text.scale      = c(1.3, 1.1, 1, 1, 1.2, 1),
  point.size      = 3,
  line.size       = 1,
  mainbar.y.label = "Number of divergent genes",
  sets.x.label    = "Total divergent genes",
  main.bar.color  = "#2171B5",
  matrix.color    = "#2171B5"
)
dev.off()

# Also save PNG
png(file.path(OUT_DIR, "plestiodon_FST_upset.png"),
    width = 12, height = 6, units = "in", res = 300)
upset(
  upset_matrix,
  sets            = col_order,
  sets.bar.color  = set_colours[col_order],
  order.by        = "freq",
  decreasing      = TRUE,
  mb.ratio        = c(0.6, 0.4),
  text.scale      = c(1.3, 1.1, 1, 1, 1.2, 1),
  point.size      = 3,
  line.size       = 1,
  mainbar.y.label = "Number of divergent genes",
  sets.x.label    = "Total divergent genes",
  main.bar.color  = "#2171B5",
  matrix.color    = "#2171B5"
)
dev.off()
cat("Saved: plestiodon_FST_upset.pdf/.png\n")


# =============================================================================
# PLOT 3: PER-PAIR MANHATTAN PLOTS
# =============================================================================

cat("Generating per-pair Manhattan plots...\n")

# Prepare FST data with chromosome coordinates for plotting
fst_plot_data <- fst_chr %>%
  dplyr::mutate(CHROM = factor(CHROM, levels = chr_order)) %>%
  dplyr::filter(!is.na(CHROM)) %>%
  dplyr::left_join(
    scaffold_lengths %>% dplyr::select(CHROM, cum_start),
    by = "CHROM"
  ) %>%
  dplyr::mutate(cum_pos = POS + cum_start)

# Add threshold per comparison
fst_plot_data <- fst_plot_data %>%
  dplyr::group_by(comparison) %>%
  dplyr::mutate(threshold = quantile(FST, 0.99, na.rm = TRUE)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    outlier      = FST >= threshold,
    colour_group = factor(match(CHROM, chr_order) %% 2),
    comparison   = gsub("_vs_", " vs ", comparison)
  )

# Clean comparison order
comp_order_plot <- c("CAI vs NS", "CAI vs SI", "NS vs SI",
                     "CAI vs SB", "NS vs SB", "SB vs SI")
fst_plot_data$comparison <- factor(fst_plot_data$comparison, levels = comp_order_plot)

# Build one plot per comparison then combine with patchwork
plot_list <- lapply(comp_order_plot, function(comp) {
  
  df      <- dplyr::filter(fst_plot_data, comparison == comp)
  thresh  <- unique(df$threshold)
  is_mainland <- grepl("SB", comp)
  
  ggplot(df, aes(x = cum_pos, y = FST, colour = colour_group)) +
    geom_point(size = 0.3, alpha = 0.4, show.legend = FALSE) +
    geom_point(
      data    = dplyr::filter(df, outlier),
      colour  = ifelse(is_mainland, "#CC0000", "#2171B5"),
      size    = 0.9, alpha = 0.9
    ) +
    geom_hline(yintercept = thresh, linetype = "dashed",
               colour = "black", linewidth = 0.4) +
    scale_x_continuous(
      breaks = scaffold_mids$mid,
      labels = gsub("_RagTag.*", "", scaffold_mids$CHROM),
      expand = c(0.01, 0)
    ) +
    scale_colour_manual(values = c("0" = "#AAAAAA", "1" = "#CCCCCC")) +
    scale_y_continuous(limits = c(0, max(df$FST, na.rm = TRUE) * 1.05)) +
    labs(title = comp, x = NULL, y = expression(F[ST])) +
    theme_classic() +
    theme(
      axis.text.x      = element_text(angle = 45, hjust = 1, size = 5),
      axis.text.y      = element_text(size = 7),
      axis.title.y     = element_text(size = 8),
      plot.title       = element_text(size = 9, face = "bold",
                                      colour = ifelse(is_mainland, "#CC0000", "#2171B5")),
      panel.grid.major.y = element_line(colour = "grey95", linewidth = 0.3)
    )
})

# Combine into 2-column layout with patchwork
p_manhattan <- patchwork::wrap_plots(plot_list, ncol = 2) +
  patchwork::plot_annotation(
    title    = "Genome-wide FST Manhattan plots — Plestiodon longirostris population pairs",
    subtitle = "Outliers = 99th percentile FST per comparison | Blue = within islet cluster | Red = islet vs SB mainland",
    theme    = theme(
      plot.title    = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(size = 9, colour = "grey40")
    )
  )

print(p_manhattan)
ggsave(file.path(OUT_DIR, "plestiodon_FST_manhattan_pairs.pdf"),
       plot = p_manhattan, width = 16, height = 14, dpi = 300)
ggsave(file.path(OUT_DIR, "plestiodon_FST_manhattan_pairs.png"),
       plot = p_manhattan, width = 16, height = 14, dpi = 300)
cat("Saved: plestiodon_FST_manhattan_pairs.pdf/.png\n")

cat("\nAll three visualisations complete.\n")

#comparing our outlier genes to other reptile species
# Get full product names for all FST outlier genes
outlier_gene_products <- gff_df %>%
  dplyr::filter(type == "mRNA") %>%
  dplyr::mutate(
    parent_id   = sapply(Parent, function(x) if (length(x) > 0) x[[1]] else NA_character_),
    product_chr = sapply(product, function(x) if (length(x) > 0) x[[1]] else NA_character_),
    dbxref_chr  = sapply(Dbxref, function(x) if (length(x) > 0) paste(x, collapse = "; ") else NA_character_)
  ) %>%
  dplyr::filter(parent_id %in% unique(fst_overlapping_genes$gene_id)) %>%
  dplyr::select(gene_id = parent_id, transcript_id = ID, product = product_chr, Dbxref = dbxref_chr) %>%
  dplyr::distinct(gene_id, .keep_all = TRUE) %>%
  dplyr::left_join(
    fst_overlapping_genes %>% 
      dplyr::select(gene_id, comparison) %>% 
      dplyr::distinct() %>%
      dplyr::group_by(gene_id) %>%
      dplyr::summarise(comparisons = paste(sort(unique(gsub("_vs_", " vs ", comparison))), 
                                           collapse = "; "),
                       n_comparisons = dplyr::n_distinct(comparison),
                       .groups = "drop"),
    by = "gene_id"
  ) %>%
  dplyr::arrange(desc(n_comparisons), gene_id)

# Also do the same for Tajima's D outlier genes
tajima_gene_products <- gff_df %>%
  dplyr::filter(type == "mRNA") %>%
  dplyr::mutate(
    parent_id   = sapply(Parent, function(x) if (length(x) > 0) x[[1]] else NA_character_),
    product_chr = sapply(product, function(x) if (length(x) > 0) x[[1]] else NA_character_),
    dbxref_chr  = sapply(Dbxref, function(x) if (length(x) > 0) paste(x, collapse = "; ") else NA_character_)
  ) %>%
  dplyr::filter(parent_id %in% unique(overlapping_genes$gene_id)) %>%
  dplyr::select(gene_id = parent_id, transcript_id = ID, 
                product = product_chr, Dbxref = dbxref_chr) %>%
  dplyr::distinct(gene_id, .keep_all = TRUE) %>%
  dplyr::left_join(
    overlapping_genes %>%
      dplyr::select(gene_id, outlier_type) %>%
      dplyr::distinct(),
    by = "gene_id"
  ) %>%
  dplyr::arrange(outlier_type, gene_id)

# Save both
write.csv(outlier_gene_products,
          file.path(OUT_DIR, "plestiodon_FST_outlier_gene_products.csv"),
          row.names = FALSE)
write.csv(tajima_gene_products,
          file.path(OUT_DIR, "plestiodon_TajimaD_outlier_gene_products.csv"),
          row.names = FALSE)

cat("Genes with recognisable product names (FST outliers):\n")
print(outlier_gene_products %>% 
        dplyr::filter(!grepl("hypothet", product, ignore.case = TRUE),
                      !is.na(product)) %>%
        dplyr::select(gene_id, product, n_comparisons, comparisons) %>%
        head(30))
