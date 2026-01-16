# --- Load libraries ---
library(Biostrings)
library(phyloseq)
library(pwalign)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(patchwork)
library(rstudioapi)
library(qiime2R)

# --- Set working directory ---
setwd(dirname(getActiveDocumentContext()$path))
###import env data
save.image(file = "final-clamp-percentage.RData")
# --- Load Qiime2 data ---
ASVs <- read_qza("feature-table-16S-Plant-all.qza")      # ASV table
metadata <- read_q2metadata("metadata_Plant.txt")        # sample metadata
taxonomy <- read_qza("taxonomy_SILVA-16S-Plant-all.qza") # taxonomy
taxonomy <- parse_taxonomy(taxonomy$data)

physeq <- qza_to_phyloseq(
  features = "feature-table-16S-Plant-all.qza",
  tree = "rooted-tree-16S-Plant-all.qza",
  taxonomy = "taxonomy_SILVA-16S-Plant-all.qza",
  metadata = "metadata_Plant.txt"
)

# --- Reproducibility ---
set.seed(1)  
physeq <- rarefy_even_depth(physeq, sample.size=40000, replace=FALSE, rngseed = 1) 

# --- Representative sequences ---
asvs <- readDNAStringSet("rep_seqs.fasta")
names(asvs) <- sub(" .*", "", names(asvs))  # Qiime2 hash IDs

# --- OTU matrix (ASVs in rows, samples in columns) ---
otu <- as(ASVs$data, "matrix")

# --- Ensure SampleID column in metadata ---
if (!"SampleID" %in% colnames(metadata)) {
  metadata <- metadata %>% rownames_to_column("SampleID")
} else {
  metadata$SampleID <- as.character(metadata$SampleID)
}

# --- Plants of interest (fixed order) ---
plants_of_interest <- c("Buttercup", "Holcus", "White clover")

# --- Function to align PNA probe and generate summary + plot ---
analyze_pna <- function(probe, probe_name, otu, asvs, metadata, plants_of_interest,
                        colors, plot_title, output_prefix) {
  
  rc_probe <- reverseComplement(probe)
  mat <- nucleotideSubstitutionMatrix(match = 1, mismatch = -1, baseOnly = TRUE)
  
  align_no_gaps <- function(query, subject) {
    pairwiseAlignment(query, subject,
                      type = "local",
                      substitutionMatrix = mat,
                      gapOpening = 1e6,
                      gapExtension = 1e6)
  }
  
  # Align both directions
  res_fwd <- lapply(asvs, function(subject) align_no_gaps(probe, subject))
  res_rev <- lapply(asvs, function(subject) align_no_gaps(rc_probe, subject))
  
  extract_info <- function(aln, id) {
    if (nchar(as.character(alignedPattern(aln))) == 0) return(NULL)
    data.frame(ASV = id, Mismatches = nrow(mismatchTable(aln)), stringsAsFactors = FALSE)
  }
  
  df_fwd <- do.call(rbind, Map(extract_info, res_fwd, names(asvs)))
  df_rev <- do.call(rbind, Map(extract_info, res_rev, names(asvs)))
  hits <- rbind(df_fwd, df_rev)
  
  # Keep best hit per ASV
  best_hits <- aggregate(Mismatches ~ ASV, data = hits, min)
  best_hits <- best_hits[best_hits$ASV %in% rownames(otu), ]
  best_hits$TotalReads <- rowSums(otu[best_hits$ASV, , drop = FALSE])
  
  # OTU long format with metadata
  otu_long <- as.data.frame(otu) %>%
    rownames_to_column("ASV") %>%
    pivot_longer(-ASV, names_to = "SampleID", values_to = "Reads") %>%
    left_join(best_hits, by = "ASV") %>%
    left_join(metadata %>% select(SampleID, Plant, Clamp, sequencing), by = "SampleID") %>%
    filter(Plant %in% plants_of_interest) %>%
    mutate(
      Plant = factor(Plant, levels = plants_of_interest),  # fixed plant order
      MismatchCategory = case_when(
        Mismatches == 0 ~ "Perfect Match",
        Mismatches == 1 ~ "1 Mismatch",
        Mismatches == 2 ~ "2 Mismatches",
        Mismatches >= 3 ~ "3+ Mismatches"
      ),
      MismatchCategory = factor(
        MismatchCategory,
        levels = c("Perfect Match", "1 Mismatch", "2 Mismatches", "3+ Mismatches")  # fixed legend order
      )
    )
  
  # Summarize by plant / condition
  summary_df <- otu_long %>%
    group_by(Plant, Clamp, sequencing, MismatchCategory) %>%
    summarise(TotalReads = sum(Reads, na.rm = TRUE), .groups = "drop") %>%
    group_by(Plant, Clamp, sequencing) %>%
    mutate(Percent = 100 * TotalReads / sum(TotalReads))
  
  # Save tables
  write.table(best_hits, paste0(output_prefix, "_perASV.txt"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(summary_df, paste0(output_prefix, "_summary.txt"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  
  # Plot
  p <- ggplot(summary_df, aes(x = Plant, y = Percent, fill = MismatchCategory)) +
    geom_bar(stat = "identity", position = "stack", width = 0.75) +  # wider bars
    geom_text(aes(label = round(Percent, 1)),
              position = position_stack(vjust = 0.5),
              size = 3.0) +  # larger text inside bars
    facet_grid(Clamp ~ sequencing, scales = "free_x", space = "free_x") +
    scale_fill_manual(values = colors) +
    scale_y_continuous(limits = c(0, 100), oob = scales::squish)+  # unified y-axis 0–100%
    theme_minimal(base_size = 16) +
    theme(
      axis.text.x = element_text(size = 14, angle = 45, hjust = 1, vjust = 1),
      axis.text.y = element_text(size = 14),
      strip.text = element_text(size = 15, face = "bold"),
      legend.text = element_text(size = 14),
      legend.title = element_text(size = 15, face = "bold"),
      plot.title = element_text(size = 17, face = "bold"),
      axis.title = element_text(size = 15, face = "bold"),
      panel.spacing = unit(1.2, "lines")
    ) +
    labs(title = plot_title, x = "Plant", y = "Percent of Reads")
  
  return(list(summary = summary_df, plot = p))
}

# --- Run analysis for both PNAs ---
chloroplast_results <- analyze_pna(
  probe = DNAString("GGCTCAACCCTGGACAG"),
  probe_name = "Chloroplast",
  otu = otu,
  asvs = asvs,
  metadata = metadata,
  plants_of_interest = plants_of_interest,
  colors = c("Perfect Match" = "#009E73",
             "1 Mismatch" = "#E69F00",
             "2 Mismatches" = "#56B4E9",
             "3+ Mismatches" = "#5D5D5D"),
  plot_title = "a) PNA Clamp Chloroplast",
  output_prefix = "chloroplast_pna"
)

mitochondria_results <- analyze_pna(
  probe = DNAString("GGCAAGTGTTCTTCGGA"),
  probe_name = "Mitochondria",
  otu = otu,
  asvs = asvs,
  metadata = metadata,
  plants_of_interest = plants_of_interest,
  colors = c("Perfect Match" = "#882255",
             "1 Mismatch" = "#E69F00",
             "2 Mismatches" = "#56B4E9",
             "3+ Mismatches" = "#5D5D5D"),
  plot_title = "b) PNA Clamp Mitochondria",
  output_prefix = "mitochondria_pna"
)

# --- Combine plots ---
combined_plot <- chloroplast_results$plot + mitochondria_results$plot

# --- Display combined plot ---
combined_plot

# --- Save combined plot as PDF and PNG ---
ggsave("PNA_clamps_comparison.pdf", combined_plot,
       width = 12, height = 6, units = "in")

ggsave("PNA_clamps_comparison.png", combined_plot,
       width = 12, height = 6, units = "in", dpi = 300)
