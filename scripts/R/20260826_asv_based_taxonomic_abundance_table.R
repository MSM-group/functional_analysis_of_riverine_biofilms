library(qiime2R)
library(tidyverse)

# --- STEP 1: Load Counts ---
message("Loading counts...")
counts_raw <- read_qza("amplicon_analyses/dada_results/filtered_no_singleton_table.qza")$data

counts <- counts_raw %>% 
  as.data.frame() %>% 
  rownames_to_column("FeatureID")

# --- STEP 2: Load Taxonomy ---
message("Loading taxonomy...")
tax_qza <- read_qza("amplicon_analyses/silva_tax_classification/taxonomy.qza")

# Check what structure taxonomy data has
tax_raw <- tax_qza$data

# If parse_taxonomy fails, fallback to direct column extraction
if ("Taxon" %in% colnames(tax_raw)) {
  message("Parsing SILVA taxonomy string...")
  taxa_parsed <- tax_raw %>% 
    select(FeatureID = Feature.ID, Taxon) %>% 
    separate(
      Taxon, 
      into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), 
      sep = ";", 
      fill = "right"
    )
} else {
  taxa_parsed <- parse_taxonomy(tax_raw) %>% 
    as.data.frame() %>% 
    rownames_to_column("FeatureID")
}

# --- STEP 3: Merge ---
message("Merging tables...")
annotated_table <- counts %>% 
  left_join(taxa_parsed, by = "FeatureID")

# Dynamic column selection to prevent crashes if 'Species' is missing
tax_cols <- intersect(c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), colnames(annotated_table))
annotated_table <- annotated_table %>% 
  select(FeatureID, all_of(tax_cols), everything())

# --- STEP 4: Save CSV ---
message("Saving CSV...")
write_csv(annotated_table, "output/amplicon_analyses/20260826_annotated_ASV_abundance_table.csv")
message("Success! Table saved.")
