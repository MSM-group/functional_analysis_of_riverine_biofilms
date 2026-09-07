library(qiime2R)
library(tidyverse)

# 1. Read count matrix & raw SILVA taxonomy
counts <- read_qza("amplicon_analyses/dada_results/filtered_no_singleton_table.qza")$data %>% 
  as.data.frame() %>% 
  rownames_to_column("FeatureID")

taxa <- read_qza("amplicon_analyses/silva_tax_classification/taxonomy.qza")$data %>% 
  select(FeatureID = Feature.ID, Taxon)

# 2. Join taxonomy string to counts
faprotax_table <- taxa %>% 
  inner_join(counts, by = "FeatureID") %>% 
  rename(`#OTU ID` = FeatureID)

# 3. Save file for python execution
write_tsv(faprotax_table, "output/amplicon_analyses/20260826_faprotax_input.tsv")