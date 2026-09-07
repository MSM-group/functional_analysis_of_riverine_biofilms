library(tidyverse)

# 1. Load FAPROTAX functional table and Metadata
fap_res <- read_tsv("output/amplicon_analyses/20260826_functional_table_faprotax1.2.12.tsv") %>%
  rename(Functional_Group = group)

metadata <- read_csv("output/amplicon_analyses/metadata_necessary.csv")

# 2. Pivot long & filter for River Biofilms
plot_data_raw <- fap_res %>% 
  pivot_longer(-Functional_Group, names_to = "SampleID", values_to = "Reads") %>% 
  left_join(metadata, by = "SampleID") %>% 
  filter(str_detect(`Sample-type`, regex("RiverBiof[li]{2}m", ignore_case = TRUE))) %>% 
  filter(!is.na(Location))

# Filter out functional groups that have 0 reads across ALL river biofilm samples
active_groups <- plot_data_raw %>% 
  group_by(Functional_Group) %>% 
  summarise(total = sum(Reads), .groups = "drop") %>% 
  filter(total > 0) %>% 
  pull(Functional_Group)

filtered_df <- plot_data_raw %>% 
  filter(Functional_Group %in% active_groups) %>% 
  arrange(Location, `short-identifier`) %>% 
  mutate(
    `short-identifier` = factor(`short-identifier`, levels = unique(`short-identifier`))
  )

# 3. Format Functional_Group factor levels (put 'other' first so it stacks at the bottom)
assigned_groups <- sort(setdiff(unique(filtered_df$Functional_Group), "other"))

plot_data <- filtered_df %>% 
  mutate(
    Category = ifelse(Functional_Group == "other", "Missing Annotation", Functional_Group),
    Category = factor(Category, levels = c("Missing Annotation", assigned_groups))
  )

# 4. Generate High-Contrast Color Palette (Spectral)
num_assigned <- length(assigned_groups)

assigned_colors <- hcl.colors(n = num_assigned, palette = "Spectral")
names(assigned_colors) <- assigned_groups

category_colors <- c("Missing Annotation" = "gray85", assigned_colors)

# 5. Compute x-axis elements for Location separators, Labels, and Fontface
x_levels <- levels(plot_data$`short-identifier`)
x_fontface <- ifelse(str_detect(x_levels, regex("dw", ignore_case = TRUE)), "bold", "plain")

location_boundaries <- plot_data %>% 
  distinct(`short-identifier`, Location) %>% 
  mutate(x_pos = as.numeric(`short-identifier`)) %>% 
  mutate(location_change = Location != lead(Location)) %>% 
  filter(location_change) %>% 
  pull(x_pos) + 0.5

location_labels <- plot_data %>% 
  distinct(`short-identifier`, Location) %>% 
  mutate(x_pos = as.numeric(`short-identifier`)) %>% 
  group_by(Location) %>% 
  summarise(x_mid = mean(x_pos), .groups = "drop")

# 6. Render Plot
ggplot(plot_data, aes(x = `short-identifier`, y = Reads, fill = Category)) +
  geom_bar(stat = "identity", position = "fill", width = 0.8) +
  geom_vline(xintercept = location_boundaries, color = "gray30", linewidth = 0.7) +
  geom_text(
    data = location_labels,
    aes(x = x_mid, y = 1.03, label = Location),
    inherit.aes = FALSE,
    fontface = "bold",
    size = 3,
    vjust = 0,
    color = "gray20"
  ) +
  scale_fill_manual(values = category_colors) +
  scale_y_continuous(
    labels = scales::percent,
    expand = expansion(mult = c(0, 0.08))
  ) +
  coord_cartesian(clip = "off") +
  labs(
    x = "Sample",
    y = "Relative Abundance (%)",
    fill = "FAPROTAX Group"
  ) +
  theme_light(base_size = 11) +
  theme(
    axis.text.x = element_text(
      angle = 45, 
      hjust = 1, 
      vjust = 1, 
      size = 8, 
      face = x_fontface
    ),
    panel.grid.major.x = element_blank(),
    legend.position = "bottom",
    legend.key.size = unit(0.2, "cm"),
    legend.text = element_text(size = 6),
    plot.margin = margin(t = 15, r = 10, b = 10, l = 10)
  )

# Save plot
ggsave("output/amplicon_analyses/20260902_RiverBiofilm_all_faprotax_groups_barplot_wo_title.pdf", width = 13, height = 7)
ggsave("output/amplicon_analyses/20260827_RiverBiofilm_all_faprotax_groups_barplot.png", width = 13, height = 7, dpi = 300)