library(tidyverse)

# 1. Load data
ratio_df <- read_csv("amplicon_analyses/20260827_faprotax_1.2.12_auto_hetero_ratios_summary.csv")

# 2. Filter and order samples grouped by Location
filtered_df <- ratio_df %>% 
  filter(str_detect(`Sample-type`, regex("RiverBiof[li]{2}m", ignore_case = TRUE))) %>% 
  filter(!is.na(Location)) %>% 
  arrange(Location, `short-identifier`) %>% 
  mutate(
    `short-identifier` = factor(`short-identifier`, levels = unique(`short-identifier`))
  )

# 3. Compute vector specifying fontface for each sample tick label
x_levels <- levels(filtered_df$`short-identifier`)
x_fontface <- ifelse(str_detect(x_levels, regex("dw", ignore_case = TRUE)), "bold", "plain")

# 4. Calculate x-positions for vertical dashed lines (boundaries)
location_boundaries <- filtered_df %>% 
  distinct(`short-identifier`, Location) %>% 
  mutate(x_pos = as.numeric(`short-identifier`)) %>% 
  mutate(location_change = Location != lead(Location)) %>% 
  filter(location_change) %>% 
  pull(x_pos) + 0.5

# 5. Calculate x-axis midpoints for location text labels
location_labels <- filtered_df %>% 
  distinct(`short-identifier`, Location) %>% 
  mutate(x_pos = as.numeric(`short-identifier`)) %>% 
  group_by(Location) %>% 
  summarise(
    x_mid = mean(x_pos),
    .groups = "drop"
  )

# 6. Pivot into long format for plotting (including 'Other_Reads')
plot_data <- filtered_df %>% 
  select(Location, `short-identifier`, Autotroph_Reads, Heterotroph_Reads, Other_Reads, Unannotated_Reads) %>% 
  rename(
    autotroph = Autotroph_Reads,
    heterotroph = Heterotroph_Reads,
    other = Other_Reads,
    missing_annotation = Unannotated_Reads
  ) %>% 
  pivot_longer(
    cols = c(autotroph, heterotroph, other, missing_annotation),
    names_to = "Category",
    values_to = "Reads"
  ) %>% 
  mutate(
    Category = as.character(Category),
    Category = factor(Category, levels = c("missing_annotation", "other", "heterotroph", "autotroph"))
  )

# 7. Color palette including 'other'
category_colors <- c(
  "autotroph" = "forestgreen",
  "heterotroph" = "coral2",
  "other" = "khaki3",
  "missing_annotation" = "gray70"
)

# 8. Render plot
ggplot(plot_data, aes(x = `short-identifier`, y = Reads, fill = Category)) +
  geom_bar(stat = "identity", position = "fill", width = 0.8) +
  geom_vline(xintercept = location_boundaries, color = "gray30", linewidth = 0.7) +
  geom_text(
    data = location_labels,
    aes(x = x_mid, y = 1.03, label = Location),
    inherit.aes = FALSE,
    fontface = "bold",
    size = 4,
    vjust = 0,
    color = "gray20"
  ) +
  scale_fill_manual(
    values = category_colors,
    labels = c(
      "autotroph" = "Autotroph", 
      "heterotroph" = "Heterotroph", 
      "other" = "Other",
      "missing_annotation" = "Missing Annotation"
    )
  ) +
  scale_y_continuous(
    labels = scales::percent,
    expand = expansion(mult = c(0, 0.08))
  ) +
  coord_cartesian(clip = "off") +
  labs(
    title = "Prokaryotes Autotroph Heterotroph Ratio Estimation Based on Taxonomical Annotation",
    subtitle = "Grouped by River Location",
    x = "Sample",
    y = "Relative Abundance (%)",
    fill = "Functional Guild"
  ) +
  theme_light(base_size = 11) +
  theme(
    # Conditional bolding applied via x_fontface vector
    axis.text.x = element_text(
      angle = 45, 
      hjust = 1, 
      vjust = 1, 
      size = 5, 
      face = x_fontface
    ),
    panel.grid.major.x = element_blank(),
    legend.position = "top",
    plot.margin = margin(t = 15, r = 10, b = 10, l = 10)
  )

# Save plot
ggsave("output/amplicon_analyses/20260827_RiverBiofilm_functional_composition_barplot.pdf", width = 11, height = 6)
ggsave("output/amplicon_analyses/20260827_RiverBiofilm_functional_composition_barplot.png", width = 11, height = 6, dpi = 300)