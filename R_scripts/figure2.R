# Packages ----
install.packages(c("tidyverse", "readxl"))

library(tidyverse)  # for ggplot2, dplyr, etc.
library(readxl)     # for reading Excel files
library(ggrepel)
library(viridis)
library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)   # for color functions
library(tibble)

# Load data from excel sheet (Source Data) ----
# Path to your Excel file
excel_file <- "41467_2023_42762_MOESM6_ESM.xlsx"

df_2b <- read_excel(excel_file, sheet = "Fig 1c, 2b, Supp Fig 1a")
df_2ef <- read_excel(excel_file, sheet = "Fig 2e-f")
df_2g <- read_excel(excel_file, sheet = "Fig 2g")

# Figure 2b: chr22q loss UMAP -----
## Color by chr22q loss ----
# Examine data structure
head(df_2b)
names(df_2b)
#see counts of 0 and 1 per cluster.
table(df_2b$final_label, df_2b$chr22q_loss)
#                  0     1
#  myeSC        1086   209
#  nmSC        13888  5178

# Make chr22q_loss a factor for coloring
df_2b <- df_2b %>% mutate(chr22q_loss = factor(chr22q_loss))

# Create UMAP plot
fig2b_1 <- ggplot(df_2b, aes(x = UMAP_1, y = UMAP_2, color = chr22q_loss)) +
  geom_point(size = 0.5, alpha = 0.8) +
  scale_color_manual(
    values = c("0" = "grey80", "1" = "#4B0082"),  # adjust purple-blue hex as needed
    name = "chr22q_loss",
    labels = c("0" = "Neutral", "1" = "Loss")
  ) +
  geom_text_repel(
    data = labels_df,
    aes(x = UMAP_1, y = UMAP_2, label = final_label),
    size = 4,
    color = "black",
    fontface = "bold",
    segment.color = "black",
    segment.size = 0.3
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "right"
  ) +
  labs(
    title = "Fig.2B: chr22q loss UMAP",
    x = "UMAP 1",
    y = "UMAP 2"
  ) 


# Print the plot
print(fig2b_1)

## Colored nmSC and myeSC as shown in the paper -----
# Add a new column for the color group
df_2b <- df_2b %>%
  mutate(chr22q_status = ifelse(final_label %in% c("myeSC", "nmSC"), "Loss", "Neutral"))

# Define colors for the 2 groups
color_map_2 <- c(
  "Neutral" = "grey90",
  "Loss" = "#2166AC"
)

# Plot using the new column
fig2b_2 <- ggplot(df_2b, aes(x = UMAP_1, y = UMAP_2, color = chr22q_status)) +
  geom_point(size = 0.5, alpha = 0.8) +
  scale_color_manual(
    values = color_map_2,
    name = "chr22q status"   # legend title
  ) +
  geom_text_repel(
    data = labels_df,
    aes(x = UMAP_1, y = UMAP_2, label = final_label),
    size = 4,
    color = "black",
    fontface = "bold",
    segment.color = "black",
    segment.size = 0.3
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "right"
  ) +
  labs(
    title = "Fig.2B: chr22q loss UMAP",
    x = "UMAP 1",
    y = "UMAP 2"
  )


# Print the plot
print(fig2b_2)


# Save fig 1c to file
ggsave(
  "Figure2b_chr22q_loss_UMAP.png",
  plot = fig2b_2,
  width = 6,    # in inches
  height = 5,   # in inches
  dpi = 300     # good for publication
)

# Figure 2e: UMAP for VS Schwann Cell Subclusters -----
#df_2ef <- read_excel(excel_file, sheet = "Fig 2e-f") %>% mutate(cluster = factor(cluster))

#check what final_label associate with what cluster number
table(df_2ef$final_label, df_2ef$cluster)

# Custom blue shades for each cluster
n_clusters <- length(levels(df_2ef$cluster))
base_blues <- brewer.pal(9, "Blues") 
my_blues <- colorRampPalette(base_blues[4:9])(n_clusters)

# Label positions
labels_df <- df_2ef %>%
  group_by(cluster) %>%
  summarize(
    UMAP_1 = median(UMAP_1),
    UMAP_2 = median(UMAP_2)
  ) %>%
  mutate(final_label = cluster)

# Optional: Move the label positions to the right and up
labels_df <- df_2ef %>%
  group_by(final_label) %>%
  summarize(
    UMAP_1 = median(UMAP_1, na.rm = TRUE),
    UMAP_2 = median(UMAP_2, na.rm = TRUE)
  ) %>%
  mutate(
    UMAP_1_label = UMAP_1 + 2,
    UMAP_2_label = UMAP_2 + 2
  )

## geom_label_repel ----
# Plot
fig2e <- ggplot(df_2ef, aes(x = UMAP_1, y = UMAP_2, color = cluster)) +
  geom_point(size = 0.5) +
  scale_color_manual(
    values = my_blues,
    name = "Cluster #"
  ) +
  geom_label_repel(
    data = labels_df,
    aes(x = UMAP_1, y = UMAP_2, label = final_label),
    size = 4,
    fontface = "bold",
    color = "black",   # text color
    fill = "white",    # box background color
    segment.color = "black",
    segment.size = 0.7
  ) +
  theme_classic() +
  labs(
    title = "Fig.2E: VS Schwann Cell Subclusters",
    x = "UMAP_1",
    y = "UMAP_2"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.title = element_text(face = "bold")
  )

print(fig2e)

ggsave(
  "Figure2e_VSSchwann_UMAP.png",
  plot = fig2e,
  width = 6,    # in inches
  height = 5,   # in inches
  dpi = 300     # good for publication
)

## trying geom_segment -----
fig2e <- ggplot(df_2ef, aes(x = UMAP_1, y = UMAP_2, color = cluster)) +
  geom_point(size = 0.5) +
  scale_color_manual(values = my_blues, name = "Cluster #") +
  geom_segment(
    data = labels_df,
    aes(x = UMAP_1, y = UMAP_2, xend = UMAP_1_label, yend = UMAP_2_label),
    color = "black",
    linewidth = 0.8
  ) +
  geom_label(
    data = labels_df,
    aes(x = UMAP_1_label, y = UMAP_2_label, label = final_label),
    size = 4,
    fontface = "bold",
    fill = "white",
    color = "black"
  ) +
  theme_classic() +
  labs(
    title = "Fig.2E: VS Schwann Cell Subclusters",
    x = "UMAP_1",
    y = "UMAP_2"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.title = element_text(face = "bold")
  )

print(fig2e)


# Figure 2f: UMAP for Ch22q loss - Schwann Cell -----
#see counts of 0 and 1 per final_label
table(df_2ef$final_label, df_2ef$chr22q_loss)

# Make chr22q_loss a factor for coloring
df_2ef <- df_2ef %>% mutate(chr22q_loss = factor(chr22q_loss))

# Create UMAP plot
fig2f <- ggplot(df_2ef, aes(x = UMAP_1, y = UMAP_2, color = chr22q_loss)) +
  geom_point(size = 0.5, alpha = 0.8) +
  scale_color_manual(
    values = c("0" = "grey90", "1" = "#08306b"),
    name = "chr22q_loss",
    labels = c("0" = "Neutral", "1" = "Loss")
  ) +
  geom_label_repel(
    data = labels_df,
    aes(x = UMAP_1, y = UMAP_2, label = final_label),
    size = 4,
    fontface = "bold",
    color = "black",   # text color
    fill = "white",    # box background color
    segment.color = "black",
    segment.size = 0.7
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "right"
  ) +
  labs(
    title = "Fig.2f: Chr22q Loss - Schwann Cells",
    x = "UMAP 1",
    y = "UMAP 2"
  ) 

print(fig2f)

ggsave(
  "Figure2f_Chr22q_Schwann_UMAP.png",
  plot = fig2f,
  width = 6,    # in inches
  height = 5,   # in inches
  dpi = 300     # good for publication
)

# Figure 2g: Heatmap VS-SC cluster in diiferent nerve type ---------
row_order <- c("Stress", "IFN response", "MHC II", "Hypoxia", "Repair-like", "Myelinating", "Core")
column_order <- c("mSC cluster 1", "nmSC", "mSC cluster 2", "Myelinating SC", "nm(R)SC", "mSC cluster 3",
                  "SC1", "Schwann cells", "SC3", "Dividing SC", "SC2",
                  "tSC", "iSC", "pmSC", "prol. SC")
# reorder row and column
df_2g_ordered <- df_2g %>%
  filter(`VS-SC subtype` %in% row_order) %>%
  mutate(`VS-SC subtype` = factor(`VS-SC subtype`, levels = row_order)) %>%
  arrange(`VS-SC subtype`) %>%
  select(`VS-SC subtype`, all_of(column_order))

mat_2g <- df_2g_ordered %>%
  column_to_rownames("VS-SC subtype") %>%
  as.matrix()

# Vector for each col
nerve_type <- factor(
  c(
    rep("Adult", 6),
    rep("Injured", 5),
    rep("Developing", 4)
  ),
  levels = c("Adult", "Injured", "Developing")  # Force desired order!
)

study <- c("Gerber", "Yim", "Gerber", "Yim", "Gerber", "Gerber",
           "Kalinski", "Carr", "Kalinski", "Carr", "Kalinski", 
           "Gerber", "Gerber", "Gerber", "Gerber")
# Colors
nerve_colors <- c("Adult" = "#1f78b4", "Injured" = "#cb181d", "Developing" = "#B2DF8A")
study_colors <- c(
  "Carr" = "#A6CEE3",
  "Gerber" = "#1F78B4",
  "Kalinski" = "#33A02C",
  "Yim" = "#B2DF8A"
)

# Annotation
col_annotation <- HeatmapAnnotation(
  df = data.frame(Nerve = nerve_type, Study = study),
  col = list(Nerve = nerve_colors, Study = study_colors),
  annotation_name_side = "left"
)

# Color ramp
col_fun <- colorRamp2(c(-2, 0, 2), c("#1F78B4", "white", "#cb181d"))

# Heatmap
fig2g <- Heatmap(
  mat_2g,
  name = "Z-score",
  col = col_fun,
  top_annotation = col_annotation,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_names_side = "left",
  column_names_side = "bottom",
  column_title = "Murine Nerve Schwann Cell Signatures",
  row_title = "VS Associated Schwann Subtypes",
  rect_gp = gpar(col = "white", lwd = 1),  # order color + width
  column_split = nerve_type,         # split by nerve group
  gap = unit(2, "mm"),               # adds clear ghttp://127.0.0.1:31595/graphics/1402fe7d-65ad-487a-a701-14f1a504eac2.pngap between splits
  cluster_column_slices = FALSE      # don’t cluster splits
)

fig2g

# Draw heatmap and save as png
png("Figure2g_Heatmap.png",
    width = 8,    # inches
    height = 6,   # inches
    units = "in", # important!
    res = 300     # dots per inch
)
draw(fig2g)
dev.off()
