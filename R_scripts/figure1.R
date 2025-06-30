# Packages ----
install.packages("tidyverse")     # For ggplot2, dplyr, tidyr
install.packages("readxl")        # For reading Excel files
install.packages("ggrepel")       # makes text labels automatically repel each other so they don’t overlap
install.packages("concaveman")    # to outline clusters
# install.packages("pheatmap")      # For simple heatmaps
# OR for more control:
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ComplexHeatmap")

library(tidyverse)
library(readxl)
library(ggrepel)
library(concaveman)
library(dplyr)
#library(pheatmap)
library(ComplexHeatmap)
library(circlize)  # for colorRamp2

# Load data from excel sheet (Source Data) ----
# Path to your Excel file
excel_file <- "41467_2023_42762_MOESM6_ESM.xlsx"

# Fig 1c UMAP data
df_1c <- read_excel(excel_file, sheet = "Fig 1c, 2b, Supp Fig 1a")

# Fig 1d cell type labels
df_1d <- read_excel(excel_file, sheet = "Fig 1d, Supp Fig 2d")

# Fig 1e marker genes
df_1e <- read_excel(excel_file, sheet = "Fig 1e")

# Fig 1f dot plot data
df_1f <- read_excel(excel_file, sheet = "Fig 1f")

# Figure 1c: scRNA-seq UMAP ----
glimpse(df_1c) #alternative way to see its data structure

## UMAP without dashed outlines ----
# Note: "stat_density_2d" didn't work

### Solid color for each cluter -----
labels_to_show <- c(
  "Myeloid", "BC", "Mast", "NKC", "TC",
  "Mucosa", "Fibroblast", "PC_VSMC", "Endothelial",
  "Cycling", "myeSC", "nmSC"
)

labels_df <- df_1c %>%
  filter(final_label %in% labels_to_show) %>%
  group_by(final_label) %>%
  summarise(
    UMAP_1 = median(UMAP_1),
    UMAP_2 = median(UMAP_2)
  )

color_map <- c(
  "Myeloid" = "#cb181d",      # strong red
  "BC" = "#FC8D59",           # light red/orange
  "Mast" = "#F46D43",         # reddish orange
  "NKC" = "#8C510A",          # brown
  "TC" = "#FDAE61",           # orange
  "Mucosa" = "#66C2A5",       # teal
  "Fibroblast" = "#A6D854",   # light green
  "PC_VSMC" = "#1B7837",      # dark green
  "Endothelial" = "#4DAF4A",  # medium green
  "Cycling" = "#E6A3C8",      # pink
  "myeSC" = "#2166AC",        # dark blue
  "nmSC" = "#67A9CF"          # sky blue
)

groups_to_outline <- c("Myeloid", "NKC", "TC", "Fibroblast", "PC_VSMC", "Endothelial", "myeSC", "nmSC")

fig1c_plot <- ggplot(df_1c, aes(x = UMAP_1, y = UMAP_2, color = final_label)) +
  geom_point(size = 0.5, alpha = 0.8) +
  stat_density_2d(
    data = df_1c %>% filter(final_label %in% groups_to_outline),
    aes(x = UMAP_1, y = UMAP_2, group = final_label),
    color = "black",
    linetype = "dashed",
    bins = 1,
    size = 0.5
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
    legend.position = "none"
  ) +
  labs(
    title = "Fig.1C: scRNA-seq UMAP",
    x = "UMAP 1",
    y = "UMAP 2"
  ) +
  scale_color_manual(values = color_map)


# Show plot
print(fig1c_plot)

# Save fig 1c to file
ggsave(
  "Figure1c_scRNAseq_UMAP.png",
  plot = fig1c_plot,
  width = 6,    # in inches
  height = 5,   # in inches
  dpi = 300     # good for publication
)

#ggsave("Figure1c_scRNAseq_UMAP.pdf", plot = fig1c_plot, width = 6, height = 5)

### Gradient color for each cluster (needs fine tune) -----
# Cluster ID to shade mapping — approximate, visually close

my_cluster_colors <- c(
  # Myeloid shades — red family
  "0" = "#67000d",
  "1" = "#a50f15",
  "2" = "#cb181d",
  "3" = "#ef3b2c",
  "4" = "#fb6a4a",
  "5" = "#fc9272",
  "10" = "#fcbba1",
  "11" = "#fee0d2",
  "12" = "#99000d",
  "13" = "#ef6548",
  "14" = "#f16913",
  
  # Lymphoid — T/NK/B — orange/brown
  "6" = "#7f2704",
  "7" = "#d94801",
  "8" = "#f16913",
  "9" = "#fd8d3c",
  "15" = "#fdae6b",
  "16" = "#fdd0a2",
  "17" = "#fff5eb",
  
  # SC (myeSC/nmSC) — blue shades
  "18" = "#08306b",
  "19" = "#08519c",
  "20" = "#2171b5",
  "21" = "#4292c6",
  
  # Stroma — Fibroblast, PC/SMC, Endothelial — green family
  "22" = "#00441b",
  "23" = "#238b45",
  "24" = "#41ae76",
  
  # Cycling, Mucosa, Mast — pink/purple
  "25" = "#810f7c",
  "26" = "#8856a7",
  "27" = "#8c96c6",
  "28" = "#b3cde3"
)

# Use this in ggplot:
fig1c <- ggplot(df_1c, aes(UMAP_1, UMAP_2, color = factor(seurat_clusters))) +
  geom_point(size = 0.5, alpha = 0.8) +
  scale_color_manual(values = my_cluster_colors) +
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
    legend.position = "none"
  ) +
  labs(
    title = "Fig.1C: scRNA-seq UMAP",
    x = "UMAP 1",
    y = "UMAP 2"
  )

print(fig1c)

## UMAP with dashed outlines using convex hulls (not working yet)----
# define groups for outlining
group1 <- df_1c %>% filter(final_label == "Myeloid")
group2 <- df_1c %>% filter(final_label %in% c("NKC", "TC"))
group3 <- df_1c %>% filter(final_label %in% c("Fibroblast", "PC_VSMC", "Endothelial"))
group4 <- df_1c %>% filter(final_label %in% c("myeSC", "nmSC"))

# Compute convex hulls
hull1 <- concaveman(as.matrix(group1[, c("UMAP_1", "UMAP_2")])) %>%
  as.data.frame() %>%
  mutate(final_label = "Myeloid")

hull2 <- concaveman(as.matrix(group2[, c("UMAP_1", "UMAP_2")])) %>%
  as.data.frame() %>%
  mutate(final_label = "NKC_TC")
hull3 <- concaveman(as.matrix(group3[, c("UMAP_1", "UMAP_2")])) %>%
  as.data.frame() %>%
  mutate(final_label = "Fibro_PC_Endo")
hull4 <- concaveman(as.matrix(group4[, c("UMAP_1", "UMAP_2")])) %>%
  as.data.frame() %>%
  mutate(final_label = "myeSC_nmSC")

# Combine them
hulls_combined <- bind_rows(hull1, hull2, hull3, hull4)

# Plot with convex hull outlines
fig1c_plot_hulls <- ggplot(df_1c, aes(x = UMAP_1, y = UMAP_2, color = final_label)) +
  geom_point(size = 0.5, alpha = 0.8) +
  geom_polygon(
    data = hulls_combined,
    aes(x = V1, y = V2, group = final_label),
    fill = NA,
    color = "black",
    linetype = "dashed",
    size = 0.5
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
    legend.position = "none"
  ) +
  labs(
    title = "Fig.1C: scRNA-seq UMAP",
    x = "UMAP 1",
    y = "UMAP 2"
  ) +
  scale_color_manual(values = color_map)

print(fig1c_plot_hulls)

# Figure 1d: scATAC-seq UMAP ----
# Note: "stat_density_2d" didn't work
Labels_df <- df_1d %>%
  filter(Final_Label %in% labels_to_show) %>%
  group_by(Final_Label) %>%
  summarise(
    umap_1 = median(umap_1),
    umap_2 = median(umap_2) 
  )

fig1d_plot <- ggplot(df_1d, aes(x = umap_1, y = umap_2, color = Final_Label)) +
  geom_point(size = 0.5, alpha = 0.8) +
  stat_density_2d(
    data = df_1d %>% filter(Final_Label %in% groups_to_outline),
    aes(x = umap_1, y = umap_2, group = Final_Label),
    color = "black",
    linetype = "dashed",
    bins = 1,
    size = 0.5
  ) +
  geom_text_repel(
    data = Labels_df,
    aes(x = umap_1, y = umap_2, label = Final_Label),
    size = 4,
    color = "black",
    fontface = "bold",
    segment.color = "black",
    segment.size = 0.3
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "none"
  ) +
  labs(
    title = "Fig.1D: scATAC-seq UMAP",
    x = "UMAP 1",
    y = "UMAP 2"
  ) +
  scale_color_manual(values = color_map)


# Show plot
print(fig1d_plot)

# Save fig 1c to file
ggsave(
  "Figure1d_scATACseq_UMAP.png",
  plot = fig1d_plot,
  width = 6,    # in inches
  height = 5,   # in inches
  dpi = 300     # good for publication
)

# Figure 1e: Dot Plot — Marker Genes per Cell Type -----
# Make the dot plot 
# Factor: Cluster order (top to bottom)
df_1e$Cluster <- factor(df_1e$Cluster, levels = rev(c(
  "myeSC", "nmSC",
  "Fibroblast", "PC_VSMC", "Endothelial",
  "Myeloid", "Mast",
  "TC", "NKC", "BC",
  "Cycling", "Mucosa"
)))

# Factor: Group order
df_1e$Group <- factor(df_1e$Group, levels = c("Schwann", "Stroma", "Myeloid", "Lymphoid", "Other"))

# Factor: Gene order (matches paper figure)
df_1e$Gene <- factor(df_1e$Gene, levels = c(
  "PRX", "ERBB4", "MLIP", "MPZ", "SOX10", "S100B", "CHL1", "NRXN1", "NCAM1", 
  "CADM2", "ITGB8", "SCN7A", "DCN", "IGFBP5","SPARCL1", "ACTA2", "MYL9", "NOTCH3", 
  "TAGLN", "VWF", "FLT1", "EGFL7", "PTPRB", "PTPRC", "MS4A7", "CD68", "CD163", "CSF1R", "MSR1",
  "CD79A", "BANK1", "TPSAB1", "TPSB2", "CD3E", "CD8A", "IL7R", "KLRD1", "NKG7", "GNLY",
  "PRF1", "MKI67", "TOP2A", "KRT14", "SLPI"
))

# Dot plot with facet_grid (split y-axis by Group)
fig1e_plot <- ggplot(df_1e %>% filter(`% Expressing` > 0),   # 🚫 filter zeros out!
                     aes(x = Gene, y = Cluster)) +
  geom_point(
    aes(size = `% Expressing`, color = Scaled_Exp)
  ) +
  scale_size(
    range = c(0, 4)
    #breaks = c(25, 50, 75)
  ) +
  scale_color_gradient2(
    low = "lightgrey",
    mid = "orange",
    high = "darkred",
    midpoint = 0
  ) +
  theme_bw() +
  labs(
    title = "Fig.1E: Dot Plot",
    x = NULL,
    y = NULL,
    size = "% Expressing",
    color = "Avg Scaled Exp"
  ) +
  facet_grid(Group ~ ., scales = "free_y", space = "free_y") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.background = element_rect(fill = "grey90", color = NA),
    strip.text.y = element_text(face = "bold")
  )

print(fig1e_plot)


# Save as PNG
ggsave(
  "Figure1e_DotPlot.png",
  plot = fig1e_plot,
  width = 8,
  height = 5,
  dpi = 300
)

# Figure 1f: Heatmap — Mouse Nerve Cell Meta-Signatures -----
## Heatmap with ComplexHeatmap -----
df_1f$Group <- c(
  "Schwann",
  rep("Fibroblast", 4),
  rep("Vascular", 3),
  rep("Immune", 2),
  "Cycling"
)

# Order rows and clusters
row_order <- c("Schwann",
                "Fibroblast", "Endoneurial", "Perineurial", "Epineurial",
                "Pericyte/VSMC", "Endothelial", "Lymphatic",
                "Myeloid", "Lymphoid",
                "Cycling")

column_order <- c(
  "nmSC", "myeSC",
  "Fibroblast", "PC_VSMC", "Endothelial",
  "Myeloid", "Mast", "TC", "NKC", "BC",
  "Cycling", "Mucosa")

# Make sure row order is a factor
df_1f$`Mouse Peripheral Nerve Cell Type` <- factor(df_1f$`Mouse Peripheral Nerve Cell Type`, levels = row_order)

# Reorder data frame by row
df_1f <- df_1f %>% arrange(`Mouse Peripheral Nerve Cell Type`)

# build matrix
mat <- df_1f %>%
  select(all_of(column_order)) %>%
  as.matrix()

rownames(mat) <- df_1f$`Mouse Peripheral Nerve Cell Type`

# Color scale
col_fun <- circlize::colorRamp2(
  c(-2, 0, 2),
  c("blue", "white", "red")
)

# Row split for group labels
row_split <- factor(df_1f$Group, levels = c("Schwann", "Fibroblast", "Vascular", "Immune", "Cycling"))

# Heatmap
ht <- Heatmap(
  mat,
  name = "Mean\nModule\nScore",
  col = col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_split = row_split,
  row_gap = unit(2, "mm"),
  rect_gp = gpar(col = "black"),  # grid lines
  column_names_rot = 45,
  row_title = "Mouse Peripheral Nerve Cell Type",
  column_title = "VS TME Celltype Clusters"
)

# Draw heatmap and save as png
png("Figure1f_Heatmap.png",
    width = 8,    # inches
    height = 6,   # inches
    units = "in", # important!
    res = 300     # dots per inch
)
draw(ht)
dev.off()



## Heatmap with ggplot -----
# Reshape: wide → long
df_1f_long <- df_1f %>%
  pivot_longer(
    cols = -`Mouse Peripheral Nerve Cell Type`,
    names_to = "Cluster",
    values_to = "Mean"
  )

# Order rows and clusters
df_1f_long$`Mouse Peripheral Nerve Cell Type` <- factor(df_1f_long$`Mouse Peripheral Nerve Cell Type`,
                                                        levels = rev(c(
                                                          "Schwann",
                                                          "Fibroblast", "Endoneurial", "Perineurial", "Epineurial",
                                                          "Pericyte/VSMC", "Endothelial", "Lymphatic",
                                                          "Myeloid", "Lymphoid",
                                                          "Cycling"
                                                        ))
)

df_1f_long$Cluster <- factor(df_1f_long$Cluster, levels = c(
  "nmSC", "myeSC",
  "Fibroblast", "PC_VSMC", "Endothelial",
  "Myeloid", "Mast", "TC", "NKC", "BC",
  "Cycling", "Mucosa"
))

# Heatmap
fig1f_plot <- ggplot(df_1f_long, aes(x = Cluster, y = `Mouse Peripheral Nerve Cell Type`, fill = Mean)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,
    limits = c(-2, 2),
    name = "Mean\nModule\nScore"
  ) +
  theme_minimal() +
  labs(
    title = "Fig.1F: Heatmap",
    x = "VS TME Celltype Clusters",
    y = "Mouse Peripheral Nerve Cell Type"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  )

print(fig1f_plot)

ggsave(
  "Figure1f_Heatmap.png",
  plot = fig1f_plot,
  width = 8,
  height = 6,
  dpi = 300
)
