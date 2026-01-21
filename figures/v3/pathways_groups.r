library(tidyverse)
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(grid))
suppressPackageStartupMessages(library(gridExtra))


# Number of bigpts in an inch
pt_per_in <- 72

# Font family
font_family <- "Helvetica"

# Nature suggested font size: 5-7 pt
font_size <- 6

# Empty rectangle grob
empty_rect_grob <- rectGrob(gp = gpar(fill = NA))

# Black rectangle grob
black_rect_grob <- rectGrob(gp = gpar(fill = "black"))

# Maximal dimensions of figure in bigpts
fig2_total_width <- 510
fig2_total_height <- 481

fig2_border_padding_width <- 4
fig2_border_padding_height <- fig2_border_padding_width

fig2_width <- fig2_total_width - 2 * fig2_border_padding_width


#
file <- "data/mouse/derivatives/v3/107/clusters/resolution_0.2/clusters.csv"
#
clusters <- read_csv(file)
#
#
# nk_max <- 10
# list_clusters <- vector(mode = "list", length = nk_max-2)
# names(list_clusters) <- nk_max:3
# for (nk in nk_max:3) {
#
#   cols <- paste0("nk", c(nk, nk-1))
#   clusters_nk <- clusters[,cols]
#   colnames(clusters_nk) <- c("target_k", "input_k")
#
#   list_clusters_nk <- vector(mode = "list", length = nk)
#   for (k in 1:nk) {
#     idx <- clusters_nk[["target_k"]] == k
#     list_clusters_nk[[k]] <- clusters_nk[idx,] %>%
#       group_by(input_k) %>%
#       summarise(prop = n()/nrow(.),
#                 .groups = "drop") %>%
#       mutate(target_k = k, target_nk = nk, input_nk = nk-1)
#   }
#
#   list_clusters[[as.character(nk)]] <- bind_rows(list_clusters_nk)
#
# }
#
#
# list_clusters[["10"]] %>%
#   group_by(target_k) %>%
#   filter(prop == max(prop)) %>%
#   arrange(target_k)
#
# list_clusters[["3"]] %>%
#   group_by(target_k) %>%
#   filter(prop == max(prop)) %>%
#   arrange(input_k)
#
#
# for (nk in nk_max:3) {
#   nk <- 10
#
#   df_clusters_nk <- list_clusters[[as.character(nk)]] %>%
#     group_by(target_k) %>%
#     filter(prop == max(prop)) %>%
#     arrange(target_k)
#
#   for (k in 1:nk) {
#
#
#     df_clusters_nk_k <- df_clusters_nk %>%
#       filter(target_k == k)
#
#     print(nrow(df_clusters_nk_k))
#
#   }
# }
#
# tmp %>%
#   filter(target_k == k)
#

SRCPATH <- Sys.getenv("SRCPATH")
PROJECTPATH <- Sys.getenv("PROJECTPATH")

source(file.path(SRCPATH, "utils.R"))
source(file.path(SRCPATH, "processing.R"))
source(file.path(SRCPATH, "analysis.R"))
source(file.path(SRCPATH, "enrichment.R"))


mouse_pipeline_dir <- file.path(PROJECTPATH, "data/mouse/derivatives/")
mouse_pipeline_dir <- file.path(mouse_pipeline_dir, "v3", "107")


enrichment_dir <- file.path(PROJECTPATH, "data", "enrichment")

# Path to Reactome hierarchy file
reactome_hierarchy <- file.path(enrichment_dir, "reactome_hierarchy.csv")

# Import Reactome hierarchy
df_reactome_hierarchy <- read_csv(reactome_hierarchy, show_col_types = FALSE)

# Initialize tree levels for root pathways
df_reactome_hierarchy <- df_reactome_hierarchy %>%
  mutate(level = ifelse(is.na(ParentID), 0, NA))

# Iterate lvl down the tree
nmissing <- sum(is.na(df_reactome_hierarchy$level))
lvl <- 0
while (nmissing > 0) {
  df_reactome_lvl <- df_reactome_hierarchy %>%
    filter(level == lvl)

  # Iterate over nodes at given level
  for (i in 1:nrow(df_reactome_lvl)) {
    node <- df_reactome_lvl[[i, "Name"]]

    df_reactome_hierarchy <- df_reactome_hierarchy %>%
      mutate(
        level = ifelse(is.na(Parent), 0, ifelse(Parent == node, lvl + 1, level))
      )
  }

  nmissing <- sum(is.na(df_reactome_hierarchy$level))
  lvl <- lvl + 1
}

# Clean up names
df_reactome_hierarchy <- df_reactome_hierarchy %>%
  mutate(
    Name = str_replace_all(Name, "  ", " "),
    Name = str_replace_all(Name, "/", " ")
  )

# Bader Reactome gene sets
module_file <- file.path(
  enrichment_dir,
  "Human_Reactome_June_01_2025_symbol.gmt"
)

# Import Bader sets
df_modules <- get_module_sizes(modules = module_file)

df_modules <- df_modules %>%
  rename(IDBader = ID) %>%
  mutate(
    ID = IDBader %>%
      str_split_i("%", i = -1) %>%
      str_remove("R-HSA-") %>%
      str_remove("\\.[0-9]+"),
    ID = paste0("R-HSA-", ID)
  ) %>%
  mutate(Title = str_replace_all(Title, "  ", " "))

df_reactome <- inner_join(df_reactome_hierarchy, df_modules, by = "ID")

# Pathways to keep
pathway_ids_keep <- df_reactome %>%
  filter(B > 10) %>%
  pull(IDBader)


# Database versions
stringdb_version <- "12.0"
bader_version <- "2025"
stringdb_threshold <- 950

# Base directory for pathway data files
pathways_dir <- file.path(mouse_pipeline_dir, "enrichment")
pathways_dir <- file.path(
  pathways_dir,
  paste("StringDB", stringdb_version, "Bader", bader_version, sep = "_")
)
pathways_dir <- file.path(pathways_dir, stringdb_threshold)
pathways_dir <- file.path(pathways_dir, "NeighbourhoodEnrichment")

if (length(list.files(pathways_dir)) == 0) {
  stop("No files in specified directory")
}

# Prefix for pathway data files
pathways_file_prefix <- "cluster_pathway_enrichment"

nk_max <- 10

# Iterate over cluster solutions and import mouse pathway enrichment files
list_pathways <- vector(mode = "list", length = nk_max - 1)
names(list_pathways) <- 2:nk_max
for (nk in 2:nk_max) {
  # Iterate over cluster number
  list_pathways[[nk - 1]] <- vector(mode = "list", length = nk)
  for (k in 1:nk) {
    pathways_file <- paste(
      pathways_file_prefix,
      nk,
      k,
      stringdb_threshold,
      sep = "_"
    )
    pathways_file <- paste0(pathways_file, ".csv")
    pathways_file <- file.path(pathways_dir, pathways_file)
    list_pathways[[nk - 1]][[k]] <- read_csv(
      pathways_file,
      show_col_types = FALSE
    ) %>%
      semi_join(df_reactome, by = c("ID" = "IDBader")) %>%
      left_join(
        df_reactome %>%
          select(ID = IDBader, level, Parent, Root),
        by = "ID"
      ) %>%
      filter(ID %in% pathway_ids_keep) %>%
      filter(!(Root == "Signal Transduction" & Title == "Integrin signaling"))
  }

  # Combine clusters per solution
  list_pathways[[nk - 1]] <- bind_rows(list_pathways[[nk - 1]]) %>%
    mutate(nk = nk)
}

# Reduce all pathway data frames into one
df_pathways_all <- bind_rows(list_pathways)

df_pathways_all <- df_pathways_all %>%
  rename(pathway = Title) %>%
  mutate(NLQ = -log10(adj.P.Val)) %>%
  mutate(
    pathway = ifelse(
      pathway == "Signaling Pathways",
      "Signal Transduction",
      pathway
    )
  ) %>%
  unite(col = "cluster_id", nk, k, sep = "-", remove = FALSE)


pathways_root_excl <- c(
  "Circadian clock",
  "Disease",
  "Digestion and absorption",
  "DNA Repair",
  "DNA Replication",
  "Drug ADME",
  "Muscle contraction",
  "Protein localization",
  "Reproduction",
  "Sensory Perception",
  "Transport of small molecules"
)

pathways_root_w_lvl_1 <- c(
  "Chromatin organization",
  "Gene expression (Transcription)",
  "Immune System",
  "Neuronal System",
  "Signal Transduction"
)

df_pathways_lvl_0 <- df_pathways_all %>%
  filter(level == 0) %>%
  filter(!(Root %in% pathways_root_excl))

df_pathways_lvl_1 <- df_pathways_all %>%
  filter(Root %in% pathways_root_w_lvl_1, level == 1)

lvl_1_enriched <- df_pathways_lvl_1 %>%
  mutate(significant = adj.P.Val < 0.01) %>%
  group_by(Root, pathway) %>%
  summarise(n_significant = sum(significant)) %>%
  ungroup() %>%
  filter(n_significant > 0) %>%
  pull(pathway)

pathway_roots <- df_pathways_lvl_0 %>%
  pull(Root) %>%
  unique() %>%
  sort()

pathways_enriched <- c(pathway_roots, lvl_1_enriched)

df_pathways_lvl_0_1 <- bind_rows(df_pathways_lvl_0, df_pathways_lvl_1)

df_pathways_lvl_0_1 <- df_pathways_lvl_0_1 %>%
  filter(pathway %in% pathways_enriched)

list_pathways_lvl_1 <- vector(mode = "list", length = length(pathway_roots))
names(list_pathways_lvl_1) <- pathway_roots
for (i in 1:length(list_pathways_lvl_1)) {
  lvls <- df_pathways_lvl_0_1 %>%
    filter(Root == pathway_roots[i], level == 1) %>%
    pull(pathway) %>%
    unique() %>%
    sort()
  list_pathways_lvl_1[[i]] <- c(pathway_roots[i], lvls)
}

pathway_lvls <- reduce(list_pathways_lvl_1, c)
pathway_lvls_trunc <- str_trunc(pathway_lvls, width = 40)


mat_pathways <- df_pathways_lvl_0_1 %>%
  group_by(pathway) %>%
  mutate(NLQ_norm = NLQ / max(NLQ)) %>%
  ungroup() %>%
  mutate(NLQ_norm = ifelse(is.nan(NLQ_norm), 0, NLQ_norm)) %>%
  select(pathway, cluster_id, NLQ_norm) %>%
  pivot_wider(
    id_cols = pathway,
    names_from = cluster_id,
    values_from = NLQ_norm
  ) %>%
  column_to_rownames("pathway") %>%
  as.matrix()


# Hierarchical clustering of the pathways
pathway_hc <- hclust(d = dist(mat_pathways, method = "euclidean"))

# Extract pathway order according to clustering
pathway_lvls <- pathway_hc[["labels"]]
pathway_order <- pathway_hc[["order"]]
pathway_lvls_clustered <- pathway_lvls[pathway_order]

pathway_lvls_clustered_trunc <- str_trunc(pathway_lvls_clustered, width = 40)

# Obtain pathway cluster order at selected solution
df_pathway_hclust <- cutree(pathway_hc, k = 6) %>%
  enframe(name = "pathway", value = "pathway_cluster")

# Relevel pathways to follow dendrogram order
df_pathway_cluster_lvls <- df_pathway_hclust %>%
  mutate(pathway = factor(pathway, levels = pathway_lvls_clustered)) %>%
  arrange(pathway) %>%
  select(pathway_cluster) %>%
  distinct() %>%
  mutate(pathway_cluster_new = 1:nrow(.))

df_pathway_hclust <- df_pathway_hclust %>%
  left_join(df_pathway_cluster_lvls, by = "pathway_cluster") %>%
  select(-pathway_cluster, pathway_cluster = pathway_cluster_new)


df_pathways_heatmap <- df_pathways_lvl_0_1 %>%
  left_join(df_pathway_hclust, by = "pathway") %>%
  mutate(
    pathway_trunc = str_trunc(pathway, width = 40),
    pathway_trunc = factor(pathway_trunc, levels = pathway_lvls_clustered_trunc)
  )


NLQ_max <- df_pathways_heatmap %>%
  pull(NLQ) %>%
  max()

# Clamped enrichment statistic
NLQ_threshold <- 30
df_pathways_heatmap <- df_pathways_heatmap %>%
  mutate(intensity = ifelse(NLQ > NLQ_threshold, NLQ_threshold, NLQ))

list_cluster_groups <- list(
  A = c("5-1", "6-1", "7-1", "8-1", "9-1", "10-1"),
  B = c("5-3", "6-3", "7-3", "8-5", "9-4", "10-3"),
  C = c("5-5", "6-4", "7-4", "8-2", "9-5", "10-5"),
  D = c("5-4", "6-5", "7-7", "8-6", "9-2", "10-4"),
  E = c("5-1", "6-6", "7-6", "8-7", "9-9", "10-9"),
  G = c("5-2", "6-2", "7-2", "8-8", "9-6", "10-2")
)

names(list_cluster_groups) <- c(
  "A (dark green)",
  "B (orange)",
  "C (light green)",
  "D (purple)",
  "E",
  "F"
)

# cluster_ids <- c("5-1", "6-1", "7-1", "8-1", "9-1", "10-1")
# cluster_ids <- c("5-3", "6-3", "7-3", "8-5", "9-4", "10-3")
# cluster_ids <- c("5-5", "6-4", "7-4", "8-2", "9-5", "10-5")
# cluster_ids <- c("5-4", "6-5", "7-7", "8-6", "9-2", "10-4")

list_pathways_groups <- vector(
  mode = "list",
  length = length(list_cluster_groups)
)
for (i in 1:length(list_pathways_groups)) {
  list_pathways_groups[[i]] <- df_pathways_heatmap %>%
    filter(cluster_id %in% list_cluster_groups[[i]]) %>%
    mutate(group = names(list_cluster_groups)[i])

  # cluster_lvls <- df_pathways_heatmap %>%
  #   select(cluster_id, nk) %>%
  #   distinct() %>%
  #   arrange(nk) %>%
  #   pull(cluster_id)
  #
  # df_pathways_heatmap <- df_pathways_heatmap %>%
  #   mutate(cluster_id = factor(cluster_id, levels = cluster_lvls))
}

df_pathways_groups <- bind_rows(list_pathways_groups)
cluster_lvls <- unique(reduce(list_cluster_groups, c))

df_pathways_groups <- df_pathways_groups %>%
  mutate(cluster_id = factor(cluster_id, levels = cluster_lvls))

# df_pathways_heatmap <- df_pathways_heatmap %>%
#   filter(cluster_id %in% cluster_ids)
#
# cluster_lvls <- df_pathways_heatmap %>%
#   select(cluster_id, nk) %>%
#   distinct() %>%
#   arrange(nk) %>%
#   pull(cluster_id)
#
# df_pathways_heatmap <- df_pathways_heatmap %>%
#   mutate(cluster_id = factor(cluster_id, levels = cluster_lvls))

heatmap_limits <- c(2, NLQ_threshold)
heatmap_fill_lab <- "Enrichment (-log10(q))"

# Heatmap palette
heatmap_palette_cols <- brewer.pal(n = 9, name = "OrRd")
heatmap_palette <- colorRampPalette(colors = heatmap_palette_cols)(255)

# Heatmap legend width
fig2_heatmap_pathways_legend_width_pt <- 50

plt <- ggplot(
  df_pathways_groups,
  aes(x = cluster_id, y = fct_rev(pathway_trunc), fill = NLQ)
) +
  geom_tile(col = "grey60") +
  facet_grid(
    pathway_cluster ~ group,
    scales = "free",
    space = "free",
    switch = "y"
  ) +
  scale_x_discrete(expand = expansion()) +
  scale_y_discrete(expand = expansion(), position = "right") +
  scale_fill_gradientn(
    colors = heatmap_palette,
    limits = heatmap_limits,
    na.value = "grey85",
    guide = guide_colourbar(
      title.position = "top",
      title.hjust = 0.5,
      barwidth = unit(fig2_heatmap_pathways_legend_width_pt, "bigpts"),
      barheight = unit(10, "bigpts")
    )
  ) +
  labs(x = "Mouse cluster", y = "Pathway", fill = heatmap_fill_lab) +
  theme_bw() +
  theme(
    strip.background.y = element_blank(),
    strip.text.y = element_blank(),
    panel.spacing = unit(2, "bigpts"),
    axis.title.x = element_text(size = font_size, family = font_family),
    axis.text.x = element_text(size = font_size, family = font_family),
    axis.ticks.x = element_line(size = 0.25),
    axis.title.y = element_text(size = font_size, family = font_family),
    axis.text.y = element_text(size = font_size - 1, family = font_family),
    axis.ticks.y = element_line(size = 0.25),
    legend.position = "bottom",
    legend.title = element_text(size = font_size - 1, family = font_family),
    legend.text = element_text(size = font_size - 1, family = font_family),
    legend.margin = margin(),
    plot.margin = margin(l = 0)
  )


outfile <- pdf(
  "heatmap_groups_clustered.pdf",
  width = unit(10, "in"),
  height = unit(10, "in")
)
print(plt)
dev.off()
