library(tidyverse)
library(ggalluvial)
library(grid)
library(viridis)
library(RColorBrewer)
library(rcartocolor)
library(readxl)

mouse_cluster_file <- "data/mouse/derivatives/v3/107/clusters/resolution_0.2/clusters.csv"
df_mouse_clusters <- read_csv(mouse_cluster_file, show_col_types = FALSE)

palette_length <- 1000
# palette_cols <- colorRampPalette(c("springgreen4", "sienna1", "darkorchid1"))(palette_length)
palette_cols <- colorRampPalette(brewer.pal(n = 11, name = "Spectral"))(
  palette_length
)

# palette_cols <- viridis(n = 1000, begin = 0.9, end = 0)
# palette_cols <- magma(n = 1000, begin = 0.8, end = 0.2)
# palette_cols <- colorRampPalette(rev(brewer.pal(n = 11, name = "PRGn")))(palette_length)
# palette_cols <- colorRampPalette(rev(brewer.pal(n = 11, name = "PiYG")))(palette_length)

# df_mouse_clusters %>%
#   group_by(nk7, nk2) %>%
#   count() %>%
#   filter(nk7 %in% c(1, 3, 4, 7), )

nk_subset <- 2:10

list_palette_x <- vector(mode = "list", length = length(nk_subset))
list_palette_x[[1]] <- tibble(nk = 2, k = 1:2, x = c(1, palette_length))
for (i in 2:length(nk_subset)) {
  nk <- nk_subset[i]
  nk_prev <- nk_subset[i - 1]

  cols <- paste0("nk", c(nk_prev, nk))
  df_clusters_nk <- df_mouse_clusters[, cols]
  colnames(df_clusters_nk) <- c("nk_prev", "nk_current")

  df_palette_x <- list_palette_x[[i - 1]]

  list_palette_x[[i]] <- df_clusters_nk %>%
    group_by(nk_current) %>%
    mutate(n_current = n()) %>%
    ungroup() %>%
    group_by(nk_current, nk_prev) %>%
    mutate(n_current_prev = n()) %>%
    ungroup() %>%
    mutate(prop = n_current_prev / n_current) %>%
    distinct() %>%
    left_join(df_palette_x, by = c("nk_prev" = "k")) %>%
    group_by(nk_current) %>%
    summarise(x_new = weighted.mean(x, w = prop)) %>%
    mutate(nk = nk) %>%
    select(nk, k = nk_current, x = x_new)
}

df_palette_x <- bind_rows(list_palette_x) %>%
  mutate(
    x = round(x),
    colors = palette_cols[x],
    cluster_id = paste0(nk, "-", k)
  )

palette <- df_palette_x$colors
names(palette) <- df_palette_x$cluster_id

# nk_subset <- c(2, 7, 10)

# Convert cluster information to long format
df_mouse_clusters_long <- df_mouse_clusters %>%
  pivot_longer(cols = -ID, names_to = "nk_name", values_to = "k") %>%
  mutate(
    nk = str_remove(nk_name, "nk"),
    nk = as.numeric(nk),
    cluster_id = paste0(nk, "-", k)
  ) %>%
  filter(nk %in% nk_subset)
# left_join(df_mouse_clusters %>%
#             select(ID, nk10), by = "ID")

# df_mouse_clusters_long <- df_mouse_clusters_long %>%
#   left_join(df_cluster_groups, by = c("nk", "k"))

cluster_ids <- df_mouse_clusters_long %>%
  select(nk, k, cluster_id) %>%
  distinct() %>%
  arrange(nk, k) %>%
  pull(cluster_id)

df_mouse_sankey <- df_mouse_clusters_long %>%
  mutate(
    nk = factor(nk),
    k = factor(k),
    cluster_id = factor(cluster_id, levels = cluster_ids)
  )

stratum_width <- 0.5

# palette <- c("A" = "springgreen4",
#              "B" = "sienna2",
#              "C" = "seagreen2",
#              "D" = "darkorchid1")
#
# palette <- c("2-1" = "springgreen4",
#              "2-2" = "darkorchid1",
#              "7-1" = "springgreen4",
#              "7-2" = "grey80",
#              "7-3" = "seagreen2",
#              "7-4" = "#446F87",
#              "7-5" = "grey80",
#              "7-6" = "grey80",
#              "7-7" = "darkorchid1",
#              "10-1" = "springgreen4",
#              "10-2" = "grey80",
#              "10-3" = "seagreen2",
#              "10-4" = "darkorchid1",
#              "10-5" = "#446F87",
#              "10-6" = "grey80",
#              "10-7" = "grey80",
#              "10-8" = "grey80",
#              "10-9" = "grey80",
#              "10-10" = "grey80")

df_palette <- read_excel("figures/v3/mouse_cluster_colours.xlsx")
colour_set <- 2
palette <- df_palette[[paste0("colour_", colour_set)]]
names(palette) <- df_palette$cluster_id


col2hex <- function(x, alpha = FALSE) {
  args <- as.data.frame(t(col2rgb(x, alpha = alpha)))
  args <- c(args, list(names = x, maxColorValue = 255))
  do.call(rgb, args)
}

col2hex("springgreen4")
col2hex("darkorchid1")
col2hex("seagreen2")


# Mouse cluster Sankey plot
plt_mouse_sankey <- ggplot(
  df_mouse_sankey,
  aes(
    x = nk,
    stratum = k,
    alluvium = ID,
    # fill = group,
    fill = cluster_id,
    label = k
  )
) +
  geom_flow(
    stat = "alluvium",
    aes.flow = "forward",
    width = stratum_width,
    col = "grey30",
    linewidth = 0.15
  ) +
  geom_stratum(alpha = 1, width = stratum_width) +
  # geom_stratum(mapping = aes(col = group), alpha = 1, width = stratum_width) +
  # scale_x_discrete(position = "top") +
  scale_x_discrete(
    position = "top",
    expand = expansion(
      add = c(stratum_width / 2 + 0.020, stratum_width / 2 + 0.005)
    )
  ) +
  # scale_x_discrete(position = "top", expand = expansion(mult = c(0.0, 0.0))) +
  scale_y_continuous(
    breaks = c(seq(0, 130, by = 20), nrow(df_mouse_clusters))
  ) +
  scale_fill_manual(values = palette, na.value = "grey80") +
  # scale_colour_manual(values = palette) +
  labs(x = "K", y = "Number of models") +
  # theme_void() +
  theme_bw() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none",
    panel.border = element_blank(),
    axis.line.y.left = element_line(colour = "grey50"),
    panel.ontop = FALSE
  )


outfile <- paste0("mouse_sankey_", colour_set, ".pdf")
outfile <- file.path("figures", "v3", outfile)
pdf(file = outfile, width = unit(10, "in"), height = unit(5, "in"))
print(plt_mouse_sankey)
dev.off()
