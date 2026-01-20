library(tidyverse)


pipeline_dir <- "data/cross_species/v3/"
similarity_mlp <- compute_similarity_significance(
  similarity = import_similarity(param_id = "375", pipeline_dir = pipeline_dir),
  permutations = import_similarity_permutations(
    param_id = "375",
    pipeline_dir = pipeline_dir
  )
)

similarity_vae <- compute_similarity_significance(
  similarity = import_similarity(param_id = "852", pipeline_dir = pipeline_dir),
  permutations = import_similarity_permutations(
    param_id = "852",
    pipeline_dir = pipeline_dir
  )
)

similarity_mlp
similarity_vae

df_compare <- tibble(
  mlp = similarity_mlp$similarity,
  vae = similarity_vae$similarity
)

ggplot(df_compare, aes(x = mlp, y = vae)) +
  geom_point()


df_compare_pval <- tibble(
  mlp = similarity_mlp$pval,
  vae = similarity_vae$pval
) %>%
  mutate(
    match_mlp = mlp < 0.05,
    match_vae = vae < 0.05,
    nlp_mlp = -log10(mlp),
    nlp_vae = -log10(vae)
  ) %>%
  mutate(
    nlp_mlp = ifelse(is.infinite(nlp_mlp), 8, nlp_mlp),
    nlp_vae = ifelse(is.infinite(nlp_vae), 8, nlp_vae)
  )

ggplot(df_compare_pval, aes(x = mlp, y = vae)) +
  geom_point()

# Convert cluster information to long format
df_mouse_clusters_long <- df_mouse_clusters %>%
  pivot_longer(cols = -ID, names_to = "nk_name", values_to = "k") %>%
  mutate(nk = str_remove(nk_name, "nk"), nk = as.numeric(nk)) %>%
  # filter(nk %in% nk_subset) %>%
  left_join(
    df_mouse_clusters %>%
      select(ID, nk10),
    by = "ID"
  )


mouse_matches_mlp <- similarity_mlp %>%
  select(
    mouse_cluster_id = img2_cluster_id,
    human_cluster_id = img1_cluster_id,
    similarity,
    pval
  ) %>%
  mutate(match = pval < 0.05) %>%
  filter(match) %>%
  pull(mouse_cluster_id) %>%
  unique()

mouse_matches_vae <- similarity_vae %>%
  select(
    mouse_cluster_id = img2_cluster_id,
    human_cluster_id = img1_cluster_id,
    similarity,
    pval
  ) %>%
  mutate(match = pval < 0.05) %>%
  filter(match) %>%
  pull(mouse_cluster_id) %>%
  unique()

mouse_sankey_mlp <- df_mouse_clusters_long %>%
  mutate(
    cluster_id = paste0(nk, "-", k),
    match = cluster_id %in% mouse_matches_mlp
  )

mouse_sankey_vae <- df_mouse_clusters_long %>%
  mutate(
    cluster_id = paste0(nk, "-", k),
    match = cluster_id %in% mouse_matches_vae
  )


similarity_mlp %>%
  select(
    mouse_cluster_id = img2_cluster_id,
    human_cluster_id = img1_cluster_id,
    similarity,
    pval
  ) %>%
  mutate(match = pval < 0.05) %>%
  filter(mouse_cluster_id == "8-1", human_cluster_id %in% c("9-5", "9-8"))

similarity_vae %>%
  select(
    mouse_cluster_id = img2_cluster_id,
    human_cluster_id = img1_cluster_id,
    similarity,
    pval
  ) %>%
  mutate(match = pval < 0.05) %>%
  filter(mouse_cluster_id == "8-1", human_cluster_id %in% c("9-5", "9-8"))


pdf(file = "mlp_HBN.pdf")
ggplot(
  mouse_sankey_mlp,
  aes(x = nk, stratum = k, alluvium = ID, fill = match, label = k)
) +
  geom_flow(stat = "alluvium", aes.flow = "forward", width = stratum_width) +
  geom_stratum(alpha = 1, width = stratum_width) +
  # geom_stratum(mapping = aes(col = group), alpha = 1, width = stratum_width) +
  # scale_x_discrete(position = "top") +
  scale_x_discrete(position = "top") +
  # scale_x_discrete(position = "top", expand = expansion(mult = c(0.0, 0.0))) +
  scale_y_continuous(
    breaks = c(seq(0, 130, by = 20), nrow(df_mouse_clusters))
  ) +
  # scale_fill_manual(values = palette, na.value = "grey80") +
  # scale_colour_manual(values = palette) +
  labs(x = "K", y = "Number of models") +
  # theme_void() +
  theme_bw()
dev.off()

pdf(file = "vae_HBN.pdf")
ggplot(
  mouse_sankey_vae,
  aes(x = nk, stratum = k, alluvium = ID, fill = match, label = k)
) +
  geom_flow(stat = "alluvium", aes.flow = "forward", width = stratum_width) +
  geom_stratum(alpha = 1, width = stratum_width) +
  # geom_stratum(mapping = aes(col = group), alpha = 1, width = stratum_width) +
  # scale_x_discrete(position = "top") +
  scale_x_discrete(position = "top") +
  # scale_x_discrete(position = "top", expand = expansion(mult = c(0.0, 0.0))) +
  scale_y_continuous(
    breaks = c(seq(0, 130, by = 20), nrow(df_mouse_clusters))
  ) +
  # scale_fill_manual(values = palette, na.value = "grey80") +
  # scale_colour_manual(values = palette) +
  labs(x = "K", y = "Number of models") +
  # theme_void() +
  theme_bw()
dev.off()


similarity_vae %>%
  filter(img2_nk == 6)


perms_mlp <- import_similarity_permutations(
  param_id = "375",
  pipeline_dir = pipeline_dir
)

perms_vae <- import_similarity_permutations(
  param_id = "852",
  pipeline_dir = pipeline_dir
)

tmp <- bind_rows(
  perms_mlp %>%
    filter(img2_cluster_id == "6-1") %>%
    mutate(model = "MLP"),
  perms_vae %>%
    filter(img2_cluster_id == "6-1") %>%
    mutate(model = "VAE")
)

ggplot(tmp, aes(x = similarity, fill = model)) +
  geom_histogram(position = "dodge")

# df_nk9 <- bind_rows(

df_nk9 <- inner_join(
  similarity_mlp %>%
    filter(img2_nk == 8) %>%
    select(
      img1_cluster_id,
      img2_cluster_id,
      mlp_similarity = similarity,
      mlp_pval = pval
    ),
  similarity_vae %>%
    filter(img2_nk == 8) %>%
    select(
      img1_cluster_id,
      img2_cluster_id,
      vae_similarity = similarity,
      vae_pval = pval
    ),
  by = c("img1_cluster_id", "img2_cluster_id")
)

ggplot(df_nk9, aes(x = mlp_pval, y = vae_pval)) +
  geom_point() +
  geom_hline(yintercept = 0.05) +
  geom_vline(xintercept = 0.05)


df_nk9 %>%
  filter(mlp_pval < 0.05) %>%
  arrange(img2_cluster_id)

ggplot(similarity_vae, aes(x = img2_nk, y = similarity)) +
  geom_jitter()

similarity_mlp %>%
  filter(img2_cluster_id == "8-1") %>%
  arrange(pval)
