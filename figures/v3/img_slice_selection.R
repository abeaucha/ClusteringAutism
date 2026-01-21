pond_img <- import_cluster_map(
  imgdir = human_centroid_dirs[1],
  mask = human_mask_file,
  nk = 7,
  k = 1
)
# ,
#                           threshold = threshold,
#                           threshold_value = threshold_value,
#                           threshold_symmetric = threshold_symmetric)

# Compute overlay thresholds
overlay_threshold <- numeric(2)
overlay_threshold[1] <- min(abs(pond_img[pond_img != 0]))
overlay_threshold[2] <- 0.8 * max(abs(pond_img[pond_img != 0]))
overlay_threshold <- round(overlay_threshold, 2)


pond_img <- mincArray(pond_img)
pond_img <- pond_img[slices_dim_1, slices_dim_2, slices_dim_3]

nk_vec <- c(5, 5, 7, 7)
k_vec <- c(3, 4, 2, 6)

list_hbn_imgs <- vector(mode = "list", length = 4)
for (i in 1:4) {
  img <- import_cluster_map(
    imgdir = human_centroid_dirs[2],
    mask = human_mask_file,
    nk = nk_vec[i],
    k = k_vec[i]
  )
  # ,
  #      threshold = threshold,
  #      threshold_value = threshold_value,
  #      threshold_symmetric = threshold_symmetric)

  # Compute overlay thresholds
  overlay_threshold <- numeric(2)
  overlay_threshold[1] <- min(abs(img[img != 0]))
  overlay_threshold[2] <- 0.8 * max(abs(img[img != 0]))
  overlay_threshold <- round(overlay_threshold, 2)

  img <- mincArray(img)
  img <- img[slices_dim_1, slices_dim_2, slices_dim_3]

  list_hbn_imgs[[i]] <- img
}

overlay_threshold <- c(0.1, 1.0)

slices <- seq(30, 160, length.out = 10)
sliceSeries(nrow = 10, ncol = 1, dimension = human_slc_dim, slices = slices) %>%
  anatomy(
    human_anat_vol_cropped,
    low = human_anat_low,
    high = human_anat_high
  ) %>%
  overlay(
    pond_img,
    low = overlay_threshold[1],
    high = overlay_threshold[2],
    symmetric = TRUE
  ) %>%
  sliceSeries(
    nrow = 10,
    ncol = 1,
    dimension = human_slc_dim,
    slices = slices
  ) %>%
  anatomy() %>%
  overlay(
    list_hbn_imgs[[1]],
    low = overlay_threshold[1],
    high = overlay_threshold[2],
    symmetric = TRUE
  ) %>%
  sliceSeries(
    nrow = 10,
    ncol = 1,
    dimension = human_slc_dim,
    slices = slices
  ) %>%
  anatomy() %>%
  overlay(
    list_hbn_imgs[[2]],
    low = overlay_threshold[1],
    high = overlay_threshold[2],
    symmetric = TRUE
  ) %>%
  sliceSeries(
    nrow = 10,
    ncol = 1,
    dimension = human_slc_dim,
    slices = slices
  ) %>%
  anatomy() %>%
  overlay(
    list_hbn_imgs[[3]],
    low = overlay_threshold[1],
    high = overlay_threshold[2],
    symmetric = TRUE
  ) %>%
  sliceSeries(
    nrow = 10,
    ncol = 1,
    dimension = human_slc_dim,
    slices = slices
  ) %>%
  anatomy() %>%
  overlay(
    list_hbn_imgs[[4]],
    low = overlay_threshold[1],
    high = overlay_threshold[2],
    symmetric = TRUE
  ) %>%
  draw()


overlay_threshold <- c(0.3, 1.0)

slices <- 66
sliceSeries(nrow = 1, ncol = 1, dimension = human_slc_dim, slices = slices) %>%
  anatomy(
    human_anat_vol_cropped,
    low = human_anat_low,
    high = human_anat_high
  ) %>%
  overlay(
    pond_img,
    low = overlay_threshold[1],
    high = overlay_threshold[2],
    symmetric = TRUE
  ) %>%
  sliceSeries(
    nrow = 1,
    ncol = 1,
    dimension = human_slc_dim,
    slices = slices
  ) %>%
  anatomy() %>%
  overlay(
    list_hbn_imgs[[1]],
    low = overlay_threshold[1],
    high = overlay_threshold[2],
    symmetric = TRUE
  ) %>%
  sliceSeries(
    nrow = 1,
    ncol = 1,
    dimension = human_slc_dim,
    slices = slices
  ) %>%
  anatomy() %>%
  overlay(
    list_hbn_imgs[[2]],
    low = overlay_threshold[1],
    high = overlay_threshold[2],
    symmetric = TRUE
  ) %>%
  sliceSeries(
    nrow = 1,
    ncol = 1,
    dimension = human_slc_dim,
    slices = slices
  ) %>%
  anatomy() %>%
  overlay(
    list_hbn_imgs[[3]],
    low = overlay_threshold[1],
    high = overlay_threshold[2],
    symmetric = TRUE
  ) %>%
  sliceSeries(
    nrow = 1,
    ncol = 1,
    dimension = human_slc_dim,
    slices = slices
  ) %>%
  anatomy() %>%
  overlay(
    list_hbn_imgs[[4]],
    low = overlay_threshold[1],
    high = overlay_threshold[2],
    symmetric = TRUE
  ) %>%
  draw()
