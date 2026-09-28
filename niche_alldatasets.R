# ============================================================================
# STATIC VS CURRENTLY SUPPORTED CLIMATIC NICHES: FOUR COMPARISON DATASETS
#
# Run this script after multi_dataset_temporal_revised.R.
#
# Required objects:
#   dataset_settings
#   grid_list
#   all_species_cell_models
#   all_coverage
#   all_calibration_thresholds
#
# Niche definitions:
#   Static niche:
#     every cell occupied at least once; weight = 1.
#
#   Currently supported niche:
#     Increasing cells; weight = abs(DP) * McFadden R2.
#     Flat high cells; weight = 1.
#
#   Decreasing, Flat low, Noisy, and unclassified cells do not contribute to
#   the currently supported niche.
#
# Climatic space and temporal thresholds are dataset-specific.
# ============================================================================

library(dplyr)
library(tidyr)
library(purrr)
library(tibble)
library(readr)
library(sf)
library(terra)
library(geodata)
library(usdm)
library(MASS)
library(ggplot2)
library(patchwork)


# ============================================================================
# 1. SETTINGS
# ============================================================================

minimum_kde_cells <- 5
niche_mass <- 0.90
kde_grid_size <- 250
maximum_replication <- 100

# Retain the same WorldClim resolution used in the Swedish analysis.
worldclim_resolution_minutes <- 0.5

output_directory <- "all_datasets_climatic_niche_outputs"
dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)

climate_data_path <- file.path(getwd(), "geodata")
dir.create(climate_data_path, recursive = TRUE, showWarnings = FALSE)

dataset_order <- c(
  "Birds", "Fish", "Invertebrates", "Phytoplankton"
)

dataset_labels <- c(
  "Birds" = "South African birds",
  "Fish" = "Finnish fishes",
  "Invertebrates" = "UK invertebrates",
  "Phytoplankton" = "Black Sea phytoplankton"
)

dataset_colours <- c(
  "Birds" = "#440154FF",
  "Fish" = "#31688EFF",
  "Invertebrates" = "#35B779FF",
  "Phytoplankton" = "#FDE725FF"
)


# ============================================================================
# 2. CHECK REQUIRED OBJECTS
# ============================================================================

required_objects <- c(
  "dataset_settings",
  "grid_list",
  "all_species_cell_models",
  "all_coverage",
  "all_calibration_thresholds"
)

missing_objects <- required_objects[
  !vapply(required_objects, exists, logical(1), inherits = TRUE)
]

if (length(missing_objects) > 0) {
  stop(
    paste(
      "Run multi_dataset_temporal_revised.R first. Missing objects:",
      paste(missing_objects, collapse = ", ")
    )
  )
}


# ============================================================================
# 3. RECLASSIFY USING EACH DATASET'S UPDATED MONTE CARLO THRESHOLDS
# ============================================================================

temporal_for_niches <- all_species_cell_models %>%
  dplyr::select(
    dataset,
    valid_name,
    cell_id,
    R2_McFadden,
    DP,
    mean_probability
  ) %>%
  dplyr::distinct(dataset, valid_name, cell_id, .keep_all = TRUE) %>%
  dplyr::left_join(
    all_calibration_thresholds %>%
      dplyr::select(
        dataset,
        R2_threshold,
        low_probability_threshold,
        high_probability_threshold
      ),
    by = "dataset"
  ) %>%
  dplyr::mutate(
    temporal_class = dplyr::case_when(
      !is.finite(R2_McFadden) |
        !is.finite(DP) |
        !is.finite(mean_probability) ~ NA_character_,
      R2_McFadden > R2_threshold & DP > 0 ~ "Increasing",
      R2_McFadden > R2_threshold & DP < 0 ~ "Decreasing",
      R2_McFadden <= R2_threshold &
        mean_probability < low_probability_threshold ~ "Flat low",
      R2_McFadden <= R2_threshold &
        mean_probability > high_probability_threshold ~ "Flat high",
      R2_McFadden <= R2_threshold ~ "Noisy",
      TRUE ~ NA_character_
    ),
    current_weight = dplyr::case_when(
      temporal_class == "Increasing" ~ abs(DP) * R2_McFadden,
      temporal_class == "Flat high" ~ 1,
      TRUE ~ 0
    )
  )


# ============================================================================
# 4. RECONSTRUCT EVER-OCCUPIED CELLS FOR EVERY SPECIES
# ============================================================================

occurrence_list <- list()

for (dataset_index in seq_len(nrow(dataset_settings))) {
  dataset_name <- dataset_settings$dataset[dataset_index]
  dataset_file <- dataset_settings$file[dataset_index]
  
  raw_data <- readr::read_csv(
    dataset_file,
    show_col_types = FALSE
  ) %>%
    dplyr::filter(
      !is.na(valid_name),
      !is.na(LONGITUDE),
      !is.na(LATITUDE),
      !is.na(YEAR)
    )
  
  species_entering <- all_coverage %>%
    dplyr::filter(dataset == dataset_name) %>%
    dplyr::distinct(valid_name) %>%
    dplyr::pull(valid_name)
  
  raw_data <- raw_data %>%
    dplyr::filter(valid_name %in% species_entering)
  
  dataset_grid <- grid_list[[dataset_name]]
  
  raw_sf <- sf::st_as_sf(
    raw_data,
    coords = c("LONGITUDE", "LATITUDE"),
    crs = 4326,
    remove = FALSE
  ) %>%
    sf::st_transform(sf::st_crs(dataset_grid))
  
  occurrence_list[[dataset_name]] <- sf::st_join(
    raw_sf,
    dataset_grid %>% dplyr::select(cell_id),
    join = sf::st_intersects,
    left = FALSE
  ) %>%
    sf::st_drop_geometry() %>%
    dplyr::transmute(
      dataset = dataset_name,
      valid_name,
      cell_id = as.integer(cell_id)
    ) %>%
    dplyr::distinct()
}

all_occupied_cells <- dplyr::bind_rows(occurrence_list)


# ============================================================================
# 5. DOWNLOAD WORLDCLIM AND PREPARE A DATASET-SPECIFIC PCA SPACE
# ============================================================================

bioclim_global <- geodata::worldclim_global(
  var = "bio",
  res = worldclim_resolution_minutes,
  version = "2.1",
  path = climate_data_path
)

prepare_climatic_space <- function(dataset_name) {
  message("Preparing climatic space: ", dataset_name)
  
  # st_intersection() can leave clipped grid cells as GEOMETRYCOLLECTIONS.
  # Convert them to valid polygonal geometries and dissolve any split pieces
  # belonging to the same cell before conversion to a terra SpatVector.
  dataset_grid <- grid_list[[dataset_name]] %>%
    dplyr::select(cell_id) %>%
    sf::st_make_valid()
  
  dataset_grid <- suppressWarnings(
    sf::st_collection_extract(dataset_grid, "POLYGON")
  )
  
  # Do not assume that the active sf geometry column is named "geometry".
  dataset_grid <- dataset_grid[
    !sf::st_is_empty(dataset_grid),
    ,
    drop = FALSE
  ]
  
  dataset_grid <- dataset_grid %>%
    dplyr::group_by(cell_id) %>%
    dplyr::summarise(.groups = "drop", do_union = TRUE) %>%
    sf::st_cast("MULTIPOLYGON", warn = FALSE) %>%
    dplyr::mutate(
      cell_id = as.integer(cell_id),
      extraction_id = dplyr::row_number()
    )
  
  grid_lookup <- dataset_grid %>%
    sf::st_drop_geometry() %>%
    dplyr::select(extraction_id, cell_id)
  
  grid_vector <- terra::vect(dataset_grid)
  
  # Crop in the original geographic CRS before projection to reduce memory.
  grid_geographic <- sf::st_transform(dataset_grid, 4326)
  geographic_extent <- terra::ext(terra::vect(grid_geographic))
  
  climate_crop <- terra::crop(
    bioclim_global,
    geographic_extent,
    snap = "out"
  )
  
  climate_projected <- terra::project(
    climate_crop,
    terra::crs(grid_vector)
  )
  
  set.seed(103)
  
  climate_sample <- terra::spatSample(
    climate_projected,
    size = min(10000, terra::ncell(climate_projected)),
    method = "random",
    na.rm = TRUE,
    as.df = TRUE
  ) %>%
    tidyr::drop_na()
  
  vif_result <- usdm::vifstep(climate_sample, th = 5)
  selected_variables <- vif_result@results$Variables
  
  climate_by_cell <- terra::extract(
    climate_projected[[selected_variables]],
    grid_vector,
    fun = mean,
    na.rm = TRUE,
    ID = TRUE
  ) %>%
    as.data.frame() %>%
    dplyr::rename(extraction_id = ID) %>%
    dplyr::left_join(grid_lookup, by = "extraction_id") %>%
    dplyr::select(cell_id, dplyr::all_of(selected_variables)) %>%
    tidyr::drop_na(dplyr::all_of(selected_variables))
  
  pca_data <- climate_by_cell %>%
    dplyr::select(dplyr::all_of(selected_variables))
  
  climate_pca <- stats::prcomp(
    pca_data,
    center = TRUE,
    scale. = TRUE
  )
  
  pca_scores <- as.data.frame(
    climate_pca$x[, 1:2, drop = FALSE]
  )
  names(pca_scores) <- c("PC1", "PC2")
  
  pca_scores <- dplyr::bind_cols(
    climate_by_cell %>% dplyr::select(cell_id),
    pca_scores
  )
  
  variance_explained <-
    summary(climate_pca)$importance[2, 1:2] * 100
  
  list(
    dataset = dataset_name,
    selected_variables = selected_variables,
    pca = climate_pca,
    pca_scores = pca_scores,
    variance_PC1 = unname(variance_explained[1]),
    variance_PC2 = unname(variance_explained[2])
  )
}

climatic_space_list <- purrr::map(
  dataset_order,
  prepare_climatic_space
)
names(climatic_space_list) <- dataset_order


# ============================================================================
# 6. KDE AND POLYGON HELPERS
# ============================================================================

pad_range <- function(x, padding = 0.50) {
  observed_range <- range(x, na.rm = TRUE)
  width <- diff(observed_range)
  if (!is.finite(width) || width == 0) width <- 1
  
  c(
    observed_range[1] - padding * width,
    observed_range[2] + padding * width
  )
}

make_weighted_kde <- function(
    data,
    weight_column,
    limits,
    maximum_replication = 100,
    n_grid = 250,
    minimum_cells = 5) {
  
  kde_data <- data %>%
    dplyr::filter(
      is.finite(PC1),
      is.finite(PC2),
      is.finite(.data[[weight_column]]),
      .data[[weight_column]] > 0
    )
  
  if (nrow(kde_data) < minimum_cells) return(NULL)
  
  maximum_weight <- max(kde_data[[weight_column]], na.rm = TRUE)
  if (!is.finite(maximum_weight) || maximum_weight <= 0) return(NULL)
  
  kde_data <- kde_data %>%
    dplyr::mutate(
      replication = pmax(
        1L,
        as.integer(round(
          maximum_replication *
            .data[[weight_column]] / maximum_weight
        ))
      )
    )
  
  expanded <- kde_data[
    rep(seq_len(nrow(kde_data)), kde_data$replication),
    c("PC1", "PC2")
  ]
  
  if (
    length(unique(expanded$PC1)) < 2 ||
    length(unique(expanded$PC2)) < 2
  ) return(NULL)
  
  tryCatch(
    MASS::kde2d(
      x = expanded$PC1,
      y = expanded$PC2,
      n = n_grid,
      lims = limits
    ),
    error = function(e) NULL
  )
}

normalize_kde <- function(kde) {
  density <- kde$z
  density[!is.finite(density) | density < 0] <- 0
  total <- sum(density)
  if (!is.finite(total) || total <= 0) return(NULL)
  density / total
}

extract_mass_polygon <- function(kde, probability = 0.90) {
  dx <- median(diff(kde$x))
  dy <- median(diff(kde$y))
  density <- as.vector(kde$z)
  density[!is.finite(density)] <- 0
  mass <- density * dx * dy
  
  density_order <- order(density, decreasing = TRUE)
  cumulative_mass <- cumsum(mass[density_order]) / sum(mass)
  threshold <- density[
    density_order[which(cumulative_mass >= probability)[1]]
  ]
  
  contours <- grDevices::contourLines(
    x = kde$x,
    y = kde$y,
    z = kde$z,
    levels = threshold
  )
  
  polygons <- lapply(contours, function(contour) {
    coordinates <- cbind(contour$x, contour$y)
    if (nrow(coordinates) < 3) return(NULL)
    
    if (!all(coordinates[1, ] == coordinates[nrow(coordinates), ])) {
      coordinates <- rbind(coordinates, coordinates[1, ])
    }
    
    tryCatch(
      sf::st_polygon(list(coordinates)),
      error = function(e) NULL
    )
  })
  
  polygons <- Filter(Negate(is.null), polygons)
  if (length(polygons) == 0) return(NULL)
  
  sf::st_sf(
    geometry = sf::st_union(
      sf::st_make_valid(sf::st_sfc(polygons))
    )
  )
}

containment <- function(A, B) {
  if (is.null(A) || is.null(B)) return(NA_real_)
  
  A <- sf::st_make_valid(A)
  B <- sf::st_make_valid(B)
  area_A <- sum(as.numeric(sf::st_area(A)))
  
  if (!is.finite(area_A) || area_A <= 0) return(NA_real_)
  
  intersection <- tryCatch(
    suppressWarnings(sf::st_intersection(A, B)),
    error = function(e) NULL
  )
  
  if (is.null(intersection) || nrow(intersection) == 0) return(0)
  
  sum(as.numeric(sf::st_area(intersection))) / area_A
}


# ============================================================================
# 7. RUN ONE SPECIES
# ============================================================================

excluded_result <- function(
    dataset_name,
    species_name,
    reason,
    n_static_cells = NA_integer_,
    n_current_cells = NA_integer_) {
  
  tibble::tibble(
    dataset = dataset_name,
    species = species_name,
    retained = FALSE,
    exclusion_reason = reason,
    n_static_cells = n_static_cells,
    n_current_cells = n_current_cells,
    n_increasing = NA_integer_,
    n_flat_high = NA_integer_,
    D = NA_real_,
    ND = NA_real_,
    PI = NA_real_,
    OE = NA_real_
  )
}

run_one_species_niche <- function(dataset_name, species_name) {
  pca_scores <- climatic_space_list[[dataset_name]]$pca_scores
  
  niche_data <- all_occupied_cells %>%
    dplyr::filter(
      dataset == dataset_name,
      valid_name == species_name
    ) %>%
    dplyr::select(cell_id) %>%
    dplyr::inner_join(pca_scores, by = "cell_id") %>%
    dplyr::left_join(
      temporal_for_niches %>%
        dplyr::filter(
          dataset == dataset_name,
          valid_name == species_name
        ) %>%
        dplyr::select(
          cell_id,
          temporal_class,
          current_weight
        ),
      by = "cell_id"
    ) %>%
    dplyr::mutate(static_weight = 1)
  
  n_static_cells <- nrow(niche_data)
  n_current_cells <- sum(
    niche_data$current_weight > 0,
    na.rm = TRUE
  )
  
  if (n_static_cells < minimum_kde_cells) {
    return(excluded_result(
      dataset_name,
      species_name,
      "Fewer than five static cells with climate data",
      n_static_cells,
      n_current_cells
    ))
  }
  
  if (n_current_cells < minimum_kde_cells) {
    return(excluded_result(
      dataset_name,
      species_name,
      "Fewer than five currently supported cells",
      n_static_cells,
      n_current_cells
    ))
  }
  
  common_limits <- c(
    pad_range(pca_scores$PC1),
    pad_range(pca_scores$PC2)
  )
  
  KDE_static <- make_weighted_kde(
    niche_data,
    "static_weight",
    common_limits,
    maximum_replication = 1,
    n_grid = kde_grid_size,
    minimum_cells = minimum_kde_cells
  )
  
  KDE_current <- make_weighted_kde(
    niche_data,
    "current_weight",
    common_limits,
    maximum_replication = maximum_replication,
    n_grid = kde_grid_size,
    minimum_cells = minimum_kde_cells
  )
  
  if (is.null(KDE_static) || is.null(KDE_current)) {
    return(excluded_result(
      dataset_name,
      species_name,
      "At least one KDE surface could not be estimated",
      n_static_cells,
      n_current_cells
    ))
  }
  
  polygon_static <- extract_mass_polygon(KDE_static, niche_mass)
  polygon_current <- extract_mass_polygon(KDE_current, niche_mass)
  
  if (is.null(polygon_static) || is.null(polygon_current)) {
    return(excluded_result(
      dataset_name,
      species_name,
      "At least one 90% niche envelope could not be extracted",
      n_static_cells,
      n_current_cells
    ))
  }
  
  static_density <- normalize_kde(KDE_static)
  current_density <- normalize_kde(KDE_current)
  
  D <- 1 - 0.5 * sum(abs(static_density - current_density))
  ND <- 1 - D
  PI <- containment(polygon_current, polygon_static)
  OE <- 1 - containment(polygon_static, polygon_current)
  
  tibble::tibble(
    dataset = dataset_name,
    species = species_name,
    retained = TRUE,
    exclusion_reason = NA_character_,
    n_static_cells = n_static_cells,
    n_current_cells = n_current_cells,
    n_increasing = sum(
      niche_data$temporal_class == "Increasing",
      na.rm = TRUE
    ),
    n_flat_high = sum(
      niche_data$temporal_class == "Flat high",
      na.rm = TRUE
    ),
    D = D,
    ND = ND,
    PI = PI,
    OE = OE
  )
}


# ============================================================================
# 8. RUN ALL SPECIES IN ALL DATASETS
# ============================================================================

species_to_analyse <- all_coverage %>%
  dplyr::distinct(dataset, valid_name) %>%
  dplyr::arrange(
    factor(dataset, levels = dataset_order),
    valid_name
  )

all_dataset_niche_metrics <- purrr::pmap_dfr(
  species_to_analyse,
  function(dataset, valid_name) {
    message("Processing ", dataset, ": ", valid_name)
    
    tryCatch(
      run_one_species_niche(dataset, valid_name),
      error = function(e) excluded_result(
        dataset,
        valid_name,
        paste0("Unexpected error: ", conditionMessage(e))
      )
    )
  }
)

retained_dataset_niche_metrics <- all_dataset_niche_metrics %>%
  dplyr::filter(
    retained,
    is.finite(D),
    is.finite(ND),
    is.finite(PI),
    is.finite(OE)
  )


# ============================================================================
# 9. ELIGIBILITY AND METRIC SUMMARIES
# ============================================================================

niche_eligibility_summary <- all_dataset_niche_metrics %>%
  dplyr::group_by(dataset) %>%
  dplyr::summarise(
    n_species_entering = dplyr::n(),
    n_species_retained = sum(retained),
    n_species_excluded = sum(!retained),
    percentage_retained = 100 * n_species_retained / n_species_entering,
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    dataset = factor(dataset, levels = dataset_order)
  ) %>%
  dplyr::arrange(dataset)

niche_exclusion_summary <- all_dataset_niche_metrics %>%
  dplyr::filter(!retained) %>%
  dplyr::count(dataset, exclusion_reason, name = "n_species") %>%
  dplyr::left_join(
    niche_eligibility_summary %>%
      dplyr::select(dataset, n_species_entering),
    by = "dataset"
  ) %>%
  dplyr::mutate(
    percentage_of_entering_species =
      100 * n_species / n_species_entering,
    dataset = factor(dataset, levels = dataset_order)
  ) %>%
  dplyr::arrange(dataset, dplyr::desc(n_species))

niche_metric_summary <- retained_dataset_niche_metrics %>%
  dplyr::select(dataset, species, D, ND, PI, OE) %>%
  tidyr::pivot_longer(
    cols = c(D, ND, PI, OE),
    names_to = "metric",
    values_to = "value"
  ) %>%
  dplyr::group_by(dataset, metric) %>%
  dplyr::summarise(
    n_species = sum(is.finite(value)),
    mean = mean(value, na.rm = TRUE),
    median = median(value, na.rm = TRUE),
    q25 = quantile(value, 0.25, na.rm = TRUE),
    q75 = quantile(value, 0.75, na.rm = TRUE),
    minimum = min(value, na.rm = TRUE),
    maximum = max(value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    dataset = factor(dataset, levels = dataset_order),
    metric = factor(metric, levels = c("D", "ND", "PI", "OE"))
  ) %>%
  dplyr::arrange(dataset, metric)

niche_summary_table <- niche_metric_summary %>%
  dplyr::select(dataset, metric, median) %>%
  tidyr::pivot_wider(
    names_from = metric,
    values_from = median
  ) %>%
  dplyr::left_join(
    niche_eligibility_summary,
    by = "dataset"
  ) %>%
  dplyr::select(
    dataset,
    n_species_entering,
    n_species_retained,
    percentage_retained,
    D,
    ND,
    PI,
    OE
  )

niche_weighted_summary <- retained_dataset_niche_metrics %>%
  dplyr::group_by(dataset) %>%
  dplyr::summarise(
    PI_weighted_mean = weighted.mean(
      PI,
      w = n_current_cells,
      na.rm = TRUE
    ),
    OE_weighted_mean = weighted.mean(
      OE,
      w = n_static_cells,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  dplyr::mutate(dataset = factor(dataset, levels = dataset_order)) %>%
  dplyr::arrange(dataset)


# ============================================================================
# 10. DENSITY PLOTS WITH DATASET-SPECIFIC MEDIANS
# ============================================================================

density_plot_data <- retained_dataset_niche_metrics %>%
  dplyr::select(dataset, species, D, ND, PI, OE) %>%
  tidyr::pivot_longer(
    cols = c(D, ND, PI, OE),
    names_to = "metric",
    values_to = "value"
  ) %>%
  dplyr::filter(is.finite(value)) %>%
  dplyr::mutate(
    dataset = factor(dataset, levels = dataset_order)
  )

density_medians <- density_plot_data %>%
  dplyr::group_by(dataset, metric) %>%
  dplyr::summarise(
    median = median(value),
    .groups = "drop"
  )

make_dataset_metric_density <- function(metric_name, x_label) {
  plot_data <- density_plot_data %>%
    dplyr::filter(metric == metric_name)
  
  median_data <- density_medians %>%
    dplyr::filter(metric == metric_name)
  
  ggplot(
    plot_data,
    aes(x = value, colour = dataset, fill = dataset)
  ) +
    geom_density(
      alpha = 0.16,
      linewidth = 0.9,
      trim = FALSE,
      bounds = c(0, 1)
    ) +
    geom_vline(
      data = median_data,
      aes(xintercept = median, colour = dataset),
      linetype = "dashed",
      linewidth = 0.65,
      show.legend = FALSE
    ) +
    scale_colour_manual(
      values = dataset_colours,
      labels = dataset_labels,
      name = "Dataset"
    ) +
    scale_fill_manual(
      values = dataset_colours,
      labels = dataset_labels,
      name = "Dataset"
    ) +
    scale_x_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, 0.2),
      labels = c("0", "0.2", "0.4", "0.6", "0.8", "1"),
      expand = c(0, 0)
    ) +
    scale_y_continuous(
      expand = expansion(mult = c(0, 0.06))
    ) +
    labs(x = x_label, y = "Density") +
    theme_classic(base_size = 14) +
    theme(
      legend.position = "bottom",
      legend.title = element_text(face = "bold"),
      legend.text = element_text(size = 11)
    )
}

plot_D <- make_dataset_metric_density(
  "D",
  "Schoener's D"
)

plot_ND <- make_dataset_metric_density(
  "ND",
  "Niche Divergence (ND)"
)

plot_PI <- make_dataset_metric_density(
  "PI",
  "Proportion Inside (PI)"
)

plot_OE <- make_dataset_metric_density(
  "OE",
  "Overestimation (OE)"
)

all_dataset_density_plot <-
  (plot_D + theme(legend.position = "none") |
     plot_ND + theme(legend.position = "none")) /
  (plot_PI + theme(legend.position = "none") |
     plot_OE + theme(legend.position = "none")) +
  patchwork::plot_layout(guides = "collect") +
  patchwork::plot_annotation(tag_levels = "A") &
  theme(
    legend.position = "bottom",
    plot.tag = element_text(size = 17, face = "bold")
  )

