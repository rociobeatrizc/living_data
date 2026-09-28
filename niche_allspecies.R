# ============================================================================
# STATIC VS CURRENTLY SUPPORTED CLIMATIC NICHES: ALL SWEDISH SPECIES
#
# Required objects:
#   occ_table
#   final_cell_table_all_species
#   grid25
#
# Entry population:
#   species already retained for temporal analysis (>20 occurrence records),
#   as represented in final_cell_table_all_species.
#
# Additional niche-analysis eligibility:
#   at least 5 static cells with climate data and at least 5 currently
#   supported cells with positive weights.
# ============================================================================

library(dplyr)
library(tidyr)
library(purrr)
library(tibble)
library(terra)
library(geodata)
library(usdm)
library(MASS)
library(sf)
library(readr)
library(ggplot2)
library(patchwork)

# ----------------------------------------------------------------------------
# 1. Settings
# ----------------------------------------------------------------------------

R2_threshold <- 0.1304739
low_probability_threshold <- 0.3333333
high_probability_threshold <- 0.6666667

niche_mass <- 0.90
minimum_kde_cells <- 5
kde_grid_size <- 250
maximum_replication <- 100

output_directory <- "all_species_climatic_niche_outputs"
dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)

# ----------------------------------------------------------------------------
# 2. Reclassify species x cell histories with the updated thresholds
# ----------------------------------------------------------------------------

required_columns <- c(
  "valid_name", "cell_id", "R2_McFadden", "DP", "mean_probability"
)

missing_columns <- setdiff(
  required_columns,
  names(final_cell_table_all_species)
)

if (length(missing_columns) > 0) {
  stop(
    paste(
      "Missing columns in final_cell_table_all_species:",
      paste(missing_columns, collapse = ", ")
    )
  )
}

temporal_all_species <- final_cell_table_all_species %>%
  dplyr::select(dplyr::all_of(required_columns)) %>%
  dplyr::distinct(valid_name, cell_id, .keep_all = TRUE) %>%
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

# Only species that passed the temporal-analysis occurrence filter enter here.
species_to_analyse <- temporal_all_species %>%
  dplyr::distinct(valid_name) %>%
  dplyr::arrange(valid_name) %>%
  dplyr::pull(valid_name)

# ----------------------------------------------------------------------------
# 3. Download climate data, select variables, and aggregate to the grid
# ----------------------------------------------------------------------------

climate_data_path <- file.path(getwd(), "geodata")
dir.create(climate_data_path, recursive = TRUE, showWarnings = FALSE)

bioclim <- geodata::worldclim_country(
  country = "SWE",
  var = "bio",
  res = 0.5,
  version = "2.1",
  path = climate_data_path
)

set.seed(103)

climate_sample <- terra::spatSample(
  bioclim,
  size = 10000,
  method = "random",
  na.rm = TRUE,
  as.df = TRUE
) %>%
  tidyr::drop_na()

vif_result <- usdm::vifstep(climate_sample, th = 5)
selected_variables <- vif_result@results$Variables

grid_for_extraction <- grid25 %>%
  dplyr::mutate(extraction_id = dplyr::row_number())

grid_lookup <- grid_for_extraction %>%
  sf::st_drop_geometry() %>%
  dplyr::select(extraction_id, cell_id)

grid_vector <- terra::vect(grid_for_extraction)

climate_projected <- bioclim[[selected_variables]] %>%
  terra::project(terra::crs(grid_vector)) %>%
  terra::crop(grid_vector)

climate_by_cell <- terra::extract(
  climate_projected,
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

# ----------------------------------------------------------------------------
# 4. Define one common Swedish climatic PCA space
# ----------------------------------------------------------------------------

short_variable_names <- paste0(
  "bio",
  sprintf(
    "%02d",
    as.integer(sub(".*bio_", "", selected_variables))
  )
)

if (anyNA(short_variable_names) ||
    length(unique(short_variable_names)) != length(selected_variables)) {
  short_variable_names <- selected_variables
}

pca_variables <- climate_by_cell %>%
  dplyr::select(dplyr::all_of(selected_variables))

names(pca_variables) <- short_variable_names

climate_pca <- stats::prcomp(
  pca_variables,
  center = TRUE,
  scale. = TRUE
)

variance_explained <- summary(climate_pca)$importance[2, 1:2] * 100

pca_score_values <- as.data.frame(
  climate_pca$x[, 1:2, drop = FALSE]
)
names(pca_score_values) <- c("PC1", "PC2")

pca_scores <- dplyr::bind_cols(
  climate_by_cell %>% dplyr::select(cell_id),
  pca_score_values
)

pad_range <- function(x, padding = 0.50) {
  observed_range <- range(x, na.rm = TRUE)
  range_width <- diff(observed_range)
  if (!is.finite(range_width) || range_width == 0) range_width <- 1
  
  c(
    observed_range[1] - padding * range_width,
    observed_range[2] + padding * range_width
  )
}

common_limits <- c(
  pad_range(pca_scores$PC1),
  pad_range(pca_scores$PC2)
)

# ----------------------------------------------------------------------------
# 5. KDE and polygon helpers
# ----------------------------------------------------------------------------

# MASS::kde2d has no weights argument. Weights are approximated by
# proportional row replication after scaling by the largest niche weight.
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
          maximum_replication * .data[[weight_column]] / maximum_weight
        ))
      )
    )
  
  expanded_data <- kde_data[
    rep(seq_len(nrow(kde_data)), kde_data$replication),
    c("PC1", "PC2")
  ]
  
  if (length(unique(expanded_data$PC1)) < 2 ||
      length(unique(expanded_data$PC2)) < 2) return(NULL)
  
  tryCatch(
    MASS::kde2d(
      x = expanded_data$PC1,
      y = expanded_data$PC2,
      n = n_grid,
      lims = limits
    ),
    error = function(e) NULL
  )
}

extract_mass_polygon <- function(kde, probability = 0.90) {
  dx <- median(diff(kde$x))
  dy <- median(diff(kde$y))
  density <- as.vector(kde$z)
  density[!is.finite(density)] <- 0
  mass <- density * dx * dy
  
  density_order <- order(density, decreasing = TRUE)
  cumulative_mass <- cumsum(mass[density_order]) / sum(mass)
  threshold <- density[density_order[which(cumulative_mass >= probability)[1]]]
  
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
    tryCatch(sf::st_polygon(list(coordinates)), error = function(e) NULL)
  })
  
  polygons <- Filter(Negate(is.null), polygons)
  if (length(polygons) == 0) return(NULL)
  
  sf::st_sf(
    geometry = sf::st_union(
      sf::st_make_valid(sf::st_sfc(polygons))
    )
  )
}

normalize_kde <- function(kde) {
  density <- kde$z
  density[!is.finite(density) | density < 0] <- 0
  total <- sum(density)
  if (!is.finite(total) || total <= 0) return(NULL)
  density / total
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

# ----------------------------------------------------------------------------
# 6. Result template for excluded species
# ----------------------------------------------------------------------------

excluded_result <- function(
    species_name,
    reason,
    n_static_cells = NA_integer_,
    n_current_cells = NA_integer_) {
  
  tibble::tibble(
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

# ----------------------------------------------------------------------------
# 7. Run the niche comparison for one species
# ----------------------------------------------------------------------------

run_niche_analysis_one_species <- function(species_name) {
  ever_cells <- occ_table %>%
    dplyr::filter(
      species == species_name,
      !is.na(cell_id),
      year %in% 1996:2019
    ) %>%
    dplyr::distinct(cell_id)
  
  niche_data <- ever_cells %>%
    dplyr::inner_join(pca_scores, by = "cell_id") %>%
    dplyr::left_join(
      temporal_all_species %>%
        dplyr::filter(valid_name == species_name) %>%
        dplyr::select(
          cell_id, temporal_class, current_weight
        ),
      by = "cell_id"
    ) %>%
    dplyr::mutate(static_weight = 1)
  
  n_static_cells <- nrow(niche_data)
  n_current_cells <- sum(niche_data$current_weight > 0, na.rm = TRUE)
  
  if (n_static_cells < minimum_kde_cells) {
    return(excluded_result(
      species_name,
      "Fewer than five static cells with climate data",
      n_static_cells,
      n_current_cells
    ))
  }
  
  if (n_current_cells < minimum_kde_cells) {
    return(excluded_result(
      species_name,
      "Fewer than five currently supported cells",
      n_static_cells,
      n_current_cells
    ))
  }
  
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

# ----------------------------------------------------------------------------
# 8. Run all retained temporal-analysis species
# ----------------------------------------------------------------------------

all_species_niche_metrics <- purrr::map_dfr(
  species_to_analyse,
  function(species_name) {
    message("Processing: ", species_name)
    
    tryCatch(
      run_niche_analysis_one_species(species_name),
      error = function(e) excluded_result(
        species_name,
        paste0("Unexpected error: ", conditionMessage(e))
      )
    )
  }
)

retained_niche_metrics <- all_species_niche_metrics %>%
  dplyr::filter(
    retained,
    is.finite(D),
    is.finite(ND),
    is.finite(PI),
    is.finite(OE)
  )

# ----------------------------------------------------------------------------
# 9. Eligibility and exclusion summaries
# ----------------------------------------------------------------------------

eligibility_summary <- tibble::tibble(
  n_species_entering_niche_analysis = length(species_to_analyse),
  n_species_retained = nrow(retained_niche_metrics),
  n_species_excluded = length(species_to_analyse) -
    nrow(retained_niche_metrics),
  percentage_retained = 100 * nrow(retained_niche_metrics) /
    length(species_to_analyse)
)

exclusion_summary <- all_species_niche_metrics %>%
  dplyr::filter(!retained) %>%
  dplyr::count(exclusion_reason, name = "n_species", sort = TRUE) %>%
  dplyr::mutate(
    percentage_of_entering_species =
      100 * n_species / length(species_to_analyse)
  )

# ----------------------------------------------------------------------------
# 10. Metric summaries
# ----------------------------------------------------------------------------

metric_summary <- retained_niche_metrics %>%
  dplyr::select(species, D, ND, PI, OE) %>%
  tidyr::pivot_longer(
    cols = c(D, ND, PI, OE),
    names_to = "metric",
    values_to = "value"
  ) %>%
  dplyr::group_by(metric) %>%
  dplyr::summarise(
    n_species = sum(is.finite(value)),
    mean = mean(value, na.rm = TRUE),
    median = median(value, na.rm = TRUE),
    q25 = quantile(value, 0.25, na.rm = TRUE),
    q75 = quantile(value, 0.75, na.rm = TRUE),
    minimum = min(value, na.rm = TRUE),
    maximum = max(value, na.rm = TRUE),
    .groups = "drop"
  )

weighted_summary <- retained_niche_metrics %>%
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
    )
  )

# ----------------------------------------------------------------------------
# 11. Density plots with metric medians
# ----------------------------------------------------------------------------

make_metric_density <- function(data, metric, x_label) {
  plot_data <- data %>%
    dplyr::transmute(value = .data[[metric]]) %>%
    dplyr::filter(is.finite(value))
  
  metric_median <- median(plot_data$value)
  
  ggplot(plot_data, aes(x = value)) +
    geom_density(
      fill = "#BFD7EA",
      colour = "grey20",
      alpha = 0.82,
      linewidth = 0.8,
      trim = FALSE,
      bounds = c(0, 1)
    ) +
    geom_vline(
      xintercept = metric_median,
      colour = "#D55E00",
      linetype = "dashed",
      linewidth = 0.9
    ) +
    scale_x_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, by = 0.2),
      labels = c("0", "0.2", "0.4", "0.6", "0.8", "1"),
      expand = c(0, 0)
    ) +
    scale_y_continuous(
      expand = expansion(mult = c(0, 0.08))
    ) +
    labs(x = x_label, y = "Density") +
    theme_classic(base_size = 14)
}

plot_D <- make_metric_density(
  retained_niche_metrics,
  "D",
  "Schoener's D"
)

plot_ND <- make_metric_density(
  retained_niche_metrics,
  "ND",
  "Niche Divergence (ND)"
)

plot_PI <- make_metric_density(
  retained_niche_metrics,
  "PI",
  "Proportion Inside (PI)"
)

plot_OE <- make_metric_density(
  retained_niche_metrics,
  "OE",
  "Overestimation (OE)"
)

density_metrics_plot <-
  (plot_D | plot_ND) /
  (plot_PI | plot_OE) +
  patchwork::plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 17, face = "bold"))

print(density_metrics_plot)
