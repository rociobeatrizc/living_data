# ============================================================================
# CLIMATIC NICHE COMPARISON FOR ONE SPECIES
#
# Static niche:
#   all ever-occupied cells, each with weight 1.
#
# Currently supported niche:
#   Increasing cells: abs(DP) * McFadden R2
#   Flat high cells:  weight 1
#
# Required objects:
#   species_name, occ_table, final_cell_table, grid25
# ============================================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(terra)
library(geodata)
library(usdm)
library(MASS)
library(sf)
library(vegan)
library(rnaturalearth)

# Updated Monte Carlo thresholds
species_name <- "Turdus iliacus"

R2_threshold <- 0.1304739
low_probability_threshold <- 0.3333333
high_probability_threshold <- 0.6666667

niche_mass <- 0.90
minimum_kde_cells <- 5

# ----------------------------------------------------------------------------
# 1. Temporal classes and niche weights
# ----------------------------------------------------------------------------

ever_cells <- occ_table %>%
  dplyr::filter(
    species == species_name,
    !is.na(cell_id),
    year %in% 1996:2019
  ) %>%
  dplyr::distinct(cell_id)

temporal_table <- final_cell_table %>%
  dplyr::mutate(
    DP = pred_2019 - pred_1996,
    
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
      temporal_class == "Increasing" ~
        abs(DP) * R2_McFadden,
      
      temporal_class == "Flat high" ~ 1,
      
      TRUE ~ 0
    )
  )
# ----------------------------------------------------------------------------
# 2. Climate variables and 25-km cell means
# ----------------------------------------------------------------------------

climate_data_path <- file.path(getwd(), "geodata")
dir.create(
  climate_data_path,
  recursive = TRUE,
  showWarnings = FALSE
)

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
  drop_na()

vif_result <- usdm::vifstep(climate_sample, th = 5)
selected_variables <- vif_result@results$Variables

grid_for_extraction <- grid25 %>%
  mutate(extraction_id = row_number())

grid_lookup <- grid_for_extraction %>%
  st_drop_geometry() %>%
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
  rename(extraction_id = ID) %>%
  left_join(grid_lookup, by = "extraction_id") %>%
  dplyr::select(cell_id, all_of(selected_variables)) %>%
  drop_na(all_of(selected_variables))

# ----------------------------------------------------------------------------
# 3. PCA: fit on all grid cells and project niche cells
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
  dplyr::select(all_of(selected_variables))

names(pca_variables) <- short_variable_names

pca <- stats::prcomp(
  pca_variables,
  center = TRUE,
  scale. = TRUE
)

variance_explained <- summary(pca)$importance[2, 1:2] * 100

pca_scores <- as.data.frame(pca$x[, 1:2, drop = FALSE]) %>%
  setNames(c("PC1", "PC2")) %>%
  bind_cols(climate_by_cell %>% dplyr::select(cell_id), .)

niche_data <- ever_cells %>%
  dplyr::inner_join(
    pca_scores,
    by = "cell_id"
  ) %>%
  dplyr::left_join(
    temporal_table %>%
      dplyr::select(
        cell_id,
        R2_McFadden,
        DP,
        mean_probability,
        temporal_class,
        current_weight
      ),
    by = "cell_id"
  ) %>%
  dplyr::mutate(static_weight = 1)
# ----------------------------------------------------------------------------
# 4. Common KDE extent and weighted KDE function
# ----------------------------------------------------------------------------

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
    filter(
      is.finite(PC1),
      is.finite(PC2),
      is.finite(.data[[weight_column]]),
      .data[[weight_column]] > 0
    )
  
  if (nrow(kde_data) < minimum_cells) return(NULL)
  
  maximum_weight <- max(kde_data[[weight_column]], na.rm = TRUE)
  
  kde_data <- kde_data %>%
    mutate(
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
  
  kde <- MASS::kde2d(
    x = expanded_data$PC1,
    y = expanded_data$PC2,
    n = n_grid,
    lims = limits
  )
  
  grid <- expand.grid(PC1 = kde$x, PC2 = kde$y) %>%
    mutate(
      density = as.vector(kde$z),
      display_density = sqrt(density)
    )
  
  list(kde = kde, grid = grid, data = kde_data)
}

KDE_static <- make_weighted_kde(
  niche_data,
  "static_weight",
  common_limits,
  maximum_replication = 1,
  minimum_cells = minimum_kde_cells
)

KDE_current <- make_weighted_kde(
  niche_data,
  "current_weight",
  common_limits,
  maximum_replication = 100,
  minimum_cells = minimum_kde_cells
)

if (is.null(KDE_static) || is.null(KDE_current)) {
  stop("At least one niche had too few cells for KDE estimation.")
}

# ----------------------------------------------------------------------------
# 5. Extract the 90% KDE envelopes
# ----------------------------------------------------------------------------

extract_mass_polygon <- function(kde, probability = 0.85) {
  dx <- median(diff(kde$x))
  dy <- median(diff(kde$y))
  density <- as.vector(kde$z)
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
    sf::st_polygon(list(coordinates))
  })
  
  polygons <- Filter(Negate(is.null), polygons)
  if (length(polygons) == 0) return(NULL)
  
  st_sf(
    geometry = st_union(st_make_valid(st_sfc(polygons)))
  )
}

polygon_static <- extract_mass_polygon(KDE_static$kde, niche_mass)
polygon_current <- extract_mass_polygon(KDE_current$kde, niche_mass)

# ----------------------------------------------------------------------------
# 6. Environmental vectors
# ----------------------------------------------------------------------------

set.seed(103)

environment_fit <- vegan::envfit(
  pca$x[, 1:2, drop = FALSE],
  pca_variables,
  permutations = 999,
  na.rm = TRUE
)

environment_arrows <- as.data.frame(
  vegan::scores(environment_fit, display = "vectors")
) %>%
  rename(PC1_vector = PC1, PC2_vector = PC2) %>%
  mutate(variable = rownames(.))

arrow_scale <- 0.35 * min(
  diff(common_limits[1:2]),
  diff(common_limits[3:4])
) / max(abs(as.matrix(
  environment_arrows[, c("PC1_vector", "PC2_vector")]
)))

environment_arrows <- environment_arrows %>%
  mutate(
    PC1_end = PC1_vector * arrow_scale,
    PC2_end = PC2_vector * arrow_scale,
    label_PC1 = PC1_end * 1.08,
    label_PC2 = PC2_end * 1.08
  )

# ============================================================================
# NICHE PLOT
# ============================================================================

niche_plot <- ggplot() +
  
  # --------------------------------------------------------------------------
# KDE density surfaces
# --------------------------------------------------------------------------

geom_raster(
  data = KDE_static$grid,
  aes(
    x = PC1,
    y = PC2,
    alpha = display_density,
    fill = "Static"
  ),
  interpolate = TRUE
) +
  
  geom_raster(
    data = KDE_current$grid,
    aes(
      x = PC1,
      y = PC2,
      alpha = display_density,
      fill = "Currently supported"
    ),
    interpolate = TRUE
  ) +
  
  # --------------------------------------------------------------------------
# 90% niche envelopes
# --------------------------------------------------------------------------

geom_sf(
  data = polygon_static,
  aes(colour = "Static"),
  fill = NA,
  linewidth = 1.2,
  inherit.aes = FALSE
) +
  
  geom_sf(
    data = polygon_current,
    aes(colour = "Currently supported"),
    fill = NA,
    linewidth = 1.2,
    inherit.aes = FALSE
  ) +
  
  # --------------------------------------------------------------------------
# Environmental vectors
# --------------------------------------------------------------------------

geom_segment(
  data = environment_arrows,
  aes(
    x = 0,
    y = 0,
    xend = PC1_end,
    yend = PC2_end
  ),
  arrow = grid::arrow(
    length = grid::unit(0.22, "cm")
  ),
  colour = "black",
  linewidth = 0.5,
  alpha = 0.65,
  inherit.aes = FALSE
) +
  
  geom_text(
    data = environment_arrows,
    aes(
      x = label_PC1,
      y = label_PC2,
      label = variable
    ),
    size = 5,
    colour = "black",
    alpha = 0.75,
    inherit.aes = FALSE
  ) +
  
  # --------------------------------------------------------------------------
# Niche colours
# --------------------------------------------------------------------------

scale_fill_manual(
  values = c(
    "Static" = "#BDBDBD",
    "Currently supported" = "#CC79A7"
  ),
  breaks = c(
    "Currently supported",
    "Static"
  ),
  name = "Niche representation"
) +
  
  scale_colour_manual(
    values = c(
      "Static" = "#4D4D4D",
      "Currently supported" = "#8E3B70"
    ),
    breaks = c(
      "Currently supported",
      "Static"
    ),
    name = "Niche representation"
  ) +
  
  scale_alpha(
    range = c(0, 0.45),
    guide = "none"
  ) +
  
  # --------------------------------------------------------------------------
# Axis scales
# --------------------------------------------------------------------------

scale_x_continuous(
  breaks = scales::pretty_breaks(n = 5)
) +
  
  scale_y_continuous(
    breaks = scales::pretty_breaks(n = 5)
  ) +
  
  # --------------------------------------------------------------------------
# Draw axes along the bottom and left boundaries
# --------------------------------------------------------------------------

annotate(
  "segment",
  x = common_limits[1],
  xend = common_limits[2],
  y = common_limits[3],
  yend = common_limits[3],
  colour = "black",
  linewidth = 0.8
) +
  
  annotate(
    "segment",
    x = common_limits[1],
    xend = common_limits[1],
    y = common_limits[3],
    yend = common_limits[4],
    colour = "black",
    linewidth = 0.8
  ) +
  
  # --------------------------------------------------------------------------
# Plot limits
# --------------------------------------------------------------------------

coord_sf(
  xlim = common_limits[1:2],
  ylim = common_limits[3:4],
  expand = FALSE,
  datum = NA
) +
  
  # --------------------------------------------------------------------------
# Axis labels
# --------------------------------------------------------------------------

labs(
  x = paste0(
    "PC1 (",
    round(variance_explained[1], 1),
    "%)"
  ),
  y = paste0(
    "PC2 (",
    round(variance_explained[2], 1),
    "%)"
  )
) +
  
  # --------------------------------------------------------------------------
# Theme
# --------------------------------------------------------------------------

theme_classic(base_size = 17) +
  
  theme(
    axis.line = element_blank(),
    
    axis.ticks = element_line(
      colour = "black",
      linewidth = 0.6
    ),
    
    axis.ticks.length = grid::unit(
      0.18,
      "cm"
    ),
    
    axis.title = element_text(
      size = 19,
      face = "bold",
      colour = "black"
    ),
    
    axis.text = element_text(
      size = 16,
      colour = "black"
    ),
    
    legend.position = "right",
    
    legend.title = element_text(
      size = 17,
      face = "bold"
    ),
    
    legend.text = element_text(
      size = 16
    ),
    
    legend.key.height = grid::unit(
      0.8,
      "cm"
    ),
    
    legend.key.width = grid::unit(
      0.9,
      "cm"
    ),
    
    panel.border = element_blank(),
    
    plot.margin = margin(
      15,
      20,
      15,
      15
    )
  )


# Display figure
print(niche_plot)


# Save enlarged figure
ggsave(
  filename = paste0(
    gsub(" ", "_", species_name),
    "_static_vs_current_niche.png"
  ),
  plot = niche_plot,
  width = 11,
  height = 8.5,
  units = "in",
  dpi = 400,
  bg = "white"
)



# ----------------------------------------------------------------------------
# 8. Niche-comparison metrics
# ----------------------------------------------------------------------------

normalize_kde <- function(kde) {
  density <- kde$z
  density[!is.finite(density) | density < 0] <- 0
  density / sum(density)
}

static_density <- normalize_kde(KDE_static$kde)
current_density <- normalize_kde(KDE_current$kde)

schoener_D <- 1 - 0.5 * sum(abs(static_density - current_density))
niche_divergence <- 1 - schoener_D

containment <- function(A, B) {
  if (is.null(A) || is.null(B)) return(NA_real_)
  
  A <- st_make_valid(A)
  B <- st_make_valid(B)
  area_A <- sum(as.numeric(st_area(A)))
  intersection <- suppressWarnings(st_intersection(A, B))
  
  if (!is.finite(area_A) || area_A <= 0) return(NA_real_)
  if (is.null(intersection) || nrow(intersection) == 0) return(0)
  
  sum(as.numeric(st_area(intersection))) / area_A
}

proportion_inside <- containment(polygon_current, polygon_static)
overestimation <- 1 - containment(polygon_static, polygon_current)

results <- tibble(
  species = species_name,
  n_static_cells = nrow(niche_data),
  n_current_cells = sum(niche_data$current_weight > 0, na.rm = TRUE),
  n_increasing = sum(
    niche_data$temporal_class == "Increasing",
    na.rm = TRUE
  ),
  n_flat_high = sum(
    niche_data$temporal_class == "Flat high",
    na.rm = TRUE
  ),
  PCA_variance_PC1 = variance_explained[1],
  PCA_variance_PC2 = variance_explained[2],
  Schoener_D = schoener_D,
  ND = niche_divergence,
  PI = proportion_inside,
  OE = overestimation
)

print(results, width = Inf)

