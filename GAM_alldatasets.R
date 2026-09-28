# ============================================================================
# REVISED MULTI-DATASET TEMPORAL-SIGNAL WORKFLOW
#
# Analytical unit:
#   one species x spatial cell presence-absence history
#
# Main logic:
#   1. Create the dataset-specific spatial grid.
#   2. Calibrate the decision tree for the dataset's time-series length.
#   3. Fit and classify every eligible species x cell history.
#   4. Summarise the resulting classes within species and within cells.
#   5. Identify dominant classes only after species x cell classification.
#   6. Retain median metrics as descriptive summaries only.
#
# IMPORTANT:
#   Median R2, DP and mean probability are never passed through the decision
#   tree. The decision tree is calibrated and applied only to individual
#   species x cell histories.
# ============================================================================

library(tidyverse)
library(sf)
library(mgcv)
library(viridis)


# ============================================================================
# DATASETS AND ANALYSIS SETTINGS
# ============================================================================

dataset_settings <- tribble(
  ~dataset,          ~file,                    ~cell_size_m,
  "Birds",           "bird_sa.csv",                   25000,
  "Fish",            "fish.csv",                      10000,
  "Invertebrates",   "invertebrates.csv",              5000,
  "Phytoplankton",   "phytoplankton.csv",              5000
)

minimum_species_occurrences <- 20
minimum_presences_per_cell <- 3
minimum_absences_per_cell <- 3

gam_k <- 4
gam_method <- "REML"

# Number of simulations for each of the seven trajectories.
n_mc_replicates <- 1000
mc_seed <- 49

temporal_class_levels <- c(
  "Increasing",
  "Decreasing",
  "Flat low",
  "Flat high",
  "Noisy"
)

summary_class_levels <- c(
  temporal_class_levels,
  "Mixed/tied",
  "Unclassified"
)

class_colours <- c(
  "Increasing" = "#CC79A7",
  "Decreasing" = "#0072B2",
  "Flat low" = "#08306B",
  "Flat high" = "#4A1486",
  "Noisy" = "#BDBDBD",
  "Mixed/tied" = "#666666",
  "Unclassified" = "#F0F0F0"
)


# ============================================================================
# SAFE SUMMARY FUNCTIONS
# ============================================================================

safe_mean <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  mean(x)
}

safe_median <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  median(x)
}

safe_sd <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 2) return(NA_real_)
  sd(x)
}


# ============================================================================
# FIT ONE GAM AND EXTRACT THE THREE CLASSIFICATION METRICS
# ============================================================================

fit_temporal_metrics <- function(presence, year_values) {
  model_data <- tibble(
    presence = as.integer(presence),
    YEAR_rel = year_values - min(year_values)
  )
  
  if (
    nrow(model_data) == 0 ||
    length(unique(model_data$presence)) < 2
  ) {
    return(tibble(
      pred_start = NA_real_,
      pred_end = NA_real_,
      DP = NA_real_,
      mean_probability = NA_real_,
      R2_McFadden = NA_real_
    ))
  }
  
  full_fit <- tryCatch(
    mgcv::gam(
      presence ~ s(YEAR_rel, bs = "cs", k = gam_k),
      family = binomial(link = "logit"),
      method = gam_method,
      data = model_data
    ),
    error = function(e) NULL
  )
  
  null_fit <- tryCatch(
    mgcv::gam(
      presence ~ 1,
      family = binomial(link = "logit"),
      method = gam_method,
      data = model_data
    ),
    error = function(e) NULL
  )
  
  if (is.null(full_fit) || is.null(null_fit)) {
    return(tibble(
      pred_start = NA_real_,
      pred_end = NA_real_,
      DP = NA_real_,
      mean_probability = NA_real_,
      R2_McFadden = NA_real_
    ))
  }
  
  predicted_probability <- tryCatch(
    as.numeric(predict(full_fit, type = "response")),
    error = function(e) rep(NA_real_, nrow(model_data))
  )
  
  loglik_full <- tryCatch(
    as.numeric(logLik(full_fit)),
    error = function(e) NA_real_
  )
  
  loglik_null <- tryCatch(
    as.numeric(logLik(null_fit)),
    error = function(e) NA_real_
  )
  
  R2_McFadden <- if (
    is.finite(loglik_full) &&
    is.finite(loglik_null) &&
    loglik_null != 0
  ) {
    1 - loglik_full / loglik_null
  } else {
    NA_real_
  }
  
  tibble(
    pred_start = predicted_probability[1],
    pred_end = predicted_probability[length(predicted_probability)],
    DP = predicted_probability[length(predicted_probability)] -
      predicted_probability[1],
    mean_probability = safe_mean(predicted_probability),
    R2_McFadden = R2_McFadden
  )
}


# ============================================================================
# MONTE CARLO CALIBRATION
# ============================================================================

make_probability_profiles <- function(n_years) {
  first_half <- floor(n_years / 2)
  second_half <- n_years - first_half
  
  list(
    increasing_abrupt = c(
      rep(0.05, first_half),
      rep(0.95, second_half)
    ),
    decreasing_abrupt = c(
      rep(0.95, first_half),
      rep(0.05, second_half)
    ),
    increasing_gradual = seq(0.05, 0.95, length.out = n_years),
    decreasing_gradual = seq(0.95, 0.05, length.out = n_years),
    noisy = rep(0.50, n_years)
  )
}

simulate_eligible_series <- function(probabilities) {
  attempts <- 0L
  
  repeat {
    attempts <- attempts + 1L
    
    presence <- rbinom(
      n = length(probabilities),
      size = 1,
      prob = probabilities
    )
    
    n_presences <- sum(presence == 1L)
    n_absences <- sum(presence == 0L)
    
    if (
      n_presences >= minimum_presences_per_cell &&
      n_absences >= minimum_absences_per_cell
    ) {
      return(list(
        presence = presence,
        attempts = attempts,
        n_presences = n_presences,
        n_absences = n_absences
      ))
    }
  }
}

classify_with_thresholds <- function(
    R2,
    DP,
    mean_probability,
    R2_threshold,
    low_probability_threshold,
    high_probability_threshold
) {
  case_when(
    !is.finite(R2) |
      !is.finite(DP) |
      !is.finite(mean_probability) ~ "Unclassified",
    R2 > R2_threshold & DP > 0 ~ "Increasing",
    R2 > R2_threshold & DP < 0 ~ "Decreasing",
    R2 <= R2_threshold &
      mean_probability < low_probability_threshold ~ "Flat low",
    R2 <= R2_threshold &
      mean_probability > high_probability_threshold ~ "Flat high",
    R2 <= R2_threshold ~ "Noisy",
    TRUE ~ "Unclassified"
  )
}

calibrate_decision_tree <- function(year_values, dataset_name) {
  n_years <- length(year_values)
  
  if (
    n_years < minimum_presences_per_cell +
    minimum_absences_per_cell
  ) {
    stop(
      paste(
        dataset_name,
        "does not contain enough years for the eligibility criterion."
      )
    )
  }
  
  set.seed(mc_seed + match(dataset_name, dataset_settings$dataset))
  
  fixed_profiles <- make_probability_profiles(n_years)
  
  simulation_definitions <- tribble(
    ~profile,                ~temporal_class,
    "increasing_abrupt",     "Increasing",
    "increasing_gradual",    "Increasing",
    "decreasing_abrupt",     "Decreasing",
    "decreasing_gradual",    "Decreasing",
    "flat_low",              "Flat low",
    "flat_high",             "Flat high",
    "noisy",                 "Noisy"
  )
  
  mc_results <- map_dfr(
    seq_len(nrow(simulation_definitions)),
    function(profile_index) {
      profile_name <- simulation_definitions$profile[profile_index]
      known_class <- simulation_definitions$temporal_class[profile_index]
      
      map_dfr(
        seq_len(n_mc_replicates),
        function(replicate_number) {
          # Flat-low and flat-high probabilities vary among simulations.
          # Their ranges allow eligible histories with at least three
          # observations of the less frequent state.
          probabilities <- if (profile_name == "flat_low") {
            rep(runif(1, 0.10, 0.20), n_years)
          } else if (profile_name == "flat_high") {
            rep(runif(1, 0.80, 0.90), n_years)
          } else {
            fixed_profiles[[profile_name]]
          }
          
          simulation <- simulate_eligible_series(probabilities)
          
          metrics <- fit_temporal_metrics(
            presence = simulation$presence,
            year_values = year_values
          )
          
          metrics %>%
            mutate(
              dataset = dataset_name,
              profile = profile_name,
              temporal_class = known_class,
              replicate = replicate_number,
              generating_probability = mean(probabilities),
              n_presences = simulation$n_presences,
              n_absences = simulation$n_absences,
              attempts = simulation$attempts
            )
        }
      )
    }
  )
  
  usable_mc <- mc_results %>%
    filter(
      is.finite(R2_McFadden),
      is.finite(DP),
      is.finite(mean_probability)
    )
  
  directional_truth <- usable_mc$temporal_class %in%
    c("Increasing", "Decreasing")
  
  R2_candidates <- sort(unique(usable_mc$R2_McFadden))
  
  R2_error <- map_dbl(
    R2_candidates,
    function(candidate) {
      directional_prediction <- usable_mc$R2_McFadden > candidate
      mean(directional_prediction != directional_truth)
    }
  )
  
  R2_threshold <- R2_candidates[which.min(R2_error)]
  
  non_directional_mc <- usable_mc %>%
    filter(
      temporal_class %in% c("Flat low", "Flat high", "Noisy"),
      R2_McFadden <= R2_threshold
    )
  
  probability_candidates <- sort(
    unique(non_directional_mc$mean_probability)
  )
  
  low_error <- map_dbl(
    probability_candidates,
    function(candidate) {
      prediction <- non_directional_mc$mean_probability < candidate
      truth <- non_directional_mc$temporal_class == "Flat low"
      mean(prediction != truth)
    }
  )
  
  high_error <- map_dbl(
    probability_candidates,
    function(candidate) {
      prediction <- non_directional_mc$mean_probability > candidate
      truth <- non_directional_mc$temporal_class == "Flat high"
      mean(prediction != truth)
    }
  )
  
  low_probability_threshold <-
    probability_candidates[which.min(low_error)]
  
  high_probability_threshold <-
    probability_candidates[which.min(high_error)]
  
  classified_mc <- usable_mc %>%
    mutate(
      predicted_class = classify_with_thresholds(
        R2 = R2_McFadden,
        DP = DP,
        mean_probability = mean_probability,
        R2_threshold = R2_threshold,
        low_probability_threshold = low_probability_threshold,
        high_probability_threshold = high_probability_threshold
      )
    )
  
  confusion <- classified_mc %>%
    count(
      known_class = temporal_class,
      predicted_class,
      name = "n"
    ) %>%
    complete(
      known_class = temporal_class_levels,
      predicted_class = c(temporal_class_levels, "Unclassified"),
      fill = list(n = 0L)
    )
  
  calibration_accuracy <- mean(
    classified_mc$predicted_class == classified_mc$temporal_class
  )
  
  thresholds <- tibble(
    dataset = dataset_name,
    start_year = min(year_values),
    end_year = max(year_values),
    n_years = n_years,
    R2_threshold = R2_threshold,
    low_probability_threshold = low_probability_threshold,
    high_probability_threshold = high_probability_threshold,
    calibration_accuracy = calibration_accuracy
  )
  
  diagnostics <- mc_results %>%
    group_by(dataset, profile, temporal_class) %>%
    summarise(
      n_successful = sum(is.finite(R2_McFadden)),
      mean_generating_probability = mean(generating_probability),
      min_presences = min(n_presences),
      min_absences = min(n_absences),
      mean_attempts = mean(attempts),
      max_attempts = max(attempts),
      .groups = "drop"
    )
  
  list(
    thresholds = thresholds,
    simulations = mc_results,
    confusion = confusion,
    diagnostics = diagnostics
  )
}


# ============================================================================
# STORAGE
# ============================================================================

species_cell_results_list <- list()
coverage_list <- list()
grid_list <- list()
calibration_threshold_list <- list()
calibration_confusion_list <- list()
calibration_diagnostics_list <- list()


# ============================================================================
# DATASET LOOP
# ============================================================================

for (dataset_index in seq_len(nrow(dataset_settings))) {
  dataset_name <- dataset_settings$dataset[dataset_index]
  dataset_file <- dataset_settings$file[dataset_index]
  grid_cell_size <- dataset_settings$cell_size_m[dataset_index]
  
  cat(
    "\n====================================================\n",
    "Processing dataset: ", dataset_name, "\n",
    "Grid-cell size: ", grid_cell_size / 1000, " km\n",
    "====================================================\n",
    sep = ""
  )
  
  raw_data <- readr::read_csv(
    dataset_file,
    show_col_types = FALSE
  )
  
  required_columns <- c(
    "valid_name",
    "LONGITUDE",
    "LATITUDE",
    "YEAR"
  )
  
  missing_columns <- setdiff(required_columns, names(raw_data))
  
  if (length(missing_columns) > 0) {
    stop(
      paste0(
        "Dataset ", dataset_name,
        " is missing: ",
        paste(missing_columns, collapse = ", ")
      )
    )
  }
  
  raw_data <- raw_data %>%
    filter(
      !is.na(valid_name),
      !is.na(LONGITUDE),
      !is.na(LATITUDE),
      !is.na(YEAR)
    ) %>%
    mutate(YEAR = as.integer(YEAR))
  
  species_counts <- raw_data %>%
    count(valid_name, name = "n_occurrences")
  
  valid_species <- species_counts %>%
    filter(n_occurrences >= minimum_species_occurrences) %>%
    pull(valid_name)
  
  data_filtered <- raw_data %>%
    filter(valid_name %in% valid_species)
  
  if (nrow(data_filtered) == 0) {
    warning(paste("No species retained for", dataset_name))
    next
  }
  
  data_sf <- st_as_sf(
    data_filtered,
    coords = c("LONGITUDE", "LATITUDE"),
    crs = 4326,
    remove = FALSE
  )
  
  dataset_centroid <- data_sf %>%
    st_union() %>%
    st_centroid()
  
  centroid_coordinates <- st_coordinates(dataset_centroid)
  centroid_longitude <- centroid_coordinates[1, 1]
  centroid_latitude <- centroid_coordinates[1, 2]
  
  utm_zone <- floor((centroid_longitude + 180) / 6) + 1
  
  utm_epsg <- as.numeric(
    paste0(
      ifelse(centroid_latitude >= 0, "326", "327"),
      sprintf("%02d", utm_zone)
    )
  )
  
  data_projected <- st_transform(data_sf, crs = utm_epsg)
  
  study_area <- data_projected %>%
    st_union() %>%
    st_buffer(dist = grid_cell_size / 2) %>%
    st_union()
  
  dataset_grid <- st_make_grid(
    study_area,
    cellsize = grid_cell_size,
    square = TRUE
  ) %>%
    st_as_sf() %>%
    st_intersection(study_area) %>%
    st_make_valid() %>%
    mutate(
      cell_id = row_number(),
      dataset = dataset_name
    )
  
  grid_list[[dataset_name]] <- dataset_grid
  
  data_with_cell <- st_join(
    data_projected,
    dataset_grid %>% dplyr::select(cell_id),
    join = st_intersects,
    left = TRUE
  )
  
  occurrence_table <- data_with_cell %>%
    st_drop_geometry() %>%
    transmute(
      year = YEAR,
      species = valid_name,
      cell_id = cell_id
    ) %>%
    filter(
      !is.na(year),
      !is.na(species),
      !is.na(cell_id)
    )
  
  if (nrow(occurrence_table) == 0) {
    warning(paste("No records assigned to cells for", dataset_name))
    next
  }
  
  # Only years represented in the dataset are treated as surveyed years.
  # Missing calendar years are not automatically converted into absences.
  years <- sort(unique(occurrence_table$year))
  
  cat(
    "Calibrating ", length(years),
    "-year decision tree...\n",
    sep = ""
  )
  
  calibration <- calibrate_decision_tree(
    year_values = years,
    dataset_name = dataset_name
  )
  
  thresholds <- calibration$thresholds
  
  calibration_threshold_list[[dataset_name]] <- thresholds
  calibration_confusion_list[[dataset_name]] <-
    calibration$confusion %>% mutate(dataset = dataset_name, .before = 1)
  calibration_diagnostics_list[[dataset_name]] <- calibration$diagnostics
  
  cat(
    "R2 threshold: ", thresholds$R2_threshold, "\n",
    "Low probability threshold: ",
    thresholds$low_probability_threshold, "\n",
    "High probability threshold: ",
    thresholds$high_probability_threshold, "\n",
    "Calibration accuracy: ",
    round(100 * thresholds$calibration_accuracy, 1), "%\n",
    sep = ""
  )
  
  dataset_results <- list()
  dataset_coverage <- list()
  species_names <- sort(unique(occurrence_table$species))
  
  for (species_index in seq_along(species_names)) {
    species_name <- species_names[species_index]
    
    cat(
      "[", species_index, "/", length(species_names), "] ",
      species_name, "\n",
      sep = ""
    )
    
    species_presence <- occurrence_table %>%
      filter(species == species_name) %>%
      distinct(cell_id, year) %>%
      mutate(presence = 1L)
    
    occupied_cells <- species_presence %>%
      distinct(cell_id)
    
    species_cell_year <- species_presence %>%
      group_by(cell_id) %>%
      complete(
        year = years,
        fill = list(presence = 0L)
      ) %>%
      ungroup()
    
    eligible_cells <- species_cell_year %>%
      group_by(cell_id) %>%
      summarise(
        n_presences = sum(presence == 1L),
        n_absences = sum(presence == 0L),
        .groups = "drop"
      ) %>%
      filter(
        n_presences >= minimum_presences_per_cell,
        n_absences >= minimum_absences_per_cell
      )
    
    dataset_coverage[[species_name]] <- tibble(
      dataset = dataset_name,
      valid_name = species_name,
      n_occurrences = species_counts$n_occurrences[
        match(species_name, species_counts$valid_name)
      ],
      occupied_cells = nrow(occupied_cells),
      eligible_cells = nrow(eligible_cells),
      eligible_ratio = if_else(
        nrow(occupied_cells) > 0,
        nrow(eligible_cells) / nrow(occupied_cells),
        NA_real_
      )
    )
    
    if (nrow(eligible_cells) == 0) next
    
    model_data <- species_cell_year %>%
      semi_join(eligible_cells, by = "cell_id")
    
    species_results <- model_data %>%
      group_by(cell_id) %>%
      group_modify(
        ~ fit_temporal_metrics(
          presence = .x$presence,
          year_values = .x$year
        )
      ) %>%
      ungroup() %>%
      mutate(
        dataset = dataset_name,
        valid_name = species_name,
        category = classify_with_thresholds(
          R2 = R2_McFadden,
          DP = DP,
          mean_probability = mean_probability,
          R2_threshold = thresholds$R2_threshold,
          low_probability_threshold =
            thresholds$low_probability_threshold,
          high_probability_threshold =
            thresholds$high_probability_threshold
        )
      ) %>%
      dplyr::select(
        dataset,
        valid_name,
        cell_id,
        pred_start,
        pred_end,
        DP,
        mean_probability,
        R2_McFadden,
        category
      )
    
    dataset_results[[species_name]] <- species_results
  }
  
  species_cell_results_list[[dataset_name]] <- bind_rows(dataset_results)
  coverage_list[[dataset_name]] <- bind_rows(dataset_coverage)
}


# ============================================================================
# COMBINE DATASETS
# ============================================================================

all_species_cell_models <- bind_rows(species_cell_results_list)
all_coverage <- bind_rows(coverage_list)
all_calibration_thresholds <- bind_rows(calibration_threshold_list)
all_calibration_confusions <- bind_rows(calibration_confusion_list)
all_calibration_diagnostics <- bind_rows(calibration_diagnostics_list)


# ============================================================================
# OVERALL SPECIES x CELL COMPOSITION WITHIN EACH DATASET
# ============================================================================

overall_class_summary <- all_species_cell_models %>%
  mutate(
    category = factor(
      category,
      levels = c(temporal_class_levels, "Unclassified")
    )
  ) %>%
  count(dataset, category, name = "n_histories", .drop = FALSE) %>%
  complete(
    dataset,
    category = factor(
      c(temporal_class_levels, "Unclassified"),
      levels = c(temporal_class_levels, "Unclassified")
    ),
    fill = list(n_histories = 0L)
  ) %>%
  group_by(dataset) %>%
  mutate(
    total_histories = sum(n_histories),
    percentage = 100 * n_histories / total_histories
  ) %>%
  ungroup()


# ============================================================================
# SPECIES-LEVEL DESCRIPTIVE METRICS AND CLASS COMPOSITION
# ============================================================================

species_metric_summary <- all_species_cell_models %>%
  group_by(dataset, valid_name) %>%
  summarise(
    n_eligible_cells = n(),
    n_classified_cells = sum(category != "Unclassified"),
    n_unclassified_cells = sum(category == "Unclassified"),
    median_R2 = safe_median(R2_McFadden),
    median_DP = safe_median(DP),
    median_probability = safe_median(mean_probability),
    .groups = "drop"
  ) %>%
  left_join(
    all_coverage %>%
      dplyr::select(
        dataset,
        valid_name,
        occupied_cells,
        eligible_cells,
        eligible_ratio
      ),
    by = c("dataset", "valid_name")
  )

species_class_counts <- all_species_cell_models %>%
  filter(category != "Unclassified") %>%
  count(dataset, valid_name, category, name = "n_cells")

species_class_composition <- species_metric_summary %>%
  dplyr::select(dataset, valid_name) %>%
  crossing(category = temporal_class_levels) %>%
  left_join(
    species_class_counts,
    by = c("dataset", "valid_name", "category")
  ) %>%
  mutate(n_cells = replace_na(n_cells, 0L)) %>%
  left_join(
    species_metric_summary %>%
      dplyr::select(
        dataset,
        valid_name,
        n_eligible_cells,
        n_classified_cells,
        n_unclassified_cells
      ),
    by = c("dataset", "valid_name")
  ) %>%
  mutate(
    percentage = 100 * n_cells / n_eligible_cells,
    category = factor(category, levels = temporal_class_levels)
  )

species_dominant_class <- species_class_composition %>%
  group_by(dataset, valid_name) %>%
  filter(
    n_classified_cells > 0,
    n_cells == max(n_cells)
  ) %>%
  summarise(
    dominant_category = if (n() == 1) {
      as.character(first(category))
    } else {
      "Mixed/tied"
    },
    dominant_n_cells = max(n_cells),
    dominant_percentage = max(percentage),
    n_eligible_cells = first(n_eligible_cells),
    n_classified_cells = first(n_classified_cells),
    n_unclassified_cells = first(n_unclassified_cells),
    .groups = "drop"
  ) %>%
  right_join(
    species_metric_summary %>%
      dplyr::select(
        dataset,
        valid_name,
        n_eligible_cells,
        n_classified_cells,
        n_unclassified_cells
      ),
    by = c("dataset", "valid_name"),
    suffix = c("", ".all")
  ) %>%
  transmute(
    dataset,
    valid_name,
    dominant_category = replace_na(
      dominant_category,
      "Unclassified"
    ),
    dominant_n_cells = replace_na(dominant_n_cells, 0L),
    dominant_percentage = replace_na(dominant_percentage, 0),
    n_eligible_cells = coalesce(
      n_eligible_cells,
      n_eligible_cells.all
    ),
    n_classified_cells = coalesce(
      n_classified_cells,
      n_classified_cells.all
    ),
    n_unclassified_cells = coalesce(
      n_unclassified_cells,
      n_unclassified_cells.all
    )
  )

species_class_summary <- species_dominant_class %>%
  count(dataset, dominant_category, name = "n_species") %>%
  complete(
    dataset,
    dominant_category = summary_class_levels,
    fill = list(n_species = 0L)
  ) %>%
  group_by(dataset) %>%
  mutate(
    total_species = sum(n_species),
    percentage = 100 * n_species / total_species,
    dominant_category = factor(
      dominant_category,
      levels = summary_class_levels
    )
  ) %>%
  ungroup()

species_composition_stats <- species_class_composition %>%
  group_by(dataset, category) %>%
  summarise(
    minimum = min(percentage, na.rm = TRUE),
    q25 = quantile(percentage, 0.25, na.rm = TRUE),
    median = median(percentage, na.rm = TRUE),
    q75 = quantile(percentage, 0.75, na.rm = TRUE),
    maximum = max(percentage, na.rm = TRUE),
    .groups = "drop"
  )


# ============================================================================
# CELL-LEVEL DESCRIPTIVE METRICS AND CLASS COMPOSITION
# ============================================================================

cell_metric_summary <- all_species_cell_models %>%
  group_by(dataset, cell_id) %>%
  summarise(
    n_eligible_species = n_distinct(valid_name),
    n_classified_species = n_distinct(
      valid_name[category != "Unclassified"]
    ),
    n_unclassified_species = n_distinct(
      valid_name[category == "Unclassified"]
    ),
    median_R2 = safe_median(R2_McFadden),
    median_DP = safe_median(DP),
    median_probability = safe_median(mean_probability),
    .groups = "drop"
  )

cell_class_counts <- all_species_cell_models %>%
  filter(category != "Unclassified") %>%
  count(dataset, cell_id, category, name = "n_species")

cell_class_composition <- cell_metric_summary %>%
  dplyr::select(dataset, cell_id) %>%
  crossing(category = temporal_class_levels) %>%
  left_join(
    cell_class_counts,
    by = c("dataset", "cell_id", "category")
  ) %>%
  mutate(n_species = replace_na(n_species, 0L)) %>%
  left_join(
    cell_metric_summary %>%
      dplyr::select(
        dataset,
        cell_id,
        n_eligible_species,
        n_classified_species,
        n_unclassified_species
      ),
    by = c("dataset", "cell_id")
  ) %>%
  mutate(
    percentage = 100 * n_species / n_eligible_species,
    category = factor(category, levels = temporal_class_levels)
  )

cell_dominant_class <- cell_class_composition %>%
  group_by(dataset, cell_id) %>%
  filter(
    n_classified_species > 0,
    n_species == max(n_species)
  ) %>%
  summarise(
    dominant_category = if (n() == 1) {
      as.character(first(category))
    } else {
      "Mixed/tied"
    },
    dominant_n_species = max(n_species),
    dominant_percentage = max(percentage),
    n_eligible_species = first(n_eligible_species),
    n_classified_species = first(n_classified_species),
    n_unclassified_species = first(n_unclassified_species),
    .groups = "drop"
  ) %>%
  right_join(
    cell_metric_summary %>%
      dplyr::select(
        dataset,
        cell_id,
        n_eligible_species,
        n_classified_species,
        n_unclassified_species
      ),
    by = c("dataset", "cell_id"),
    suffix = c("", ".all")
  ) %>%
  transmute(
    dataset,
    cell_id,
    dominant_category = replace_na(
      dominant_category,
      "Unclassified"
    ),
    dominant_n_species = replace_na(dominant_n_species, 0L),
    dominant_percentage = replace_na(dominant_percentage, 0),
    n_eligible_species = coalesce(
      n_eligible_species,
      n_eligible_species.all
    ),
    n_classified_species = coalesce(
      n_classified_species,
      n_classified_species.all
    ),
    n_unclassified_species = coalesce(
      n_unclassified_species,
      n_unclassified_species.all
    )
  )

cell_class_summary <- cell_dominant_class %>%
  count(dataset, dominant_category, name = "n_cells") %>%
  complete(
    dataset,
    dominant_category = summary_class_levels,
    fill = list(n_cells = 0L)
  ) %>%
  group_by(dataset) %>%
  mutate(
    total_cells = sum(n_cells),
    percentage = 100 * n_cells / total_cells,
    dominant_category = factor(
      dominant_category,
      levels = summary_class_levels
    )
  ) %>%
  ungroup()

cell_composition_stats <- cell_class_composition %>%
  group_by(dataset, category) %>%
  summarise(
    minimum = min(percentage, na.rm = TRUE),
    q25 = quantile(percentage, 0.25, na.rm = TRUE),
    median = median(percentage, na.rm = TRUE),
    q75 = quantile(percentage, 0.75, na.rm = TRUE),
    maximum = max(percentage, na.rm = TRUE),
    .groups = "drop"
  )


# ============================================================================
# DESCRIPTIVE DISTRIBUTIONS OF MEDIAN METRICS
# ============================================================================

species_metric_distribution <- species_metric_summary %>%
  dplyr::select(
    dataset,
    valid_name,
    eligible_ratio,
    median_R2,
    median_DP,
    median_probability
  ) %>%
  pivot_longer(
    cols = c(
      eligible_ratio,
      median_R2,
      median_DP,
      median_probability
    ),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(
    metric = factor(
      metric,
      levels = c(
        "eligible_ratio",
        "median_R2",
        "median_DP",
        "median_probability"
      ),
      labels = c(
        "Eligible-cell ratio",
        "Median McFadden R2",
        "Median DP",
        "Median P-bar"
      )
    )
  )

cell_metric_distribution <- cell_metric_summary %>%
  dplyr::select(
    dataset,
    cell_id,
    median_R2,
    median_DP,
    median_probability
  ) %>%
  pivot_longer(
    cols = c(median_R2, median_DP, median_probability),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(
    metric = factor(
      metric,
      levels = c("median_R2", "median_DP", "median_probability"),
      labels = c(
        "Median McFadden R2",
        "Median DP",
        "Median P-bar"
      )
    )
  )

species_distribution_summary <- species_metric_distribution %>%
  group_by(dataset, metric) %>%
  summarise(
    n = sum(is.finite(value)),
    mean = safe_mean(value),
    sd = safe_sd(value),
    median = safe_median(value),
    q25 = quantile(value[is.finite(value)], 0.25, names = FALSE),
    q75 = quantile(value[is.finite(value)], 0.75, names = FALSE),
    minimum = min(value[is.finite(value)]),
    maximum = max(value[is.finite(value)]),
    .groups = "drop"
  )

cell_distribution_summary <- cell_metric_distribution %>%
  group_by(dataset, metric) %>%
  summarise(
    n = sum(is.finite(value)),
    mean = safe_mean(value),
    sd = safe_sd(value),
    median = safe_median(value),
    q25 = quantile(value[is.finite(value)], 0.25, names = FALSE),
    q75 = quantile(value[is.finite(value)], 0.75, names = FALSE),
    minimum = min(value[is.finite(value)]),
    maximum = max(value[is.finite(value)]),
    .groups = "drop"
  )

eligible_ratio_summary <- all_coverage %>%
  group_by(dataset) %>%
  summarise(
    n_species_after_occurrence_filter = n(),
    n_species_with_eligible_cells = sum(eligible_cells > 0),
    percentage_with_eligible_cells =
      100 * n_species_with_eligible_cells /
      n_species_after_occurrence_filter,
    mean_eligible_ratio = safe_mean(eligible_ratio),
    median_eligible_ratio = safe_median(eligible_ratio),
    sd_eligible_ratio = safe_sd(eligible_ratio),
    .groups = "drop"
  )
# COMPOSITE FIGURES: COMPARISON ACROSS DATASETS
#
# Run after the multi-dataset analysis has created:
#   all_species_cell_models
#   all_coverage
#   overall_class_summary
#   species_metric_summary
#   species_class_summary
#   cell_metric_summary
#   cell_class_summary
#
# Outputs:
#   Figure 1  — dominant temporal classes (main text)
#   Figure S1 — eligible ratios and species x cell histories
#   Figure S2 — R2McF and DP distributions at species and cell levels
# ============================================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

# ----------------------------------------------------------------------------
# Shared settings and aesthetics
# ----------------------------------------------------------------------------

figure_directory <- "all_datasets_figure_outputs"
dir.create(figure_directory, showWarnings = FALSE, recursive = TRUE)

dataset_levels <- c(
  "Birds", "Fish", "Invertebrates", "Phytoplankton"
)

temporal_class_levels <- c(
  "Increasing", "Decreasing", "Flat low", "Flat high", "Noisy"
)

summary_class_levels <- c(
  temporal_class_levels, "Mixed/tied", "Unclassified"
)

# Same temporal-class palette used for the Swedish dataset.
class_colours <- c(
  "Increasing"   = "#CC79A7",
  "Decreasing"   = "#0072B2",
  "Flat low"     = "#08306B",
  "Flat high"    = "#4A1486",
  "Noisy"        = "#BDBDBD",
  "Mixed/tied"   = "#666666",
  "Unclassified" = "#F2F2F2"
)

# Colourblind-friendly dataset palette, used for every density plot.
dataset_colours <- c(
  "Birds"         = "#6A3D9A",
  "Fish"          = "#1F78B4",
  "Invertebrates" = "#009E73",
  "Phytoplankton" = "#D55E00"
)

theme_temporal <- function(base_size = 15) {
  theme_classic(base_size = base_size) +
    theme(
      axis.title = element_text(size = base_size + 1),
      axis.text = element_text(size = base_size),
      plot.title = element_text(size = base_size + 2, face = "bold"),
      legend.title = element_text(size = base_size),
      legend.text = element_text(size = base_size - 1),
      legend.key.height = grid::unit(0.55, "cm"),
      legend.key.width = grid::unit(0.85, "cm"),
      plot.tag = element_text(size = base_size + 4, face = "bold"),
      plot.margin = margin(10, 12, 10, 10)
    )
}

probability_breaks <- seq(0, 1, by = 0.2)
probability_labels <- c("0", "0.2", "0.4", "0.6", "0.8", "1")
DP_breaks <- seq(-1, 1, by = 0.5)
DP_labels <- c("-1", "-0.5", "0", "0.5", "1")

# Keep ordering identical in all panels.
overall_class_summary <- overall_class_summary %>%
  mutate(
    dataset = factor(dataset, levels = dataset_levels),
    category = factor(category, levels = c(
      temporal_class_levels, "Unclassified"
    ))
  )

species_class_summary <- species_class_summary %>%
  mutate(
    dataset = factor(dataset, levels = dataset_levels),
    dominant_category = factor(
      dominant_category,
      levels = summary_class_levels
    )
  )

cell_class_summary <- cell_class_summary %>%
  mutate(
    dataset = factor(dataset, levels = dataset_levels),
    dominant_category = factor(
      dominant_category,
      levels = summary_class_levels
    )
  )

all_coverage <- all_coverage %>%
  mutate(dataset = factor(dataset, levels = dataset_levels))

species_metric_summary <- species_metric_summary %>%
  mutate(dataset = factor(dataset, levels = dataset_levels))

cell_metric_summary <- cell_metric_summary %>%
  mutate(dataset = factor(dataset, levels = dataset_levels))

# ============================================================================
# FIGURE 1 — MAIN TEXT: DOMINANT TEMPORAL CLASSES
# ============================================================================
# ============================================================================
# FIGURE 1 — MAIN TEXT: DOMINANT TEMPORAL CLASSES
# ============================================================================

species_dominant_plot <- ggplot(
  species_class_summary,
  aes(x = dataset, y = percentage, fill = dominant_category)
) +
  geom_col(width = 0.68, colour = "grey20", linewidth = 0.30) +
  scale_fill_manual(
    values = class_colours,
    breaks = summary_class_levels,
    drop = FALSE,
    name = "Dominant temporal class"
  ) +
  scale_y_continuous(
    breaks = seq(0, 100, by = 20),
    expand = c(0, 0)
  ) +
  coord_cartesian(ylim = c(0, 100)) +
  labs(
    title = "Species level",
    x = NULL,
    y = "Species (%)"
  ) +
  theme_temporal() +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 20, hjust = 1)
  )

cell_dominant_plot <- ggplot(
  cell_class_summary,
  aes(x = dataset, y = percentage, fill = dominant_category)
) +
  geom_col(width = 0.68, colour = "grey20", linewidth = 0.30) +
  scale_fill_manual(
    values = class_colours,
    breaks = summary_class_levels,
    drop = FALSE,
    name = "Dominant temporal class"
  ) +
  scale_y_continuous(
    breaks = seq(0, 100, by = 20),
    expand = c(0, 0)
  ) +
  coord_cartesian(ylim = c(0, 100)) +
  labs(
    title = "Spatial-cell level",
    x = NULL,
    y = "Spatial cells (%)"
  ) +
  theme_temporal() +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 20, hjust = 1)
  )

figure_main_dominant <-
  (species_dominant_plot | cell_dominant_plot) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "A") &
  theme(
    legend.position = "bottom",
    plot.tag = element_text(size = 19, face = "bold")
  )

# ============================================================================
# FIGURE S1 — ELIGIBILITY AND SPECIES x CELL HISTORIES
# ============================================================================

eligible_ratio_plot <- all_coverage %>%
  filter(is.finite(eligible_ratio)) %>%
  ggplot(
    aes(
      x = eligible_ratio,
      colour = dataset,
      fill = dataset
    )
  ) +
  geom_density(
    alpha = 0.14,
    linewidth = 1.0,
    adjust = 1,
    trim = FALSE,
    bounds = c(0, 1)
  ) +
  scale_colour_manual(
    values = dataset_colours,
    breaks = dataset_levels,
    drop = FALSE,
    name = "Dataset"
  ) +
  scale_fill_manual(
    values = dataset_colours,
    breaks = dataset_levels,
    drop = FALSE,
    name = "Dataset"
  ) +
  scale_x_continuous(
    limits = c(0, 1),
    breaks = probability_breaks,
    labels = probability_labels,
    expand = c(0, 0)
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(
    x = "Eligible cells / occupied cells",
    y = "Density"
  ) +
  theme_temporal() +
  theme(legend.position = "bottom")

history_class_plot <- ggplot(
  overall_class_summary,
  aes(x = dataset, y = percentage, fill = category)
) +
  geom_col(width = 0.68, colour = "grey20", linewidth = 0.30) +
  scale_fill_manual(
    values = class_colours,
    breaks = c(temporal_class_levels, "Unclassified"),
    drop = FALSE,
    name = "Temporal class"
  ) +
  scale_y_continuous(
    limits = c(0, 100),
    breaks = seq(0, 100, by = 20),
    expand = c(0, 0)
  ) +
  labs(
    x = NULL,
    y = "Species-cell histories (%)"
  ) +
  theme_temporal() +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 20, hjust = 1)
  )

# The panels use different legends, so they remain next to their own plots.
figure_supp_eligibility_histories <-
  (eligible_ratio_plot | history_class_plot) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 19, face = "bold"))

print(figure_supp_eligibility_histories)


# ============================================================================
# FIGURE S2 — R2McF AND DP ACROSS DATASETS
#
# Rows distinguish aggregation levels. Columns distinguish metrics.
# These medians are descriptive and are not classified by the decision tree.
# ============================================================================

make_dataset_density <- function(
    data, variable, x_label, panel_title,
    x_limits, x_breaks, x_labels,
    show_legend = FALSE) {
  
  density_data <- data %>%
    filter(is.finite(.data[[variable]]))
  
  ggplot(
    density_data,
    aes(
      x = .data[[variable]],
      colour = dataset,
      fill = dataset
    )
  ) +
    geom_density(
      alpha = 0.14,
      linewidth = 1.0,
      adjust = 1,
      trim = FALSE,
      bounds = x_limits
    ) +
    scale_colour_manual(
      values = dataset_colours,
      breaks = dataset_levels,
      drop = FALSE,
      name = "Dataset"
    ) +
    scale_fill_manual(
      values = dataset_colours,
      breaks = dataset_levels,
      drop = FALSE,
      name = "Dataset"
    ) +
    scale_x_continuous(
      limits = x_limits,
      breaks = x_breaks,
      labels = x_labels,
      expand = c(0, 0)
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    labs(
      title = panel_title,
      x = x_label,
      y = "Density"
    ) +
    theme_temporal() +
    theme(
      legend.position = if (show_legend) "bottom" else "none"
    )
}

species_R2_plot <- make_dataset_density(
  species_metric_summary,
  "median_R2",
  expression("Median McFadden " * R^2 * " across cells"),
  "Species level: temporal structure",
  c(0, 1),
  probability_breaks,
  probability_labels
)

cell_R2_plot <- make_dataset_density(
  cell_metric_summary,
  "median_R2",
  expression("Median McFadden " * R^2 * " across species"),
  "Spatial-cell level: temporal structure",
  c(0, 1),
  probability_breaks,
  probability_labels
)

species_DP_plot <- make_dataset_density(
  species_metric_summary,
  "median_DP",
  expression("Median " * Delta * P * " across cells"),
  "Species level: direction of change",
  c(-1, 1),
  DP_breaks,
  DP_labels
)

cell_DP_plot <- make_dataset_density(
  cell_metric_summary,
  "median_DP",
  expression("Median " * Delta * P * " across species"),
  "Spatial-cell level: direction of change",
  c(-1, 1),
  DP_breaks,
  DP_labels,
  show_legend = TRUE
)

# Column A/C = species summaries; column B/D = spatial-cell summaries.
figure_supp_metric_densities <-
  (species_R2_plot | cell_R2_plot) /
  (species_DP_plot | cell_DP_plot) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "A") &
  theme(
    legend.position = "bottom",
    plot.tag = element_text(size = 19, face = "bold")
  )

print(figure_supp_metric_densities)
