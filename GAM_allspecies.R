# =====================================================================
# TEMPORAL-SIGNAL ANALYSIS ACROSS ALL SPECIES
#
# Analytical logic
# ----------------
# 1. Reconstruct one annual presence-absence history per species x cell.
# 2. Retain cells with at least 3 presences and 3 absences.
# 3. Fit the same binomial GAM used in the Monte Carlo calibration.
# 4. Classify every species x cell trajectory with its own GAM metrics.
# 5. Summarise the resulting classes within species and within cells.
#
# IMPORTANT:
# The decision tree is applied only to individual species x cell histories.
# Median metrics across cells or species are descriptive and are not passed
# through the calibrated decision tree.
# =====================================================================

library(mgcv)
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)

# ---------------------------------------------------------------------
# Settings and recalibrated thresholds
# ---------------------------------------------------------------------

years <- 1996:2019
n_years <- length(years)

minimum_occurrences <- 20
minimum_presences <- 3
minimum_absences <- 3

R2_threshold <- 0.1304739
low_occ_threshold <- 0.3333333
high_occ_threshold <- 0.6666667

category_levels <- c(
  "Increasing",
  "Decreasing",
  "Flat low",
  "Flat high",
  "Noisy"
)

category_colours <- c(
  "Increasing" = "#CC79A7",
  "Decreasing" = "#0072B2",
  "Flat low" = "#08306B",
  "Flat high" = "#4A1486",
  "Noisy" = "#BDBDBD",
  "Mixed/tied" = "#666666"
)

# ---------------------------------------------------------------------
# Species inclusion criterion
# ---------------------------------------------------------------------

species_occurrence_counts <- occ_table %>%
  filter(
    !is.na(species),
    !is.na(cell_id),
    year %in% years
  ) %>%
  count(species, name = "n_occurrence_records")

species_to_process <- species_occurrence_counts %>%
  filter(n_occurrence_records > minimum_occurrences) %>%
  pull(species)

cat("Species retained:", length(species_to_process), "\n")

# ---------------------------------------------------------------------
# Fit one cell-level GAM and extract calibration-compatible metrics
# ---------------------------------------------------------------------

summarise_one_cell <- function(cell_data) {
  model_data <- cell_data %>%
    arrange(YEAR) %>%
    mutate(YEAR_rel = YEAR - min(years))
  
  fitted_model <- tryCatch(
    gam(
      presence ~ s(YEAR_rel, bs = "cs", k = 4),
      family = binomial,
      method = "REML",
      data = model_data
    ),
    error = function(e) NULL
  )
  
  null_model <- tryCatch(
    gam(
      presence ~ 1,
      family = binomial,
      method = "REML",
      data = model_data
    ),
    error = function(e) NULL
  )
  
  if (is.null(fitted_model) || is.null(null_model)) {
    return(tibble(
      pred_1996 = NA_real_,
      pred_2019 = NA_real_,
      DP = NA_real_,
      mean_probability = NA_real_,
      R2_McFadden = NA_real_
    ))
  }
  
  predicted_probability <- predict(
    fitted_model,
    newdata = model_data,
    type = "response"
  )
  
  ll_model <- as.numeric(logLik(fitted_model))
  ll_null <- as.numeric(logLik(null_model))
  
  tibble(
    pred_1996 = predicted_probability[1],
    pred_2019 = predicted_probability[n_years],
    DP = predicted_probability[n_years] - predicted_probability[1],
    # Mean across all 24 annual predictions, matching the calibration
    mean_probability = mean(predicted_probability),
    R2_McFadden = 1 - (ll_model / ll_null)
  )
}

# ---------------------------------------------------------------------
# Apply the calibrated decision tree to one trajectory
# ---------------------------------------------------------------------

classify_trajectory <- function(R2_McFadden, DP, mean_probability) {
  case_when(
    is.na(R2_McFadden) |
      is.na(DP) |
      is.na(mean_probability) ~ NA_character_,
    
    R2_McFadden > R2_threshold & DP > 0 ~ "Increasing",
    R2_McFadden > R2_threshold & DP < 0 ~ "Decreasing",
    
    R2_McFadden <= R2_threshold &
      mean_probability < low_occ_threshold ~ "Flat low",
    
    R2_McFadden <= R2_threshold &
      mean_probability > high_occ_threshold ~ "Flat high",
    
    R2_McFadden <= R2_threshold ~ "Noisy",
    
    # Includes the rare case R2 > threshold and DP == 0
    TRUE ~ NA_character_
  )
}

# =====================================================================
# FIT ALL ELIGIBLE SPECIES x CELL HISTORIES
# =====================================================================

results_list <- vector("list", length(species_to_process))
names(results_list) <- species_to_process

eligibility_list <- vector("list", length(species_to_process))
names(eligibility_list) <- species_to_process

for (species_name in species_to_process) {
  cat("Processing:", species_name, "\n")
  
  # Presence records for the focal species
  species_presence <- occ_table %>%
    filter(
      species == species_name,
      !is.na(cell_id),
      year %in% years
    ) %>%
    distinct(cell_id, year) %>%
    mutate(presence = 1L)
  
  occupied_cells <- species_presence %>%
    distinct(cell_id)
  
  # Complete every occupied cell over all study years
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
      n_pres = sum(presence == 1L),
      n_abs = sum(presence == 0L),
      .groups = "drop"
    ) %>%
    filter(
      n_pres >= minimum_presences,
      n_abs >= minimum_absences
    )
  
  eligibility_list[[species_name]] <- tibble(
    valid_name = species_name,
    n_occupied_cells = nrow(occupied_cells),
    n_eligible_cells = nrow(eligible_cells),
    eligible_ratio = if_else(
      nrow(occupied_cells) > 0,
      nrow(eligible_cells) / nrow(occupied_cells),
      NA_real_
    )
  )
  
  if (nrow(eligible_cells) == 0) {
    next
  }
  
  model_data <- species_cell_year %>%
    semi_join(eligible_cells, by = "cell_id") %>%
    transmute(
      cell_id,
      YEAR = year,
      presence
    )
  
  species_results <- model_data %>%
    group_by(cell_id) %>%
    group_modify(~ summarise_one_cell(.x)) %>%
    ungroup() %>%
    mutate(
      valid_name = species_name,
      category = classify_trajectory(
        R2_McFadden,
        DP,
        mean_probability
      )
    ) %>%
    dplyr::select(
      valid_name,
      cell_id,
      pred_1996,
      pred_2019,
      DP,
      mean_probability,
      R2_McFadden,
      category
    )
  
  results_list[[species_name]] <- species_results
}

final_cell_table_all_species <- bind_rows(results_list)
species_eligibility_summary <- bind_rows(eligibility_list)

# ----------------------------------------------------------------------------
# Output directory and shared aesthetics
# ----------------------------------------------------------------------------

figure_directory <- "all_species_figure_outputs"
dir.create(figure_directory, showWarnings = FALSE, recursive = TRUE)

temporal_classes <- c(
  "Increasing", "Decreasing", "Flat low", "Flat high", "Noisy"
)

summary_classes <- c(
  temporal_classes, "Mixed/tied", "Unclassified"
)

class_colours <- c(
  "Increasing"   = "#CC79A7",
  "Decreasing"   = "#0072B2",
  "Flat low"     = "#08306B",
  "Flat high"    = "#4A1486",
  "Noisy"        = "#BDBDBD",
  "Mixed/tied"   = "#666666",
  "Unclassified" = "#F2F2F2"
)

density_fill <- "#BFD7EA"
density_line <- "#333333"
median_colour <- "#D55E00"
mean_colour <- "#0072B2"

theme_temporal <- function(base_size = 15) {
  theme_classic(base_size = base_size) +
    theme(
      axis.title = element_text(size = base_size + 1),
      axis.text = element_text(size = base_size),
      plot.title = element_text(size = base_size + 2, face = "bold"),
      plot.subtitle = element_text(size = base_size - 1),
      legend.title = element_text(size = base_size),
      legend.text = element_text(size = base_size - 1),
      plot.tag = element_text(size = base_size + 4, face = "bold"),
      plot.margin = margin(10, 12, 10, 10)
    )
}

# Use 0 and 1 at the endpoints rather than 0.0 and 1.0.
probability_breaks <- seq(0, 1, by = 0.2)
probability_labels <- c("0", "0.2", "0.4", "0.6", "0.8", "1")

# ----------------------------------------------------------------------------
# Rebuild the summaries needed by the figures
# ----------------------------------------------------------------------------

# All eligible species x cell histories, including failed/unresolved models.
overall_class_summary <- final_cell_table_all_species %>%
  transmute(
    category = replace_na(category, "Unclassified")
  ) %>%
  count(category, name = "n_histories") %>%
  complete(
    category = summary_classes,
    fill = list(n_histories = 0L)
  ) %>%
  mutate(
    total_histories = sum(n_histories),
    percentage = 100 * n_histories / total_histories,
    category = factor(category, levels = summary_classes)
  )

# Species-level descriptive metrics and totals.
species_metric_summary <- final_cell_table_all_species %>%
  group_by(valid_name) %>%
  summarise(
    median_R2 = median(R2_McFadden, na.rm = TRUE),
    median_DP = median(DP, na.rm = TRUE),
    median_probability = median(mean_probability, na.rm = TRUE),
    n_eligible_cells = n(),
    n_classified_cells = sum(!is.na(category)),
    n_unclassified_cells = sum(is.na(category)),
    .groups = "drop"
  ) %>%
  left_join(
    species_eligibility_summary %>%
      dplyr::select(valid_name, n_occupied_cells, eligible_ratio),
    by = "valid_name"
  )

species_class_counts <- final_cell_table_all_species %>%
  filter(!is.na(category)) %>%
  count(valid_name, category, name = "n_cells")

species_class_composition <- crossing(
  valid_name = species_metric_summary$valid_name,
  category = temporal_classes
) %>%
  left_join(species_class_counts, by = c("valid_name", "category")) %>%
  mutate(n_cells = replace_na(n_cells, 0L)) %>%
  left_join(
    species_metric_summary %>%
      dplyr::select(
        valid_name, n_eligible_cells,
        n_classified_cells, n_unclassified_cells
      ),
    by = "valid_name"
  ) %>%
  mutate(
    percentage = if_else(
      n_eligible_cells > 0,
      100 * n_cells / n_eligible_cells,
      NA_real_
    )
  )

# Assign each species its most frequent cell-level class.
# Equal maxima are retained as Mixed/tied; no classified cells = Unclassified.
species_dominant_class <- species_class_composition %>%
  group_by(valid_name) %>%
  summarise(
    dominant_category = case_when(
      first(n_classified_cells) == 0 ~ "Unclassified",
      sum(n_cells == max(n_cells)) > 1 ~ "Mixed/tied",
      TRUE ~ category[which.max(n_cells)]
    ),
    .groups = "drop"
  )

species_category_summary <- species_dominant_class %>%
  count(dominant_category, name = "n_species") %>%
  complete(
    dominant_category = summary_classes,
    fill = list(n_species = 0L)
  ) %>%
  mutate(
    total_species = sum(n_species),
    percentage = 100 * n_species / total_species,
    dominant_category = factor(dominant_category, levels = summary_classes)
  )

# Cell-level descriptive metrics and totals.
cell_metric_summary <- final_cell_table_all_species %>%
  group_by(cell_id) %>%
  summarise(
    median_R2 = median(R2_McFadden, na.rm = TRUE),
    median_DP = median(DP, na.rm = TRUE),
    median_probability = median(mean_probability, na.rm = TRUE),
    n_eligible_species = n_distinct(valid_name),
    n_classified_species = n_distinct(valid_name[!is.na(category)]),
    n_unclassified_species = n_distinct(valid_name[is.na(category)]),
    .groups = "drop"
  )

cell_class_counts <- final_cell_table_all_species %>%
  filter(!is.na(category)) %>%
  count(cell_id, category, name = "n_species")

cell_class_composition <- crossing(
  cell_id = cell_metric_summary$cell_id,
  category = temporal_classes
) %>%
  left_join(cell_class_counts, by = c("cell_id", "category")) %>%
  mutate(n_species = replace_na(n_species, 0L)) %>%
  left_join(
    cell_metric_summary %>%
      dplyr::select(
        cell_id, n_eligible_species,
        n_classified_species, n_unclassified_species
      ),
    by = "cell_id"
  ) %>%
  mutate(
    percentage = if_else(
      n_eligible_species > 0,
      100 * n_species / n_eligible_species,
      NA_real_
    )
  )

# Assign each cell its most frequent species-level class.
cell_dominant_class <- cell_class_composition %>%
  group_by(cell_id) %>%
  summarise(
    dominant_category = case_when(
      first(n_classified_species) == 0 ~ "Unclassified",
      sum(n_species == max(n_species)) > 1 ~ "Mixed/tied",
      TRUE ~ category[which.max(n_species)]
    ),
    .groups = "drop"
  )

cell_category_summary <- cell_dominant_class %>%
  count(dominant_category, name = "n_cells") %>%
  complete(
    dominant_category = summary_classes,
    fill = list(n_cells = 0L)
  ) %>%
  mutate(
    total_cells = sum(n_cells),
    percentage = 100 * n_cells / total_cells,
    dominant_category = factor(dominant_category, levels = summary_classes)
  )

# ----------------------------------------------------------------------------
# Helper functions
# ----------------------------------------------------------------------------

make_dominant_plot <- function(data, count_column, y_label, panel_title) {
  count_name <- rlang::as_name(rlang::ensym(count_column))
  
  plot_data <- data %>%
    filter(.data[[count_name]] > 0) %>%
    mutate(
      label = sprintf(
        "%s\n(%.1f%%)",
        format(.data[[count_name]], big.mark = ","),
        percentage
      )
    )
  
  ggplot(
    plot_data,
    aes(x = dominant_category, y = percentage, fill = dominant_category)
  ) +
    geom_col(width = 0.68, colour = "grey25", linewidth = 0.35) +
    geom_text(aes(label = label), vjust = -0.35, size = 4.4) +
    scale_fill_manual(values = class_colours, drop = FALSE) +
    scale_y_continuous(
      limits = c(0, 100),
      breaks = seq(0, 100, by = 20),
      expand = expansion(mult = c(0, 0.02))
    ) +
    labs(title = panel_title, x = NULL, y = y_label) +
    theme_temporal() +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 25, hjust = 1)
    )
}

make_metric_density <- function(
    data, variable, median_value, x_label, panel_title,
    x_limits, x_breaks, x_labels) {
  
  ggplot(
    data %>% filter(is.finite(.data[[variable]])),
    aes(x = .data[[variable]])
  ) +
    geom_density(
      fill = density_fill,
      colour = density_line,
      alpha = 0.82,
      linewidth = 0.8,
      trim = FALSE,
      bounds = x_limits
    ) +
    geom_vline(
      xintercept = median_value,
      colour = median_colour,
      linetype = "dashed",
      linewidth = 0.9
    ) +
    annotate(
      "text",
      x = median_value,
      y = Inf,
      label = sprintf("Median = %.3f", median_value),
      colour = median_colour,
      angle = 90,
      vjust = 1.35,
      hjust = 1.05,
      size = 4.0
    ) +
    scale_x_continuous(
      limits = x_limits,
      breaks = x_breaks,
      labels = x_labels,
      expand = c(0, 0)
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
    labs(title = panel_title, x = x_label, y = "Density") +
    theme_temporal(base_size = 14)
}

# ============================================================================
# FIGURE 1 MAIN TEXT: DOMINANT TEMPORAL CLASSES
# ============================================================================

species_dominant_plot <- make_dominant_plot(
  species_category_summary,
  n_species,
  "Species (%)",
  "Species level"
)

cell_dominant_plot <- make_dominant_plot(
  cell_category_summary,
  n_cells,
  "Spatial cells (%)",
  "Spatial-cell level"
)

figure_main_dominant <-
  (species_dominant_plot | cell_dominant_plot) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 19, face = "bold"))


# ============================================================================
# FIGURE S1 ELIGIBILITY AND POOLED SPECIES x CELL HISTORIES
# ============================================================================

eligible_ratio_stats <- species_eligibility_summary %>%
  summarise(
    mean_ratio = mean(eligible_ratio, na.rm = TRUE),
    median_ratio = median(eligible_ratio, na.rm = TRUE)
  )

eligible_ratio_plot <- species_eligibility_summary %>%
  filter(is.finite(eligible_ratio)) %>%
  ggplot(aes(x = eligible_ratio)) +
  geom_density(
    fill = density_fill,
    colour = density_line,
    alpha = 0.82,
    linewidth = 0.8,
    trim = FALSE,
    bounds = c(0, 1)
  ) +
  geom_vline(
    xintercept = eligible_ratio_stats$mean_ratio,
    aes(linetype = "Mean"),
    colour = mean_colour,
    linewidth = 0.9
  ) +
  geom_vline(
    xintercept = eligible_ratio_stats$median_ratio,
    aes(linetype = "Median"),
    colour = median_colour,
    linewidth = 0.9
  ) +
  scale_linetype_manual(
    values = c("Mean" = "dotted", "Median" = "dashed"),
    name = NULL
  ) +
  scale_x_continuous(
    limits = c(0, 1),
    breaks = probability_breaks,
    labels = probability_labels,
    expand = c(0, 0)
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.06))) +
  labs(
    x = "Eligible cells / occupied cells",
    y = "Density"
  ) +
  theme_temporal() +
  theme(legend.position = c(0.82, 0.84))

history_class_plot <- overall_class_summary %>%
  filter(n_histories > 0) %>%
  mutate(
    label = sprintf(
      "%s\n(%.1f%%)",
      format(n_histories, big.mark = ","),
      percentage
    )
  ) %>%
  ggplot(aes(x = category, y = percentage, fill = category)) +
  geom_col(width = 0.68, colour = "grey25", linewidth = 0.35) +
  geom_text(aes(label = label), vjust = -0.35, size = 4.2) +
  scale_fill_manual(values = class_colours, drop = FALSE) +
  scale_y_continuous(
    limits = c(0, 100),
    breaks = seq(0, 100, by = 20),
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(x = NULL, y = "Species-cell histories (%)") +
  theme_temporal() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 25, hjust = 1)
  )

figure_supp_eligibility_histories <-
  (eligible_ratio_plot | history_class_plot) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 19, face = "bold"))


# ============================================================================
# FIGURE S2 â DESCRIPTIVE METRIC DISTRIBUTIONS
# ============================================================================

species_metric_medians <- species_metric_summary %>%
  dplyr::summarise(
    R2 = median(median_R2, na.rm = TRUE),
    DP = median(median_DP, na.rm = TRUE),
    probability = median(median_probability, na.rm = TRUE)
  )

cell_metric_medians <- cell_metric_summary %>%
  dplyr::summarise(
    R2 = median(median_R2, na.rm = TRUE),
    DP = median(median_DP, na.rm = TRUE),
    probability = median(median_probability, na.rm = TRUE)
  )

species_R2_plot <- make_metric_density(
  species_metric_summary, "median_R2", species_metric_medians$R2,
  expression("Median McFadden " * R^2 * " across cells"),
  "Species level: temporal structure",
  c(0, 1), probability_breaks, probability_labels
)

species_DP_plot <- make_metric_density(
  species_metric_summary, "median_DP", species_metric_medians$DP,
  expression("Median " * Delta * P * " across cells"),
  "Species level: direction of change",
  c(-1, 1), seq(-1, 1, by = 0.5), c("-1", "-0.5", "0", "0.5", "1")
)

species_P_plot <- make_metric_density(
  species_metric_summary, "median_probability", species_metric_medians$probability,
  expression("Median " * bar(P) * " across cells"),
  "Species level: mean probability",
  c(0, 1), probability_breaks, probability_labels
)

cell_R2_plot <- make_metric_density(
  cell_metric_summary, "median_R2", cell_metric_medians$R2,
  expression("Median McFadden " * R^2 * " across species"),
  "Spatial-cell level: temporal structure",
  c(0, 1), probability_breaks, probability_labels
)

cell_DP_plot <- make_metric_density(
  cell_metric_summary, "median_DP", cell_metric_medians$DP,
  expression("Median " * Delta * P * " across species"),
  "Spatial-cell level: direction of change",
  c(-1, 1), seq(-1, 1, by = 0.5), c("-1", "-0.5", "0", "0.5", "1")
)

cell_P_plot <- make_metric_density(
  cell_metric_summary, "median_probability", cell_metric_medians$probability,
  expression("Median " * bar(P) * " across species"),
  "Spatial-cell level: mean probability",
  c(0, 1), probability_breaks, probability_labels
)

figure_supp_metric_densities <-
  (species_R2_plot | species_DP_plot | species_P_plot) /
  (cell_R2_plot | cell_DP_plot | cell_P_plot) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 18, face = "bold"))
