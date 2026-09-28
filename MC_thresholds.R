# =====================================================================
# MONTE CARLO CALIBRATION OF TEMPORAL STRUCTURE METRICS
#
# Key consistency rule:
# Every simulated time series must satisfy the same eligibility criterion
# used for empirical species x cell histories: at least 3 presences and
# at least 3 absences.
#
# Five temporal classes are represented by seven trajectory variants:
#   Increasing: abrupt and gradual
#   Decreasing: abrupt and gradual
#   Flat low
#   Noisy
#   Flat high
# =====================================================================

library(mgcv)
library(dplyr)
library(purrr)
library(ggplot2)

set.seed(49)

# ---------------------------------------------------------------------
# Settings
# ---------------------------------------------------------------------

years <- 1996:2019
n_years <- length(years)       # 24 annual observations
n_reps <- 1000
min_presences <- 3
min_absences <- 3

# ---------------------------------------------------------------------
# Probability generators
#
# Flat low and Flat high cover ranges rather than one fixed probability.
# A single constant probability is drawn for each replicate, so these
# profiles remain temporally flat.
# ---------------------------------------------------------------------

profile_generators <- list(
  increasing_abrupt = function() {
    c(rep(0.05, n_years / 2), rep(0.95, n_years / 2))
  },
  increasing_gradual = function() {
    seq(0.05, 0.95, length.out = n_years)
  },
  decreasing_abrupt = function() {
    c(rep(0.95, n_years / 2), rep(0.05, n_years / 2))
  },
  decreasing_gradual = function() {
    seq(0.95, 0.05, length.out = n_years)
  },
  flat_low = function() {
    rep(runif(1, min = 0.10, max = 0.20), n_years)
  },
  noisy = function() {
    rep(0.50, n_years)
  },
  flat_high = function() {
    rep(runif(1, min = 0.80, max = 0.90), n_years)
  }
)

profile_key <- tibble(
  profile = names(profile_generators),
  temporal_class = c(
    "Increasing", "Increasing",
    "Decreasing", "Decreasing",
    "Flat low", "Noisy", "Flat high"
  ),
  directional = c(TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE)
)

# ---------------------------------------------------------------------
# Simulate one eligible presence-absence series
#
# Rejection sampling is used so that each profile contributes exactly
# n_reps series satisfying the empirical eligibility criterion.
# ---------------------------------------------------------------------

simulate_eligible_series <- function(
    probability_generator,
    min_pres = min_presences,
    min_abs = min_absences,
    max_attempts = 100000) {
  
  for (attempt in seq_len(max_attempts)) {
    probability <- probability_generator()
    presence <- rbinom(n_years, size = 1, prob = probability)
    
    n_pres <- sum(presence == 1)
    n_abs <- sum(presence == 0)
    
    if (n_pres >= min_pres && n_abs >= min_abs) {
      return(list(
        presence = presence,
        probability = probability,
        n_presence = n_pres,
        n_absence = n_abs,
        attempts = attempt
      ))
    }
  }
  
  stop("Unable to generate an eligible time series within max_attempts.")
}

# ---------------------------------------------------------------------
# Fit the same GAM used for the empirical time series
# ---------------------------------------------------------------------

fit_metrics <- function(presence) {
  df <- tibble(
    presence = presence,
    YEAR_rel = years - min(years)
  )
  
  fit <- tryCatch(
    gam(
      presence ~ s(YEAR_rel, bs = "cs", k = 4),
      family = binomial,
      method = "REML",
      data = df
    ),
    error = function(e) NULL
  )
  
  if (is.null(fit)) {
    return(NULL)
  }
  
  null_fit <- tryCatch(
    gam(
      presence ~ 1,
      family = binomial,
      method = "REML",
      data = df
    ),
    error = function(e) NULL
  )
  
  if (is.null(null_fit)) {
    return(NULL)
  }
  
  ll_model <- as.numeric(logLik(fit))
  ll_null <- as.numeric(logLik(null_fit))
  predicted <- predict(fit, newdata = df, type = "response")
  
  tibble(
    R2_McFadden = 1 - (ll_model / ll_null),
    DP = predicted[n_years] - predicted[1],
    mean_probability = mean(predicted)
  )
}

# ---------------------------------------------------------------------
# Run Monte Carlo calibration
# ---------------------------------------------------------------------

mc_results <- pmap_dfr(
  profile_key,
  function(profile, temporal_class, directional) {
    map_dfr(
      seq_len(n_reps),
      function(replicate_id) {
        simulated <- simulate_eligible_series(
          profile_generators[[profile]]
        )
        
        metrics <- fit_metrics(simulated$presence)
        
        if (is.null(metrics)) {
          return(NULL)
        }
        
        metrics %>%
          mutate(
            profile = profile,
            temporal_class = temporal_class,
            directional = directional,
            replicate = replicate_id,
            generating_probability = mean(simulated$probability),
            n_presence = simulated$n_presence,
            n_absence = simulated$n_absence,
            simulation_attempts = simulated$attempts
          )
      }
    )
  }
)

# Confirm that the simulation and eligibility requirements were met
stopifnot(
  all(mc_results$n_presence >= min_presences),
  all(mc_results$n_absence >= min_absences)
)

simulation_diagnostics <- mc_results %>%
  group_by(profile, temporal_class) %>%
  summarise(
    n_successful = n(),
    mean_generating_probability = mean(generating_probability),
    min_presences_observed = min(n_presence),
    min_absences_observed = min(n_absence),
    mean_attempts = mean(simulation_attempts),
    max_attempts = max(simulation_attempts),
    .groups = "drop"
  )

print(simulation_diagnostics)

# ---------------------------------------------------------------------
# Calibrate R2 threshold: directional versus non-directional
# ---------------------------------------------------------------------

candidate_R2 <- sort(unique(mc_results$R2_McFadden))

R2_error <- map_dbl(
  candidate_R2,
  function(threshold) {
    predicted_directional <- mc_results$R2_McFadden > threshold
    mean(predicted_directional != mc_results$directional)
  }
)

R2_threshold <- candidate_R2[which.min(R2_error)]

# ---------------------------------------------------------------------
# Calibrate mean-probability thresholds sequentially
#
# Only non-directional simulations that reach the non-directional branch
# of the calibrated tree are used here.
# ---------------------------------------------------------------------

non_directional_results <- mc_results %>%
  filter(
    !directional,
    R2_McFadden <= R2_threshold
  )

candidate_probability <- sort(
  unique(non_directional_results$mean_probability)
)

low_error <- map_dbl(
  candidate_probability,
  function(threshold) {
    predicted_flat_low <-
      non_directional_results$mean_probability < threshold
    true_flat_low <-
      non_directional_results$temporal_class == "Flat low"
    
    mean(predicted_flat_low != true_flat_low)
  }
)

low_occ_threshold <- candidate_probability[which.min(low_error)]

high_error <- map_dbl(
  candidate_probability,
  function(threshold) {
    predicted_flat_high <-
      non_directional_results$mean_probability > threshold
    true_flat_high <-
      non_directional_results$temporal_class == "Flat high"
    
    mean(predicted_flat_high != true_flat_high)
  }
)

high_occ_threshold <- candidate_probability[which.min(high_error)]

# ---------------------------------------------------------------------
# Apply the calibrated tree back to simulations for diagnostics
# ---------------------------------------------------------------------

mc_classified <- mc_results %>%
  mutate(
    predicted_class = case_when(
      R2_McFadden > R2_threshold & DP > 0 ~ "Increasing",
      R2_McFadden > R2_threshold & DP < 0 ~ "Decreasing",
      R2_McFadden <= R2_threshold &
        mean_probability < low_occ_threshold ~ "Flat low",
      R2_McFadden <= R2_threshold &
        mean_probability > high_occ_threshold ~ "Flat high",
      R2_McFadden <= R2_threshold ~ "Noisy",
      TRUE ~ NA_character_
    )
  )

calibration_confusion <- table(
  Known = mc_classified$temporal_class,
  Predicted = mc_classified$predicted_class,
  useNA = "ifany"
)

calibration_accuracy <- mean(
  mc_classified$temporal_class == mc_classified$predicted_class,
  na.rm = TRUE
)

print(calibration_confusion)

cat("\n---------------------------------------\n")
cat("CALIBRATED THRESHOLDS\n")
cat("---------------------------------------\n")
cat("McFadden R2 threshold:", R2_threshold, "\n")
cat("Low mean-probability threshold:", low_occ_threshold, "\n")
cat("High mean-probability threshold:", high_occ_threshold, "\n")
cat("Five-class calibration accuracy:", calibration_accuracy, "\n")

# ---------------------------------------------------------------------
# Metric distributions
# ---------------------------------------------------------------------

ggplot(mc_results, aes(profile, R2_McFadden, fill = temporal_class)) +
  geom_boxplot() +
  geom_hline(
    yintercept = R2_threshold,
    linetype = "dashed",
    colour = "red"
  ) +
  coord_flip() +
  theme_bw() +
  labs(x = NULL, y = "McFadden's pseudo-R2", fill = "Class")

ggplot(
  non_directional_results,
  aes(mean_probability, colour = temporal_class)
) +
  geom_density(linewidth = 1) +
  geom_vline(
    xintercept = c(low_occ_threshold, high_occ_threshold),
    linetype = "dashed"
  ) +
  theme_bw() +
  labs(
    x = "Mean predicted probability",
    y = "Density",
    colour = "Class"
  )

ggplot(mc_results, aes(profile, DP, fill = temporal_class)) +
  geom_boxplot() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  coord_flip() +
  theme_bw() +
  labs(x = NULL, y = "DP", fill = "Class")

# ---------------------------------------------------------------------
# Representative eligible time series for the calibration figure
# ---------------------------------------------------------------------

set.seed(99)

profile_labels <- c(
  increasing_abrupt = "Increasing - abrupt",
  increasing_gradual = "Increasing - gradual",
  decreasing_abrupt = "Decreasing - abrupt",
  decreasing_gradual = "Decreasing - gradual",
  flat_low = "Flat low",
  noisy = "Noisy",
  flat_high = "Flat high"
)

profile_examples <- map_dfr(
  profile_key$profile,
  function(profile_name) {
    simulated <- simulate_eligible_series(
      profile_generators[[profile_name]]
    )
    
    example_data <- tibble(
      year = years,
      YEAR_rel = years - min(years),
      probability = simulated$probability,
      presence = simulated$presence
    )
    
    example_fit <- gam(
      presence ~ s(YEAR_rel, bs = "cs", k = 4),
      family = binomial,
      method = "REML",
      data = example_data
    )
    
    example_data %>%
      mutate(
        fitted_probability = predict(
          example_fit,
          newdata = example_data,
          type = "response"
        ),
        profile = profile_name
      )
  }
) %>%
  mutate(
    profile = factor(
      profile,
      levels = profile_key$profile,
      labels = profile_labels[profile_key$profile]
    )
  )

profile_figure <- ggplot(profile_examples, aes(x = year)) +
  geom_step(
    aes(y = presence),
    direction = "hv",
    colour = "black",
    linewidth = 0.6
  ) +
  geom_point(
    aes(y = presence),
    colour = "black",
    size = 1.8
  ) +
  geom_line(
    aes(y = fitted_probability),
    colour = "red",
    linewidth = 1.1
  ) +
  facet_wrap(~ profile, ncol = 2) +
  scale_x_continuous(breaks = seq(1996, 2019, by = 5)) +
  scale_y_continuous(
    limits = c(-0.05, 1.05),
    breaks = c(0, 1)
  ) +
  labs(x = "Year", y = "Presence-absence") +
  theme_bw(base_size = 12) +
  theme(
    strip.background = element_rect(fill = "grey90", colour = "grey30"),
    strip.text = element_text(size = 12, face = "bold"),
    panel.grid.minor = element_blank(),
    legend.position = "none"
  )

profile_figure

ggsave(
  filename = "canonical_temporal_profiles.png",
  plot = profile_figure,
  width = 10,
  height = 10,
  dpi = 300
)

mc_classified %>%
  filter(is.na(predicted_class)) %>%
  select(
    profile,
    temporal_class,
    R2_McFadden,
    DP,
    mean_probability
  )

mc_classified <- mc_results %>%
  mutate(
    predicted_class = case_when(
      R2_McFadden > R2_threshold & DP > 0 ~ "Increasing",
      R2_McFadden > R2_threshold & DP < 0 ~ "Decreasing",
      R2_McFadden > R2_threshold & DP == 0 ~ "Noisy",
      
      R2_McFadden <= R2_threshold &
        mean_probability < low_occ_threshold ~ "Flat low",
      
      R2_McFadden <= R2_threshold &
        mean_probability > high_occ_threshold ~ "Flat high",
      
      R2_McFadden <= R2_threshold ~ "Noisy"
    )
  )

calibration_accuracy <- mean(
  mc_classified$temporal_class ==
    mc_classified$predicted_class
)
