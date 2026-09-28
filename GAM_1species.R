# ================================================================
# APPLY CALIBRATED GAM CLASSIFICATION TO REAL SPECIES DATA
# ================================================================
library(dplyr)
library(tidyr)
library(purrr)
library(mgcv)
library(sf)
library(ggplot2)
library(viridis)
library(geodata)
library(rnaturalearth)

# ---------------------- Data ------------------------------------------------------------------------------- 
# weight according to p(last year), sign DP, r^2= low as threshold
# SINGLE SPECIES
birds <- readr::read_csv("raw_data_634 (3).csv")
unique(birds$valid_name)

# occurrences per species
species_counts <- birds %>%
  count(valid_name, sort = TRUE)

# filter: no species with less than 20 occurrences
birds_filtered <- birds %>%
  semi_join(
    species_counts %>% filter(n >= 20),
    by = "valid_name"
  )

valid_species <- species_counts %>%
  filter(n >= 20) %>%
  pull(valid_name)

birds <- birds %>%
  filter(valid_name %in% valid_species)

length(unique(birds$valid_name))

# ---------------------- Country polygon --------------------------------------------------------------------
swe_wgs  <- ne_countries(scale = "large", country = "Sweden", returnclass = "sf")
swe_3006 <- st_transform(swe_wgs, 3006)

# ---------------------- Grid -------------------------------------------------------------------------------
grid25 <- st_make_grid(swe_3006, cellsize = 25000, square = TRUE) |>
  st_as_sf() |>
  st_intersection(swe_3006) |>
  st_make_valid() |>
  mutate(cell_id = row_number())

# ---------------------- Data -> sf points -> EPSG ----------------------------------------------------------
birds_sf <- birds |>
  filter(!is.na(LONGITUDE), !is.na(LATITUDE)) |>
  st_as_sf(coords = c("LONGITUDE", "LATITUDE"), crs = 4326, remove = FALSE) |>
  st_transform(3006)


# ---------------------- spatial join: attach cell_id -------------------------------------------------------
birds_with_cell <- st_join(
  birds_sf,
  grid25[, "cell_id"],          # keep only this column + geometry
  join = st_intersects,
  left = TRUE
)

# ---------------------- Final dataset: lon, lat, year, species name, cell_id -------------------------------
occ_table <- birds_with_cell |>
  st_drop_geometry() |>
  transmute(
    lon      = LONGITUDE,
    lat      = LATITUDE,
    year     = YEAR,
    species  = valid_name,
    cell_id  = cell_id
  )

species_name <- "Turdus iliacus"

years <- 1996:2019


# ---------------------------------------------------------------
# Build presence/absence cell-year table
# ---------------------------------------------------------------

sp_cell_year_long <- occ_table %>%
  filter(
    species == species_name,
    !is.na(cell_id),
    year %in% years
  ) %>%
  distinct(cell_id, year) %>%
  mutate(
    presence = 1L
  ) %>%
  group_by(cell_id) %>%
  complete(
    year = years,
    fill = list(presence = 0L)
  ) %>%
  ungroup()



# ---------------------------------------------------------------
# Filter cells
#
# Keep only cells where the GAM has both presences and absences
# ---------------------------------------------------------------

good_cells <- sp_cell_year_long %>%
  group_by(cell_id) %>%
  summarise(
    n_pres = sum(presence == 1),
    n_abs  = sum(presence == 0),
    .groups = "drop"
  ) %>%
  filter(
    n_pres >= 3,
    n_abs >= 3
  )


dat_sp <- sp_cell_year_long %>%
  semi_join(
    good_cells,
    by="cell_id"
  ) %>%
  rename(
    YEAR = year
  )



# ---------------------------------------------------------------
# GAM per cell
#
# Outputs:
# R2_McFadden
# DP
# mean_probability
# ---------------------------------------------------------------

summarise_one_cell <- function(df_cell){
  
  
  cell <- unique(df_cell$cell_id)
  
  
  if(length(unique(df_cell$presence)) == 1){
    
    return(
      tibble(
        cell_id = cell,
        pred_1996 = NA,
        pred_2019 = NA,
        DP = NA,
        mean_probability = NA,
        R2_McFadden = NA
      )
    )
    
  }
  
  
  df_cell <- df_cell %>%
    mutate(
      YEAR_rel = YEAR - min(YEAR)
    )
  
  
  fit <- tryCatch(
    
    gam(
      presence ~ s(YEAR_rel,
                   bs="cs",
                   k=4),
      family=binomial,
      data=df_cell,
      method="REML"
    ),
    
    error=function(e) NULL
    
  )
  
  
  if(is.null(fit)){
    
    return(
      tibble(
        cell_id = cell,
        pred_1996 = NA,
        pred_2019 = NA,
        DP = NA,
        mean_probability = NA,
        R2_McFadden = NA
      )
    )
    
  }
  
  
  
  # prediction for every year
  
  nd <- tibble(
    
    YEAR = years,
    
    YEAR_rel = years-min(years)
    
  )
  
  
  pred <- predict(
    fit,
    newdata=nd,
    type="response"
  )
  
  
  
  # null model
  
  null_fit <- gam(
    
    presence ~ 1,
    
    family=binomial,
    
    data=df_cell,
    
    method="REML"
    
  )
  
  
  ll_mod <- as.numeric(logLik(fit))
  ll_null <- as.numeric(logLik(null_fit))
  
  
  r2 <- 1 - ll_mod/ll_null
  
  
  
  tibble(
    
    cell_id = cell,
    
    pred_1996 = pred[1],
    
    pred_2019 = pred[length(pred)],
    
    DP = pred[length(pred)] - pred[1],
    
    mean_probability = mean(pred),
    
    R2_McFadden = r2
    
  )
  
}



# ---------------------------------------------------------------
# Run all cells
# ---------------------------------------------------------------

cell_metrics <- dat_sp %>%
  group_by(cell_id) %>%
  group_modify(~summarise_one_cell(.x)) %>%
  ungroup()



# ================================================================
# CLASSIFICATION TREE
# ================================================================


R2_threshold <- 0.13

low_occ_threshold <- 0.333

high_occ_threshold <- 0.666



classified_cells <- cell_metrics %>%
  
  mutate(
    
    category = case_when(
      
      
      R2_McFadden > R2_threshold &
        DP > 0 ~
        
        "Increasing",
      
      
      
      R2_McFadden > R2_threshold &
        DP < 0 ~
        
        "Decreasing",
      
      
      
      R2_McFadden <= R2_threshold &
        mean_probability < low_occ_threshold ~
        
        "Flat low",
      
      
      
      R2_McFadden <= R2_threshold &
        mean_probability > high_occ_threshold ~
        
        "Flat high",
      
      
      
      R2_McFadden <= R2_threshold ~
        
        "Noisy",
      
      
      
      TRUE ~
        
        NA_character_
      
    )
    
  )



# ================================================================
# SUMMARY OF CATEGORIES
# ================================================================

category_summary <- classified_cells %>%
  
  count(category) %>%
  
  arrange(desc(n))


print(category_summary)



# ================================================================
# FINAL TABLE
#
# One row per cell with all properties
# ================================================================

final_cell_table <- classified_cells %>%
  
  dplyr::select(
    cell_id,
    pred_1996,
    pred_2019,
    DP,
    mean_probability,
    R2_McFadden,
    category
  )


head(final_cell_table)

# ================================================================
# MAP
#
# Assumes:
# grid25_sf = your Sweden 25 km grid
# containing cell_id and geometry
# ================================================================


map_data <- grid25 %>%
  
  left_join(
    final_cell_table,
    by="cell_id"
  )



# -----------------------------
# Increasing cells
# -----------------------------

ggplot() +
  
  geom_sf(
    data=swe_3006,
    fill="grey95",
    colour="grey50"
  )+
  
  geom_sf(
    data=
      map_data %>%
      filter(category=="Increasing"),
    
    aes(
      fill=R2_McFadden
    ),
    
    colour=NA
  )+
  
  scale_fill_viridis_c(
    option="C"
  )+
  
  theme_bw()+
  
  labs(
    title="Increasing cells",
    fill="McFadden R²"
  )



# -----------------------------
# Decreasing cells
# -----------------------------

ggplot() +
  
  geom_sf(
    data=swe_3006,
    fill="grey95",
    colour="grey50"
  )+
  
  geom_sf(
    data=
      map_data %>%
      filter(category=="Decreasing"),
    
    aes(
      fill=R2_McFadden
    ),
    
    colour=NA
  )+
  
  scale_fill_viridis_c(
    option="B"
  )+
  
  theme_bw()+
  
  labs(
    title="Decreasing cells",
    fill="McFadden R²"
  )



# -----------------------------
# No temporal structure
#
# split into occupancy classes
# -----------------------------

ggplot() +
  
  geom_sf(
    data=swe_3006,
    fill="grey95",
    colour="grey50"
  )+
  
  geom_sf(
    data=
      map_data %>%
      filter(
        category %in%
          c(
            "Flat low",
            "No detectable structure",
            "Flat high"
          )
      ),
    
    aes(
      fill=category
    ),
    
    colour=NA
  )+
  
  scale_fill_manual(
    values=c(
      "Flat low"="blue",
      "No detectable structure"="grey60",
      "Flat high"="red"
    )
  )+
  
  theme_bw()+
  
  labs(
    title="Non-directional cells"
  )



# ================================================================
# FINAL TEMPORAL DYNAMICS MAP
#
# Directional cells:
#   fill   = DP
#   blue   = decreasing
#   purple = increasing
#   alpha  = McFadden R²
#
# Non-directional cells:
#   Flat low  = dark blue
#   Noisy     = grey
#   Flat high = dark purple
# ================================================================

library(ggplot2)
library(sf)
library(dplyr)
library(scales)
library(ggnewscale)

# ---------------------------------------------------------------
# Calculate visual scaling ranges
# ---------------------------------------------------------------

directional_data <- classified_cells %>%
  filter(
    category %in% c("Increasing", "Decreasing"),
    is.finite(R2_McFadden),
    is.finite(DP)
  )

directional_R2_range <- range(
  directional_data$R2_McFadden,
  na.rm = TRUE
)

# Use symmetrical limits so zero remains at the centre
DP_limit <- max(
  abs(directional_data$DP),
  na.rm = TRUE
)

# ---------------------------------------------------------------
# Prepare map data
# ---------------------------------------------------------------

map_data <- grid25 %>%
  left_join(
    classified_cells,
    by = "cell_id"
  ) %>%
  mutate(
    
    # DP is displayed only for directional cells
    DP_plot = if_else(
      category %in% c("Increasing", "Decreasing"),
      DP,
      NA_real_
    ),
    
    # R² controls transparency only for directional cells
    alpha_R2 = if_else(
      category %in% c("Increasing", "Decreasing"),
      scales::rescale(
        R2_McFadden,
        to = c(0.35, 1),
        from = directional_R2_range
      ),
      NA_real_
    ),
    
    # Categorical variable for non-directional cells
    nondirectional_class = if_else(
      category %in% c(
        "Flat low",
        "Noisy",
        "Flat high"
      ),
      category,
      NA_character_
    )
  )

# ---------------------------------------------------------------
# Create map
# ---------------------------------------------------------------

temporal_map <- ggplot() +
  
  # Sweden background
  geom_sf(
    data = swe_3006,
    fill = "grey97",
    colour = "grey65",
    linewidth = 0.3
  ) +
  
  # -------------------------------------------------------------
# Non-directional cells
# -------------------------------------------------------------

geom_sf(
  data = map_data %>%
    filter(!is.na(nondirectional_class)),
  aes(fill = nondirectional_class),
  colour = "grey10",
  linewidth = 0.18
) +
  
  scale_fill_manual(
    values = c(
      "Flat low"  = "#08306B",
      "Noisy"     = "#BDBDBD",
      "Flat high" = "#4A1486"
    ),
    breaks = c(
      "Flat high",
      "Noisy",
      "Flat low"
    ),
    labels = c(
      "Flat high",
      "Noisy",
      "Flat low"
    ),
    name = "Non-directional class",    
    na.translate = FALSE,
    guide = guide_legend(
      order = 1,
      title.position = "top",
      title.hjust = 0,
      keyheight = grid::unit(6, "mm"),
      keywidth = grid::unit(6, "mm")
    )
  ) +
  
  # Start a new fill scale for directional cells
  ggnewscale::new_scale_fill() +
  
  # -------------------------------------------------------------
# Directional cells
# -------------------------------------------------------------

geom_sf(
  data = map_data %>%
    filter(category %in% c("Increasing", "Decreasing")),
  aes(
    fill = DP_plot,
    alpha = alpha_R2
  ),
  colour = "grey45",
  linewidth = 0.08
) +
  
  scale_fill_gradient2(
    low = "#0072B2",
    mid = "#F7F7F7",
    high = "#CC79A7",
    midpoint = 0,
    limits = c(-DP_limit, DP_limit),
    oob = scales::squish,
    name = expression(Delta * P),
    guide = guide_colourbar(
      order = 2,
      title.position = "top",
      title.hjust = 0.5,
      barheight = grid::unit(35, "mm"),
      barwidth = grid::unit(6, "mm"),
      ticks = TRUE,
      frame.colour = "grey30"
    )
  ) +
  
  # Apply R² transparency without displaying an alpha legend
  scale_alpha_identity(
    guide = "none"
  ) +
  
  # -------------------------------------------------------------
# Map layout
# -------------------------------------------------------------

coord_sf(datum = NA) +
  
  theme_void() +
  
  theme(
    legend.position = "right",
    
    # Arrange the two legends next to each other
    legend.box = "horizontal",
    legend.box.just = "top",
    
    legend.title = element_text(
      size = 19,
      face = "bold"
    ),
    
    legend.text = element_text(
      size = 17
    ),
    
    legend.key.height = grid::unit(6, "mm"),
    legend.key.width = grid::unit(6, "mm"),
    
    # Space between the two legends
    legend.spacing.x = grid::unit(6, "mm"),
    legend.spacing.y = grid::unit(2, "mm"),
    legend.box.spacing = grid::unit(5, "mm"),
    
    plot.margin = margin(
      t = 5,
      r = 8,
      b = 5,
      l = 5
    )
  )

# Display map
temporal_map

# stats
library(dplyr)
library(ggplot2)

#------------------------------------------------------------
# Summary statistics
#------------------------------------------------------------

threshold <- 0.13

summary_stats <- final_cell_table %>%
  summarise(
    Mean   = mean(R2_McFadden, na.rm = TRUE),
    Median = median(R2_McFadden, na.rm = TRUE),
    SD     = sd(R2_McFadden, na.rm = TRUE),
    Min    = min(R2_McFadden, na.rm = TRUE),
    Max    = max(R2_McFadden, na.rm = TRUE),
    N      = n()
  )

print(summary_stats)

#------------------------------------------------------------
# Density plot
#------------------------------------------------------------

ggplot(
  final_cell_table,
  aes(x = R2_McFadden)
) +
  
  geom_density(
    fill = "#BFD7EA",
    colour = "black",
    linewidth = 0.8,
    alpha = 0.8,
    trim = TRUE
  ) +
  
  # Calibrated threshold
  geom_vline(
    xintercept = threshold,
    colour = "black",
    linewidth = 1
  ) +
  
  # Mean
  geom_vline(
    xintercept = summary_stats$Mean,
    colour = "red3",
    linetype = "dashed",
    linewidth = 0.8
  ) +
  
  # Median
  geom_vline(
    xintercept = summary_stats$Median,
    colour = "blue3",
    linetype = "dotted",
    linewidth = 0.8
  ) +
  
  scale_x_continuous(
    expand = c(0, 0)
  ) +
  
  scale_y_continuous(
    expand = c(0, 0)
  ) +
  
  labs(
    x = expression("McFadden " * R^2),
    y = "Density"
  ) +
  
  theme_classic(base_size = 19)
# number of cells per category
category_summary <- final_cell_table %>%
  count(category) %>%
  mutate(
    percentage = 100 * n / sum(n),
    percentage = round(percentage, 1)
  ) %>%
  arrange(desc(n))

category_summary

# --------------------------------------------------------------------------------------
# ---------------------- Supplementary Information -------------------------------------
# --------------------------------------------------------------------------------------

## Sensitivity analysis on eligible cells: how much does it change when changing the filter? -------------
thresholds <- 2:5

# Store the model outputs for each threshold
model_results <- list()

sensitivity_summary <- map_dfr(thresholds, function(th){
  
  # ---------------- Filter cells ----------------
  good_cells <- sp_cell_year_long %>%
    group_by(cell_id) %>%
    summarise(
      n_pres = sum(presence == 1L),
      n_abs  = sum(presence == 0L),
      .groups = "drop"
    ) %>%
    filter(n_pres >= th, n_abs >= th)
  
  dat_sp <- sp_cell_year_long %>%
    semi_join(good_cells, by = "cell_id") %>%
    transmute(
      cell_id,
      YEAR = year,
      presence
    )
  
# ---------------- Fit models ----------------
summary_table_sensitivity <- dat_sp %>%
    group_by(cell_id) %>%
    group_modify(~ summarise_one_cell(.x)) %>%
    ungroup()
  
# Save for later
model_results[[paste0("th", th)]] <<- summary_table_sensitivity
  
# ---------------- Summaries ----------------
summary_table_sensitivity %>%
    mutate(delta_p = pred_2019 - pred_1996) %>%
    summarise(
      threshold      = th,
      retained_cells = n(),
      median_R2      = median(R2_McFadden, na.rm = TRUE),
      mean_R2        = mean(R2_McFadden, na.rm = TRUE),
      median_P1996   = median(pred_1996, na.rm = TRUE),
      median_P2019   = median(pred_2019, na.rm = TRUE),
      median_deltaP  = median(delta_p, na.rm = TRUE),
      mean_deltaP    = mean(delta_p, na.rm = TRUE)
    )
  
})

# Correlations relative to threshold = 3
ref <- model_results$th3

cor_summary <- map_dfr(thresholds, function(th){
  
  dat <- inner_join(
    ref,
    model_results[[paste0("th", th)]],
    by = "cell_id",
    suffix = c("_3", "_x")
  )
  
  tibble(
    threshold = th,
    cor_P1996 = cor(dat$pred_1996_3, dat$pred_1996_x, use = "complete.obs"),
    cor_P2019 = cor(dat$pred_2019_3, dat$pred_2019_x, use = "complete.obs"),
    cor_R2    = cor(dat$R2_McFadden_3, dat$R2_McFadden_x, use = "complete.obs")
  )
})

# ---------------- Final table ----------------

sensitivity_summary <- sensitivity_summary %>%
  left_join(cor_summary, by = "threshold")

sensitivity_summary
print(sensitivity_summary, width = Inf)

## P/A ratio
cell_ratio <- sp_cell_year_long %>%
  group_by(cell_id) %>%
  summarise(
    n_presence = sum(presence),
    n_absence  = sum(presence == 0),
    pa_ratio   = n_presence / n_absence,
    .groups = "drop"
  )

cell_ratio

# Join the ratio to the grid
grid_ratio <- grid25 %>%
  left_join(cell_ratio, by = "cell_id")


# Compute the number of presences/absences and log-odds for each cell
cell_ratio <- sp_cell_year_long %>%
  group_by(cell_id) %>%
  summarise(
    n_presence = sum(presence),
    n_absence  = sum(presence == 0),
    log_odds   = log2((n_presence + 0.5) / (n_absence + 0.5)),
    .groups = "drop"
  )

# Join to the grid
grid_ratio <- grid25 %>%
  left_join(cell_ratio, by = "cell_id")

# Plot
ggplot(grid_ratio) +
  geom_sf(aes(fill = log_odds), color = NA) +
  scale_fill_viridis_c(
    option = "viridis",
    name = "log2(P/A)",
    na.value = "grey90"
  ) +
  coord_sf() +
  labs(
    title = species_name,
    subtitle = "Presence vs. absence balance (1996–2019)",
    x = NULL,
    y = NULL
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )
