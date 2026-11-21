############################################################################
# Pareus WP4 - shared utilities
# Single home for helpers and constants used by 00-04 function files and the
# case scripts (TRD_case.R, SK021_case.R, FR_case.R). Source this file first;
# the 01-04 function files do so themselves. Packages: sf, terra, dplyr
############################################################################

if (!exists(".wp4_utils_loaded")) {
.wp4_utils_loaded <- TRUE

## ---- 1. Constants ------------------------------------------------------------
default_iucn_rank <- c(
  "Not Applicable" = 1, "Not Assigned" = 1, "Not Reported" = 1,
  "VI" = 2, "V" = 3, "IV" = 4, "III" = 5, "II" = 6, "Ib" = 7, "Ia" = 8
)
default_lulc_names <- c("1" = "built-up", "2" = "agricultural", "3" = "forest",
                        "4" = "wetland", "5" = "water")
default_iucn_colours <- c(
  "Ia" = "#00A300", "Ib" = "#007A00", "II" = "#00FF00",
  "III" = "#804FB3", "IV" = "#9969C7", "V" = "#B589D6", "VI" = "#CDB4E6",
  "Not Applicable" = "#F1EE8D", "Not Assigned" = "#F1EE8D", "Not Reported" = "#F1EE8D"
)
pa_levels <- c("existing core pa", "new core PA", "proposed upgrade existing PA",
               "other protected areas", "not protected")
pa_cols <- c("existing core pa" = "#00A300", "new core PA" = "#bf8bff",
             "proposed upgrade existing PA" = "#e5d0ff",
             "other protected areas" = "#ADD8E6", "not protected" = "#1A1A1A")

## ---- 2. Helpers --------------------------------------------------------------
# projected CRS (EPSG) per study site
site_crs <- function(site, lookup = c(TRD = 25833, FRA_BAR2 = 2154, FRA = 2154, FRL04 = 2154, SVK = 5514), default = 3035)  {
  if (site %in% names(lookup)) unname(lookup[[site]]) else default
}

# fill "{site}" placeholder in a path template
fill_path <- function(template, site) gsub("{site}", site, template, fixed = TRUE)

# rescale a numeric vector to 0-1
scale01 <- function(x) {
  r <- range(x, na.rm = TRUE)
  if (diff(r) == 0) return(rep(0, length(x)))
  (x - r[1]) / diff(r)
}

# rescale columns of an sf/data.frame to 0-1 (adds "<col>_scaled")
zero_one_scale <- function(sf_df, cols = NULL, na.rm = TRUE) {
  
  # If no columns specified: use all numeric columns
  if (is.null(cols)) {
    cols <- sf_df %>%
      st_drop_geometry() %>%
      select(where(is.numeric)) %>%
      names()
  }
  
  sf_df %>%
    mutate(
      across(
        all_of(cols),
        .fns = function(x) {
          
          min_x <- min(x, na.rm = na.rm)
          max_x <- max(x, na.rm = na.rm)
          rng   <- max_x - min_x
          
          if (rng == 0) {
            return(rep(0, length(x)))
          } else {
            return((x - min_x) / rng)
          }
        },
        .names = "{.col}_scaled"
      )
    )
}

# rescale each raster layer to 0-1
min_max_normalize <- function(r) {
  # Normalize each layer individually
  norm_r <- lapp(r, function(x) {
    r_min <- min(x, na.rm = TRUE)
    r_max <- max(x, na.rm = TRUE)
    (x - r_min) / (r_max - r_min)
  })
  return(norm_r)
  
}

}
