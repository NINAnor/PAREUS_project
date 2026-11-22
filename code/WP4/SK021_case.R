############################################################################
# Pareus WP4 - Case study SK021 (Slovakia)
# Workflow (uses the streamlined function files 01-04 + wp4_functions_utils.R):
#   01 PA gap analysis -> 02 prepare optimisation grid -> 03 PA optimisation
#   (agriculture allowed vs. excluded, global target) -> 04 OECM / PCA corridors
# Inputs: data/shared/SK021.gpkg, data/WP2/SK021, data/WP3, data/WP4
# Outputs: outputs/WP4/01_PA_analysis, 02_optim, 03_pca_landscape
############################################################################
library(prioritizr)
library(ggplot2)
library(sf)
library(terra)
library(dplyr)
library(tidyverse)
library(ggnewscale)
stud_site<-"SK021"
shared_dir = "data/shared"
stud_area <- sf::read_sf(file.path(shared_dir, paste0(stud_site, ".gpkg")))
source("code/WP4/wp4_functions_utils.R")
source("code/WP4/01_pa_analysis_funtions.R")

out <- "outputs/WP4/01_PA_analysis"

x2 <- build_pa_grid(stud_site, lulc_method = "majority", out_dir = out)
res <- analyse_pa_gap(x2, global_target = 0.1,
                      lulc_targets = c(forest = 0.1, water = 0.1, wetland = 0.1, agricultural = 0),
                      out_dir = out)

#### 2. prep optimization
source("code/WP4/02_prep_optim_funtions.R")
op <- prepare_optim_grid(stud_site, na_action = "neighbour_mean", skip_missing = TRUE,     layer_files  = c(
  eco_cond    = "data/WP2/{site}_ECSI.tif",
  cost_policy = "data/WP3/output/pol_resistance_{site}.tif"
))

### 3. optim PA
source("code/WP4/03_pa_optim_functions.R")

feat_cols <- c("sampled_eco_cond_scaled", "sampled_es_reg_scaled", "sampled_cost_policy_scaled",
               "sampled_es_prov_scaled", "sampled_es_cult_scaled",
               "min_distance_scaled")

SK_feature_sets <- list(
  nature = data.frame(column = feat_cols, weight = 1,
                      direction = c("max", "max", "max", "min", "min", "min")),
  multifunctional = data.frame(column = feat_cols, weight = 1,
                               direction = c("max", "max", "max", "max", "max", "min"))
)


res_SK_glob <- run_pa_optimization(
  stud_site,
  scenario         = "nature",
  feature_sets     = SK_feature_sets,
  pu = NULL, gap = NULL, stud_area = NULL,    # objects, or read from the dirs below
  gap_dir          = "outputs/WP4/01_PA_analysis",
  pu_dir           = "outputs/WP4/02_optim",
  shared_dir       = "data/shared",
  out_dir          = "outputs/WP4/02_optim",
  budget_col       = "gap_area",           # or "gap_area"
  dist_transform   = "log",
  min_area_km2     = 0.075, keep_core = TRUE,
  area_weighted    = T,
  lulcs            = NULL, run_global = TRUE,
  global_lockout   = c("built-up", "agricultural"),
  boundary_penalty = 5e-4, edge_factor = 0.5, share_boundary = TRUE,
  mip_gap = 0.01, time_limit = 600,
  write            = T
)

res_SK_agri <- run_pa_optimization(
  stud_site,
  scenario         = "nature",
  feature_sets     = SK_feature_sets,
  pu = NULL, gap = NULL, stud_area = NULL,    # objects, or read from the dirs below
  gap_dir          = "outputs/WP4/01_PA_analysis",
  pu_dir           = "outputs/WP4/02_optim",
  shared_dir       = "data/shared",
  out_dir          = "outputs/WP4/02_optim",
  budget_col       = "gap_area",           # or "gap_area"
  dist_transform   = "log",
  min_area_km2     = 0.075, keep_core = TRUE,
  area_weighted    = T,
  lulcs            = NULL, run_global = TRUE,
  global_lockout   = c("built-up"),
  boundary_penalty = 5e-4, edge_factor = 0.5, share_boundary = TRUE,
  mip_gap = 0.01, time_limit = 600,
  write            = T
)

plot(res_SK_agri$plots$map_global)
res_SK_NOagri$log

## PA area stats:

stats<-rbind(res_SK_agri$stats%>%filter(scenario == "GLOBAL")%>%mutate(scenario = case_when(scenario == "GLOBAL" ~"AGRI_ALLOWED")),
             res_SK_glob$stats%>%filter(scenario == "GLOBAL"))

exclude_lulc = c("built-up")
d <- stats[stats$core_pa %in% pa_levels[1:3] & !stats$lulc_name %in% exclude_lulc, ]
ggplot2::ggplot(d, ggplot2::aes(x = scenario, y = km2, fill = core_pa)) +
  ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.8), width = 0.7) +
  ggplot2::facet_wrap(~ lulc_name, scales = "free_y") +
  ggplot2::scale_fill_manual(values = pa_cols, breaks = pa_levels[1:3]) +
  ggplot2::labs(x = NULL, y = expression("Area (km"^2*")"), fill = NULL) +
  ggplot2::theme_bw() +
  ggplot2::theme(legend.position = "top", strip.background = ggplot2::element_blank(),
                 axis.text.x = ggplot2::element_text(angle = 30, hjust = 1))




source("code/WP4/04_OECM_functions.R")
scenario<-"nature"# "nature" #"intermediate" "multifunctional"

# PUs containing the optimized core PA
pu_agri<-st_read(paste0("outputs/WP4/02_optim/PA_optim_grid_",stud_site,"_agri.geojson"))
#pu_glob<-st_read(paste0("outputs/WP4/02_optim/PA_optim_grid_",stud_site,"_global.geojson"))
#agri <- run_pca_corridors("SK021", "agri")
agri<-run_pca_corridors_grid("SK021", scenario = "agri",
                                         levels = c(low = 0.1, medium = 0.5, high = 0.8),
                                         n_rep = 3, seed = 1, k_nearest = 3,
                                         optim_dir = "outputs/WP4/02_optim",
                                         shared_dir = "data/shared",
                                         out_dir = "outputs/WP4/03_pca_landscape",
                                         write_maps = FALSE)

glob<-run_pca_corridors_grid("SK021", scenario = "global",
                             levels = c(low = 0.1, medium = 0.5, high = 0.8),
                             n_rep = 3, seed = 1, k_nearest = 3,
                             optim_dir = "outputs/WP4/02_optim",
                             shared_dir = "data/shared",
                             out_dir = "outputs/WP4/03_pca_landscape",
                             write_maps = FALSE)


human_c<-list(
  w_cult = 0.45,
  w_prov =0.45,
  w_reg = 0.1
)

oecm_lin_sk<-function(cult_es,prov_es,ec,w_cult,w_prov,w_ec){
  weighted_sum <- cult_es * w_cult +
    prov_es * w_prov +
    ec * w_ec 
  
  weight_total <- w_cult + w_prov + w_ec 
  
  return(weighted_sum / weight_total)
}
glob$pu$oecm_suit <- oecm_lin_sk(glob$pu$sampled_es_cult_scaled,
                                 glob$pu$sampled_es_prov_scaled,
                                 glob$pu$sampled_es_reg_scaled,
                                 human_c$w_cult, human_c$w_prov, human_c$w_reg)  # check third weight

agri$pu$oecm_suit <- oecm_lin_sk(agri$pu$sampled_es_cult_scaled,
                                 agri$pu$sampled_es_prov_scaled,
                                 agri$pu$sampled_es_reg_scaled,
                                 human_c$w_cult, human_c$w_prov, human_c$w_reg)

glob$pu <- glob$pu %>%
  dplyr::mutate(oecm_suit = dplyr::case_when(
    n_paths > 0 ~ 1,          # cell is crossed by at least one corridor path
    .default    = oecm_suit   # otherwise keep the calculated suitability
  ))

agri$pu <- agri$pu %>%
  dplyr::mutate(oecm_suit = dplyr::case_when(
    n_paths > 0 ~ 1,          # cell is crossed by at least one corridor path
    .default    = oecm_suit   # otherwise keep the calculated suitability
  ))


## calculate oecm_suit for both (paths will be set as max suitability later)
oecm_glob<-select_oecm(pu=glob$pu, 
                       mode= "global", 
                       coverage = 0.2,
                       suitability_col = "oecm_suit",
                       corePA_col = "core_pa_global",
                       search_oecm_in = c("not protected","other protected areas"),
                       w= 0.9, order = 2)

oecm_agri<-select_oecm(pu=agri$pu, 
                       mode= "global", 
                       coverage = 0.2,
                       suitability_col = "oecm_suit",
                       corePA_col = "core_pa_global",
                       search_oecm_in = c("not protected","other protected areas"),
                       w= 0.9, order = 2)

base_glob<-plot_pca_map(pu=glob$pu,
                        corePA_col = "core_pa_global",
                        oecm_df=oecm_glob,
                        stud_area,
                        scen = "Global",
                        save_output = F)

base_agri<-plot_pca_map(pu=agri$pu,
                        corePA_col = "core_pa_global",
                        oecm_df=oecm_agri,
                        stud_area,
                        scen = "AGRI",
                        save_output = T)



# pies and comparison table (summarise_area() returns the area as `km2`)
pca_labels <- c(
  "Other PA (IUCN III-VI) high suitability for OECM"     = "Other PA (IUCN III-VI)",
  "other_pa"                                             = "Other PA (IUCN III-VI)",
  "core_PA"                                              = "Optimized IUCN Ia or II",
  "not_protected"                                        = "No protection",
  "High OECM suitability not protected - potential OECM" = "Potential OECM"
)
base_agri$stats$scenario<-"AGRI_ALLOWED"

stats_all <- rbind(base_agri$stats, base_glob$stats) 
stats_all$AREA_SQKM<-stats_all$area/10^6

fill_values = c(  "existing core pa" = "#00A300",
                  "new core PA" = "#bf8bff",
                  "proposed upgrade existing PA" = "#e5d0ff",
                  "Other PA (IUCN III-VI) high suitability for OECM" = "#44a6c6",
                  "High OECM suitability not protected - potential OECM" = "#194553",
                  "Other PA (IUCN III-VI)" = "#ADD8E6",
                  "not protected" = "#999999"
)

plot_area_pies(stats_all%>%filter(lulc_name == "agricultural"), slice_col = "pca", label = "pct", keep_slices = NULL, 
               facet_col= NULL, fill_values = fill_values)

