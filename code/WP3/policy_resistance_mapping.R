############################################################################
# Pareus WP3 - mapping policy resistance
# Date: 06.08.2026
# Author: Reto Spielhofer (Norwegian Institute for Nature Research)
############################################################################

# ------------------------------------------------------------------
# FRA: from ES impact values + LULC mean impact -> final_resistance_norm
# Inputs (CSV):
#   xx_es_impact.csv : kpt, kpt_name, 10 ES columns (averaged participant ES impact values)
#   xx_lulc_mean.csv : kpt, 5 LULC columns          (averaged participant LULC impact values)
# ------------------------------------------------------------------

library(terra)
library(sf)
library(dplyr)
## map policy according to influence on landscape type
proj_id<-"FRA"

## CORINE land use land cover with 24 classes
lulc_rast<-rast(paste0("data/shared/",proj_id,"_lulc.tif"))
lulc_rast[lulc_rast == 0] <- NA

es_file   <- paste0("data/WP3/curated_input_mapping",proj_id,"_es_impact.csv")
lulc_file <- paste0("data/WP3/curated_input_mapping",proj_id,"_lulc_impact.csv")

es   <- read.csv(es_file,   sep = ";", dec=",")
lulc <- read.csv(lulc_file, sep=";", dec = ",")

# ES matrix: rows = KPT, cols = 10 ES (drop id columns)
es_mat <- as.matrix(es[, !(names(es) %in% c("KPT", "name"))])
rownames(es_mat) <- es$KPT

# LULC mean matrix: rows = KPT, cols = LULC classes
lulc_mat <- as.matrix(lulc[, names(lulc) != "KPT"])
rownames(lulc_mat) <- lulc$KPT
lulc_mat <- lulc_mat[rownames(es_mat), , drop = FALSE]   # same KPT order

n_kpt <- nrow(es_mat)

# 1) Correlation between KPTs across the 10 ES  (Excel CORREL)  -> corr_es [KPT x KPT]
corr_es <- cor(t(es_mat))                      # Pearson, 7x7

# 2) Rescale correlation from [-1,1] to [0,1]                   -> norm_corr
norm_corr <- (corr_es + 1) / 2

# 3) Normalise LULC mean impact (scale -3..+3) to [0,1]         -> LULC_IMP "norm"
#    (Excel block J3:N9 is (mean + 3) / 6; verified against the sheet)
lulc_norm <- (lulc_mat + 3) / 6

# 4) Correlation-weighted LULC impact                           -> mult_norm_mat
#    Excel: MMULT(norm_corr, lulc_norm) / COUNT(kpt)
mult_norm_mat <- (norm_corr %*% lulc_norm) / n_kpt

# 5) AHP weights per KPT (equal weights, 1/n)                   -> AHP KPT
ahp_w <- rep(1 / n_kpt, n_kpt)

# 6) Weighted sum over KPTs -> one value per LULC class         -> final_resistance_norm
#    Excel: MMULT(AHP_row, mult_norm_mat)
lulc_resistance<-data_frame(CLC = c(1,2,3,4,5), policy_resistance = as.vector(ahp_w %*% mult_norm_mat))
#lulc_resistance$policy_resistance <- as.vector(ahp_w %*% mult_norm_mat)


#recl lulc to 5 classes
lulc_rast <- lulc_rast %>%
  (\(x) floor(x / 100))() %>%
  disagg(fact = 5, method = "near")


rcl <- as.matrix(lulc_resistance[, c("CLC", "policy_resistance")])



# Map LULC classes to lookup values
pol_res <- classify(
  lulc_rast,
  rcl,
  others = NA
)

writeRaster(pol_res,paste0("data/WP3/output/pol_resistance_",proj_id,".tif"))

write.csv(data.frame(lulc_resistance),
          paste0("data/WP3/output/pol_resistance_norm_",proj_id,".csv"), row.names = FALSE)
