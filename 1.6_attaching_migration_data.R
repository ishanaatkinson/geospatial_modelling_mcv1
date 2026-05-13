# 13_attach_migration.R — attach Admin1 net migration rate (Sorichetta et al. 2016)
library(dplyr); library(sf); library(geodata); library(readr); library(writexl)
sf::sf_use_s2(FALSE)

input_rds  <- "C:/Users/Ishana/OneDrive - University College London/Documents/YEAR 1/PhD/MEASLES PROJECT/Data/DHS_final/all_KR_surveys_with_gps_HR_IR_maps_GC_ACLED.rds"
output_rds  <- "C:/Users/Ishana/OneDrive - University College London/Documents/YEAR 1/PhD/MEASLES PROJECT/Data/DHS_final/all_KR_surveys_with_gps_HR_IR_maps_GC_ACLED_migration.rds"
mig_dir    <- "C:/.../Data/WorldPop_Migration/SexDisaggregated_Migration"  # unzipped
gadm_dir   <- "C:/.../Data/GADM"

all_KR <- readRDS(input_rds)

dhs_to_iso3 <- c(AO="AGO",BJ="BEN",BF="BFA",BU="BDI",CM="CMR",CF="CAF",TD="TCD",
                 KM="COM",CG="COG",CD="COD",CI="CIV",ER="ERI",SZ="SWZ",ET="ETH",GA="GAB",
                 GM="GMB",GH="GHA",GN="GIN",KE="KEN",LS="LSO",LB="LBR",MD="MDG",MW="MWI",
                 ML="MLI",MR="MRT",MZ="MOZ",NM="NAM",NI="NER",NG="NGA",RW="RWA",SN="SEN",
                 SL="SLE",ZA="ZAF",SD="SDN",TZ="TZA",TG="TGO",UG="UGA",ZM="ZMB",ZW="ZWE")

# Load flow file — WorldPop ships per-country CSVs of O-D flows by Admin1
load_country_netmig <- function(iso3) {
  f <- list.files(mig_dir, pattern = paste0(iso3, ".*\\.csv$"),
                  full.names = TRUE, recursive = TRUE, ignore.case = TRUE)
  if (!length(f)) return(NULL)
  flows <- read_csv(f[1], show_col_types = FALSE)
  # Expected columns: orig_GID_1, dest_GID_1, flow_total (or sex-specific)
  flow_col <- grep("flow|total|migrants", names(flows), value = TRUE, ignore.case = TRUE)[1]
  orig_col <- grep("orig",  names(flows), value = TRUE, ignore.case = TRUE)[1]
  dest_col <- grep("dest",  names(flows), value = TRUE, ignore.case = TRUE)[1]
  out <- flows %>% group_by(GID_1 = .data[[orig_col]]) %>%
    summarise(out_flow = sum(.data[[flow_col]], na.rm = TRUE), .groups = "drop")
  ins <- flows %>% group_by(GID_1 = .data[[dest_col]]) %>%
    summarise(in_flow  = sum(.data[[flow_col]], na.rm = TRUE), .groups = "drop")
  full_join(out, ins, by = "GID_1") %>%
    mutate(net_migration = (coalesce(in_flow, 0) - coalesce(out_flow, 0)),
           net_migration_rate = net_migration / pmax(coalesce(in_flow,0)+coalesce(out_flow,0), 1))
}

gadm_cache <- list(); mig_cache <- list()
out_list <- list()

for (sname in names(all_KR)) {
  df <- all_KR[[sname]]
  iso3 <- dhs_to_iso3[toupper(substr(sname,1,2))]
  if (is.na(iso3) || !all(c("v001","LATNUM","LONGNUM") %in% names(df))) {
    df$net_migration_rate_adm1 <- NA_real_; out_list[[sname]] <- df; next
  }
  if (is.null(gadm_cache[[iso3]])) gadm_cache[[iso3]] <-
    sf::st_make_valid(sf::st_as_sf(geodata::gadm(iso3, level=1, path=gadm_dir)))
  if (is.null(mig_cache[[iso3]]))  mig_cache[[iso3]]  <- load_country_netmig(iso3)
  g1 <- gadm_cache[[iso3]]; mig <- mig_cache[[iso3]]
  if (is.null(mig)) { df$net_migration_rate_adm1 <- NA_real_; out_list[[sname]] <- df; next }
  
  cl <- df %>% distinct(v001, LATNUM, LONGNUM) %>%
    filter(!is.na(LATNUM), !(LATNUM==0 & LONGNUM==0)) %>%
    st_as_sf(coords=c("LONGNUM","LATNUM"), crs=4326) %>%
    st_join(g1[,"GID_1"], join=st_intersects) %>%
    st_drop_geometry() %>% left_join(mig %>% select(GID_1, net_migration_rate), by="GID_1")
  
  df <- df %>% mutate(v001 = as.numeric(v001)) %>%
    left_join(cl %>% mutate(v001 = as.numeric(v001)) %>%
                select(v001, net_migration_rate_adm1 = net_migration_rate), by="v001")
  out_list[[sname]] <- df
  message(sname, " ✅")
}

saveRDS(out_list, output_rds)