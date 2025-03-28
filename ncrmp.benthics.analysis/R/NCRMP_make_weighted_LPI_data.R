## Function to calculate weighted percent cover by strata and protected area and then at the regional level

# Purpose:
# support function to calculate weighted percent cover data


## Tag: data analysis


# outputs created in this file --------------
# unwh_cover_strata
# domain estimates

# Current weighting scheme:
# STRAT + PROT


# CallS:
# analysis ready data

# output gets called by:
# NCRMP_calculate_cover.R
#

# NCRMP Caribbean Benthic analytics team: Groves, viehman, Williams
# Last update: Jan 2025


##############################################################################################################################

#' Creates weighted benthic cover data
#'
#' Calcualtes weighted benthic cover data. NCRMP utilizes a stratified random
#' sampling design. Regional estimates of benthic cover are weighted by the
#' number of grid cells of a stratum in the sample frame.
#' Function produces strata means, weighted strata means,
#' and weighted regional estimates for benthic cover data.
#' Support function called by [NCRMP_calculate_cover()].
#'
#'
#'
#'
#'
#' @param inputdata A dataframe of benthic cover data summarized by cover group at each site in a single region.
#' @param region A string indicating the region. Options are: "SEFCRI", "FLK", "Tortugas", "STX", "STTSTJ", "PRICO", and "FGB".
#' @param project A string indicating the project. "NCRMP" is the only option.
#' @return A list of dataframes, including a dataframe of strata means of cover groups
#' and a dataframe of weighted regional estimates of cover groups, for specified region.
#' @importFrom magrittr "%>%"
#' @export
#'
#'
#'
NCRMP_make_weighted_LPI_data <- function(inputdata, region, project = "NULL") {
  
  # Define regional groups
  FL <- c("SEFCRI", "FLK", "Tortugas")
  GOM <- "GOM"
  Carib <- c("STTSTJ", "STX", "PRICO")
  
  #####Load NTOT data####
  ntot <- load_NTOT(region = region, inputdata = inputdata, project = project)
  
  #####processing cover data (used by both FL and non-FL regions)####
  process_cover_data <- function(data, region_is_FL = TRUE) {
    # Group data by YEAR, ANALYSIS_STRATUM, STRAT, PROT (if FL) or cover_group
    cover_est <- data %>%
      dplyr::group_by(YEAR, ANALYSIS_STRATUM, STRAT, cover_group) %>%
      dplyr::summarise(
        avcvr = mean(Percent_Cvr, na.rm = TRUE),
        svar = var(Percent_Cvr, na.rm = TRUE),
        n = length(unique(PRIMARY_SAMPLE_UNIT)),
        MIN_DEPTH = mean(MIN_DEPTH, na.rm = TRUE),
        MAX_DEPTH = mean(MAX_DEPTH, na.rm = TRUE),
        DEPTH_M = (MIN_DEPTH + MAX_DEPTH) / 2
      ) %>%
      # Handle 0 variance values
      dplyr::mutate(
        svar = ifelse(svar == 0, 0.00000001, svar),
        Var = svar / n,
        std = sqrt(svar),
        SE = sqrt(Var),
        CV_perc = (SE / avcvr) * 100
      )
    
    # Merge with NTOT data
    cover_est <- cover_est %>%
      dplyr::full_join(ntot) %>%
      # Calculate weighted estimates and handle missing values
      dplyr::mutate(
        whavcvr = wh * avcvr,
        whsvar = wh^2 * Var,
        n = tidyr::replace_na(n, 0),
        PROT = ifelse(region_is_FL, PROT, NA)  # Set PROT to NA for non-FL regions
      ) %>%
      dplyr::filter(cover_group != "NA")
    
    return(cover_est)
  }
  
  #####Process Cover Data (region dependent)####
  
  # Process cover data for FL regions
  if (region %in% FL) {
    cover_est <- process_cover_data(inputdata, region_is_FL = TRUE)
  }
  #Process cover data for GOM and Caribbean regions
  else if (region %in% GOM | region %in% Carib) {
    cover_est <- process_cover_data(inputdata, region_is_FL = FALSE)
  }
  
  #####Create strata means####
  cover_strata <- cover_est %>%
    dplyr::select(REGION, YEAR, ANALYSIS_STRATUM, STRAT, PROT, DEPTH_M, cover_group, n, avcvr, Var, SE, CV_perc) %>%
    dplyr::mutate(n = tidyr::replace_na(n, 0)) %>%
    dplyr::mutate(CV_perc = ifelse(CV_perc == Inf, NA_real_, CV_perc))
  
  ####Create domain estimates####
  Domain_est <- cover_est %>%
    dplyr::group_by(REGION, YEAR, cover_group) %>%
    dplyr::summarise(
      avCvr = sum(whavcvr, na.rm = TRUE),
      Var = sum(whsvar, na.rm = TRUE),
      SE = sqrt(Var),
      CV_perc = (SE / avCvr) * 100,
      n_sites = sum(n),
      n_strat = length(unique(ANALYSIS_STRATUM)),
      ngrtot = sum(NTOT)
    ) %>%
    dplyr::ungroup()
  
  #####Export results####
  output <- list(
    "cover_strata" = cover_strata,
    "Domain_est" = Domain_est
  )
  
  return(output)
}
