## Function to calculate weighted invert density by strata

# Purpose:
# support function to calculate weighted invert density data


## Tag: data analysis


# outputs created in this file --------------
# invert_density_site
# unwh_invert_strata
# Domain est

# Current weighting scheme:
# STRAT + PROT (FL)
# STRAT (Carib/GOM)

# CallS:
# analysis ready data

# output gets called by:
# NCRMP_calculate_invert_density.R
#

# NCRMP Caribbean Benthic analytics team: Groves, Viehman
# Last update: Feb 2025


##############################################################################################################################

#' Creates weighted invertebrate density data
#'
#' Calculates weighted invertebrate density data collected by benthic assessment.
#' NCRMP utilizes a stratified random sampling design.
#' Regional estimates of invertebrate density are weighted by the
#' number of grid cells of a stratum in the sample frame.
#' Function produces strata means, weighted strata means,
#' and weighted regional estimates for invertebrate data.
#' SUpport function called by [NCRMP_calculate_invert_density()].
#'
#'
#'
#'
#'
#' @param inputdata A dataframe of site level invertebrate density.
#' @param region A string indicating the region. Options are: "SEFCRI", "FLK", "Tortugas", "STX", "STTSTJ", "PRICO", and "GOM".
#' @param project A string indicating the project: "NCRMP" or "MIR". Default is NCRMP.
#' @return A list of dataframes, including a dataframe of strata means of invertebrate
#' and a dataframe of weighted regional estimates of invertebrate density, for specified region.
#' @importFrom magrittr "%>%"
#' @export
#'
#'

NCRMP_make_weighted_invert_density_data <- function(inputdata, region, project = "NULL") {

  #### Define Regional Groups ####
  FL <- c("SEFCRI", "FLK", "Tortugas")
  GOM <- "GOM"
  Carib <- c("STTSTJ", "STX", "PRICO")


  #### NTOT Data ####
  ntot <- load_NTOT(region = region, inputdata = inputdata, project = project) %>%
    dplyr::mutate(
      wh = NTOT / ngrtot,
      PROT = as.factor(PROT))


  #### Density Calc Helper Function ####
  density_calculator <- function(data){
    data <- data %>%
      dplyr::summarise(
        avden = mean(Diadema_dens), # calculate mean density
        svar = var(Diadema_dens), # calculate stratum variance
        n = length(Diadema_dens)) %>%
      # convert 0 for stratum variance so that the sqrt is a small num but not a 0
      dplyr::mutate(
        svar = dplyr::if_else(svar == 0, 1e-8, svar),
        std = sqrt(svar)
      )
  }

  #### Stratum Est Helper Function ####
  stratum_estimates <- function(data){
    data <- data %>%
      dplyr::mutate(whavden = wh * avden,
                    whsvar = wh^2 * svar,
                    whstd = wh * std,
                    n = tidyr::replace_na(n, 0),
                    # Add the following to match FL format for non FL regions
                    PROT = NA,
                    RUG_CD = NA)
  }


  if(region %in% FL) {

    #### Calculate avdns, svar, n and std at the strata + PROT level ####
    dens_est <- inputdata %>%
      dplyr::mutate(ANALYSIS_STRATUM = paste(STRAT, "/ PROT =", PROT, sep = " ")) %>%
      dplyr::group_by(YEAR, ANALYSIS_STRATUM, STRAT, PROT) %>% # group by analysis level strata
      density_calculator() %>%
      dplyr::full_join(., ntot) %>% # Merge ntot with coral_est_spp
      # stratum estimates
      stratum_estimates
  }


  if(region %in% GOM | region %in% Carib) {

    #### Calculate avdns, svar, n and std at the strata ####
    dens_est <- inputdata %>%
      dplyr::mutate(ANALYSIS_STRATUM = STRAT) %>%
      # group by analysis level strata
      dplyr::group_by(YEAR, ANALYSIS_STRATUM, STRAT) %>% # Modify this line to change analysis stratum
      density_calculator() %>%
      # Merge ntot with coral_est_spp
      dplyr::full_join(., ntot) %>%
      # stratum estimates
      stratum_estimates()
  }

  # Reformat output
  dens_est <- dens_est %>%
    dplyr::select(REGION, YEAR, ANALYSIS_STRATUM, STRAT, PROT, NTOT, ngrtot, wh, n, avden, svar, std, whavden, whsvar, whstd)

  # cover, unweighted by strata
  invert_strata <- dens_est %>%
    dplyr::select(REGION, YEAR, ANALYSIS_STRATUM, STRAT, PROT, NTOT, ngrtot, wh, n, avden, svar, std)


  ## Domain Estimates
  Domain_est <- dens_est %>%
    dplyr::group_by(REGION, YEAR) %>%
    dplyr::summarise(avDen = sum(whavden, na.rm = T),
                     var = sum(whsvar, na.rm = T),
                     std = sqrt(var),
                     ngrtot = sum(NTOT) )

  ####Export####
  output <- list(
    "invert_strata" = invert_strata,
    "Domain_est" = Domain_est

  )
}
