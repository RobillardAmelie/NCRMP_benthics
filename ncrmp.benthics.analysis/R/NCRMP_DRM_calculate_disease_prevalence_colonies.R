## Function to calculate disease prevalence & bleaching prevalence for NCRMP, MIR, and NCRMP + DRM data (FL only) by calculating species then colony prevalence (%) at the site level,
## taking the mean of all sites within each strata, strata area weighting each strata and summing all strata means to reach the domain estimate.

# Purpose:
# creates csv files with disease/bleaching prevalence by region


## Tag: data analysis


# outputs created in this file --------------
# disease_prev_site
# disease_prev_strata
# Domain estimates

# CallS:
# analysis ready data

# output gets called by:
# NCRMP_DRM_calculate_dis_ble_prevalence_species_domain.R
#

# NCRMP Caribbean Benthic analytics team: Groves, Viehman, Williams
# Last update: Feb 2025


##############################################################################################################################

#' Calculate disease prevalence & bleaching prevalence at the species/site, site, strata and domain levels
#'
#' Calculates disease and bleaching prevalence for each species at each site,
#' at each site across all species, at each strata across all species, and
#' regional estimates for each year of a given region.
#' NCRMP utilizes a stratified random sampling design.
#' Regional estimates of disease and bleaching prevalence are weighted by the number of
#' grid cells of a stratum in the sample frame. Species-level outputs from
#' this function are utilized by [NCRMP_DRM_calculate_dis_ble_prevalence_species_domain()].
#'
#'
#'
#'
#' @param project A string indicating the project, NCRMP, MIR, or NCRMP and DRM combined ("NCRMP_DRM").
#' @param region A string indicating the region. Options are: "SEFCRI", "FLK", "Tortugas", "STX", "STTSTJ", "PRICO", and "FGB".
#' @param species_filter An optional concatenated string indicating whether to filter to a subset of species
#' @return A list of dataframes including 1) bleaching and disease prevalence by species
#' and site, 2) bleaching and disease prevalence by site, 3) disease prevalence by
#' strata, 4) bleaching prevalence by strata, and 5) regional estimates for disease
#' and bleaching prevalence.
#' @importFrom magrittr "%>%"
#' @export
#'
#'


NCRMP_DRM_calculate_disease_prevalence_colonies <- function(project, region, species_filter = "NULL"){

  ####prep data####
  # Load data
  tmp <- load_NCRMP_DRM_demo_data(project = project, region = region, species_filter = species_filter)
  #unpack list
  list2env(tmp, envir = environment())

  ####Helper Function: clean data, no disease na####
  FL_dis_NA <- function(data){
    data %>% dplyr::filter(N == 1,
                           DISEASE != "N/A",
                           JUV == 0,
                           SUB_REGION_NAME != "Marquesas",
                           SUB_REGION_NAME != "Marquesas-Tortugas Trans")
  }

 ####Helper Function: Mutate Formatting
  mutate_formatting <- function(data){
    data %>%
      dplyr::mutate(PRIMARY_SAMPLE_UNIT = as.factor(PRIMARY_SAMPLE_UNIT),
                    LAT_DEGREES = sprintf("%0.4f", as.numeric(LAT_DEGREES)), 
                    LON_DEGREES = sprintf("%0.4f", as.numeric(LON_DEGREES)))
  }


  ####Helper Function: process data####
  process_data <- function(data, include_PL = FALSE) {
    
    bleach_conditions <- if (include_PL == TRUE) {
      c("P", "T", "B", "PB", "PL")  
    } else {
      c("P", "T", "B", "PB")  
    }
    
    data %>%
      dplyr::mutate(PROT = as.factor(PROT),
                    DISEASE = as.integer(DISEASE %in% c("P", "F", "S")),
                    BLEACH = as.integer(BLEACH_CONDITION %in% bleach_conditions)) %>%
      #call mutate formatting helper function
      mutate_formatting()
  }

  
  ####Helper Function: Sum ble/dis####
  sum_disease_bleaching <- function(data){
    data%>% dplyr::summarise(Total_dis = sum(DISEASE),
                       Total_ble = sum(BLEACH),
                       Total_col = sum(N),
                       DIS_PREV = (Total_dis/Total_col)*100,
                       BLE_PREV = (Total_ble/Total_col)*100, .groups = "keep")
  }

  ####Helper Function: mean ble/dis####
  mean_disease_bleaching <- function(data){
    data %>% dplyr::summarise(Total_dis = mean(Total_dis),
                       Total_ble = mean(Total_ble),
                       Total_col = mean(Total_col),
                       DIS_PREV = mean(DIS_PREV),
                       BLE_PREV = mean(BLE_PREV), .groups = "keep")
  }
  ####Helper Function: Format ble/dis####
  format_dis_ble_prev <- function(data) {
    data %>% dplyr::mutate(DIS_PREV = as.numeric(sprintf("%0.1f", DIS_PREV)),
                          BLE_PREV = as.numeric(sprintf("%0.1f", BLE_PREV)))
  }


  # Calulate site level disease prevalence

  if (project == "NCRMP_DRM" || (project == "NCRMP" && (region == "SEFCRI" || region == "Tortugas"))) {

    dat1_1stage <- dat_1stage %>%
      FL_dis_NA() %>%
      dplyr::filter(!(is.na(DISEASE))) %>%
      mutate_formatting() %>%
      mutate(DATE = paste(MONTH, DAY, YEAR, sep = "/")) %>%
      process_data(include_PL = FALSE) %>%
      dplyr::group_by(SURVEY, REGION, YEAR, SUB_REGION_NAME, PRIMARY_SAMPLE_UNIT, LAT_DEGREES, LON_DEGREES, STRAT, HABITAT_CD, PROT) %>%
      sum_disease_bleaching() %>%
      format_dis_ble_prev()


    dis_species_1stage <- dat_1stage %>%
      FL_dis_NA() %>%
      dplyr::filter(!(is.na(DISEASE))) %>%
      mutate_formatting() %>%
      mutate(DATE = paste(MONTH, DAY, YEAR, sep = "/")) %>%
      process_data(include_PL = FALSE) %>%
      dplyr::group_by(SURVEY, REGION, YEAR, SUB_REGION_NAME, PRIMARY_SAMPLE_UNIT, LAT_DEGREES, LON_DEGREES, STRAT, HABITAT_CD, PROT, SPECIES_CD) %>%
      sum_disease_bleaching() %>%
      format_dis_ble_prev()


    dat1_2stage <- dat_2stage %>%
      FL_dis_NA() %>%
      process_data(include_PL = TRUE) %>%
      dplyr::group_by(SURVEY, REGION, YEAR, SUB_REGION_NAME, PRIMARY_SAMPLE_UNIT, STATION_NR, LAT_DEGREES, LON_DEGREES, STRAT, HABITAT_CD, PROT) %>%
      sum_disease_bleaching() %>%
      dplyr::group_by(SURVEY, REGION, YEAR, SUB_REGION_NAME, PRIMARY_SAMPLE_UNIT, LAT_DEGREES, LON_DEGREES, STRAT, HABITAT_CD, PROT) %>%
      mean_disease_bleaching() %>%
      format_dis_ble_prev()

    dis_species_2stage <- dat_2stage %>%
      FL_dis_NA() %>%
      mutate_formatting() %>%
      process_data(include_PL = FALSE) %>%
      dplyr::group_by(SURVEY, REGION, YEAR, SUB_REGION_NAME, PRIMARY_SAMPLE_UNIT, STATION_NR, LAT_DEGREES, LON_DEGREES, STRAT, HABITAT_CD, PROT, SPECIES_CD) %>%
      sum_disease_bleaching() %>%
      dplyr::group_by(SURVEY, REGION, YEAR, SUB_REGION_NAME, PRIMARY_SAMPLE_UNIT, LAT_DEGREES, LON_DEGREES, STRAT, HABITAT_CD, PROT, SPECIES_CD) %>%
      mean_disease_bleaching() %>%
      format_dis_ble_prev()

    disease_prev_species <-dplyr::bind_rows(dis_species_1stage, dis_species_2stage)

    disease_prev_site <-dplyr::bind_rows(dat1_1stage, dat1_2stage)

  } else {



    disease_prev_site <- dat_1stage %>%
      dplyr::filter(N == 1,
                    #DISEASE != "N/A", # keep NA's because disease info wasn't consistently collected in 2013
                    JUV == 0) %>%
      process_data(include_PL = FALSE) %>%
      dplyr::group_by(SURVEY, REGION, YEAR, SUB_REGION_NAME, PRIMARY_SAMPLE_UNIT, LAT_DEGREES, LON_DEGREES, STRAT, HABITAT_CD, PROT) %>%
      sum_disease_bleaching() %>%
      format_dis_ble_prev()

    disease_prev_species <- dat_1stage %>%
      dplyr::filter(N == 1,
                    #DISEASE != "N/A", # keep NA's because disease info wasn't consistently collected in 2013
                    JUV == 0) %>%
      process_data(include_PL = FALSE) %>%
      dplyr::group_by(SURVEY, REGION, YEAR, SUB_REGION_NAME, PRIMARY_SAMPLE_UNIT, LAT_DEGREES, LON_DEGREES, STRAT, HABITAT_CD, PROT, SPECIES_CD) %>%
      sum_disease_bleaching() %>%
      format_dis_ble_prev()

  }

  # Run through the weighting function
  tmp  <- NCRMP_make_weighted_demo_data(project, inputdata = disease_prev_site, region, datatype = "disease")
  # unpack list
  list2env(tmp, envir = environment())

  ####Export####
  output <- list(
    'dis_ble_prev_species' = disease_prev_species,
    "dis_ble_prev_site" = disease_prev_site,
    "dis_prev_strata" = dis_prev_strata,
    'ble_prev_strata' = ble_prev_strata,
    "Domain_est" = Domain_est)

  return(output)
}
