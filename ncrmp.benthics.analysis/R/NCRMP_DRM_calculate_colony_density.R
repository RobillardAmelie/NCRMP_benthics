## Function to calculate colony density for combined NCRMP and DRM  data

# Purpose:
# creates csv files with colony density.


## Tag: data analysis


# outputs created in this file --------------
# density_site
# density_strata,
# Domain_est


# CallS:
# analysis ready data

# output gets called by:
# Analysis Rmarkdown, etc.
#

# NCRMP Caribbean Benthic analytics team: Groves, Viehman, Williams
# Last update: Jul 2025


##############################################################################################################################

#' Creates colony density summary dataframes
#'
#' Calculates coral density from NCRMP coral demographic data at the site level,
#' strata mean level, and regional weighted mean level (all across species), as well
#' site at the species level.  NCRMP utilizes a stratified random sampling design.
#' Regional estimates of coral density are weighted by the number of
#' grid cells of a stratum in the sample frame.
#'
#'
#'
#'
#'
#' @param project A string indicating the project, NCRMP or NCRMP and DRM combined ("NCRMP_DRM").
#' @param region A string indicating the region.  Options are: "SEFCRI", "FLK", "Tortugas", "STX", "STTSTJ", "PRICO", and "FGB".
#' @param species_filter A string indicating whether to filter to a subset of species.
#' @return A list of dataframes includeind 1) coral density by species at each site,
#' 2) total coral density at each site, 3) mean total coral density for each strata,
#' and 4) weighted mean regional total coral density.
#' @importFrom magrittr "%>%"
#' @export
#'
#'
#'

NCRMP_DRM_calculate_colony_density <- function(project = "NULL", region, species_filter = "NULL") {
  
  ####Load and preprocess data ####
  tmp <- load_NCRMP_DRM_demo_data(project = project, region = region, species_filter = species_filter) 
  list2env(tmp, envir = environment())

  ####Helper function to make the dataframe format####
  
  #vars: this is a string of col names that serve as the id cols 
  pivot_function <- function( data, vars, spp) {
    data %>%
      #format the data wider, so each unique specied_cd is a column, with the density as the value for each PSU
      #ex: one col for fake PSU 1 would include all possible species, and the density values would be just for the species obs
      tidyr::pivot_wider(id_cols = all_of(vars),
                         names_from = SPECIES_CD, values_from = DENSITY) %>%
      #we turn it back to the longer format, which introduces NA values were a species
      tidyr::pivot_longer(cols = spp$SPECIES_CD, names_to = "SPECIES_CD", values_to = "DENSITY") %>%
      #turn the NA values into 0, representing no density (the coral wasn't seen at the PSU in that year etc)
      dplyr::mutate(DENSITY = case_when(is.na(DENSITY) ~ 0, TRUE ~ DENSITY)) %>%
      #grab the species full names
      dplyr::left_join(spp, by = "SPECIES_CD")
  }

   
  #### Helper Function: calculate abundance and density ####
  calc_dens <- function(data, transect_2 = FALSE){
    abund <- data %>% dplyr::summarise(ABUNDANCE = sum(N), .groups = "keep")
    
    if(transect_2 == TRUE){
      dens <- abund %>% dplyr::mutate(DENSITY_transect = ABUNDANCE/METERS_COMPLETED)
    } else{
      dens <- abund %>%
      dplyr::mutate(DENSITY = ABUNDANCE/METERS_COMPLETED) %>%
      dplyr::select(-ABUNDANCE, -METERS_COMPLETED)
    }
    return(dens)
  }
  
  #### to get distinct species and species names ####
  spp_helper <- function(spp_data){
    spp_data %>%
      dplyr::select(SPECIES_CD, SPECIES_NAME) %>%
      dplyr::distinct(.)
  }
  
  ####Main Helper Function: Process dens data####
  process_density_data <- function(data, groups, transect_2 = FALSE){
    data %>%
      dplyr::filter(!(SUB_REGION_NAME %in% c("Marquesas", "Marquesas-Tortugas Trans")), JUV == 0) %>%
      dplyr::mutate(PROT = as.factor(PROT),
                    LAT_DEGREES = sprintf("%0.4f", LAT_DEGREES),
                    LON_DEGREES = sprintf("%0.4f", LON_DEGREES)) %>% 
      dplyr::mutate(PRIMARY_SAMPLE_UNIT = as.factor(PRIMARY_SAMPLE_UNIT)) %>%
      # dplyr::group_by(REGION, SURVEY, YEAR, SUB_REGION_NAME, ADMIN, PRIMARY_SAMPLE_UNIT, LAT_DEGREES, LON_DEGREES, STRAT, HABITAT_CD, PROT, METERS_COMPLETED, SPECIES_CD, SPECIES_NAME) %>%
      dplyr::group_by(across(all_of(groups)))%>%
      calc_dens(transect_2) %>%
      dplyr::ungroup()
  }
  
  if (project == "NCRMP_DRM" || (project == "NCRMP" && (region == "SEFCRI" || region == "Tortugas"))) {
      

      dat1_1stage <-process_density_data(dat_1stage, c("REGION", "SURVEY", "YEAR", "SUB_REGION_NAME", "ADMIN",   "PRIMARY_SAMPLE_UNIT", "LAT_DEGREES", "LON_DEGREES", 
                                                       "STRAT", "HABITAT_CD", "MIN_DEPTH", "MAX_DEPTH", "PROT", "METERS_COMPLETED"))

      density_species_1stage <-process_density_data(dat_1stage, groups =   c("REGION", "SURVEY", "YEAR", "SUB_REGION_NAME", "ADMIN", "PRIMARY_SAMPLE_UNIT", "LAT_DEGREES", 
                                                         "LON_DEGREES", "STRAT", "HABITAT_CD",  "PROT", "METERS_COMPLETED", "SPECIES_CD", "SPECIES_NAME"))
      
      dat1_2stage <- process_density_data(dat_2stage, groups =  c("REGION", "SURVEY", "YEAR", "SUB_REGION_NAME", "ADMIN", "PRIMARY_SAMPLE_UNIT", "STATION_NR", "LAT_DEGREES", 
                                                          "LON_DEGREES", "STRAT", "HABITAT_CD", "MIN_DEPTH", "MAX_DEPTH", "PROT", "METERS_COMPLETED"), transect_2 = TRUE) %>%
        dplyr::group_by(REGION, SURVEY, YEAR, SUB_REGION_NAME, ADMIN, PRIMARY_SAMPLE_UNIT, LAT_DEGREES, LON_DEGREES, STRAT, HABITAT_CD, PROT) %>%
        dplyr::summarise(DENSITY = mean(DENSITY_transect),
                         MIN_DEPTH = mean(MIN_DEPTH),
                         MAX_DEPTH = mean(MAX_DEPTH), .groups = "keep") %>% dplyr::ungroup() 
      
      density_species_2stage <- process_density_data(dat_2stage, groups =    c("REGION", "SURVEY", "YEAR", "SUB_REGION_NAME", "ADMIN", "PRIMARY_SAMPLE_UNIT", "STATION_NR", "LAT_DEGREES", 
                                                                     "LON_DEGREES", "STRAT", "HABITAT_CD",  "PROT", "METERS_COMPLETED", "SPECIES_CD", "SPECIES_NAME"), transect_2 = TRUE) %>%
        dplyr::group_by(REGION, SURVEY, YEAR, SUB_REGION_NAME, ADMIN, PRIMARY_SAMPLE_UNIT, LAT_DEGREES, LON_DEGREES, STRAT, HABITAT_CD, PROT, SPECIES_CD, SPECIES_NAME) %>%
        dplyr::summarise(DENSITY = mean(DENSITY_transect),
                         ABUNDANCE = sum(ABUNDANCE), .groups = "keep") %>% dplyr::ungroup() 
      
      density_site <- dplyr::bind_rows(dat1_1stage, dat1_2stage)

      density_species <- dplyr::bind_rows(density_species_1stage, density_species_2stage)

      density_species <- pivot_function(
        data = density_species,
        vars = c("REGION", "SURVEY", "YEAR", "SUB_REGION_NAME", "ADMIN", "PRIMARY_SAMPLE_UNIT", "LAT_DEGREES", "LON_DEGREES", "STRAT", "HABITAT_CD", "PROT"),
        spp = spp_helper(density_species)) %>%
        dplyr::select(REGION, SURVEY, YEAR, SUB_REGION_NAME, ADMIN, PRIMARY_SAMPLE_UNIT, LAT_DEGREES, LON_DEGREES, STRAT, HABITAT_CD, PROT, SPECIES_CD, SPECIES_NAME, DENSITY)

    } else {
      
      density_species <-process_density_data(dat_1stage, groups =    c("REGION", "SURVEY", "YEAR", "SUB_REGION_NAME", "ADMIN", "PRIMARY_SAMPLE_UNIT", "LAT_DEGREES", 
                                                         "LON_DEGREES", "STRAT", "HABITAT_CD", "MIN_DEPTH", "MAX_DEPTH", "PROT", "METERS_COMPLETED", "SPECIES_CD", "SPECIES_NAME"))
    
      #call the pivot helper function to ensure that there is '0 density' placed for corals that were not observed
      density_species <- pivot_function(
        data = density_species,
        vars = c("REGION", "SURVEY", "YEAR", "SUB_REGION_NAME", "ADMIN", "PRIMARY_SAMPLE_UNIT", "LAT_DEGREES", "LON_DEGREES", "STRAT", "HABITAT_CD", "PROT", "MIN_DEPTH", "MAX_DEPTH"),
        spp = spp_helper(density_species)) %>%
        dplyr::select(REGION, SURVEY, YEAR, SUB_REGION_NAME, ADMIN, PRIMARY_SAMPLE_UNIT, LAT_DEGREES, LON_DEGREES, STRAT, HABITAT_CD, MIN_DEPTH, MAX_DEPTH,PROT, SPECIES_CD, SPECIES_NAME, DENSITY)

      #calc density by site
      density_site <-process_density_data(dat_1stage, groups =    c("REGION", "SURVEY", "YEAR", "SUB_REGION_NAME", "ADMIN", "PRIMARY_SAMPLE_UNIT", "LAT_DEGREES", 
                                                             "LON_DEGREES", "STRAT", "HABITAT_CD",  "PROT", "METERS_COMPLETED"))
    }

  #### Run weighting function ####
  tmp <- NCRMP_make_weighted_demo_data(project, inputdata = density_site, datatype = 'density', region, species_filter = species_filter)
  list2env(tmp, envir = environment())


  ####Export####
  output <- list(
    "density_species" = density_species,
    "density_site" = density_site,
    "density_strata" = density_strata,
    "Domain_est" = Domain_est
  )
  
  if (project == "MIR") {
    output$Domain_est_PROT <- Domain_est_PROT
  }
  return(output)
}
