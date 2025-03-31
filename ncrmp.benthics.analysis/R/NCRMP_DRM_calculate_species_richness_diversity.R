## Function to create a complete species list, calculate species richness and species diversity for combined NCRMP and DRM data (FL only)

# Purpose:
# creates csv files with species list, species richness and species diversity


## Tag: data analysis


# outputs created in this file --------------
# species_list
# richness_site
# unwh_richness_strata
# Domain_est
# species_diversity_site
# species_diversity_strata


# CallS:
# analysis ready data

# output gets called by:
# Analysis Rmarkdown, etc.
#

# NCRMP Caribbean Benthic analytics team: Viehman and Groves
# Last update: Jan 2025


##############################################################################################################################

#' Creates species list, species richness and species diversity dataframes from NCRMP and DRM data
#'
#' Creates data summaries of species richness and diversity, based on NCRMP
#' coral demographic data.
#' Species richness includes juveniles. Diversity is based only on adults.
#' Note that there is no accounting for transect length at this point, for example,
#' we do not calculate # of species/m^2. However, transect length may vary.
#' For 2-stage data (Florida, NCRMP+DRM), richness and diversity are averaged
#' between the 2 stations. No coral sites are not included in these estimates of
#' richness and no adult coral sites are not included in the estimates of diversity.
#'
#'
#'
#'
#' @param project A string indicating the project, NCRMP ("NCRMP") or NCRMP and DRM combined ("NCRMP_DRM").
#' @param region A string indicating the region. Options are: "SEFCRI", "FLK", "Tortugas", "STTSTJ", "STX, "FGB" or "PRICO".
#' @return A list dataframes
#' @importFrom magrittr "%>%"
#' @importFrom vegan "diversity"
#' @export
#'
#'
NCRMP_DRM_calculate_species_richness_diversity <- function(project, region){
  
  #### Prep Data ####
  
  # Define regional groups
  FL <- c("SEFCRI", "FLK", "Tortugas")
  FGB <- "FGB"
  Carib <- c("STTSTJ", "STX", "PRICO")
  
  # Load demo data
  demo_data <- load_NCRMP_DRM_demo_data(project = project, region = region)
  
  # Unpack data
  list2env(demo_data, envir = environment())
  
  #### Recode and Clean Species ####
  recode_and_clean_species <- function(data){
    data %>%
      dplyr::mutate(SPECIES_NAME = dplyr::case_when(
        SPECIES_CD == "MEAN JACK" ~ "Meandrina jacksoni",
        SPECIES_CD == "DIP STRI" ~ "Pseudodiploria strigosa",
        SPECIES_CD == "CLA ARBU" ~ "Cladacora arbuscula",
        SPECIES_CD == "CLA ABRU" ~ "Cladacora arbuscula",
        SPECIES_CD == "DIP CLIV" ~ "Pseudodiploria clivosa",
        SPECIES_CD == "PSE CLIV" ~ "Pseudodiploria clivosa",
        SPECIES_CD == "PSE STRI" ~ "Pseudodiploria clivosa",
        TRUE ~ as.character(SPECIES_NAME))) %>%
      dplyr::mutate(SPECIES_CD = dplyr::case_when(
        SPECIES_NAME == "Pseudodiploria strigosa" ~ "PSE STRI",
        SPECIES_NAME == "Pseudodiploria clivosa" ~ "PSE CLIV",
        SPECIES_NAME == "Meandrina jacksoni" ~ "MEA JACK",
        TRUE ~ as.character(SPECIES_CD)))
  }
  
  #### Format Geospatial Data ####
  format_geospatial_data <- function(data){
    data %>%
      dplyr::mutate(
        LAT_DEGREES = sprintf("%0.4f", LAT_DEGREES),
        LON_DEGREES = sprintf("%0.4f", LON_DEGREES),
        PROTECTION_STATUS = as.factor(PROT)
      )
  }
  
  # Apply formatting
  dat_1stage <- format_geospatial_data(dat_1stage)
  if(length(demo_data) > 1){
    dat_2stage <- format_geospatial_data(dat_2stage)
  }
  
  #### Bind Data Rows ####
  combined_data <- switch(region,
                          "SEFCRI" = dplyr::bind_rows(dat_1stage, dat_2stage),
                          "FLK" = if (project == "NCRMP_DRM") dplyr::bind_rows(dat_1stage, dat_2stage) else dat_1stage,
                          "Tortugas" = dplyr::bind_rows(dat_1stage, dat_2stage),
                          "STTSTJ" = dat_1stage,
                          "STX"= dat_1stage,
                          "PRICO" = dat_1stage,
                          "FGB" = dat_1stage,
                          NULL
  )
  
  # Recode species
  combined_data <- recode_and_clean_species(combined_data)
  dat_1stage <- recode_and_clean_species(dat_1stage)
  if(project %in% c("NCRMP_DRM", "NCRMP") && region %in% c("SEFCRI", "Tortugas")){
    dat_2stage <- recode_and_clean_species(dat_2stage)
  }
  
  #### Species Richness by Site ####
  filter_select_species <- function(data){
    data %>%
      dplyr::filter(N == 1, !grepl('SPE.', SPECIES_CD), !grepl('ANCX', SPECIES_CD), SPECIES_CD != "OTH CORA") %>%
      dplyr::select(SPECIES_NAME, SPECIES_CD) %>%
      dplyr::distinct()
  }
  
  summarize_species_richness <- function(data){
    data %>%
      dplyr::mutate(PROTECTION_STATUS = as.factor(PROT)) %>%
      dplyr::group_by(REGION, YEAR, PRIMARY_SAMPLE_UNIT, LAT_DEGREES, LON_DEGREES, STRAT, HABITAT_CD, PROTECTION_STATUS, SPECIES_NAME) %>%
      dplyr::summarise(Species_Count = sum(N)) %>%
      dplyr::mutate(Present = 1) %>%
      dplyr::group_by(REGION, YEAR, PRIMARY_SAMPLE_UNIT, LAT_DEGREES, LON_DEGREES, STRAT, HABITAT_CD, PROTECTION_STATUS) %>%
      dplyr::summarise(Species_Richness = sum(Present)) %>%
      dplyr::ungroup()
  }
  
  # Process species data
  if(project %in% c("NCRMP_DRM", "NCRMP") && region %in% c("SEFCRI", "Tortugas")){
    species_1stage <- filter_select_species(dat_1stage)
    species_2stage <- filter_select_species(dat_2stage)
    species_list <- dplyr::bind_rows(species_1stage, species_2stage) %>% unique()
  } else {
    species_list <- filter_select_species(dat_1stage)
  }
  
  #### Calculate Coral Diversity ####
  clean_data <- function(data){
    data %>%
      dplyr::filter(N == 1, JUV == 0, !grepl('SPE.', SPECIES_CD), !grepl('ANCX', SPECIES_CD), SPECIES_CD != "OTH CORA") %>%
      dplyr::mutate(PRIMARY_SAMPLE_UNIT = as.character(PRIMARY_SAMPLE_UNIT))
  }
  
  # Compute species diversity indices
  calculate_species_diversity <- function(species_data){
    species_data %>%
      dplyr::mutate(
        Simpson = vegan::diversity(species_data[, -c(1:2)], index = "simpson"),
        Inverse_Simpson = vegan::diversity(species_data[, -c(1:2)], index = "invsimpson"),
        Shannon = vegan::diversity(species_data[, -c(1:2)], index = "shannon")
      ) %>%
      dplyr::select(YEAR, PRIMARY_SAMPLE_UNIT, Simpson, Inverse_Simpson, Shannon)
  }
  
  return(list(
    species_richness = summarize_species_richness(combined_data),
    species_diversity = calculate_species_diversity(combined_data)
  ))
}
