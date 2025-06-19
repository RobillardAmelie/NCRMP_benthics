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

NCRMP_DRM_calculate_species_richness_diversity <- function(project, region) {
  project <- "NCRMP_DRM"
  region <- "FLK"
  
  #### Helper Functions ####
  
  # Clean and filter valid coral species
  filter_valid_species <- function(data, include_juveniles = TRUE) {
    data <- data %>%
      dplyr::filter(N == 1, !grepl("SPE\\.", SPECIES_CD), !grepl("ANCX", SPECIES_CD), SPECIES_CD != "OTH CORA")
    if (!include_juveniles) {
      data <- data %>% dplyr::filter(JUV == 0)
    }
    return(data)
  }
  
  # Recode species names/codes
  recode_and_clean_species <- function(data) {
    data %>%
      dplyr::mutate(SPECIES_NAME = dplyr::case_when(
        SPECIES_CD == "MEAN JACK" ~ "Meandrina jacksoni",
        SPECIES_CD %in% c("DIP STRI", "PSE STRI") ~ "Pseudodiploria strigosa",
        SPECIES_CD %in% c("DIP CLIV", "PSE CLIV") ~ "Pseudodiploria clivosa",
        SPECIES_CD %in% c("CLA ARBU", "CLA ABRU") ~ "Cladacora arbuscula",
        TRUE ~ as.character(SPECIES_NAME)
      )) %>%
      dplyr::mutate(SPECIES_CD = dplyr::case_when(
        SPECIES_NAME == "Pseudodiploria strigosa" ~ "PSE STRI",
        SPECIES_NAME == "Pseudodiploria clivosa" ~ "PSE CLIV",
        SPECIES_NAME == "Meandrina jacksoni" ~ "MEA JACK",
        TRUE ~ as.character(SPECIES_CD)
      ))
  }
  
  # Format spatial fields
  format_geospatial_data <- function(data) {
    data %>%
      dplyr::mutate(
        LAT_DEGREES = sprintf("%0.4f", LAT_DEGREES),
        LON_DEGREES = sprintf("%0.4f", LON_DEGREES),
        PROTECTION_STATUS = as.factor(PROT)
      )
  }
  
  # Species richness
  summarize_species_richness <- function(data) {
    data %>%
      dplyr::mutate(PROTECTION_STATUS = as.factor(PROT)) %>%
      dplyr::group_by(REGION, YEAR, PRIMARY_SAMPLE_UNIT, LAT_DEGREES, LON_DEGREES, STRAT, HABITAT_CD, PROTECTION_STATUS, SPECIES_NAME) %>%
      dplyr::summarise(Species_Count = sum(N), .groups = "drop_last") %>%
      dplyr::mutate(Present = 1) %>%
      dplyr::group_by(REGION, YEAR, PRIMARY_SAMPLE_UNIT, LAT_DEGREES, LON_DEGREES, STRAT, HABITAT_CD, PROTECTION_STATUS) %>%
      dplyr::summarise(Species_Richness = sum(Present), .groups = "drop")
  }
  
  # Diversity index calculation
  calculate_species_diversity <- function(data) {
    
    diversity_ready <- data %>%
      filter_valid_species(include_juveniles = FALSE) %>%
      dplyr::group_by(YEAR, PRIMARY_SAMPLE_UNIT, SPECIES_CD) %>%
      dplyr::summarise(Abundance = sum(N), .groups = "drop") %>%
      tidyr::pivot_wider(names_from = SPECIES_CD, values_from = Abundance, values_fill = 0)
    
    diversity_matrix <- diversity_ready %>% dplyr::select(-YEAR, -PRIMARY_SAMPLE_UNIT)
    
    diversity_ready %>%
      dplyr::mutate(
        Simpson = vegan::diversity(diversity_matrix, index = "simpson"),
        Inverse_Simpson = vegan::diversity(diversity_matrix, index = "invsimpson"),
        Shannon = vegan::diversity(diversity_matrix, index = "shannon")
      ) %>%
      dplyr::select(YEAR, PRIMARY_SAMPLE_UNIT, Simpson, Inverse_Simpson, Shannon)
  }
  
  #### Load and Format Data ####
  
  FL <- c("SEFCRI", "FLK", "Tortugas")
  
  demo_data <- load_NCRMP_DRM_demo_data(project = project, region = region)
  list2env(demo_data, envir = environment())
  
  dat_1stage <- format_geospatial_data(dat_1stage)
  if (length(demo_data) > 1) dat_2stage <- format_geospatial_data(dat_2stage)
  
  combined_data <- switch(region,
                          "SEFCRI" = dplyr::bind_rows(dat_1stage, dat_2stage),
                          "FLK" = if (project == "NCRMP_DRM") dplyr::bind_rows(dat_1stage, dat_2stage) else dat_1stage,
                          "Tortugas" = dplyr::bind_rows(dat_1stage, dat_2stage),
                          "STTSTJ" = dat_1stage,
                          "STX" = dat_1stage,
                          "PRICO" = dat_1stage,
                          "FGB" = dat_1stage
  )
  
  dat_1stage <- recode_and_clean_species(dat_1stage)
  if (exists("dat_2stage")) dat_2stage <- recode_and_clean_species(dat_2stage)
  combined_data <- recode_and_clean_species(combined_data)
  
  #### Final Outputs ####
  
  species_list <- if (exists("dat_2stage")) {
    dplyr::bind_rows(filter_valid_species(dat_1stage), filter_valid_species(dat_2stage)) %>% dplyr::distinct()
  } else {
    filter_valid_species(dat_1stage)
  }
  
  species_richness <- summarize_species_richness(combined_data)
  species_diversity <- calculate_species_diversity(combined_data)
  
  #return list
  return(list(
    species_list = species_list,
    species_richness = species_richness,
    species_diversity = species_diversity
  ))
}
