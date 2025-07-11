


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
# Last update: Jan 2024


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
#' @param region A string indicating the region. Options are: "SEFCRI", "FLK", "Tortugas", "STTSTJ", "STX, "GOM" or "PRICO".
#' @return A list dataframes
#' @importFrom magrittr "%>%"
#' @importFrom vegan "diversity"
#' @export
#'
#'

NCRMP_DRM_calculate_species_richness_diversity <- function(project, region){
  #call vegan library
  library(vegan)
  
  ####  Define regional groups   #### 
  FL <- c("SEFCRI", "FLK", "Tortugas")
  FGB <- "FGB"
  Carib <- c("STTSTJ", "STX", "PRICO")
  
  #### Load data   #### 
  tmp <- load_NCRMP_DRM_demo_data(project = project, region = region)
  list2env(tmp, envir = environment())
  
  #### Helper Function: Format Geospatial/Prot Data   #### 
  reformat_geospatial_data <- function(data){
    data %>%    dplyr::mutate(LAT_DEGREES = sprintf("%0.4f", LAT_DEGREES),
                              LON_DEGREES = sprintf("%0.4f", LON_DEGREES),
                              PROT = as.factor(PROT))
  }

  # apply helper function to format data
  dat_1stage <- dat_1stage %>%
    reformat_geospatial_data() 
  
  if(length(tmp) > 1){ dat_2stage <- dat_2stage %>%
      reformat_geospatial_data() 
  }
  
  #### Bind the Data ####
  dat <- switch(region,"SEFCRI" = dplyr::bind_rows(dat_1stage, dat_2stage),
                          "FLK" = if (project == "NCRMP_DRM") dplyr::bind_rows(dat_1stage, dat_2stage) else dat_1stage,
                          "Tortugas" = dplyr::bind_rows(dat_1stage, dat_2stage),
                          "STTSTJ" = dat_1stage,
                          "STX" = dat_1stage,
                          "PRICO" = dat_1stage,
                          "FGB" = dat_1stage
  )
  
  #### Helper Function: Recode species names/codes   #### 
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
  
  ####  Helper Function: Species richness   #### 
  summarize_SPP_RICHNESS <- function(data) {
    data %>%
      dplyr::mutate(PROT = as.factor(PROT)) %>%
      dplyr::group_by(REGION, YEAR, PRIMARY_SAMPLE_UNIT, LAT_DEGREES, LON_DEGREES, STRAT, HABITAT_CD, PROT, SPECIES_NAME) %>%
      dplyr::summarise(Species_Count = sum(N), .groups = "drop_last") %>%
      dplyr::mutate(Present = 1) %>%
      dplyr::group_by(REGION, YEAR, PRIMARY_SAMPLE_UNIT, LAT_DEGREES, LON_DEGREES, STRAT, HABITAT_CD, PROT) %>%
      dplyr::summarise(SPP_RICHNESS = sum(Present), .groups = "drop")
  }
  
  
  
  ####  Clean up species names   #### 
  dat <- dat %>%
    dplyr::filter(SPECIES_CD != "NO. SCLE") %>%
    recode_and_clean_species()
  
  dat_1stage <- dat_1stage %>%
    dplyr::filter(SPECIES_CD != "NO. SCLE") %>%
    recode_and_clean_species()
  
  
  if(project == "NCRMP_DRM" || project == "NCRMP" && region == "SEFCRI" ||   project == "NCRMP" && region == "Tortugas") {
    
    dat_2stage <- dat_2stage %>%
      dplyr::filter(SPECIES_CD != "NO. SCLE") %>%
      recode_and_clean_species()
  }
  
  #### Helper Function: Get Distinct Species ####
  get_distinct_species <- function(data){
    data %>%
      dplyr::filter(N == 1,
                    !grepl('SPE.', SPECIES_CD),
                    !grepl('ANCX', SPECIES_CD),
                    SPECIES_CD != "OTH CORA") %>% # Filter out SPE
      dplyr::select(SPECIES_NAME, SPECIES_CD) %>%
      dplyr::distinct(SPECIES_NAME, SPECIES_CD)
  }
  
  
  ####  Calculate species richness by site   #### 
  if(project == "NCRMP_DRM" || project == "NCRMP" && region == "SEFCRI" ||   project == "NCRMP" && region == "Tortugas") {
    # Clean up species
    species_1stage <- dat_1stage %>%
      get_distinct_species() 
    
    species_2stage <- dat_2stage %>%
      get_distinct_species() 
    
    species_list <- species_1stage %>%
      rbind(., species_2stage) %>%
      unique(.)
    
    dat1_1stage <- dat_1stage %>%
      dplyr::filter(N == 1,
                    SUB_REGION_NAME != "Marquesas",
                    SUB_REGION_NAME != "Marquesas-Tortugas Trans",
                    !grepl('SPE.', SPECIES_CD),
                    !grepl('ANCX', SPECIES_CD)) %>% # Filter out SPE
      dplyr::mutate(PROT = as.factor(PROT)) %>% # Change PROT to factor for ggplot will recognize it as a grouping variable
      dplyr::group_by(REGION, SURVEY, YEAR, SUB_REGION_NAME, PRIMARY_SAMPLE_UNIT, LAT_DEGREES, LON_DEGREES, STRAT, HABITAT_CD, PROT, SPECIES_NAME) %>% #No need to include region, will be added from ntot in wh. function
      dplyr::summarise(SppSumSite = sum(N)) %>%
      dplyr::mutate(present = 1) %>%
      dplyr::group_by(REGION, SURVEY, YEAR, SUB_REGION_NAME, PRIMARY_SAMPLE_UNIT, LAT_DEGREES, LON_DEGREES, STRAT, HABITAT_CD, PROT) %>%
      dplyr::summarise(SPP_RICHNESS = sum(present)) %>%
      dplyr::ungroup()
    
    
    dat1_2stage <- dat_2stage %>%
      dplyr::filter(N == 1,
                    !grepl('SPE.', SPECIES_CD),
                    !grepl('ANCX', SPECIES_CD)) %>% # Filter out SPE
      dplyr::mutate(PROT = as.factor(PROT)) %>% # Change PROT to factor for ggplot will recognize it as a grouping variable
      # calculate number of species at each station first
      dplyr::group_by(REGION, SURVEY, YEAR, SUB_REGION_NAME, PRIMARY_SAMPLE_UNIT, STATION_NR, LAT_DEGREES, LON_DEGREES, STRAT, HABITAT_CD, PROT) %>% #No need to include region, will be added from ntot in wh. function
      dplyr::summarise(SppSumSite = length(unique(SPECIES_NAME))) %>%
      # average for 2 stage data
      dplyr::group_by(REGION, SURVEY, YEAR, SUB_REGION_NAME, PRIMARY_SAMPLE_UNIT, LAT_DEGREES, LON_DEGREES, STRAT, HABITAT_CD, PROT) %>%
      dplyr::summarise(SPP_RICHNESS = mean(SppSumSite)) %>%
      dplyr::ungroup()
    
    dat1 <- rbind(dat1_1stage, dat1_2stage)
    
    
  } else {
    
    # Clean up species
    species_list <- dat_1stage %>%
      get_distinct_species() %>%
      unique(.)
    
    dat1_1stage <- dat_1stage %>%
      dplyr::filter(N == 1,
                    SUB_REGION_NAME != "Marquesas",
                    SUB_REGION_NAME != "Marquesas-Tortugas Trans",
                    !grepl('SPE.', SPECIES_CD),
                    !grepl('ANCX', SPECIES_CD),
                    SPECIES_CD != "OTH CORA") %>% # Filter out SPE
      dplyr::mutate(PROT = as.factor(PROT)) %>% # Change PROT to factor for ggplot will recognize it as a grouping variable
      dplyr::group_by(REGION, SURVEY, YEAR, SUB_REGION_NAME, PRIMARY_SAMPLE_UNIT, LAT_DEGREES, LON_DEGREES, STRAT, HABITAT_CD, PROT, SPECIES_NAME) %>% #No need to include region, will be added from ntot in wh. function
      dplyr::summarise(SppSumSite = sum(N)) %>%
      dplyr::mutate(present = 1) %>%
      dplyr::group_by(REGION, SURVEY, YEAR, SUB_REGION_NAME, PRIMARY_SAMPLE_UNIT, LAT_DEGREES, LON_DEGREES, STRAT, HABITAT_CD, PROT) %>%
      dplyr::summarise(SPP_RICHNESS = sum(present)) %>%
      dplyr::ungroup()
    
    dat1 <- dat1_1stage
    
  }
  
  richness_site <- dat1
  
  #### Run through the weighting function   #### 
  tmp  <- NCRMP_make_weighted_demo_data(project, inputdata = richness_site,  region,    datatype = "richness",  species_filter = "FALSE")
  list2env(tmp, envir = environment())

  
  ###### Calculate coral diversity with Simpson, Inv. Simpson and Shannon Indices #####
  
  
  if(project == "NCRMP_DRM" || project == "NCRMP" && region == "SEFCRI"||  project == "NCRMP" && region == "Tortugas"){
    
    # Site level -- when there is 2 stage data, this is station level
    sites <- dat %>%
      # Update in Jan. 2024 - exclude juveniles - we don't count these (only once/site), so including them in diversity metrics is improper
      dplyr::filter(N == 1,
                    JUV == 0,
                    SUB_REGION_NAME != "Marquesas",
                    SUB_REGION_NAME != "Marquesas-Tortugas Trans",
                    !grepl('SPE.', SPECIES_CD),
                    !grepl('ANCX', SPECIES_CD))  %>%
      dplyr::filter(SPECIES_CD != "OTH CORA") %>%
      dplyr::mutate(PRIMARY_SAMPLE_UNIT = as.character(PRIMARY_SAMPLE_UNIT)) %>%
      dplyr::select(REGION, YEAR, SUB_REGION_NAME, PRIMARY_SAMPLE_UNIT, STATION_NR, LAT_DEGREES, LON_DEGREES, STRAT, PROT) %>%
      dplyr::distinct(.)
    
    species_site <- dat %>%
      # Update in Jan. 2024 - exclude juveniles - we don't count these (only once/site), so including them in diversity metrics is improper
      dplyr::filter(N == 1,
                    JUV == 0,
                    SUB_REGION_NAME != "Marquesas",
                    SUB_REGION_NAME != "Marquesas-Tortugas Trans",
                    !grepl('SPE.', SPECIES_CD),
                    !grepl('ANCX', SPECIES_CD)) %>%
      dplyr::filter(SPECIES_CD != "OTH CORA") %>%
      dplyr::mutate(PRIMARY_SAMPLE_UNIT = as.character(PRIMARY_SAMPLE_UNIT)) %>% #Convert PSU to character to make sure vegan::diversity does not use it as counts
      dplyr::group_by(REGION, YEAR, SUB_REGION_NAME, PRIMARY_SAMPLE_UNIT, STATION_NR, LAT_DEGREES, LON_DEGREES, STRAT, PROT, SPECIES_NAME) %>%
      dplyr::summarise(SPP_Count = sum(N))
    
    species_wide_site <- species_site %>%
      tidyr::spread(., key = SPECIES_NAME, value = SPP_Count, fill = 0) %>%
      dplyr::ungroup()
    
    spp <- unique(species_site$SPECIES_NAME)
    
    species_only_site <- species_wide_site %>% dplyr::select(YEAR, PRIMARY_SAMPLE_UNIT, STATION_NR, spp)
    
    species_diversity_site <- species_only_site %>%
      dplyr::mutate(Simpson =  vegan::diversity(species_only_site[, -c(1:2)], index = "simpson")) %>%
      dplyr::mutate(Inv_Simpson =  vegan::diversity(species_only_site[, -c(1:2)], index = "invsimpson")) %>%
      dplyr::mutate(Shannon =  vegan::diversity(species_only_site[, -c(1:2)], index = "shannon")) %>%
      dplyr::select(YEAR, PRIMARY_SAMPLE_UNIT, STATION_NR, Simpson, Inv_Simpson, Shannon ) %>%
      # this originally didn't account for repeated PSU's across years.
      dplyr::full_join(., sites, by = c("YEAR", "PRIMARY_SAMPLE_UNIT", "STATION_NR")) %>%
      # average 2 stage data
      dplyr::group_by(REGION, YEAR, SUB_REGION_NAME, PRIMARY_SAMPLE_UNIT, LAT_DEGREES, LON_DEGREES, STRAT, PROT) %>%
      dplyr::summarize(Simpson = mean(Simpson),
                       Inv_Simpson = mean(Inv_Simpson),
                       Shannon = mean(Shannon)) %>%
      dplyr::select(REGION, YEAR, SUB_REGION_NAME, PRIMARY_SAMPLE_UNIT, LAT_DEGREES, LON_DEGREES, STRAT, PROT, Simpson, Inv_Simpson, Shannon) %>%
      dplyr::ungroup()
    
  } else{
    # Site level
    sites <- dat %>%
      # Update in Jan. 2024 - exclude juveniles - we don't count these (only once/site), so including them in diversity metrics is improper
      dplyr::filter(N == 1,
                    JUV == 0,
                    SUB_REGION_NAME != "Marquesas",
                    SUB_REGION_NAME != "Marquesas-Tortugas Trans",
                    !grepl('SPE.', SPECIES_CD),
                    !grepl('ANCX', SPECIES_CD))  %>%
      dplyr::filter(SPECIES_CD != "OTH CORA") %>%
      dplyr::mutate(PRIMARY_SAMPLE_UNIT = as.character(PRIMARY_SAMPLE_UNIT)) %>%
      dplyr::select(REGION, YEAR, SUB_REGION_NAME, PRIMARY_SAMPLE_UNIT, LAT_DEGREES, LON_DEGREES, STRAT, PROT) %>%
      dplyr::distinct(.)
    
    species_site <- dat %>%
      # Update in Jan. 2024 - exclude juveniles - we don't count these (only once/site), so including them in diversity metrics is improper
      dplyr::filter(N == 1,
                    JUV == 0,
                    SUB_REGION_NAME != "Marquesas",
                    SUB_REGION_NAME != "Marquesas-Tortugas Trans",
                    !grepl('SPE.', SPECIES_CD),
                    !grepl('ANCX', SPECIES_CD)) %>%
      dplyr::filter(SPECIES_CD != "OTH CORA") %>%
      dplyr::mutate(PRIMARY_SAMPLE_UNIT = as.character(PRIMARY_SAMPLE_UNIT)) %>% #Convert PSU to character to make sure vegan::diversity does not use it as counts
      dplyr::group_by(REGION, YEAR, SUB_REGION_NAME, PRIMARY_SAMPLE_UNIT, LAT_DEGREES, LON_DEGREES, STRAT, PROT, SPECIES_NAME) %>%
      dplyr::summarise(SPP_Count = sum(N))
    species_wide_site <- species_site %>%
      tidyr::spread(., key = SPECIES_NAME, value = SPP_Count, fill = 0) %>%
      dplyr::ungroup()
    
    spp <- unique(species_site$SPECIES_NAME)
    
    species_only_site <- species_wide_site %>% dplyr::select(YEAR, PRIMARY_SAMPLE_UNIT, spp)
    
    species_diversity_site <- species_only_site %>%
      dplyr::mutate(Simpson =  vegan::diversity(species_only_site[, -c(1:2)], index = "simpson")) %>%
      dplyr::mutate(Inv_Simpson =  vegan::diversity(species_only_site[, -c(1:2)], index = "invsimpson")) %>%
      dplyr::mutate(Shannon =  vegan::diversity(species_only_site[, -c(1:2)], index = "shannon")) %>%
      dplyr::select(YEAR, PRIMARY_SAMPLE_UNIT, Simpson, Inv_Simpson, Shannon ) %>%
      # this originally didn't account for repeated PSU's across years.
      dplyr::full_join(., sites, by = c("YEAR", "PRIMARY_SAMPLE_UNIT")) %>%
      dplyr::select(REGION, YEAR, SUB_REGION_NAME, PRIMARY_SAMPLE_UNIT, LAT_DEGREES, LON_DEGREES, STRAT, PROT, Simpson, Inv_Simpson, Shannon) %>%
      dplyr::ungroup()
  }
  
  #### Calculate regional diversity using weighting function   #### 
  tmp  <- NCRMP_make_weighted_demo_data(project, inputdata = species_diversity_site, region, datatype = "diversity",species_filter = "FALSE")
  list2env(tmp, envir = environment())
  
  
  ####Export ####
  output <- list(
    "species_list" = species_list,
    "richness_site" = richness_site,
    "richness_strata" = richness_strata,
    "Domain_est" = Domain_est,
    "species_diversity_site" = species_diversity_site,
    "diversity_strata" = diversity_strata,
    "Domain_est_div" = Domain_est_div
  )
  return(output)
}


