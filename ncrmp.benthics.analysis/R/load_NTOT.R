## Function to load NTOT files and calculate wh

# Purpose:
# support function to load NTOT files and calculate wh


## Tag: data analysis


# outputs created in this file --------------
# ntot


# CallS:
# NTOTs

# output gets called by:
# NCRMP_make_weighted_demo_data.R

# NCRMP Caribbean Benthic analytics team: Groves, Viehman, Williams
# Last update: Mar 2025


##############################################################################################################################

#' Creates weights for coral demographic data
#'
#' Creates weighting scheme for coral demographic data from a single region.
#' NCRMP utilizes a stratified random sampling design. Generally, weighting is
#' based on the number of grid cells in the sample frame in each stratum, to provide
#' regional estimates. This function is called by broader calculation functions
#' such as [NCRMP_make_weighted_demo_data()],
#' which will provide the inputdata object.
#'
#'
#'
#' @param project A string indicating the project, "NCRMP" or NCRMP and DRM combined ("NCRMP_DRM").
#' @param inputdata A dataframe of coral demographic data. Can be in various forms.
#' @param region A string indicating the region. Options are: "SEFCRI", "FLK", "Tortugas", "STX", "STTSTJ", "PRICO", and "FGB".
#' @return A dataframe of weighting scheme for corresponding region.
#' @importFrom magrittr "%>%"
#' @export
#'
#'
#'

load_NTOT <- function(region, inputdata, project){
  
  ####Process NTOT####
  process_ntot <- function(data, year, strat_mapping = NULL){
    data <- data %>%
      dplyr::mutate(ANALYSIS_STRATUM = paste(STRAT, "/ PROT =", PROT, sep = " "),
                    ngrtot = sum(NTOT))
  }

  ####Filter By Year####
  #This replaces the original for loop
  filter_by_year <- function(NTOT_all, data, year) {
    #filter data to only include the inputed year
    strata_to_keep <- data %>%
      filter(YEAR == year) %>%
      #extract values from analysis stratum and get non redundant values
      pull(ANALYSIS_STRATUM) %>%
      unique()
    
    #do the filtering and calculate ngrtot
    NTOT_all %>%
      filter(YEAR == year, ANALYSIS_STRATUM %in% strata_to_keep) %>%
      mutate(ngrtot = sum(NTOT))
  }
  
  #### SEFCRI ####

  #Helper function to process 2016 and 2018 sefcri 
  SEFL_16_18_process <- function(data){
    data <- data %>%
      dplyr::mutate(STRAT = dplyr::case_when(STRAT == "PTSH2"~"NEAR1",
                                             STRAT == "PTDP0"~"OFFR0",
                                             STRAT == "PTDP1"~"OFFR1", TRUE ~ as.character(STRAT))) %>%
      dplyr::group_by(YEAR, REGION, STRAT, PROT, GRID_SIZE) %>%
      dplyr::summarise(NTOT = sum(NTOT)) %>%
      dplyr::ungroup() %>%
      dplyr::filter(STRAT != "RGDP1" & STRAT != "RGDP0")
  }
  

  if(region == "SEFCRI") {
    if(project == "NCRMP" || project == "NULL"){

      ntot <- bind_rows(
        SEFL_2014_NTOT %>% dplyr::filter(STRAT %in% c("MIDR1", "MIDR0")) %>% process_ntot(),
        SEFL_2016_NTOT %>% SEFL_16_18_process() %>% process_ntot(),
        SEFL_2018_NTOT %>% SEFL_16_18_process() %>% process_ntot(),
        SEFL_2020_NTOT %>% dplyr::filter(!STRAT %in% c("RGDP1", "RGDP0")) %>% process_ntot(),
        SEFL_2022_NTOT %>% dplyr::filter(!STRAT %in% c("RGDP1", "RGDP0")) %>% process_ntot()
      )
    }

    if(project == "NCRMP_DRM"){

      # Filter NTOT to only strata sampled that year
      # Make a dataframe of just the YEAR and STRAT
      tmp <- inputdata %>%
        dplyr::mutate(ANALYSIS_STRATUM = paste(STRAT, "/ PROT =", PROT, sep = " ")) %>%
        dplyr::group_by(YEAR, ANALYSIS_STRATUM) %>%
        dplyr::summarise(N = length(ANALYSIS_STRATUM))

      # Make a list of all the years
      Years <- sort(unique(tmp$YEAR))

      # create a data frame of the full NTOTs for FLK
      NTOT_all <- dplyr::bind_rows(SEFL_2014_NTOT, SEFL_2014_NTOT %>% dplyr::mutate(YEAR = 2015),
                                   SEFL_2016_NTOT, SEFL_2018_NTOT %>% dplyr::mutate(YEAR = 2017),
                                   SEFL_2018_NTOT, SEFL_2018_NTOT %>% dplyr::mutate(YEAR = 2019),
                                   SEFL_2020_NTOT, SEFL_2020_NTOT %>% dplyr::mutate(YEAR = 2021),
                                   SEFL_2022_NTOT) %>%
        dplyr::mutate(REGION = "SEFCRI", ANALYSIS_STRATUM = paste(STRAT, "/ PROT =", PROT, sep = " "))

      #this uses the map_dfr function from the purrr package (this can be found in tidyverse)
      #can be used to apply a function (in this case the helper function filter by year) to every element of a list 
      #in this case the list is Years
      #so its calling filter by year with the NTOT_alll and tmp arguments, 
      #...along with .x which represents the current year being processed
      ntot <- map_dfr(Years, ~filter_by_year(NTOT_all, tmp, .x))%>%
        mutate(PROT = as.factor(PROT))
    }
  }
  
  #### FLK ####

  #helper function for processing FLK ntot data
  process_FLK_NTOT <- function(data){
    data %>%
      mutate(PROT = 0) %>%
      dplyr::group_by(REGION, YEAR, PROT, STRAT, GRID_SIZE) %>%
      dplyr::summarise(NTOT = sum(NTOT), .groups = "drop")
  }

  if(region == "FLK") {
    if(project == "NCRMP" || project == "NULL"){

      # Use a loop to create a unique lists for each year of strata sampled
      # Filter NTOT to only strata sampled that year (previously done manually)
      tmp <- inputdata %>%
        dplyr::mutate(ANALYSIS_STRATUM = paste(STRAT, "/ PROT =", PROT, sep = " ")) %>%
        dplyr::group_by(YEAR, ANALYSIS_STRATUM) %>%
        dplyr::summarise(N = length(ANALYSIS_STRATUM), .groups = "keep")

      # Make a list of all the years
      Years <- sort(unique(tmp$YEAR))

      ### UPDATE IN DEC. 2023!!
      # PROT is re-coded here to 0 for ALL strata as fish and benthics met 12/19/23
      # to determine that it is not appropriate to keep PROT in the analysis strat
      # in FLK because the data aren't allocated that way
      # NTOTs must be re-calculated with PROT=0 here
      NTOT_all <- dplyr::bind_rows(
        process_FLK_NTOT(FLK_2014_NTOT),
        process_FLK_NTOT(FLK_2016_NTOT),
        process_FLK_NTOT(FLK_2018_NTOT),
        process_FLK_NTOT(FLK_2020_NTOT),
        process_FLK_NTOT(FLK_2022_NTOT)
      ) %>%
        dplyr::mutate(REGION = "FLK", ANALYSIS_STRATUM = paste(STRAT, "/ PROT =", PROT, sep = " "))

      #see above for explanation of what this is doing
      ntot <- map_dfr(Years, ~filter_by_year(NTOT_all, tmp, .x))%>%
        mutate(PROT = as.factor(PROT))
    }


    if(project == "NCRMP_DRM") {

      # Filter NTOT to only strata sampled that year
      # Make a dataframe of just the YEAR and STRAT
      tmp <- inputdata %>%
        dplyr::mutate(ANALYSIS_STRATUM = paste(STRAT, "/ PROT =", PROT, sep = " ")) %>%
        dplyr::group_by(YEAR, ANALYSIS_STRATUM) %>%
        dplyr::summarise(N = length(ANALYSIS_STRATUM), .groups = "keep")

      # Make a list of all the years
      Years <- sort(unique(tmp$YEAR))

      ### UPDATE IN DEC. 2023!!
      # PROT is re-coded here to 0 for ALL strata as fish and benthics met 12/19/23
      # to determine that it is not appropriate to keep PROT in the analysis strat
      # in FLK because the data aren't allocated that way
      # NTOTs must be re-calculated with PROT=0 here
      NTOT_all <- dplyr::bind_rows(
        process_FLK_NTOT(FLK_2014_NTOT),
        process_FLK_NTOT(FLK_2016_NTOT),
        process_FLK_NTOT(FLK_2018_NTOT),
        process_FLK_NTOT(FLK_2020_NTOT),
        process_FLK_NTOT(FLK_2022_NTOT)
      ) %>%
        dplyr::mutate(REGION = "FLK",
                      ANALYSIS_STRATUM = paste(STRAT, "/ PROT =", PROT, sep = " "))

      #see above for explanation of what this is doing
      ntot <- map_dfr(Years, ~filter_by_year(NTOT_all, tmp, .x))%>%
        mutate(PROT = as.factor(PROT))

    }

    #### MIR ####
    if(project == "MIR"){

      ntot22 <- FLK_MIR_2022_NTOT %>%
        dplyr::ungroup() %>%
        # Remove strata not sampled - this will be unique to each year and need to be updated. See Tortugas for examples.
        dplyr::filter(STRAT_CORA != "CS05",
                      STRAT_CORA != "CS16",
                      ANALYSIS_STRATUM != "CS02 / PROT = 1",
                      ANALYSIS_STRATUM != "CS07 / PROT = 1",
                      ANALYSIS_STRATUM != "CS12 / PROT = 0",
                      ANALYSIS_STRATUM != "CS14 / PROT = 1",
                      ANALYSIS_STRATUM != "CS15 / PROT = 1") %>%
        dplyr::mutate(YEAR = 2022,
                      ngrtot = sum(NTOT))

      ntot <- ntot22
    }
  }

  #### Tortugas ####
  if(region == "Tortugas") {
      if(project == "NCRMP" || project == "NULL") {
    
    tort_analysis_strat <- function(data){
      data %>%
        dplyr::mutate(REGION = "Tortugas",
                      ANALYSIS_STRATUM = paste(STRAT, "/ PROT =", PROT, sep = " "))
    }
    
    ntot14 <- Tort_2014_NTOT %>%
      dplyr::filter(STRAT != "SPGR_LR") %>% # Not sampled in 2014
      tort_analysis_strat() %>%
      dplyr::mutate( ngrtot = sum(NTOT))
    
    ntot16 <- Tort_2016_NTOT %>%
      tort_analysis_strat() %>%
      dplyr::filter(ANALYSIS_STRATUM != "ISOL_LR / PROT = 0") %>% #Not sampled
      dplyr::mutate(ngrtot = sum(NTOT))
    
    ntot18 <- Tort_2018_NTOT %>%
      dplyr::filter(STRAT != "SPGR_LR") %>% # Not sampled in 2018
      tort_analysis_strat() %>%
      dplyr::filter(ANALYSIS_STRATUM != "ISOL_LR / PROT = 0",
                    ANALYSIS_STRATUM != "ISOL_LR / PROT = 1") %>% # Not sampled in 2018
      dplyr::mutate(ngrtot = sum(NTOT))
    
    ntot20 <- Tort_2020_NTOT %>%
      tort_analysis_strat() %>%
      dplyr::mutate(PROT = as.integer(PROT)) %>%
      dplyr::filter(ANALYSIS_STRATUM != "T09 / PROT = 0",
                    ANALYSIS_STRATUM != "T09 / PROT = 1",
                    ANALYSIS_STRATUM != "T03 / PROT = 2",
                    ANALYSIS_STRATUM != "T06 / PROT = 2") %>% #Not sampled
      dplyr::mutate(ngrtot = sum(NTOT))
    
    ntot22 <- Tort_2022_NTOT %>%
      tort_analysis_strat() %>%
      dplyr::filter(PROT != 0) %>%
      dplyr::filter(ANALYSIS_STRATUM != "T07 / PROT = 1",
                    ANALYSIS_STRATUM != "T09 / PROT = 1",
                    ANALYSIS_STRATUM != "T04 / PROT = 2") %>% #Not sampled
      dplyr::mutate(ngrtot = sum(NTOT))
    
    ntot <- dplyr::bind_rows(ntot14, ntot16, ntot18, ntot20, ntot22)
      }


    if(project == "NCRMP_DRM") {

      # Filter NTOT to only strata sampled that year
      # Make a dataframe of just the YEAR and STRAT
      tmp <- inputdata %>%
        dplyr::mutate(ANALYSIS_STRATUM = paste(STRAT, "/ PROT =", PROT, sep = " ")) %>%
        dplyr::group_by(YEAR, ANALYSIS_STRATUM) %>%
        dplyr::summarise(N = length(ANALYSIS_STRATUM))

      # Make a list of all the years
      Years <- sort(unique(tmp$YEAR))
      # add an empty data frame to populate with the filtered NTOTs
      ntot <- data.frame()
      # create a data frame of the full NTOTs for FLK
      NTOT_all <- dplyr::bind_rows(Tort_2014_NTOT, Tort_2016_NTOT %>% dplyr::mutate(YEAR = 2015),
                                   Tort_2016_NTOT, Tort_2018_NTOT %>% dplyr::mutate(YEAR = 2017),
                                   Tort_2018_NTOT, Tort_2018_NTOT %>% dplyr::mutate(YEAR = 2019),
                                   Tort_2020_NTOT, Tort_2020_NTOT %>% dplyr::mutate(YEAR = 2021),
                                   Tort_2022_NTOT) %>%
        dplyr::mutate(REGION = "Tortugas", ANALYSIS_STRATUM = paste(STRAT, "/ PROT =", PROT, sep = " "))

      #see above for explanation of what this is doing
      ntot <- map_dfr(Years, ~filter_by_year(NTOT_all, tmp, .x))%>%
        mutate(PROT = as.factor(PROT))
    }
  }

  #### USVI Function ####

  process_USVI_NTOT <- function(data, year){
    data %>%
      dplyr::group_by(REGION, YEAR, STRAT, HABITAT_CD, DEPTH_STRAT) %>%
      dplyr::summarise(NTOT = sum(NTOT)) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(YEAR = year,
                    ANALYSIS_STRATUM = STRAT,
                    PROT = NA_character_,
                    ngrtot = sum(NTOT)) %>%
      dplyr::ungroup()

  }

  #### STTSTJ ####
  if(region == "STTSTJ"){

    ntot13 <- USVI_2021_NTOT %>%
      dplyr::filter(REGION == "STTSTJ",
                    STRAT != "HARD_SHLW") %>% # Hard shlw was not sampled in 2013
      process_USVI_NTOT(year = 2013)

    ntot15 <- USVI_2021_NTOT %>%
      dplyr::filter(REGION == "STTSTJ") %>%
      process_USVI_NTOT(year = 2015)

    ntot17 <- USVI_2021_NTOT %>%
      dplyr::filter(REGION == "STTSTJ") %>%
      process_USVI_NTOT(year = 2017)

    ntot19 <- USVI_2021_NTOT %>%
      dplyr::filter(REGION == "STTSTJ") %>%
      process_USVI_NTOT(year = 2019)

    ntot21 <- USVI_2021_NTOT %>%
      dplyr::filter(REGION == "STTSTJ") %>%
      process_USVI_NTOT(year = 2021)

    ntot23 <- USVI_2023_NTOT %>%
      dplyr::filter(REGION == "STTSTJ") %>%
      process_USVI_NTOT(year = 2023)

    ntot <- dplyr::bind_rows(ntot13, ntot15, ntot17, ntot19, ntot21, ntot23)

  }
  #### STX ####
  if(region == "STX"){

    ntot15 <- USVI_2021_NTOT %>%
      dplyr::filter(REGION == "STX",
                    STRAT != "HARD_SHLW", # Hard shlw was not sampled in 2015
                    STRAT != "HARD_DEEP") %>% # Hard deep was not sampled in 2015
      process_USVI_NTOT(year = 2015)

    ntot17 <- USVI_2021_NTOT %>%
      dplyr::filter(REGION == "STX",
                    STRAT != "HARD_SHLW") %>%
      process_USVI_NTOT(year = 2017)

    ntot19 <- USVI_2021_NTOT %>%
      dplyr::filter(REGION == "STX") %>%
      process_USVI_NTOT(year = 2019)

    ntot21 <- USVI_2021_NTOT %>%
      dplyr::filter(REGION == "STX") %>%
      dplyr::filter(!(STRAT == "BDRK_SHLW" | STRAT == "BDRK_DEEP")) %>% # 2021 we didnt sample BDRK shallow OR deep
      process_USVI_NTOT(year = 2021)

    ntot23 <- USVI_2023_NTOT %>%
      dplyr::filter(REGION == "STX") %>%
      process_USVI_NTOT(year = 2023)

    ntot <- dplyr::bind_rows(ntot15, ntot17, ntot19, ntot21, ntot23)

  }

  #### PRICO ####
  PRICO_data_clean <- function(data, year){
    
    #if 2023, filter out HARD habitat code
    if (year == 2023) {data <- dplyr::filter(data, HABITAT_CD != "HARD")}

    data <- data %>%
      dplyr::group_by(REGION, YEAR, STRAT, HABITAT_CD, DEPTH_STRAT) %>%
      dplyr::summarise(NTOT = sum(NTOT)) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(YEAR = year,
                    ANALYSIS_STRATUM = STRAT,
                    PROT = NA_character_,
                    ngrtot = sum(NTOT))
  }


  if(region == "PRICO"){

    ntot14 <- PRICO_2023_NTOT %>%
      dplyr::filter(STRAT != "HARD_DEEP", # Hard shlw was not sampled in 2014
                    STRAT != "HARD_SHLW") %>% # Hard deep was not sampled in 2014
      PRICO_data_clean(year = 2014)

    #list of all other years besides 2014 
    years <- c(2014,2019,2021, 2023)
    
    other_years <- map_dfr(years, ~PRICO_data_clean(PRICO_2023_NTOT, .x))

    ntot <- dplyr::bind_rows(ntot14, other_years)
  }

  #### FGB ####
  if(region == "FGB"){

    FGB_data_clean <- function(data, year){
      data <- data %>%
        dplyr::mutate(ANALYSIS_STRATUM = "FGBNMS",
                      PROT = NA_character_,
                      YEAR = year) %>%
        dplyr::group_by(REGION, YEAR, ANALYSIS_STRATUM, DEPTH_STRAT, PROT) %>%
        dplyr::summarise(NTOT = sum(NTOT),
                         ngrtot = sum(NTOT)) %>%
        dplyr::ungroup()
    }
    
    years <- c(2013, 2015, 2018, 2022, 2024)
    
    ntot <- map_dfr(years, ~FGB_data_clean(FGBNMS_2024_NTOT, .x))
  }

  ntot <- ntot %>%
    dplyr::mutate(wh = NTOT/ngrtot) %>%
    dplyr::mutate(PROT = as.factor(PROT))

  ####Export####
  return(ntot)
}
