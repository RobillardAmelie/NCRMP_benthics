## Function to create percent cover of macroalgae, turf algae, CCA and other inverts

# Purpose:
# create csv files with percent cover by site at the species and group levels


## Tag: data analysis


# outputs created in this file --------------
# Domain estimates
# cover_strata
# percent_cover_site
# percent_cover_species
# cover_group_key
# Domain_check


# CallS:
# analysis ready data

# output gets called by:
# Analysis Rmarkdown, etc.
#

# NCRMP Caribbean Benthic analytics team: Groves, Viehman, Williams
# Last update: Feb 2023


##############################################################################################################################

#' Creates percent cover dataframe
#'
#' Calculates percent cover at various levels and groups. Data summaries include:
#' 1) percent cover by species at each site, 2) percent cover by group at each site,
#' 3) mean percent cover by group for each stratum, and 4) weighted regional mean
#' percent cover, all for a given region. NCRMP utilizes a stratified random
#' sampling design. Regional estimates of cover are weighted by the number of
#' grid cells of a stratum in the sample frame.
#'
#'
#'
#'
#'
#' @param region A string indicating the region. Options are: "SEFCRI", "FLK", "Tortugas", "STX", "STTSTJ", "PRICO", and "FGB".
#' @param project A string indicating the project: "NCRMP" or "MIR".
#' @return A list of dataframes, including site level cover by species, site level cover by group, strata mean cover by group, and weighted regional mean cover by group.
#' @importFrom magrittr "%>%"
#' @export
#'
#'


NCRMP_calculate_cover <- function(region, project = "NULL"){

  ####Calc Cover Data####
  cover_data <- load_NCRMP_benthic_cover_data(region = region, project = project) %>%
    dplyr::mutate(LAT_DEGREES = sprintf("%0.4f", LAT_DEGREES),
                  LON_DEGREES = sprintf("%0.4f", LON_DEGREES),
                  PROT = as.factor(PROT))

  ####Calculate Percent Cover####
  calc_perc_cover <- function(data) {
    data %>%
      dplyr::mutate(Percent_Cvr = HARDBOTTOM_P+SOFTBOTTOM_P+RUBBLE_P) %>%
      dplyr::select(-HARDBOTTOM_P, -SOFTBOTTOM_P, -RUBBLE_P)
  }


  # simplify species_name for MAC FLES and MAC CALC (as were not used in NCRMP) more detailed species collected in scream survey
  ncrmp_frrp_sppcodes2 <- ncrmp_frrp_sppcodes %>%
    dplyr::select(fl_ncrmp_code, species_name, cover_group) %>%
    # remove duplicate codes
    dplyr::distinct(fl_ncrmp_code, species_name, cover_group) %>%
    dplyr::filter(!fl_ncrmp_code %in% c("ERY CARI", "POF SPE.", "ENCR GORG"), species_name != "Undaria spp") %>%
    dplyr::mutate(species_name = as.character(species_name)) %>%
    dplyr::mutate(cover_group = as.character(cover_group)) %>%
    dplyr::mutate(cover_group = dplyr::case_when(species_name == "Ramicrusta spp" ~ "RAMICRUSTA SPP.",
                                                 species_name == "Peysonnellia" ~ "PEYSONNELLIA",
                                                 species_name == "Turf Algae Free of Sediment" ~ "TURF ALGAE",
                                                 species_name == "Turf Algae with Sediment" ~ "TURF ALGAE",
                                                 cover_group == "ALGAE" ~ "MACROALGAE",
                                                 TRUE ~ cover_group))

  cover_group_key <-  ncrmp_frrp_sppcodes2


  # Create a dataframe with site level mean depths, we will add this back in later
  # do this at the site level (distinct just to site level characteristics)
  # otherwise it will average biased towards the site with more categories
  depth <- cover_data %>%
    dplyr::ungroup() %>%
    dplyr::select(REGION, YEAR, MONTH, DAY, SUB_REGION_NAME, ADMIN, PRIMARY_SAMPLE_UNIT, LAT_DEGREES, LON_DEGREES, ANALYSIS_STRATUM, STRAT, HABITAT_CD, PROT, MIN_DEPTH, MAX_DEPTH) %>%
    dplyr::distinct(.) %>%
    dplyr::group_by(REGION, YEAR, MONTH, DAY, SUB_REGION_NAME, ADMIN, PRIMARY_SAMPLE_UNIT, LAT_DEGREES, LON_DEGREES,ANALYSIS_STRATUM, STRAT, HABITAT_CD, PROT) %>%
    dplyr::summarise(MIN_DEPTH = mean(MIN_DEPTH),
                     MAX_DEPTH = mean(MAX_DEPTH))


  # Calculate percent cover of species by site

  if (project == "NCRMP" && region %in% c("SEFCRI", "FLK", "Tortugas")) {

    # UPDATED in Nov. 2023
    # previously, this only really averaged categories present in both stations
    # but if a category was present in only 1 station, it was not averaged, but
    # instead it was just assigned the % of the station it was present
    # because we don't have rows for 0% cover categories in our raw data
    n_stations <- cover_data %>%
      dplyr::ungroup() %>%
      dplyr::group_by(REGION, YEAR, MONTH, DAY, SUB_REGION_NAME, ADMIN, PRIMARY_SAMPLE_UNIT, LAT_DEGREES, LON_DEGREES, ANALYSIS_STRATUM, STRAT, HABITAT_CD, PROT) %>%
      dplyr::summarise(n_stations = length(unique(STATION_NR)))


    dat2 <- cover_data %>%
      calc_perc_cover() %>%
      dplyr::ungroup() %>%
      dplyr::group_by(REGION, YEAR, MONTH, DAY, SUB_REGION_NAME, ADMIN, PRIMARY_SAMPLE_UNIT, LAT_DEGREES, LON_DEGREES,ANALYSIS_STRATUM, STRAT, HABITAT_CD, PROT, COVER_CAT_NAME) %>%
      dplyr::summarise(Percent_Cvr_site_sum = sum(Percent_Cvr)) %>%
      dplyr::left_join(., n_stations) %>%
      dplyr::mutate(Percent_Cvr_site = Percent_Cvr_site_sum/n_stations) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(Percent_Cvr = Percent_Cvr_site) %>%
      dplyr::left_join(., depth) %>%
      dplyr::select(REGION, YEAR, MONTH, DAY, SUB_REGION_NAME, ADMIN, PRIMARY_SAMPLE_UNIT, LAT_DEGREES, LON_DEGREES, MIN_DEPTH, MAX_DEPTH, ANALYSIS_STRATUM, STRAT, HABITAT_CD, PROT, COVER_CAT_NAME, Percent_Cvr)
  }


  if((project == "MIR" && region == "FLK") || (project == "NCRMP" && region %in% c("STTSTJ", "STX", "PRICO", "FGB"))) {

    dat2 <- cover_data %>%
      calc_perc_cover()%>%
      dplyr::mutate(PROT = NA) %>%
      dplyr::select(REGION, YEAR, MONTH, DAY, SUB_REGION_NAME, ADMIN, PRIMARY_SAMPLE_UNIT, LAT_DEGREES, LON_DEGREES, MIN_DEPTH, MAX_DEPTH, ANALYSIS_STRATUM, STRAT, HABITAT_CD, PROT, COVER_CAT_NAME, Percent_Cvr)
  }

  cvr_wide <- dat2 %>%
    # add in zeros for species that didn't occur per site. Easiest to flip to wide format ( 1 row per site) for this
     tidyr::spread(., COVER_CAT_NAME, Percent_Cvr,
                  fill = 0)

  percent_cover_species <- tidyr::gather(cvr_wide, COVER_CAT_NAME, Percent_Cvr, 16:ncol(cvr_wide)) %>%
    dplyr::filter(COVER_CAT_NAME != '<NA>') %>%
    dplyr::left_join(.,ncrmp_frrp_sppcodes2 %>% dplyr::rename("COVER_CAT_NAME" = species_name),  by = "COVER_CAT_NAME")  %>%
    dplyr::select(REGION, YEAR, SUB_REGION_NAME, ADMIN,  PRIMARY_SAMPLE_UNIT, LAT_DEGREES, LON_DEGREES, MIN_DEPTH, MAX_DEPTH, ANALYSIS_STRATUM, STRAT,  HABITAT_CD, PROT, COVER_CAT_NAME, cover_group, Percent_Cvr)


  # Calculate percent cover of major biological categories by site
    dat1 <-  cover_data %>%
      calc_perc_cover()%>%
      dplyr::left_join(.,ncrmp_frrp_sppcodes2,  by = c( "COVER_CAT_CD" = "fl_ncrmp_code")) %>% # previously this was an inner join, but that's bad because it ignores things that maybe don't match with species codes
      dplyr::group_by(YEAR, REGION, SUB_REGION_NAME, ADMIN, PRIMARY_SAMPLE_UNIT, STATION_NR, LAT_DEGREES, LON_DEGREES, ANALYSIS_STRATUM, STRAT, HABITAT_CD, PROT, cover_group) %>%
      dplyr::summarise(Percent_Cvr = sum(Percent_Cvr)) %>%
      dplyr::ungroup()


  # Calculate percent cover at the STRATA level
  ### Reclassify all other categories to other or macroalgae

  dat2 <- dat1 %>%
    dplyr::mutate(cover_group = dplyr::case_when(cover_group == "SUBSTRATE" ~ "OTHER",
                                                 cover_group == "HYDROCORALS" ~ "OTHER",
                                                 cover_group == "OTHER INVERTEBRATES" ~ "OTHER",
                                                 cover_group == "CYANOBACTERIA" ~ "OTHER",
                                                 cover_group == "SEAGRASSES" ~ "OTHER",
                                                 cover_group == "PEYSONNELLIA" ~ "MACROALGAE",
                                                 cover_group == "OTHER" ~ "OTHER",
                                                 TRUE ~ cover_group))


  ### Account for SEFCRI/FLK 2014 & Tortugas 2018 2 transect data - take the transect means

  if(region %in% c("SEFCRI", "FLK", "Tortugas")){

    ### Noticed in Nov. 2023, that this doesn't work with sites where a category
    # isn't present in both stations
    # FIXED in Nov. 2023

    dat3 <- dat2 %>%
      dplyr::group_by(YEAR, REGION, SUB_REGION_NAME, ADMIN, PRIMARY_SAMPLE_UNIT, STATION_NR, LAT_DEGREES, LON_DEGREES, ANALYSIS_STRATUM, STRAT, HABITAT_CD, PROT, cover_group) %>%
      dplyr::summarise(Percent_Cvr = sum(Percent_Cvr)) %>%
      dplyr::ungroup() %>%
      dplyr::group_by(YEAR, REGION, SUB_REGION_NAME, ADMIN, PRIMARY_SAMPLE_UNIT, LAT_DEGREES, LON_DEGREES, ANALYSIS_STRATUM, STRAT, HABITAT_CD, PROT, cover_group) %>%
      dplyr::summarise(Percent_Cvr_sum = sum(Percent_Cvr)) %>%
      dplyr::left_join(., n_stations) %>%
      dplyr::mutate(Percent_Cvr = Percent_Cvr_sum/n_stations) %>%
      dplyr::select(-Percent_Cvr_sum, -n_stations) %>%
      dplyr::ungroup()

  } else{
    dat3 <- dat2 %>%
      dplyr::group_by(YEAR, REGION, SUB_REGION_NAME, ADMIN, PRIMARY_SAMPLE_UNIT, LAT_DEGREES, LON_DEGREES, ANALYSIS_STRATUM, STRAT, HABITAT_CD, PROT, cover_group) %>%
      dplyr::summarise(Percent_Cvr = sum(Percent_Cvr)) %>%
      dplyr::ungroup()
  }


  # Add in cover groups that are missing from certain sites
  groupsOfInterest <- c("CCA", "HARD CORALS", "MACROALGAE", "OTHER", "SOFT CORALS", "SPONGES", "TURF ALGAE", "RAMICRUSTA SPP.")

  # make unique sites df
  if(project == "NCRMP" & region %in% c("FLK", "SEFCRI", "Tortugas")){
    allsites <- unique(dat3[, c("YEAR", "REGION", "SUB_REGION_NAME", 'ADMIN', "PRIMARY_SAMPLE_UNIT", "MONTH", "DAY",
                                "LAT_DEGREES", "LON_DEGREES", "ANALYSIS_STRATUM", "STRAT", "HABITAT_CD", "PROT")])
  }
  else{
    allsites <- unique(dat3[, c("YEAR", "REGION", "SUB_REGION_NAME", 'ADMIN', "PRIMARY_SAMPLE_UNIT",
                                "LAT_DEGREES", "LON_DEGREES", "ANALYSIS_STRATUM", "STRAT", "HABITAT_CD", "PROT")])
  }

  # add missiong groups to data
  groupsOfInterest2 <- data.frame(cover_group = groupsOfInterest, x = 1)
  dat4 <- allsites %>%
    dplyr::mutate(x = 1) %>%
    dplyr::full_join(., groupsOfInterest2) %>%
    dplyr::full_join(., dat3) %>%
    dplyr::mutate(Percent_Cvr = case_when(is.na(Percent_Cvr) ~ 0,
                                          TRUE ~ Percent_Cvr)) %>%
    dplyr::mutate(n = dplyr::case_when(Percent_Cvr > 0 ~ 1, TRUE ~ 0)) %>%
    # add site-averaged depth (only an average for FL 2014 data)
    dplyr::left_join(., depth) %>%
    dplyr::select(-x)


  percent_cover_site <- dat4 %>%
    dplyr::select(REGION, YEAR, SUB_REGION_NAME, ADMIN,  PRIMARY_SAMPLE_UNIT, LAT_DEGREES, LON_DEGREES,MIN_DEPTH, MAX_DEPTH, ANALYSIS_STRATUM, STRAT,  HABITAT_CD, PROT, cover_group, Percent_Cvr, n)

  # site level totals are ~100% (may not be exact for 2 transect data where means are taken)
  cover_check_site <- percent_cover_site %>%
    dplyr::group_by(REGION, YEAR, PRIMARY_SAMPLE_UNIT, STRAT) %>%
    dplyr::summarise(Cover = sum(Percent_Cvr))

  lpi_data  <- NCRMP_make_weighted_LPI_data(inputdata = dat4, region = region, project = project)

  # unpack list
  list2env(lpi_data, envir = .GlobalEnv)

  # Check your Domain totals - they should equal about 100
  # If this number is not very close to 100, there is likely something wrong with your ntot. Any strata that was
  # not sampled that year needs to be removed from the ntot BEFORE ngrot is calculated. Slight variations from 100
  # could be due to using the most recent NTOT file for past years - this is OK.

  Domain_check <- Domain_est %>%
    dplyr::group_by(REGION, YEAR) %>%
    dplyr::summarise(Whole_pie = sum(avCvr))

####Export####
  output <- list(
    "percent_cover_species" = percent_cover_species,
    "percent_cover_site" = percent_cover_site,
    "cover_check_site" = cover_check_site,
    "cover_group_key" = cover_group_key,
    "cover_strata" = cover_strata,
    "Domain_est" = Domain_est,
    "Domain_check" = Domain_check
  )

  return(output)
}
