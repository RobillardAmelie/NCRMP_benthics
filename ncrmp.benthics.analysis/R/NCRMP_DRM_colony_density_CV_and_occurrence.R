## Function to calculate colony density cv by species, site level occurrence, identify species with CV less than 20%, and plot back to back figure

# Purpose:
# creates csv files with colony density, occurrence and a plot of both.


## Tag: data analysis


# outputs created in this file --------------
#
#
#


# CallS:
# analysis ready data

# output gets called by:
# Analysis Rmarkdown, etc.
#

# NCRMP Caribbean Benthic analytics team: Groves, Viehman, Williams, Sturm
# Last update: Sep 2024


##############################################################################################################################

#' Creates colony density summary dataframes
#'
#' Calculates regional estimate of coral density and coefficient of variation (CV),
#' by species, for a given region. NCRMP utilizes a stratified random
#' sampling design. Regional estimates of density are weighted by the number of grid cells of a stratum
#' in the sample frame. Function calculates weighted strata means to produce
#' regional estimates for coral density data by species.
#' Also calculates occurrence of each species in each year.
#' Additionally, function produces a figure of species occurrence and CV by species
#' for a given region and year.
#'
#'
#'
#'
#' @param region A string indicating the region. Options are: "SEFCRI", "FLK", "Tortugas", "STX", "STTSTJ", "PRICO", and "FGB".
#' @param ptitle A string indicating the plot title, usually the region.
#' @param year A numeric indicating the year of interest, which will be plotted.
#' @param path A string indicating the filepath for the figure.
#' @param project A string indicating the project, "NCRMP" or NCRMP and DRM combined ("NCRMP_DRM").
#' @return A list dataframes. First, a dataframe of regional weighted mean density, CV,
#' and occurrence, by species for a given region. Second, a dataframe of the same,
#' filtered to only species/years where CV is less than or equal to 20%.
#' @importFrom magrittr "%>%"
#' @importFrom ggplot2 "ggplot"
#' @export
#'
#'


NCRMP_DRM_colony_density_CV_and_occurrence <- function(region, ptitle, year, file_path = "NULL", species_filter = "NULL", project = "NULL"){

  #### Call Function to get Species ####
    coral_species_by_region <- species_for_CV_and_occurrence(region)

    ####Get dataset for sppdens####
  #based on region and project
  sppdens <- switch(region,
                    "FLK" = switch(project,
                                   "NCRMP" = NCRMP_FLK_2014_22_density_species,
                                   "MIR" = MIR_2022_density_species_DUMMY,
                                   "NCRMP_DRM" = NCRMP_DRM_FLK_2014_22_density_species,
                                   stop("Unknown project for FLK")),
                    "Tortugas" = switch(project,
                                    "NCRMP_DRM" = NCRMP_DRM_Tort_2014_22_density_species,
                                    "NCRMP" = NCRMP_Tort_2014_22_density_species,
                                    stop("Unknown project for Tortugas")),
                    "SEFCRI" = switch(project,
                                      "NCRMP" = NCRMP_SEFCRI_2014_22_density_species,
                                      "NCRMP_DRM" = NCRMP_DRM_SEFCRI_2014_22_density_species,
                                      stop("Unknown project for SEFCRI")),
                    "PRICO" = NCRMP_PRICO_2014_23_density_species,
                    "STTSTJ" = NCRMP_STTSTJ_2013_23_density_species,
                    "STX" = NCRMP_STX_2015_23_density_species,
                    "FGB" = NCRMP_FGBNMS_2013_24_density_species,
                    stop("Unknown region"))


  #call NCRMP_make_weighted_density_CV_data function to get region means
  region_means <- NCRMP_make_weighted_density_CV_data(region = region, sppdens = sppdens, project = project)

  #filter species if user adds a species filter
  if (species_filter == TRUE) {
    region_means <- region_means %>% dplyr::filter(SPECIES_CD %in% coral_species_by_region)
  }

  region_means <- region_means %>% filter(YEAR == year)

    plot <- function(region_means, set_occurrence, year, ptitle) {
      g.mid <- region_means %>%
        dplyr::filter(YEAR == year) %>%
        dplyr::filter(occurrence > set_occurrence) %>%
        dplyr::filter(CV < 1) %>%
        ggplot(aes(x = 1, y = reorder(SPECIES_CD, occurrence))) +
        geom_text(aes(label = SPECIES_CD,
                      fontface = ifelse(SPECIES_CD %in% c("Acropora cervicornis", "Acropora palmata",
                                                          "Dendrogyra cylindrus", "Orbicella annularis",
                                                          "Orbicella faveolata", "Orbicella franksi",
                                                          "Mycetophyllia ferox"),
                                        "bold.italic", "italic")), size = 3) +
        geom_segment(aes(x = 0.94, xend = 0.96, yend = SPECIES_CD)) +
        geom_segment(aes(x = 1.04, xend = 1.065, yend = SPECIES_CD)) +
        ggtitle(ptitle) +
        ylab(NULL) +
        scale_x_continuous(expand = c(0, 0), limits = c(0.94, 1.065)) +
        theme(axis.title = element_blank(),
              panel.grid = element_blank(),
              axis.text.y = element_blank(),
              axis.ticks.y = element_blank(),
              panel.background = element_blank(),
              axis.text.x = element_text(color = NA),
              axis.ticks.x = element_line(color = NA),
              plot.margin = unit(c(t = 1, r = -1, b = 1, l = -1), "mm"),
              plot.title = element_text(hjust = 0.5))


    g1 <-  region_means %>%
      # filter to year of interest
      dplyr::filter(YEAR == year) %>%
      # exclude occurrences of 0
      dplyr::filter(occurrence > 0.01) %>%
      # exclude CVs over 1
      dplyr::filter(CV < 1) %>%

      ggplot(.,
             aes(x = reorder(SPECIES_CD, occurrence),
                 y = occurrence,
                 fill = 'even')) +
      # geom_hline(yintercept = c(0.25, 0.5, 0.75),
      # colour = "light grey") +
      geom_bar(stat = "identity",
               fill = "deepskyblue4") +
      ggtitle(paste("Species", "Occurrence", sep = " ")) +
      # scale_fill_manual(values = c( "#0a4595")) +
      theme_light() +
      scale_y_continuous(expand = c(0,0)) +
      theme(axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            plot.margin = unit(c(t = 1, r = -0.5, b = 1, l = 1), "mm"),
            plot.title = element_text(hjust = 0.5,
                                      size = 10,
                                      face = "bold")) +
      coord_flip() +
      scale_y_reverse(expand = c(NA,0)) +
      # scale_y_reverse(breaks = c(0, 0.25, 0.5, 0.75)) +
      guides(fill = "none")


    g2 <- region_means %>%
      # filter to year of interest
      dplyr::filter(YEAR == year) %>%
      # exclude occurrences of 0
      dplyr::filter(occurrence > set_occurrence) %>%
      # exclude CVs over 1
      dplyr::filter(CV < 1) %>%
      ggplot(data = .,
             aes(x = reorder(SPECIES_CD, occurrence),
                 y = CV*100,
                 fill = 'even')) +
      xlab(NULL) +
      geom_bar(stat = "identity",
               fill = "deepskyblue4") +
      ggtitle("Coefficient of Variation (CV) of Density") +
      theme_light() +
      scale_y_continuous(expand = c(0,0)) +
      theme(axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            plot.margin = unit(c(t = 1, r = 1, b = 1, l = -2), "mm"),
            plot.title = element_text(hjust = 0.5,
                                      size = 10,
                                      face = "bold")) +
      coord_flip() +
      guides(fill = "none") +
      geom_hline(yintercept=20, linetype="dashed", color = "black")
    # scale_fill_manual(values= c( "#58babb"))

    return(list(g.mid = g.mid, g1 = g1, g2 = g2))
  }



    if (region == "STX" || region == "STTSTJ" || region == "PRICO" || region == "FGB") {
      resulting_plot_parts <- plot(region_means = region_means, set_occurrence = 0.01, year = year, ptitle)
    } else {
      resulting_plot_parts <- plot(region_means, set_occurrence = 0.02, year, ptitle)
    }

    #use cowplot to plot grid
    cowplot::plot_grid(resulting_plot_parts$g1, resulting_plot_parts$g.mid, resulting_plot_parts$g2, ncol = 3,
                       rel_widths = c(3.5 / 9, 2 / 9, 3.5 / 9))

    #use ggsave to save plots
    ggsave(filename = paste(region_means$REGION[1], "Occ_v_CV.jpeg", sep = "_"),
           plot = gridExtra::grid.arrange(resulting_plot_parts$g1, resulting_plot_parts$g.mid, resulting_plot_parts$g2,
                                          ncol = 3, widths = c(3.5 / 9, 2 / 9, 3.5 / 9)),
           width = 9.8,
           height = 6.5,
           dpi = 300,
           units = "in",
           device = "jpg")


  region_means_cv20 <- region_means %>%
    dplyr::filter(CV <= .20)

  #### Export ####

  # Create list to export
  output <- list(
    "region_means" = region_means,
    "region_means_cv20" = region_means_cv20)

  return(output)

}
