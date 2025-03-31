#' Gets allocated species for CV and occurrence functions
#'
#'
#' @param region A string indicating the region. Options are: "SEFCRI", "FLK", "Tortugas", "STX", "STTSTJ", "PRICO", "FGB".
#' @return A vector of coral species for the specified region.
#' @export
#'
species_for_CV_and_occurrence <- function(region){

  coral_species_by_region <- switch(region,
                                    "STTSTJ" = c("Colpophyllia natans", "Diploria labyrinthiformis", "Madracis decactis", "Meandrina meandrites", "Montastraea cavernosa", "Orbicella annularis", "Orbicella faveolata", "Pseudodiploria strigosa", "Siderastrea siderea"),
                                    "STX" = c("Colpophyllia natans", "Dichocoenia stokesii", "Madracis decactis", "Montastraea cavernosa", "Orbicella annularis", "Orbicella franksi", "Pseudodiploria strigosa"),
                                    "PRICO" = c("Colpophyllia natans", "Diploria labyrinthiformis", "Madracis decactis", "Meandrina meandrites", "Montastraea cavernosa", "Orbicella annularis", "Orbicella faveolata", "Orbicella franksi", "Pseudodiploria strigosa"),
                                    "FLK" = c("Colpophyllia natans", "Montastraea cavernosa", "Orbicella faveolata", "Porites astreoides", "Siderastrea siderea", "Solenastrea bournoni"),
                                    "Tortugas" = c("Colpophyllia natans", "Montastraea cavernosa", "Orbicella faveolata", "Porites astreoides", "Orbicella franksi", "Stephanocoenia intersepta"),
                                    "SEFCRI" = c("Acropora cervicornis", "Dichocoenia stokesii", "Montastraea cavernosa", "Porites astreoides", "Pseudodiploria strigosa", "Siderastrea siderea"),
                                    "FGB" = c("Montastraea cavernosa", "Orbicella faveolata", "Orbicella franksi", "Siderastrea siderea", "Stephanocoenia intersepta", "Porites porites", "Agaricia agaricites", "Colpophyllia natans", "Mussa angulosa", "Agaricia fragilis", "Madracis auretenra", "Pseudodiploria strigosa",
                                              "Orbicella annularis", "Agaricia humilis", "Scolymia cubensis", "Agaricia lamarcki", "Tubastraea coccinea", "Madracis decactis", "Porites astreoides"))
  return(coral_species_by_region)
}
