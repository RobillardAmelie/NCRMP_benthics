#' Update Protection Status in Dataset
#' 
#' ## Tag: data analysis


# output gets called by:
# 
#

# NCRMP Caribbean Benthic analytics team: Davis, Groves, Viehman, Williams
# Last update: Jan 2025
#'
#' This function updates the protection status (`PROT`) in the given dataset (`data`)
#' by joining it with a reference grid (`grid_df`). Specific hardcoded adjustments 
#' are made for known `PRIMARY_SAMPLE_UNIT` values.
#'
#' @param data A dataframe containing survey data with a `MAPGRID_NR` column.
#' @param grid_df A dataframe containing protection status (`PROT`) for each `MAPGRID_NR`.
#'
#' @return The updated dataframe with the modified `PROT` column.
#'
update_protection_status <- function(data, grid_df ) {

  new_prots <- grid_df %>%
    dplyr::select(MAPGRID_NR, PROT) %>%
    dplyr::rename("PROT_og" = PROT) %>%
    dplyr::mutate(
      MAPGRID_NR = as.numeric(MAPGRID_NR),
      PROT_og = as.numeric(PROT_og)
    )

  data <- data %>%
    dplyr::left_join(.,new_prots, by = "MAPGRID_NR") %>%
    dplyr::mutate(
      PROT_og = case_when(
        PRIMARY_SAMPLE_UNIT == 1006 ~ 0,
        PRIMARY_SAMPLE_UNIT == 1382 ~ 1,
        TRUE ~ PROT_og
      )
    ) %>%
    dplyr::select(-PROT) %>%
    dplyr::rename("PROT" = PROT_og)
}



