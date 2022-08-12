#' Pie charts of modification characteristics
#'
#' @description Distribution of characteristics of modifified
#' peptides, e.g. modified residue, localisation probability
#'
#' @param raw_modTable Output of read_phosphoSTY() or similar
#' @param modsp  Output of read_modSpPeptides()
#' @param quantification One of "LFQ" or "SILAC", depending
#' which quantification method was used for experiment
#'
#' @return three ggplot objects as list
#' @export
#'
#' @importFrom dplyr select all_of filter mutate group_by
#' @importFrom dplyr summarise ungroup n .data
#' @importFrom tidyr pivot_longer
#'
modificationCharacteristics <- function(raw_modTable,
                                        modsp,
                                        quantification = c("LFQ", "SILAC")) {

  if(!any(c("Amino acid", "Sequence window",
            "Position in peptide", "Localisation prob") %in%
          colnames(raw_modTable))){
    stop("Please use output of read_phosphoSTY()")
  }

  if(!any(c("Sequence", "Mass", "Mass Fractional Part" ) %in%
          colnames(modsp))){
    stop("Please use output of read_modSpPeptides()")
  }


  experiment_cols <- extract_experiment_cols(raw_modTable, quantification)

  locProb <-  raw_modTable %>%
    dplyr::select(.data$id,
                  .data$`Localization prob`,
                  dplyr::all_of(experiment_cols)) %>%
    tidyr::pivot_longer(-c(.data$id, .data$`Localization prob`)) %>%
    dplyr::filter({
      if (any(grepl("[Ii]ntensity", experiment_cols))) {
        .data$value > 0
      } else {
        !is.na(.data$value)
      }
    }) %>%
    dplyr::mutate(name =
                    ifelse(.data$`Localization prob` > 0.75,
                           "Class I", "Class II")) %>%
    dplyr::select(-.data$value) %>%
    unique() %>%
    dplyr::group_by(.data$name) %>%
    dplyr::summarise(n = dplyr::n(),
                     .groups = "keep") %>%
    dplyr::ungroup() %>%
    plotPie()

  residue <- raw_modTable %>%
    dplyr::select(.data$id,
                  .data$`Amino acid`,
                  dplyr::all_of(experiment_cols)) %>%
    # mutate(across(2:length(phospho), as.numeric)) %>%
    tidyr::pivot_longer(-c(.data$id, .data$`Amino acid`)) %>%
    dplyr::filter({
      if (any(grepl("[Ii]ntensity", experiment_cols))) {
        .data$value > 0
      } else {
        !is.na(.data$value)
      }

    }) %>%
    dplyr::select(-.data$value) %>%
    dplyr::mutate(name = .data$`Amino acid`) %>%
    unique() %>%
    dplyr::group_by(.data$name) %>%
    dplyr::summarise(n = dplyr::n(), .groups = "keep") %>%
    dplyr::ungroup() %>%
    plotPie()


  numberMods <- modsp %>%
    dplyr::select(.data$id,
                  .data$`Phospho (STY)`,
                  dplyr::all_of(experiment_cols)) %>%
    # mutate(across(2:length(phospho), as.numeric)) %>%
    tidyr::pivot_longer(-c(.data$id,
                           .data$`Phospho (STY)`)) %>%
    dplyr::filter({
      if (any(grepl("[Ii]ntensity", experiment_cols))) {
        .data$value > 0
      } else {
        !is.na(.data$value)
      }

    },
    .data$`Phospho (STY)` > 0) %>%
    dplyr::mutate(name =
                    ifelse(
                      .data$`Phospho (STY)` == 1,
                      "Single",
                      ifelse(
                        .data$`Phospho (STY)` == 2,
                        "Double",
                        ifelse(.data$`Phospho (STY)` > 2, ">Double", NA)
                      )
                    )) %>%
    dplyr::select(-.data$value) %>%
    unique() %>%
    dplyr::group_by(.data$name) %>%
    dplyr::summarise(n = dplyr::n(), .groups = "keep") %>%
    dplyr::ungroup() %>%
    plotPie()

  return(list(
    locProb = locProb,
    residue = residue,
    numberMods = numberMods
  ))

}

#' Calculate enrichment efficiency
#'
#' @param modsp Output of read_modSpPeptides()
#' @param quantification One of "LFQ" or "SILAC", depending
#' which quantificiation method was used for experiment
#' @param pattern Regex string to Subset of experimental columns,
#' e.g. just proteome or STY columns. Default = NULL.
#'
#' @return numeric value, percentage peptides with STY site identified
#' @export
#'
#' @importFrom dplyr filter if_all all_of .data
phosphoenrichment <- function(modsp,
                              quantification = c("LFQ", "SILAC"),
                              pattern = NULL){
  experiments <- extract_experiment_cols(modsp, quantification[1], pattern)
  experiments <- gsub("Intensity ", "Experiment ", experiments)

  modsp_flt <- modsp %>%
    dplyr::filter(is.na(.data$`Potential contaminant`),
                  is.na(.data$`Reverse`)) %>%
    dplyr::filter(
      !dplyr::if_all(
        dplyr::all_of(experiments),
        is.na
      ) |
        !dplyr::if_all(
          dplyr::all_of(experiments),
          ~ .x == 0
        )
    )

  phosphoenrichment <-
    nrow(modsp_flt[modsp_flt$`Phospho (STY)` > 0 ,]) /
    nrow(modsp_flt)

  return(phosphoenrichment)
}


