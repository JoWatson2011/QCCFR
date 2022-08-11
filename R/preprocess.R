#' Filter MaxQuant Phospho (STY).txt
#'
#' @param sty output of read_phosphoSTY
#' @param quantification one of LFQ or SILAC, depending on type
#' of quantification used in experiment.
#'
#' @return filtered STY dataframe
#' @export
#' @importFrom dplyr select filter .data all_of
#' @importFrom magrittr `%>%`
filter_sty <- function(sty,
                       quantification = c("LFQ", "SILAC")
)
{
  if(!any(c("Amino acid", "Sequence window",
            "Position in peptide", "Localisation prob") %in%
          colnames(sty))){
    stop("Please use output of read_phosphoSTY()")
  }

  experiment_cols <- extract_experiment_cols(sty, quantification)

  # Summary of variables Phospho (STY)Sites.txt will be filtered on.
  # After running the script this information will be stored in the
  # variable named report.

  report <- vector()
  report["Total Sites Identified"] <- nrow(sty)
  report["Potential Contaminants"] <- sum(!is.na(sty$`Potential contaminant`))
  report["Reverse"] <- sum(!is.na(sty$Reverse))
  report["Localisation probability < 0.75"] <- sum(sty$`Localization prob` < 0.75)

  # Filter Phospho (STY)Sites.txt to remove contaminants, reverse, and low localisation
  # probability sites. Localisation probability threshold can be changed here.
  sty_flt <- sty %>%
    dplyr::filter(is.na(.data$`Potential contaminant`),
                  is.na(.data$Reverse),
                  .data$`Localization prob` > 0.75) %>%
    dplyr::mutate(
      dplyr::across(
        dplyr::all_of(experiment_cols),
        ~ ifelse(.x == 0, NA, .x)
      )
    )

  # Add info to report
  report["Incomplete cases"] <-
    nrow(
      sty_flt[rowSums(is.na(sty_flt[experiment_cols])) == length(experiment_cols),]
    )

  # Filter rows with no quantitative info
  sty_flt_cc <-
    sty_flt[rowSums(!is.na(sty_flt[experiment_cols])) >= 1,]

  report["Sites remaining following filtering"] <- nrow(sty_flt_cc)

  print(report)
  return(sty_flt_cc)
}

#' Filter MaxQuant proteinGroups.txt
#'
#' @param proteinGroups output of read_proteinGroups
#' @param quantification one of LFQ or SILAC, depending on type
#' of quantification used in experiment.
#'
#' @return filtered proteinGroups dataframe
#' @export
#' @importFrom dplyr select filter .data all_of
#' @importFrom magrittr `%>%`
filter_proteinGroups <- function(proteinGroups,
                                 quantification = c("LFQ", "SILAC")
){

  if(!any(c("Majority protein IDs", "Razor + unique peptides",
            "Unique + razor sequence coverage [%]") %in%
          colnames(proteinGroups))){
    stop("Please use output of read_proteinGroups()")
  }

  experiment_cols <- extract_experiment_cols(proteinGroups, quantification)

  report <- vector()
  report["Total Proteins Identified"] <- nrow(proteinGroups)
  report["Potential Contaminants"] <- sum(!is.na(proteinGroups$`Potential contaminant`))
  report["Reverse"] <- sum(!is.na(proteinGroups$Reverse ))
  report["Only Identified by Site"] <- sum(!is.na(proteinGroups$`Only identified by site`))

  proteinGroups_flt <- proteinGroups %>%
    dplyr::filter(is.na(.data$`Potential contaminant`),
                  is.na(.data$`Reverse`),
                  is.na(.data$`Only identified by site`),
                  .data$`Razor + unique peptides` > 1,
                  .data$`Unique + razor sequence coverage [%]` >= 5
    )%>%
    dplyr::mutate(
      dplyr::across(
        dplyr::all_of(experiment_cols),
        ~ ifelse(.x == 0, NA, .x)
      )
    )

  report["Incomplete cases"] <- nrow(
    proteinGroups_flt[
      rowSums(!is.na(proteinGroups_flt[experiment_cols])) == length(experiment_cols),
    ]
  )

  pg_flt_cc <-
    proteinGroups_flt[rowSums(!is.na(proteinGroups_flt[experiment_cols])) >= 1,]

  report["Proteins remaining following filtering"] <- nrow(pg_flt_cc)

  print(report)
  return(pg_flt_cc)
}



#' Normalise quantitative values
#'
#' @param filteredData output of filter_sty() or filter_proteinGroups()
#' @param quantification one of LFQ or SILAC, depending on type
#' of quantification used in experiment.
#'
#' @return Normalised dataframe
#' @export
#' @importFrom limma normalizeQuantiles
#' @importFrom dplyr select mutate across everything all_of
#' @importFrom ggplot2 ggplot aes geom_density theme ggtitle
#' @importFrom tidyr pivot_longer
#' @importFrom magrittr `%$%`
#' @importFrom dplyr .data

normaliseData <- function(filteredData,
                          quantification = c("LFQ", "SILAC")){

  if(class(filteredData) != "data.frame"){
    stop("Please pass data frame output of filter_proteinGroups or filter_sty")
  }

  if(nrow(filteredData) == 0){
    stop("Data frame is empty")
  }

  experiment_cols <- extract_experiment_cols(filteredData, quantification)

  newDat <- filteredData %>%
    mutate(across(all_of(experiment_cols), log10))


  # Normailse LFQ values
  # SILAC internally normalised by MQ
  if(quantification[1] == "LFQ"){
    newDat <- cbind(
      filteredData[,!(names(filteredData) %in% c(experiment_cols))],
      limma::normalizeQuantiles(as.matrix(newDat[,experiment_cols]))
    )
  }

  # Visualise
  raw_g <- filteredData %>%
    dplyr::select(.data$id, dplyr::all_of(experiment_cols)) %>%
    dplyr::mutate(
      dplyr::across(
        -.data$id,
        ~ ifelse(.x == 0, NA, .x)
      )
    ) %>%
    dplyr::mutate(
      dplyr::across(
        dplyr::all_of(experiment_cols),
        log2)
    ) %>%
    tidyr::pivot_longer(cols = -.data$id,
                        names_to = "experiment",
                        values_to = "intensity") %>%
    dplyr::mutate(rep = substr(.data$experiment,
                               nchar(.data$experiment)-1,
                               nchar(.data$experiment))) %>%
    ggplot2::ggplot(
      ggplot2::aes(x = .data$intensity,
                   color = .data$rep,
                   group = .data$experiment)
    ) +
    ggplot2::geom_density() +
    ggplot2::theme(legend.position = 'none') +
    ggplot2::ggtitle("Raw Values") +
    ggplot2::xlab(ifelse(quantification[1] == "LFQ",
                         "Intensity",
                         "Normalised SILAC Ratio")
    )

  norm_g <- newDat %>%
    as.data.frame() %>%
    dplyr::select(.data$id,
                  dplyr::all_of(experiment_cols)) %>%
    tidyr::pivot_longer(cols = -c(.data$id),
                        names_to = "experiment",
                        values_to = "intensity") %>%
    dplyr::mutate(rep = substr(.data$experiment,
                               nchar(.data$experiment)-1,
                               nchar(.data$experiment))) %>%
    ggplot2::ggplot(
      ggplot2::aes(x = .data$intensity,
                   color = .data$rep,
                   group = .data$experiment)
    ) +
    ggplot2::geom_density() +
    ggplot2::theme(legend.position = 'none') +
    ggplot2::ggtitle(ifelse(quantification[1] == "LFQ",
                            "Normalised intensities",
                            "Transformed SILAC Ratio")
    ) +
    ggplot2::xlab(ifelse(quantification[1] == "LFQ",
                         "Intensity",
                         "Normalised SILAC Ratio")
                  )


  return(
    list(
      dat = newDat,
      plots = list(
        raw = raw_g,
        normalised = norm_g
      )
    )
  )

}

#' Impute Data by drawing randomly from a distribution
#' @description The default represents a stimulation to represent a
#' population below the limits of detection, as described in
#'  Tyanova et al., 2016
#'  ONLY suitable for label free data!!
#'
#' @param filteredData output of filter_sty() or filter_proteinGroups()
#'
#' @return imputed values
#' @importFrom msm rtnorm
#' @importFrom stats sd
#' @export

imputeTruncNorm <- function(filteredData){

  if(class(filteredData) != "data.frame"){
    stop("Please pass data frame output of filter_proteinGroups or filter_sty")
  }

  if(nrow(filteredData) == 0){
    stop("Data frame is empty")
  }

  experiment_cols <- extract_experiment_cols(filteredData)
  dat <- filteredData

  dat[,experiment_cols] <-
    apply(filteredData[,experiment_cols], 2, function(v){
      dis <- msm::rtnorm(length(v),
                         lower = -(sd(v, na.rm = T)),
                         upper = (-(sd(v, na.rm = T)) + 1.0)
      )
      v[is.na(v)] <- sample(dis, length(v[is.na(v)]), replace = F)
      return(v)
    })

  return(dat)
}

