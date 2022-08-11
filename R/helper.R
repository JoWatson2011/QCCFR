#' Extract columns containing quantitative data
#'
#' @param dat MQ data frame
#' @param quantification One of "LFQ" or "SILAC", depending
#' which quantificiation method was used for experiment
#' @param pattern Regex string to Subset of experimental columns,
#' e.g. just proteome or STY columns. Default = NULL.
#'
#' @return vector of  columns containing experimental values
#'
extract_experiment_cols <- function(dat,
                                    quantification = c("LFQ", "SILAC"),
                                    pattern = NULL) {
  if (quantification[1] == "LFQ") {
    experiment_cols <-
      grep("(LFQ )?[iI]ntensity", colnames(dat), value = T)
  } else {
    experiment_cols <-
      grep("Ratio [^m]", colnames(dat), value = T)
  }

  if (!is.null(pattern)) {
    experiment_cols <- grep(pattern, experiment_cols, value = T)
  }

  if(length(experiment_cols) == 0){
    stop("Columns containing quantitative values not found.
         Did you select correct quantification type? LFQ or SILAC")
  }

  return(experiment_cols)
}


#' plot pie charts for modification charactersitics
#'
#' @param dat 2 column data frame containing name of pie slices
#' and n of each
#'
#' @return ggplot2 object
#'
#' @importFrom ggplot2 ggplot aes geom_rect geom_text
#' @importFrom ggplot2 coord_polar xlim theme_void theme element_text
#' @importFrom dplyr mutate .data
#' @importFrom utils head
plotPie <- function(dat) {
  g <- dat %>%
    dplyr::mutate(
      Fraction = .data$n / sum(.data$n),
      # Compute percentages
      ymax = cumsum(.data$Fraction),
      # Compute the cumulative percentages (top of each rectangle)
      ymin = c(0,
               utils::head(.data$ymax,
                    n = -1)),
      # Compute the bottom of each rectangle
      labelPosition = (.data$ymax + .data$ymin) / 2,
      # Compute label position
      Percentage = .data$Fraction * 100,
      # Compute a good label
      Percentage = format(round(.data$Percentage,
                                1),
                          nsmall = 1),
      label = paste0(.data$name,
                     "\n (",
                     .data$Percentage,
                     "%)")
    ) %>%
    ggplot2::ggplot(ggplot2::aes(
      ymax = .data$ymax,
      ymin = .data$ymin,
      xmax = 1.2,
      xmin = 0.2,
      fill = .data$name,
      show.legend = T
    )) +
    ggplot2::geom_rect(alpha = 0.6) +
    ggplot2::geom_text(
      x = 3,
      ggplot2::aes(y = .data$labelPosition,
                   label = .data$label,
                   color = .data$name),
      size = 2,
      check_overlap = F,
      # x here controls label position (inner / outer)
      nudge_x = -0.2
    ) +
    ggplot2::coord_polar(theta = "y") +
    ggplot2::xlim(c(-2, 3)) +
    ggplot2::theme_void() +
    ggplot2::theme(legend.position = "none",
                   plot.title = ggplot2::element_text(hjust = 0.5))

  return(g)
}
