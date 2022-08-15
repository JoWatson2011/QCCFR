#' Plot heatmap of Pearson Correlation
#'
#' @param normDat output of either normaliseData()
#' @param quantification One of "LFQ" or "SILAC", depending
#' which quantification method was used for experiment
#'
#' @return ggplot2 heatmap
#' @export
#'
#' @importFrom dplyr select .data
#' @importFrom tibble rownames_to_column
#' @importFrom tidyr pivot_longer
#' @importFrom stats cor
#' @importFrom ggplot2 ggplot aes geom_tile geom_text theme scale_x_discrete scale_y_discrete scale_fill_gradientn
#' @importFrom ggplot2 element_text element_blank margin unit
#' @importFrom magrittr `%$%`

cor_plot <- function(normDat,
                     quantification = c("LFQ", "SILAC")) {

  if(class(normDat) == "list" &
     all(names(normDat) %in% c("dat", "plots"))){
    normDat <- normDat$dat
  } else if(class(normDat) != "data.frame" |
            !"Gene names" %in% colnames(normDat)){
    stop("Please assign normalised data output from normaliseData() function OR
         a MaxQuant output table to the normDat argument.")
  }

  experiment_cols <- extract_experiment_cols(normDat,
                                             quantification)

  cors <- normDat %>%
    dplyr::select(all_of(experiment_cols)) %>%
    stats::cor(use = "complete.obs")

  if (quantification[1] == "LFQ") {
    rownames(cors) <-
      colnames(cors) <-
      gsub(
        "(LFQ )?[iI]ntensity ", "", colnames(cors)
           )
  } else {
    rownames(cors) <-
      colnames(cors) <-
      gsub("Ratio ", "",
           gsub(" normalized", "", colnames(cors))
      )
  }

  cor_hm <- cors %>%
    as.data.frame() %>%
    tibble::rownames_to_column("Exp1") %>%
    tidyr::pivot_longer(-.data$Exp1, names_to = "Exp2") %>%
    # mutate(Exp1 = factor(Exp1, levels = rownames(cors)),
    #        Exp2, factor(Exp2, levels = colnames(cors))
    #        ) %>%
    ggplot2::ggplot(ggplot2::aes(
      x = .data$Exp1,
      y = .data$Exp2,
      fill = .data$value
    )) +
    ggplot2::geom_tile() +
    ggplot2::geom_text(aes(label = round(.data$value, 2)),
                       size = 2) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        size = 5,
        angle = 45,
        hjust = 1
      ),
      axis.text.y = ggplot2::element_text(size = 5,
                                          angle = 45),
      axis.title = ggplot2::element_blank(),
      legend.key.size = ggplot2::unit(.25,
                                      "cm"),
      legend.margin = ggplot2::margin(0, 0, 0, 0,
                                      "cm"),
      plot.margin = ggplot2::margin(0, 0, 0, 0,
                                    "cm"),
      legend.title = ggplot2::element_text(size = 5),
      legend.text = ggplot2::element_text(size = 5)
    ) +
    ggplot2::scale_x_discrete(limits = rownames(cors)) +
    ggplot2::scale_y_discrete(limits = colnames(cors)) +
    ggplot2::scale_fill_gradientn(
      name = "R",
      colours = c("#2C7BB6",
                  "#ABD9E9",
                  "#FFFFBF",
                  "#FDAE61",
                  "#D7191C"),
      #   limits = c(0,1),
      breaks = seq(0, 1, 0.2)

    )

  return(cor_hm)
}

#' Principal component analysis
#'
#' @param normDat Output of either normalise_sty() or normalise_proteinGroups()
#' @param quantification One of "LFQ" or "SILAC", depending
#' which quantification method was used for experiment
#' @param rep_pattern Regex string to group replicates of same
#' experimental conditions
#' @param colour_grp_pattern  Regex string to group experimental conditions
#' - i.e. same stimulant, same time point
#' and projected as colour aesthetic in plot
#' @param shape_grp_pattern Regex string to group experimental conditions
#' i.e. same stimulant, same time point
#' and projected as shape aesthetic in plot.
#' @param colour_scale vector (can be named based on variables names in plot) of
#' colours
#'
#' @return ggplot2 object
#' @export
#'
#' @importFrom stats prcomp na.omit
#' @importFrom dplyr .data
#' @importFrom ggplot2 ggplot aes geom_line geom_point xlab ylab theme
#' @importFrom ggplot2 scale_x_continuous scale_y_continuous scale_colour_brewer scale_colour_manual
#' @importFrom ggplot2 element_blank element_line
pca_plot <- function(normDat,
                     quantification = c("LFQ", "SILAC"),
                     rep_pattern = "_0?[1-9]$",
                     colour_grp_pattern = NULL,
                     shape_grp_pattern = NULL,
                     colour_scale = NULL) {

  if(class(normDat) == "list" &
     all(names(normDat) %in% c("dat", "plots"))){
    normDat <- normDat$dat
  } else if(class(normDat) != "data.frame" |
            !"Gene names" %in% colnames(normDat)){
    stop("Please assign normalised data output from normaliseData() function OR
         a MaxQuant output table to the normDat argument.")
  }

  experiment_cols <- extract_experiment_cols(normDat,
                                             quantification)

  pca <-
    stats::prcomp(t(as.matrix(stats::na.omit(normDat[, experiment_cols]))))
  pca_df <- as.data.frame(pca$x)

  if (quantification[1] == "LFQ") {
    pca_df$experiment <-
      gsub(rep_pattern, "",
           gsub("(LFQ )?[iI]ntensity ", "",
                rownames(pca_df)))
  } else {
    pca_df$experiment <-
      gsub(rep_pattern, "",
           gsub("Ratio ", "",
                gsub(" normalized", "",
                     rownames(pca_df))
           )
      )
  }

  eigs <- pca$sdev ^ 2

  eigs_df <- data.frame(pc = 1:length(eigs),
                        eig = sapply(eigs, function(i) {
                          signif(100 * (i / sum(eigs)), 4)
                        }))

  if(!is.null(colour_grp_pattern)){
    vec <- regexpr(colour_grp_pattern,
                   pca_df$experiment)

    pca_df$colour_group <-
      substr(
        pca_df$experiment,
        vec,
        vec + attr(vec, "match.length") - 1

      )
  }

  if(!is.null(shape_grp_pattern)){
    vec <- regexpr(shape_grp_pattern,
                   pca_df$experiment)

    pca_df$shape_group <-
      substr(
        pca_df$experiment,
        vec,
        vec + attr(vec, "match.length") - 1

      )
  }

  scree <- ggplot2::ggplot(eigs_df,
                           ggplot2::aes(x = .data$pc, y = .data$eig)) +
    ggplot2::geom_line() +
    ggplot2::geom_point() +
    ggplot2::xlab("Principal component") +
    ggplot2::ylab("Eigenvalue")

  pca_12_g <-
    ggplot2::ggplot(pca_df,
                    ggplot2::aes(.data$PC1, .data$PC2,
                                 colour = {
                                   if(is.null(colour_grp_pattern))
                                     .data$experiment
                                   else
                                     .data$colour_group
                                 }
                    )
    ) +
    # ggplot2::geom_point(size = 3,
    # alpha = 0.7) #+
    ggplot2::scale_x_continuous(paste("PC1 (", eigs_df[1, 2], "%)")) +
    ggplot2::scale_y_continuous(paste("PC2 (", eigs_df[2, 2], "%)")) +
    ggplot2::theme(
      panel.background = ggplot2::element_blank(),
      panel.grid = ggplot2::element_line("grey",
                                         linetype = "dashed"),
      legend.key = ggplot2::element_blank(),
      axis.line = ggplot2::element_line("black"),
      legend.title = ggplot2::element_blank()
    )

  pca_23_g <- ggplot2::ggplot(pca_df,
                              ggplot2::aes(.data$PC2,
                                           .data$PC3,
                                           colour = {
                                             if(is.null(colour_grp_pattern))
                                               .data$experiment
                                             else
                                               .data$colour_group
                                           })) +
    ggplot2::scale_x_continuous(paste("PC2 (", eigs_df[2, 2], "%)")) +
    ggplot2::scale_y_continuous(paste("PC3 (", eigs_df[3, 2], "%)")) +
    ggplot2::theme(
      panel.background = ggplot2::element_blank(),
      panel.grid = ggplot2::element_line("grey", linetype = "dashed"),
      legend.key = ggplot2::element_blank(),
      axis.line = ggplot2::element_line("black"),
      legend.title = ggplot2::element_blank()
    )

  if(!is.null(shape_grp_pattern)){
    pca_12_g <- pca_12_g +
      ggplot2::geom_point(
        ggplot2::aes(shape = .data$shape_group),
        size = 3,
        alpha = 0.7)
    pca_23_g <- pca_23_g +
      ggplot2::geom_point(
        ggplot2::aes(shape = .data$shape_group),
        size = 3,
        alpha = 0.7)
  } else {
    pca_12_g <- pca_12_g +
      ggplot2::geom_point(size = 3,
                          alpha = 0.7)
    pca_23_g <- pca_23_g +
      ggplot2::geom_point(size = 3,
                          alpha = 0.7)
  }

  if(!is.null(colour_scale)){
    pca_12_g <- pca_12_g +
      ggplot2::scale_colour_manual(values = colour_scale)
    pca_23_g <- pca_23_g +
      ggplot2::scale_colour_manual(values = colour_scale)

  } #else {
  #   pca_12_g <- pca_12_g +
  #     ggplot2::scale_colour_brewer(palette = "Paired")
  #   pca_23_g <- pca_23_g +
  #     ggplot2::scale_colour_brewer(palette = "Paired")
  # }

  return(list(
    scree = scree,
    PC1_PC2 = pca_12_g,
    PC2_PC3 = pca_23_g
  ))
}

#' Plot no. quantified
#'
#' @param normDat Output of either normalise_sty() or normalise_proteinGroups()
#' @param quantification One of "LFQ" or "SILAC", depending
#' which quantification method was used for experiment
#' @param rep_pattern Regex string to group experimental conditions
#' - i.e. which string defines the replicates
#'
#' @return ggplot2 object
#' @export
#'
#' @importFrom dplyr select filter group_by summarise mutate ungroup
#' @importFrom dplyr .data n row_number
#' @importFrom tidyr pivot_longer
#' @importFrom ggplot2 ggplot aes geom_col theme
#' @importFrom ggplot2 element_text element_blank
#' @importFrom magrittr `%>%`
#' @importFrom stats reorder
count_plot <- function(normDat,
                       quantification = c("LFQ", "SILAC"),
                       rep_pattern = "_0?[1-9]$") {

  if(class(normDat) == "list" &
     all(names(normDat) %in% c("dat", "plots"))){
    normDat <- normDat$dat
  } else if(class(normDat) != "data.frame" |
            !"Gene names" %in% colnames(normDat)){
    stop("Please assign normalised data output from normaliseData() function OR
         a MaxQuant output table to the normDat argument.")
  }

  experiment_cols <- extract_experiment_cols(normDat,
                                             quantification)

  count <- normDat %>%
    dplyr::select(.data$id, all_of(experiment_cols)) %>%
    tidyr::pivot_longer(cols = -.data$id,
                        names_to = "experiment",
                        values_to = "ratio") %>%

    dplyr::filter(.data$ratio != 0) %>%
    dplyr::group_by(.data$experiment) %>%
    dplyr::summarise(n = n(), .groups = "keep") %>%
    dplyr::mutate(group =   {
      if (quantification[1] == "LFQ") {
        gsub(rep_pattern,
             "",
             gsub("(LFQ )?[iI]ntensity ", "",
                  .data$experiment))
      } else {
        gsub(rep_pattern,
             "",
             gsub("Ratio [HML]/[ML] normalized ", "",
                  .data$experiment))
      }
    }) %>%
    dplyr::ungroup() %>%
    unique() %>%
    dplyr::mutate(order = row_number()) %>%
    ggplot2::ggplot(ggplot2::aes(
      x = stats::reorder(.data$experiment,
                         .data$order),
      y = n
    )) +
    ggplot2::geom_col(ggplot2::aes(fill = .data$group),
                      position = "dodge") +
    ggplot2::theme(
      legend.position = "none",
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      axis.title.x = ggplot2::element_blank()
    )

  return(count)
}

#' Plot Mass Error Density
#'
#' @param evidence Output of read_evidence()
#' @param pattern Regex pattern of any experimental columns to exclude
#'
#' @return ggplot object
#' @export
#'
#' @importFrom dplyr filter select mutate
#' @importFrom ggplot2 geom_bin2d ggplot aes scale_fill_continuous
#' @importFrom ggplot2 ylab theme_bw
massErrorDensity <- function(evidence,
                             pattern = NULL) {
  print("Note: This plot may take several minutes to render")

  if (!is.null(pattern)) {
    df <- dplyr::filter(evidence,
                        grepl(pattern, .data$Experiment))
  }

  df <- evidence %>%
    dplyr::select(.data$`Mass error [ppm]`,
                  .data$Intensity) %>%
    dplyr::mutate(Intensity = log10(.data$Intensity))

  ppm_g <- ggplot2::ggplot(df,
                           ggplot2::aes(.data$`Mass error [ppm]`,
                                        .data$Intensity)) +
    ggplot2::geom_bin2d(bins = 200) +
    ggplot2::scale_fill_continuous(type = "viridis")  +
    ggplot2::theme_bw() +
    ggplot2::ylab("Log10(Intensity)")

  return(ppm_g)
}

#' Plot Andromeda scores
#'
#' @param evidence Output of read_evidence()
#' @param pattern Regex pattern of any experimental columns to exclude
#'
#' @return ggplot object
#' @export
#' @importFrom dplyr filter .data
#' @importFrom stats median
#' @importFrom graphics hist

andromedaScorePlot <- function(evidence,
                               pattern = NULL) {
  x <- evidence %>%
    dplyr::filter(!is.na(.data$Score))

  if (!is.null(pattern)) {
    x <- dplyr::filter(x,
                       grepl(pattern,
                             .data$Experiment))
  }

  med <- stats::median(x$Score)
  y_max <- graphics::hist(x$Score,
                          breaks = 30,
                          plot = F)

  andromScore <- ggplot2::ggplot(x, ggplot2::aes(.data$Score)) +
    ggplot2::geom_histogram(bins = length(y_max$breaks)) +
    ggplot2::geom_vline(xintercept = med, colour = "red") +
    ggplot2::theme_bw() +
    ggplot2::xlab("Andromeda Score") +
    ggplot2::ylab("Count") +
    ggplot2::annotate(geom = "text",
                      y= max(y_max$counts),
                      x = med + (med/100*75),
                      colour = "red",
                      label = paste("Median = ", med)
    )

  return(andromScore)
}
