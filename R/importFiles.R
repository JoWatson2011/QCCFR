#' Import MQ output files
#'
#' @param folderPath Path of MQ txt output files


read_experiments_names <- function(folderPath){

  if(!"summary.txt" %in% list.files(folderPath)){
    stop("Please ensure MaxQuant output table summary.txt
         is in specified path location")
  }

  experiments <- readr::read_tsv(
    paste0(folderPath, "/summary.txt")
  )$Experiment

  experiments_cols <- unique(
    experiments[experiments != "" & !is.na(experiments)]
  )

  return(experiments_cols)
}

#' Import MQ output files
#'
#' @param folderPath Path of MQ txt output files
#' @param pattern Regex string to Subset of experimental columns,
#' e.g. just proteome or STY columns. Default = NULL.
#' @param quantification One of "LFQ" or "SILAC", depending
#' which quantificiation method was used for experiment
#' @param SILAC_ratios If quantification is SILAC, which ratios to include?
#' Defaults to "all".
#'
#' @return data in format compatible for further analysis
#' @export
#' @importFrom readr read_tsv
#' @importFrom dplyr all_of

read_proteinGroups <- function(
    folderPath,
    pattern = NULL,
    quantification = c("LFQ", "SILAC"),
    SILAC_ratios = c("all",
                     "H/L and M/L",
                     "H/L",
                     "M/L"
    )
){

  experiment_cols <- read_experiments_names(folderPath)

  if(!is.null(pattern)){
    experiment_cols <- experiment_cols[grepl(pattern, experiment_cols)]
  }

  if(quantification[1] == "LFQ"){
    experiment_cols <- paste("LFQ intensity", experiment_cols)
  } else {

    SILAC_ratios_val <-
      ifelse(
        SILAC_ratios[1] == "all",
        c("Ratio H/L normalized ",
          "Ratio M/L normalized ",
          "Ratio H/M normalized "),
        ifelse(
          SILAC_ratios[1] == "H/L and M/L",
          c("Ratio H/L normalized ",
            "Ratio M/L normalized "),
          ifelse(
            SILAC_ratios[1] == "H/L" | SILAC_ratios[1] == "M/L",
            paste("Ratio", SILAC_ratios[1],"normalized"),
            ""
          )
        )
      )
    experiment_cols <-
      apply(
        expand.grid(SILAC_ratios_val,
                    experiment_cols),
        1,
        function(i) paste(i, collapse = "")
      )
  }

  if(!"proteinGroups.txt" %in% list.files(folderPath)){
    stop("Please ensure MaxQuant output table proteinGroups.txt
         is in specified path location")
  }

  pg <- readr::read_tsv(
    paste0(folderPath, "/proteinGroups.txt"),
    col_select =
      dplyr::all_of(
        c("id", "Majority protein IDs", "Sequence length",
          "Gene names", "Protein names", "Reverse",
          "Potential contaminant", "Fasta headers",
          "Razor + unique peptides",
          "Unique + razor sequence coverage [%]",
          "Only identified by site",
          experiment_cols
        )
      )
  )

  return(pg)
}

#' @export
#' @rdname read_proteinGroups
read_phosphoSTY <- function(
    folderPath,
    pattern = NULL,
    quantification = c("LFQ", "SILAC"),
    SILAC_ratios = c("all",
                     "H/L and M/L",
                     "H/L",
                     "M/L"
    )
){

  experiment_cols <- read_experiments_names(folderPath)

  if(!is.null(pattern)){
    experiment_cols <- experiment_cols[grepl(pattern, experiment_cols)]
  }

  if(quantification[1] == "LFQ"){
    experiment_cols <- paste("Intensity", experiment_cols)
  } else {


    if(SILAC_ratios[1] == "all"){
      SILAC_ratios_val <- c("Ratio H/L normalized ",
                            "Ratio M/L normalized ",
                            "Ratio H/M normalized ")
    }else if(SILAC_ratios[1] == "H/L and M/L"){
      SILAC_ratios_val <- c("Ratio H/L normalized ",
                            "Ratio M/L normalized ")
    }else if(SILAC_ratios[1] == "H/L" | SILAC_ratios[1] == "M/L"){
      SILAC_ratios_val <- paste("Ratio", SILAC_ratios[1],"normalized")
    }else{
      SILAC_ratios_val <-""
    }

    experiment_cols <-
      apply(
        expand.grid(SILAC_ratios_val,
                    experiment_cols),
        1,
        function(i) paste(i, collapse = "")
      )
  }

  if(!"Phospho (STY)Sites.txt" %in% list.files(folderPath)){
    stop("Please ensure MaxQuant output table Phospho (STY)Sites.txt
         is in specified path location")
  }

  sty <- readr::read_tsv(
    paste0(folderPath, "/Phospho (STY)Sites.txt"),
    col_select =
      dplyr::all_of(
        c("id", "Protein", "Protein names",
          "Gene names",	"Amino acid", "Position",
          "Sequence window", "Reverse", "Potential contaminant",
          "Localization prob", "Mod. peptide IDs", "Fasta headers",
          experiment_cols
        )
      )
  )

  return(sty)
}

#' @export
#' @rdname read_proteinGroups
read_modSpPeptides <- function(
    folderPath,
    quantification = c("LFQ", "SILAC"),
    SILAC_ratios = c("all",
                     "H/L and M/L",
                     "H/L",
                     "M/L")
){
  experiment_names <- read_experiments_names(folderPath)

  if(quantification[1] == "LFQ"){
    experiment_cols <- paste("Intensity", experiment_names)
  } else {
    if(SILAC_ratios[1] == "all"){
      SILAC_ratios_val <- c("Ratio H/L normalized ",
                            "Ratio M/L normalized ",
                            "Ratio H/M normalized ")
    }else if(SILAC_ratios[1] == "H/L and M/L"){
      SILAC_ratios_val <- c("Ratio H/L normalized ",
                            "Ratio M/L normalized ")
    }else if(SILAC_ratios[1] == "H/L" | SILAC_ratios[1] == "M/L"){
      SILAC_ratios_val <- paste("Ratio", SILAC_ratios[1],"normalized")
    }else{
      SILAC_ratios_val <-""
    }

    experiment_cols <-
      apply(
        expand.grid(SILAC_ratios_val,
                    experiment_names),
        1,
        function(i) paste(i, collapse = "")
      )
  }

  modsp <- readr::read_tsv(
    paste0(folderPath, "/modificationSpecificPeptides.txt"),
    col_select =
      dplyr::all_of(
        c("Sequence", "Mass", "Mass Fractional Part",
          "Protein Groups", "Proteins",	"Gene Names",
          "Protein Names", "Unique (Groups)",
          "Unique (Proteins)", "Phospho (STY)",
          paste("Experiment", experiment_names),
          experiment_cols,
          "Reverse", "Potential contaminant",
          "id", "Phospho (STY) site IDs"
        )
      )
  )

  return(modsp)
}

#' @export
#' @rdname read_proteinGroups
read_peptides <- function(
    folderPath
){
  experiment_cols <- read_experiments_names(folderPath)

  peptides <- readr::read_tsv(
    paste0(folderPath, "/peptides.txt"),
    col_select =
      dplyr::all_of(
        c("id", "Sequence", "Length", "Mass", "Proteins",
          "Potential contaminant", "Reverse",
          paste("Experiment", experiment_cols)
        )
      )
  )

  return(peptides)
}

#' @export
#' @rdname read_proteinGroups
read_evidence <- function(
    folderPath
){
  experiment_cols <- read_experiments_names(folderPath)

  evidence <- readr::read_tsv(
    paste0(folderPath, "/evidence.txt"),
    col_select =
      dplyr::all_of(
        c("Mass error [ppm]", "Intensity", "Score",
          "Sequence", "Gene names",
          "Experiment",
          "Phospho (STY)",
          "Phospho (STY) site IDs")
      )
  )

  return(evidence)
}
