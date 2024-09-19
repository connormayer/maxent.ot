#' Converts an OTSoft tableaux file to a data frame
#'
#' Loads an OTSoft tableaux file and converts it to the data frame format used by
#' the maxent.ot functions.
#'
#' @param input The path to the input data file.
#'   This should contain more OT tableaux consisting of
#'   mappings between underlying and surface forms with observed frequency and
#'   violation profiles. Constraint violations must be numeric.
#'
#'   The file should be in OTSoft format.
#'   For examples of OTSoft format, see inst/extdata/sample_data_file.txt.
#' @param output_path (optional) A string specifying the path to a file to
#'   which the data frame will be saved in CSV format. If the file exists it
#'   will be overwritten. If this argument isn't provided, the output will not
#'   be written to a file.
#' @param encoding (optional) The character encoding of the input file. Defaults
#'  to "unknown".
#' @return A data frame corresponding to the input OTSoft tableau, containing
#' the columns
#' \itemize{
#'  \item `Input`: The input form.
#'  \item `Output`: The output form.
#'  \item `Frequency`: The frequency of the input/output mapping.
#'  \item One column for each constraint containing its violation counts.
#' }
#' @examples
#'   # Convert OTSoft file to data frame format
#'   otsoft_file <- system.file(
#'       "extdata", "sample_data_file_otsoft.txt", package = "maxent.ot"
#'   )
#'   df_output <- otsoft_tableaux_to_df(otsoft_file)
#'
#'   # Save data frame to a file
#'   tmp_output <- tempfile()
#'   otsoft_tableaux_to_df(otsoft_file, tmp_output)
#' @export
otsoft_tableaux_to_df <- function(input, output_path=NA, encoding='unknown') {
  results <- load_input(input, encoding)
  output_df <- data.frame(results$data)
  output_df[output_df == ''] <- 0
  if (!is.na(output_path)) {
    utils::write.csv(output_df, file=output_path, row.names=FALSE, quote=FALSE)
  }
  return(output_df)
}

#' Converts an OTSoft bias file to a data frame
#'
#' Loads an OTSoft bias file and converts it to the data frame format used by
#' the maxent.ot functions.
#'
#' @param input The path to the input bias file. This should contain more
#' OT tableaux consisting of mappings between underlying and surface forms with
#' observed frequency and violation profiles. Constraint violations must be
#' numeric.
#'
#' The file should be in OTSoft format. For examples of OTSoft format, see
#' inst/extdata/sample_bias_file_otsoft.txt.
#' @param output_path (optional) A string specifying the path to a file to
#'   which the data frame will be saved in CSV format. If the file exists it
#'   will be overwritten. If this argument isn't provided, the output will not
#'   be written to a file.
#' @return A data frame corresponding to the input OTSoft bias file, containing
#' the columns
#' \itemize{
#'  \item `Constraint`: The constraint name.
#'  \item `Mu`: The mu value for the regularization term.
#'  \item `Sigma`: The sigma value for the regularization term.
#' }
#' @examples
#'   # Convert OTSoft bias file to data frame format
#'   otsoft_file <- system.file(
#'       "extdata", "sample_bias_file_otsoft.txt", package = "maxent.ot"
#'   )
#'   df_output <- otsoft_bias_to_df(otsoft_file)
#'
#'   # Save data frame to a file
#'   tmp_output <- tempfile()
#'   otsoft_bias_to_df(otsoft_file, tmp_output)
#' @export
otsoft_bias_to_df <- function(input, output_path=NA) {
  bias_df <- load_bias_file_otsoft(input)
  if (!is.na(output_path)) {
    utils::write.csv(bias_df, file=output_path, row.names=FALSE, quote=FALSE)
  }
  return(load_bias_file_otsoft(input))
}

load_input <- function(input, encoding = 'unknown', model_name=NA) {
  if (is.data.frame(input)) {
    long_names <- colnames(input)[4:ncol(input)]
    data <- data.table::data.table(input)
    data[is.na(data)] <- 0
    data[,1] <- fill_the_blanks(data[,1])
    n <- sum(data[,3], na.rm = TRUE)
  } else {
    # Else: default -- input_file is a .txt file with the ot-soft format
    # If no model name provided, use filename sans extension
    if (is.na(model_name)) {
      model_name <- tools::file_path_sans_ext(basename(input))
    }
    input <- load_data_otsoft(input, encoding = encoding)
    long_names <- input$full_names
    data <- input$data
    colnames(data) <- c("Input", "Output", "Frequency", unlist(long_names))
    n <- input$n
  }
  result = list(
    data = data,
    long_names = long_names,
    n = n,
    model_name = model_name
  )
  return(result)
}

# Loads tableaux in OTSoft format
load_data_otsoft <- function(infile, encoding = 'unknown') {
  in.dt <- data.table::fread(
    infile, header = FALSE, sep = '\t', fill = TRUE, encoding = encoding
  )

  # Data should minimally have four columns: Input, Output, Frequency, and
  # at least one constraint.
  num_cols <- ncol(in.dt)
  if (num_cols < 4) {
    stop(sprintf("Input data only has %s columns: have you passed the
                  correct separator character into `in_sep`? For example, if
                  your data is comma separated, include the argument
                  `in_sep=',' in your call to `optimize_weights`.", num_cols))
  }

  full_names <- in.dt[1, 4:ncol(in.dt)]
  abbr_names <- in.dt[2, 4:ncol(in.dt)]

  if (any(sapply(full_names, is.numeric)) | any(sapply(abbr_names, is.numeric))) {
    stop("One or more of the constraint names is numeric. Have you forgotten
          either the long names or abbreviated names row?")
  }

  in_data <- in.dt[3:nrow(in.dt),]
  in_data[,1] <- fill_the_blanks(in_data[,1])

  n <- sum(in_data[,3], na.rm = TRUE)

  return(list(full_names = full_names, abbr_names = abbr_names,
              data = in_data, n=n))
}

# Loads bias file in OTSoft format
load_bias_file_otsoft <- function(infile) {
  in.dt <- data.table::fread(infile, header = FALSE, sep = '\t')
  names(in.dt) <- c("Constraint", "Mu", "Sigma")
  return(in.dt)
}

# Replace empty cells in a column with the closest cell above that is not empty
fill_the_blanks <- function(x, missing = ""){
  enc <- base::rle(as.character(x))
  if (is.na(missing)) {
    empty <- which(is.na(enc$value))
  }
  else {
    empty <- which(enc$value == missing)
  }
  enc$values[empty] <- enc$values[empty - 1]
  base::inverse.rle(enc)
}
