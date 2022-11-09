load_input <- function(input, sep="\t", encoding = 'unknown', model_name=NA) {
  if (is.data.frame(input)) {
    long_names <- colnames(input)[4:ncol(input)]
    data <- data.table::data.table(input)
    data[,1] <- fill_the_blanks(data[,1], missing=NA)
    n <- sum(data[,3], na.rm = TRUE)
  } else {
    # Else: default -- input_file is a .txt file with the ot-soft format
    # If no model name provided, use filename sans extension
    if (is.na(model_name)) {
      model_name <- tools::file_path_sans_ext(basename(input))
    }
    input <- load_data_otsoft(input, sep = sep, encoding = encoding)
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
load_data_otsoft <- function(infile, sep = "\t", encoding = 'unknown') {
  in.dt <- data.table::fread(
    infile, header = FALSE, sep = sep, fill = TRUE, encoding = encoding
  )

  # TODO: Some validation of format?
  full_names <- in.dt[1, 4:ncol(in.dt)]
  abbr_names <- in.dt[2, 4:ncol(in.dt)]

  in_data <- in.dt[3:nrow(in.dt),]
  in_data[,1] <- fill_the_blanks(in_data[,1])

  n <- sum(in_data[,3], na.rm = TRUE)

  return(list(full_names = full_names, abbr_names = abbr_names,
              data = in_data, n=n))
}

# Loads bias file in OTSoft format
load_bias_file_otsoft <- function(infile, sep = "\t") {
  in.dt <- data.table::fread(infile, header = FALSE, sep = sep)
  bias_params <- in.dt[,2:3]
  names(bias_params) <- c("mus", "sigmas")
  return(bias_params)
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
