# Loads tableaux in OTSoft format
load_data_otsoft <- function(infile, sep = "\t", encoding = 'unknown') {
  in.dt <- data.table::fread(
    infile, header = FALSE, sep = sep, fill = TRUE, encoding = encoding
  )

  # TODO: Some validation of format?
  full_names <- in.dt[1, 4:ncol(in.dt)]
  abbr_names <- in.dt[2, 4:ncol(in.dt)]
  data <- in.dt[3:nrow(in.dt),]
  data[,1] <- fill_the_blanks(data[,1])
  # candidate_rows <- which(in.dt$V1 != "")
  # candidate_entries <- split(
  #   in.dt,
  #   cumsum(1:nrow(in.dt) %in% which(in.dt[,1] != ""))
  # )
  # candidate_entries <- candidate_entries[2:length(candidate_entries)]
  n <- sum(data[,3], na.rm = TRUE)
  return(list(full_names = full_names, abbr_names = abbr_names,
              data = data, n=n))
}

# Loads bias file in OTSoft format
load_bias_file_otsoft <- function(infile, sep = "\t") {
  in.dt <- data.table::fread(infile, header = FALSE, sep = sep)
  bias_params <- in.dt[,2:3]
  names(bias_params) <- c("mus", "sigmas")
  return(bias_params)
}

fill_the_blanks <- function(x, missing = ""){
  rle <- rle(as.character(x))
  empty <- which(rle$value == missing)
  rle$values[empty] <- rle$value[empty - 1]
  inverse.rle(rle)
}
