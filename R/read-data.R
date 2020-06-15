#' @export
load_data_otsoft <- function(infile, sep = "\t") {
  in.dt <- data.table::fread(infile, header = FALSE, sep = sep)

  # TODO: Some validation of format?
  full_names <- in.dt[1, 4:ncol(in.dt)]
  abbr_names <- in.dt[2, 4:ncol(in.dt)]

  candidate_rows <- which(in.dt$V1 != "")
  candidate_entries <- split(
    in.dt,
    cumsum(1:nrow(in.dt) %in% which(in.dt[,1] != ""))
  )
  candidate_entries <- candidate_entries[2:length(candidate_entries)]
  list(full_names, abbr_names, candidate_entries)
}

#
# setwd("C:/Users/conno/Dropbox/ling/maxent_r/maxent.ot/R/test-files")
#
# bar <- load_data_otsoft("SampleDataFile.txt")
#
