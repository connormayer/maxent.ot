% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read-data.R
\name{otsoft_tableaux_to_df}
\alias{otsoft_tableaux_to_df}
\title{Converts an OTSoft tableaux file to a data frame}
\usage{
otsoft_tableaux_to_df(input, output_path = NA, encoding = "unknown")
}
\arguments{
\item{input}{The path to the input data file.
This should contain more OT tableaux consisting of
mappings between underlying and surface forms with observed frequency and
violation profiles. Constraint violations must be numeric.

The file should be in OTSoft format.
For examples of OTSoft format, see inst/extdata/sample_data_file.txt.}

\item{output_path}{(optional) A string specifying the path to a file to
which the data frame will be saved in CSV format. If the file exists it
will be overwritten. If this argument isn't provided, the output will not
be written to a file.}

\item{encoding}{(optional) The character encoding of the input file. Defaults
to "unknown".}
}
\value{
A data frame corresponding to the input OTSoft tableau, containing
the columns
\itemize{
\item \code{Input}: The input form.
\item \code{Output}: The output form.
\item \code{Frequency}: The frequency of the input/output mapping.
\item One column for each constraint containing its violation counts.
}
}
\description{
Loads an OTSoft tableaux file and converts it to the data frame format used by
the maxent.ot functions.
}
\examples{
  # Convert OTSoft file to data frame format
  otsoft_file <- system.file(
      "extdata", "sample_data_file_otsoft.txt", package = "maxent.ot"
  )
  df_output <- otsoft_tableaux_to_df(otsoft_file)

  # Save data frame to a file
  tmp_output <- tempfile()
  otsoft_tableaux_to_df(otsoft_file, tmp_output)
}
