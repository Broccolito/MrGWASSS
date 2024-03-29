% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/annotate_data.R
\name{annotate_data}
\alias{annotate_data}
\title{Annotate Summary Statistics}
\usage{
annotate_data(data, file_name = "data_annotated.xlsx", write_file = TRUE)
}
\arguments{
\item{data}{an data frame object loaded from load_data function}

\item{file_name}{The name of the file where the annotated results are saved}

\item{write_file}{A boolean value indicating whether the annotated results need to be saved}
}
\value{
data_annotated, a data frame of the annotated summary statistics
}
\description{
This function annotates GWAS summary statistics using the Biothings API
}
