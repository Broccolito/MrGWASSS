% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/load_data.R
\name{load_data}
\alias{load_data}
\title{Load GWAS Summary Statistics from file}
\usage{
load_data(
  file_name = "FHS_EA_MRS_5e8_snplist.txt",
  by_marker = TRUE,
  marker_name_column = "SNPID",
  chr_column = "CHR",
  pos_column = "POS",
  ref_column = "NEA",
  alt_column = "EA",
  pvalue_column = "p.value",
  delimiter = " "
)
}
\arguments{
\item{file_name}{The name of the file containing the statistics}

\item{by_marker}{A Boolean variable indicating whether the marker names will be used to specify the variant. If TRUE, marker_name_column has to be non-empty. If FALSE, chr_column, pos_column, ref_column, alt_column have to be non-empty}

\item{marker_name_column}{The column name of the 1KG markers, usually in the format of chr:pos:ref:alt}

\item{chr_column, }{The column name of the chromosome number column}

\item{pos_column, }{The column name of the chromosome position column}

\item{ref_column, }{The column name of the reference allele column}

\item{alt_column, }{The column name of the alternative allele column}

\item{pvalue_column, }{The column name of the p-value column}

\item{delimiter, }{The delimiter used in the file}
}
\value{
data, a data frame of the summary statistics
}
\description{
This function loads GWAS summary statistics from file
}
