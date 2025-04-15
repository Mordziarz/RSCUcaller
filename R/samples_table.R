#' @title Example Sample Metadata Table
#' 
#' @description 
#' Demonstration dataset containing sample identifiers and corresponding GenBank 
#' accession numbers for use with RSCUcaller functions.
#'
#' @format A data frame with 23 observations and 2 variables:
#' \describe{
#'   \item{ID}{Character vector containing unique sample identifiers. Must follow 
#'             "number_Name" format (e.g., "1_Human").}
#'   \item{GENBANK_ACCESSION}{Character vector of NCBI GenBank accession numbers 
#'                           matching sequence records.}
#' }
#'
#' @usage data(samples_table)
#'
#' @source Curated for demonstration purposes from publicly available NCBI data.
#'         Not intended for real analyses.
#'
#' @examples
#' # Load dataset
#' data(samples_table)
#' 
#' # View first 6 entries
#' head(samples_table)
#' 
#' # Use with prepare_fasta() function
#' \dontrun{
#' prepare_fasta(samples_table = samples_table, 
#'               path = "ncbi_sequences.fasta",
#'               file_out = "processed_sequences.fasta")
#' }
"samples_table"