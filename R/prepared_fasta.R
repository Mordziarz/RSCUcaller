#' @title Example Processed FASTA Sequences for RSCU Analysis
#' 
#' @description 
#' Demonstration dataset containing processed mitochondrial gene sequences in FASTA format,
#' ready for analysis with RSCUcaller functions. Sequences are properly formatted and annotated
#' for use with \code{\link{get_RSCU}} and related functions.
#'
#' @format A list object of class \code{"fasta"} containing:
#' \describe{
#'   \item{Sequence names}{Character vector following "number_Species" format (e.g. "1_Homo_sapiens")}
#'   \item{Sequences}{DNA character strings representing coding sequences}
#' }
#' 
#' @details
#' This dataset serves as a template for proper input structure requirements. Each sequence:
#' \itemize{
#'   \item Starts with a valid start codon (ATG)
#'   \item Ends with a proper stop codon
#'   \item Contains only standard nucleotide characters (A/T/C/G)
#'   \item Maintains reading frame continuity
#' }
#'
#' @usage data(prepared_fasta)
#'
#' @source Curated example data derived from public NCBI resources. Contains simulated sequences
#'         for demonstration purposes only.
#'
#' @examples
#' # Load dataset
#' data("prepared_fasta")
#' 
#' # View structure
#' str(prepared_fasta[1:2])
#' 
#' # Use with get_RSCU function
#' \dontrun{
#' rscu_results <- get_RSCU(merged_sequences = prepared_fasta)
#' }
#' 
#' @seealso
#' Related functions:
#' \code{\link{prepare_fasta}}, \code{\link{get_RSCU}}
"prepared_fasta"