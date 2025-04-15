#' Calculate Relative Synonymous Codon Usage (RSCU) for multiple sequences
#'
#' Computes RSCU values for a set of DNA sequences provided as a single FASTA file
#' (prepared using \code{\link{prepare_fasta}}) or as a list of sequences.
#'
#' @param merged_sequences Either a character string specifying the path to a FASTA file
#'   containing prepared sequences (see \code{\link{prepare_fasta}}), or a list of DNA sequences.
#'
#' @return A \code{data.frame} containing RSCU values for all codons and all samples.
#'   The columns include:
#'   \describe{
#'     \item{AA}{Amino acid}
#'     \item{codon}{Codon sequence}
#'     \item{eff}{Codon efficiency (if available)}
#'     \item{freq}{Codon frequency}
#'     \item{RSCU}{Relative Synonymous Codon Usage value}
#'     \item{Col}{Internal grouping column}
#'     \item{index}{Unique index for each codon/sample combination}
#'     \item{Species}{Sample name}
#'   }
#'
#' @details
#' The function calculates RSCU values for each codon in each sample. Input can be either:
#' \itemize{
#'   \item A FASTA file with multiple sequences, prepared using \code{\link{prepare_fasta}}
#'   \item A list of DNA sequences (as returned by \code{seqinr::read.fasta})
#' }
#' The sample names should follow the convention "number_name" (e.g., "1_Human").
#'
#' @examples
#' # Using a prepared FASTA file
#' data("prepared_fasta", package = "RSCUcaller")
#' rscu_df <- get_RSCU(merged_sequences = prepared_fasta)
#'
#' # Using a file path
#' # rscu_df <- get_RSCU(merged_sequences = "your_prepared_sequences.fasta")
#'
#' @references
#' Sharp, P. M., & Li, W. H. (1986). An evolutionary perspective on synonymous codon usage in unicellular organisms. Journal of Molecular Evolution, 24(1-2), 28-38.
#'
#' @seealso \code{\link{prepare_fasta}}, \code{\link{get_RSCU_other}}, \code{\link{get_matrix}}
#'
#' @importFrom dplyr filter group_by mutate
#' @importFrom magrittr %>%
#' @export

get_RSCU <- function(merged_sequences = ""){

  if (base::missing(merged_sequences)) {
    stop("The merged_sequences predictions are required. Please provide a valid argument.",
         call. = FALSE)
  }

  base::message(base::paste0("Loading data from "))

  if (base::is.list(merged_sequences)){

    merged_seq <- merged_sequences

    Rscu_all <- base::data.frame(row.names = 1, AA = NA, codon = NA, eff = NA, freq = NA, RSCU = NA, Col = NA, index = NA, Species = NA)
    t <- base::gsub("^\\d+_", "",names(merged_seq))

    for (i in 1:base::length(merged_seq)) {
      a <- merged_seq[[i]]
      cp1 <- seqinr::uco(a, as.data.frame = TRUE)
      cp1$Col <- 1
      cp1$index <- base::paste(cp1$AA, t[i])
      cp1$Species <- t[i]
      Rscu_all <- base::rbind(Rscu_all,cp1)
    }

    Rscu_all <- Rscu_all %>% dplyr::filter(AA != "NA") %>%
      dplyr::group_by(index) %>%
      dplyr::mutate(Col= stats::lag(base::cumsum(Col), default = 0)) %>%
      dplyr::mutate(Col=base::as.factor(Col))

  }
  else {
    if(grepl(".fasta$|.txt$", merged_sequences, ignore.case = TRUE)){

      merged_seq <- seqinr::read.fasta(file = merged_sequences ,set.attributes = TRUE, seqtype = "DNA",as.string = FALSE)

      Rscu_all <- base::data.frame(row.names = 1, AA = NA, codon = NA, eff = NA, freq = NA, RSCU = NA, Col = NA, index = NA, Species = NA)
      t <- base::gsub("^\\d+_", "",names(merged_seq))

      for (i in 1:base::length(merged_seq)) {
        a <- merged_seq[[i]]
        cp1 <- seqinr::uco(a, as.data.frame = TRUE)
        cp1$Col <- 1
        cp1$index <- base::paste(cp1$AA, t[i])
        cp1$Species <- t[i]
        Rscu_all <- base::rbind(Rscu_all,cp1)
      }

      Rscu_all <- Rscu_all %>% dplyr::filter(AA != "NA") %>%
        dplyr::group_by(index) %>%
        dplyr::mutate(Col= stats::lag(base::cumsum(Col), default = 0)) %>%
        dplyr::mutate(Col=base::as.factor(Col))
    }
  }
  base::message(base::paste0("Success"))
  return(Rscu_all)
}
