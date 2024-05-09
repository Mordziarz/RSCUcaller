#' Calculating RSCU for multiple sequences
#'
#' @param merged_sequences path to the file with prepared sequences via the prepare_fasta() function
#'
#' @return A `data.frame` object containing Rscu_all.
#' @export
#'

get_RSCU <- function(merged_sequences = ""){

  if (base::missing(merged_sequences)) {
    stop("The merged_sequences predictions are required. Please provide a valid argument.",
         call. = FALSE)
  }

  base::message(base::paste0("Loading data from ", merged_sequences))

  if (base::is.list(merged_sequences)){

    merged_seq <- merged_sequences

    Rscu_all <- base::data.frame(row.names = 1, AA = NA, codon = NA, eff = NA, freq = NA, RSCU = NA, Name = NA, Col = NA, index = NA, Species = NA)
    t <- base::gsub("^\\d+_", "",names(merged_seq))
    n <- base::gsub("_.*","",names(merged_seq))

    for (i in 1:base::length(merged_seq)) {
      a <- merged_seq[[i]]
      cp1 <- seqinr::uco(a, as.data.frame = TRUE)
      cp1$Name <- n[i]
      cp1$Col <- 1
      cp1$index <- base::paste(cp1$AA, n[i])
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

      merged_seq <- seqinr::read.fasta(file = merged_sequences ,set.attributes = T, seqtype = "DNA",as.string = F)

      Rscu_all <- base::data.frame(row.names = 1, AA = NA, codon = NA, eff = NA, freq = NA, RSCU = NA, Name = NA, Col = NA, index = NA, Species = NA)
      t <- base::gsub("^\\d+_", "",names(merged_seq))
      n <- base::gsub("_.*","",names(merged_seq))

      for (i in 1:base::length(merged_seq)) {
        a <- merged_seq[[i]]
        cp1 <- seqinr::uco(a, as.data.frame = TRUE)
        cp1$Name <- n[i]
        cp1$Col <- 1
        cp1$index <- base::paste(cp1$AA, n[i])
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
