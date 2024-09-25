#' Preparing fasta files for further processing
#'
#' @param samples_table table with sequence descriptions
#' @param path path to the original sequence file
#' @param file_out path to the file with processed sequences
#'
#' @return None (invisible)
#' @export
#'

prepare_fasta <- function(samples_table=samples_table,path="",file_out="") {

  if (base::missing(samples_table)) {
    stop("samples_table argument is missing. Please provide a valid argument.",
         call. = FALSE)
  }

  base::message(base::paste0("Loading data"))

  sequences_list <- base::vector("list", base::nrow(samples_table))

  if (base::all(colnames(samples_table)==c("sequence_path","sample_name"))){
    for (i in 1:nrow(samples_table)) {
      sequence <- seqinr::read.fasta(file = samples_table$sequence_path[i],set.attributes = T, seqtype = "DNA",as.string = T)
      for (j in 1:length(sequence)) {
        sequences_list[[i]] <- c(sequences_list[[i]], paste0(sequence[[j]]))
      }
      base::names(sequences_list)[i] <- samples_table$sample_name[i]
    }
  }

  if (base::grepl(".fasta$|.txt$", path, ignore.case = TRUE)){
    sequence <- seqinr::read.fasta(file = path, set.attributes = T, seqtype = "DNA",as.string = T)
    for (i in 1:nrow(samples_table)) {
      positions <- base::which(grepl(samples_table$GENBANK_ACCESSION[i], base::names(sequence), ignore.case = TRUE))
      a <- sequence[positions]
      for (j in 1:length(a)) {
        sequences_list[[i]] <- c(sequences_list[[i]], paste0(sequence[[j]]))
      }
      base::names(sequences_list)[i] <- samples_table$ID[i]
    }

  }
  sequences_list <- base::as.list(sequences_list)
  seqinr::write.fasta(sequences = sequences_list,names = names(sequences_list),file.out = file_out)

  base::message(base::paste0("Success"))

}
