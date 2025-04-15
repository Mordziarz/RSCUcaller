#' Calculate RSCU for Multiple Sequences with Custom Genetic Codes
#'
#' Computes Relative Synonymous Codon Usage (RSCU) values for a set of DNA sequences,
#' allowing each sample to use a different genetic code (codon table).
#'
#' @param merged_sequences Either a character string specifying the path to a FASTA file
#'   containing prepared sequences (see \code{\link{prepare_fasta}}), or a list of DNA sequences.
#' @param samples_table A data frame describing the samples. Must contain either columns
#'   \code{sample_name} or \code{ID} (matching sequence names), and a \code{codon_table_id}
#'   column specifying the genetic code for each sample (see \code{\link{get_codon_table}}).
#' @param pseudo_count A numeric value added to codon counts to avoid division by zero (default: 1).
#'
#' @return A \code{data.frame} with columns: \code{AA} (amino acid), \code{codon}, \code{eff} (codon count),
#'   \code{RSCU} (Relative Synonymous Codon Usage), and sample/grouping columns.
#'
#' @details
#' This function enables RSCU calculation for datasets where each sample may use a different
#' genetic code (e.g., mitochondrial vs. nuclear). The \code{samples_table} must specify
#' the genetic code for each sample using the \code{codon_table_id} column.
#'
#' @examples
#' prepared_fasta <- system.file("extdata/prepared_fasta.fasta", package = "RSCUcaller")
#' samples_path <- system.file("extdata/samples_table.csv", package = "RSCUcaller")
#' samples_table <- read.csv2(samples_path,sep = ";")
#' samples_table$codon_table_id <- 2
#' rscu_df <- get_RSCU_other2(
#'   merged_sequences = prepared_fasta,
#'   samples_table = samples_table,
#'   pseudo_count = 1
#' )
#'
#' @seealso \code{\link{get_RSCU}}, \code{\link{get_RSCU_other}}, \code{\link{get_codon_table}}
#' @importFrom dplyr filter group_by mutate
#' @importFrom magrittr %>%
#' @importFrom stats lag
#' @export


get_RSCU_other2 <- function(merged_sequences = "",pseudo_count=1,samples_table=samples_table){
  
  if (base::missing(merged_sequences)) {
    stop("The merged_sequences predictions are required. Please provide a valid argument.",
         call. = FALSE)
  }
  
  if (base::missing(samples_table)) {
    stop("The samples_table predictions are required. Please provide a valid argument.",
         call. = FALSE)
  }
  
  samples_table$codon_table_id <- base::as.numeric(samples_table$codon_table_id)

  if (base::is.list(merged_sequences)){
    
    merged_seq <- merged_sequences

        t <- base::gsub("^\\d+_", "",names(merged_seq))
    
    data_sequence <- RSCUcaller::seq_to_data_frame(merged_sequences)

            if ("sample_name" %in% base::colnames(samples_table)) {
      data_sequence <- base::merge(data_sequence, samples_table, by.x = "names", by.y = "sample_name")
    } else if ("ID" %in% base::colnames(samples_table)) {
      data_sequence <- base::merge(data_sequence, samples_table, by.x = "names", by.y = "ID")
    } else {
      stop("Column 'sample_name' or 'ID' not found in samples_table.")
    }
    
    Rscu_all <- base::data.frame(row.names = 1, AA = NA, codon = NA, eff = NA, RSCU = NA, Col = NA, index = NA, Species = NA)
    
    for (i in 1:base::nrow(data_sequence)) {
      
      rscu_data <- RSCUcaller::calculate_rscu(data_sequence$sequences[i], codon_table_id = data_sequence$codon_table_id[i], pseudo_count = pseudo_count)
      rscu_data$Col <- 1
      rscu_data$index <- base::paste(rscu_data$AA, t[i])
      rscu_data$Species <- t[i]
      Rscu_all <- base::rbind(Rscu_all,rscu_data)
    }
    
    Rscu_all <- Rscu_all %>% dplyr::filter(AA != "NA") %>%
      dplyr::group_by(index) %>%
      dplyr::mutate(Col= stats::lag(base::cumsum(Col), default = 0)) %>%
      dplyr::mutate(Col=base::as.factor(Col))
    
    base::message(base::paste0("Success"))
    return(Rscu_all)

  }
    if(grepl(".fasta$|.txt$", merged_sequences, ignore.case = TRUE)){
    
  base::message(base::paste0("Loading data from ", merged_sequences))
  
  data_sequence <- RSCUcaller::seq_to_data_frame(merged_sequences)

              if ("sample_name" %in% base::colnames(samples_table)) {
      data_sequence <- base::merge(data_sequence, samples_table, by.x = "names", by.y = "sample_name")
    } else if ("ID" %in% base::colnames(samples_table)) {
      data_sequence <- base::merge(data_sequence, samples_table, by.x = "names", by.y = "ID")
    } else {
      stop("Column 'sample_name' or 'ID' not found in samples_table.")
    }

           t <- base::gsub("^\\d+_", "",data_sequence$names)
  
  Rscu_all <- base::data.frame(row.names = 1, AA = NA, codon = NA, eff = NA, RSCU = NA, Col = NA, index = NA, Species = NA)
  
  for (i in 1:nrow(data_sequence)) {
    
    rscu_data <- RSCUcaller::calculate_rscu(data_sequence$sequences[i], codon_table_id = data_sequence$codon_table_id[i], pseudo_count = pseudo_count)
    rscu_data$Col <- 1
    rscu_data$index <- base::paste(rscu_data$AA, t[i])
    rscu_data$Species <- t[i]
    Rscu_all <- base::rbind(Rscu_all,rscu_data)
  }
  
  Rscu_all <- Rscu_all %>% dplyr::filter(AA != "NA") %>%
    dplyr::group_by(index) %>%
    dplyr::mutate(Col= stats::lag(base::cumsum(Col), default = 0)) %>%
    dplyr::mutate(Col=base::as.factor(Col))

  base::message(base::paste0("Success"))
  return(Rscu_all)
  }
}

#' Calculate RSCU for Multiple Sequences with a Custom Genetic Code
#'
#' Computes Relative Synonymous Codon Usage (RSCU) values for a set of DNA sequences
#' using a specified genetic code (codon table).
#'
#' @param merged_sequences Either a character string specifying the path to a FASTA file
#'   containing prepared sequences (see \code{\link{prepare_fasta}}), or a list of DNA sequences.
#' @param codon_table_id An integer specifying the codon table to use (see \code{\link{get_codon_table}}).
#'   Common values: 1 = Standard, 2 = Vertebrate Mitochondrial, etc.
#' @param pseudo_count A numeric value added to codon counts to avoid division by zero (default: 1).
#'
#' @return A \code{data.frame} with columns: \code{AA} (amino acid), \code{codon}, \code{eff} (codon count),
#'   \code{RSCU} (Relative Synonymous Codon Usage), and sample/grouping columns.
#'
#' @details
#' This function is useful for RSCU analysis of genomes using non-standard genetic codes,
#' such as mitochondrial or alternative nuclear codes.
#'
#' @examples
#' prepared_fasta <- system.file("extdata/prepared_fasta.fasta", package = "RSCUcaller")
#' rscu_df <- get_RSCU_other(
#'   merged_sequences = prepared_fasta,
#'   codon_table_id = 2,
#'   pseudo_count = 0.5
#' )
#'
#' @seealso \code{\link{get_RSCU}}, \code{\link{get_RSCU_other2}}, \code{\link{get_codon_table}}
#' @importFrom dplyr filter group_by mutate
#' @importFrom magrittr %>%
#' @importFrom stats lag
#' @importFrom seqinr read.fasta s2c
#' @export

get_RSCU_other <- function(merged_sequences = "",codon_table_id=1,pseudo_count=1){
  
  if (base::missing(merged_sequences)) {
    stop("The merged_sequences predictions are required. Please provide a valid argument.",
         call. = FALSE)
  }
  
  if (base::is.list(merged_sequences)){
    
    merged_seq <- merged_sequences

    t <- base::gsub("^\\d+_", "",names(merged_seq))
    
    data_sequence <- RSCUcaller::seq_to_data_frame(merged_seq)
    
    Rscu_all <- base::data.frame(row.names = 1, AA = NA, codon = NA, eff = NA, RSCU = NA, Col = NA, index = NA, Species = NA)
    
    for (i in 1:base::nrow(data_sequence)) {
      
      rscu_data <- RSCUcaller::calculate_rscu(data_sequence$sequences[i], codon_table_id = codon_table_id, pseudo_count = pseudo_count)
      rscu_data$Col <- 1
      rscu_data$index <- base::paste(rscu_data$AA, t[i])
      rscu_data$Species <- t[i]
      Rscu_all <- base::rbind(Rscu_all,rscu_data)
    }
    
    Rscu_all <- Rscu_all %>% dplyr::filter(AA != "NA") %>%
      dplyr::group_by(index) %>%
      dplyr::mutate(Col= stats::lag(base::cumsum(Col), default = 0)) %>%
      dplyr::mutate(Col=base::as.factor(Col))
    
    base::message(base::paste0("Success"))
    return(Rscu_all)
    
  }
      if(grepl(".fasta$|.txt$", merged_sequences, ignore.case = TRUE)){
  
  base::message(base::paste0("Loading data from ", merged_sequences))

  merged_seq <- merged_sequences

  data_sequence <- RSCUcaller::seq_to_data_frame(merged_seq)

  t <- base::gsub("^\\d+_", "",data_sequence$names)

  Rscu_all <- base::data.frame(row.names = 1, AA = NA, codon = NA, eff = NA, RSCU = NA, Col = NA, index = NA, Species = NA)
  
  for (i in 1:base::nrow(data_sequence)) {
    
    rscu_data <- RSCUcaller::calculate_rscu(data_sequence$sequences[i], codon_table_id = codon_table_id, pseudo_count = pseudo_count)
    rscu_data$Col <- 1
    rscu_data$index <- base::paste(rscu_data$AA, t[i])
    rscu_data$Species <- t[i]
    Rscu_all <- base::rbind(Rscu_all,rscu_data)
  }
  
  Rscu_all <- Rscu_all %>% dplyr::filter(AA != "NA") %>%
    dplyr::group_by(index) %>%
    dplyr::mutate(Col= stats::lag(base::cumsum(Col), default = 0)) %>%
    dplyr::mutate(Col=base::as.factor(Col))

  base::message(base::paste0("Success"))
  return(Rscu_all)
    }
}


#' Convert FASTA Sequences to a Data Frame
#'
#' Converts a set of DNA sequences (from a FASTA file or list) into a data frame
#' with sequence names and sequence strings.
#'
#' @param merged_sequences Either a character string specifying the path to a FASTA file
#'   (prepared using \code{\link{prepare_fasta}}), or a list of DNA sequences.
#'
#' @return A \code{data.frame} with columns: \code{names} (sequence names), \code{sequences} (DNA strings).
#'
#' @examples
#' prepared_fasta <- system.file("extdata/prepared_fasta.fasta", package = "RSCUcaller")
#' seq_df <- seq_to_data_frame(prepared_fasta)
#'
#' @seealso \code{\link{prepare_fasta}}
#' @importFrom dplyr filter
#' @importFrom seqinr read.fasta
#' @export

seq_to_data_frame <- function(merged_sequences = ""){
  
  if (base::missing(merged_sequences)) {
    stop("The merged_sequences predictions are required. Please provide a valid argument.", 
         call. = FALSE)
  }
  
   if (base::is.list(merged_sequences)){

  merged_seq <- merged_sequences

  for (i in base::seq_along(merged_seq)) {
    base::attr(merged_seq[[i]], "name") <- base::names(merged_seq[i])
    base::attr(merged_seq[[i]], "seq") <- merged_seq[[i]][1]
  }
  
  base::attr(merged_seq, "name") <- base::lapply(merged_seq, function(x) base::attr(x, "name", exact = FALSE))
  base::attr(merged_seq, "seq") <- base::lapply(merged_seq, function(x) base::attr(x, "seq", exact = FALSE))
  
  merged_seq <- base::as.data.frame(merged_seq)
  merged_seq <- base::t(merged_seq)
  merged_seq <- base::as.data.frame(merged_seq)
  merged_seq <- merged_seq %>% dplyr::filter(V1 != "none")
  merged_seq$names <- base::rownames(merged_seq)
  merged_seq$sequences <- merged_seq$V1
  merged_seq <- merged_seq[,c("names","sequences")]
  merged_seq$names <- base::gsub("^X", "", merged_seq$names)
  base::rownames(merged_seq) <- 1:base::nrow(merged_seq)

  return(merged_seq)
}


  if(grepl(".fasta$|.txt$", merged_sequences, ignore.case = TRUE)){
        
          merged_seq <- seqinr::read.fasta(merged_sequences,seqtype ="DNA" ,as.string = TRUE)

  
  for (i in base::seq_along(merged_seq)) {
    base::attr(merged_seq[[i]], "name") <- base::names(merged_seq[i])
    base::attr(merged_seq[[i]], "seq") <- merged_seq[[i]][1]
  }
  
  base::attr(merged_seq, "name") <- base::lapply(merged_seq, function(x) base::attr(x, "name", exact = FALSE))
  base::attr(merged_seq, "seq") <- base::lapply(merged_seq, function(x) base::attr(x, "seq", exact = FALSE))
  
  merged_seq <- base::as.data.frame(merged_seq)
  merged_seq <- base::t(merged_seq)
  merged_seq <- base::as.data.frame(merged_seq)
  merged_seq <- merged_seq %>% dplyr::filter(V1 != "none")
  merged_seq$names <- base::rownames(merged_seq)
  merged_seq$sequences <- merged_seq$V1
  merged_seq <- merged_seq[,c("names","sequences")]
  merged_seq$names <- base::gsub("^X", "", merged_seq$names)
  base::rownames(merged_seq) <- 1:base::nrow(merged_seq)
  
  return(merged_seq)
  }
}

#' Calculate RSCU for a Single Sequence
#'
#' Calculates Relative Synonymous Codon Usage (RSCU) values for a single DNA sequence
#' using a specified genetic code (codon table).
#'
#' @param nucleotide_input A character string representing the nucleotide sequence.
#' @param codon_table_id An integer specifying the codon table to use (see \code{\link{get_codon_table}}).
#' @param pseudo_count A numeric value added to codon counts to avoid division by zero (default: 1).
#'
#' @return A \code{data.frame} with columns: \code{AA} (amino acid), \code{codon}, \code{eff} (codon count),
#'   \code{RSCU} (Relative Synonymous Codon Usage).
#'
#' @examples
#' seq <- "ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG"
#' rscu_df <- calculate_rscu(seq, codon_table_id = 1, pseudo_count = 1)
#'
#' @seealso \code{\link{get_codon_table}}
#' @importFrom seqinr s2c uco
#' @importFrom stats setNames
#' @export

calculate_rscu <- function(nucleotide_input, codon_table_id = 1, pseudo_count = 1) {
  
  if (base::missing(nucleotide_input)) {
    stop("The nucleotide_input predictions are required. Please provide a valid argument.", 
         call. = FALSE)
  }
  
  nucleotide_string <- nucleotide_input
  nucleotide_string <- base::tolower(nucleotide_string)
  if (base::nchar(nucleotide_string) == 0) stop("The sequence does not contain valid nucleotides.")
  
  seq_chars <- seqinr::s2c(nucleotide_string)
  codon_counts <- seqinr::uco(seq_chars, frame = 0)
  
  codon_table <- RSCUcaller::get_codon_table(codon_table_id)
  genetic_code <- codon_table$genetic_code
  
  rscu_values <- stats::setNames(base::rep(NA, base::length(codon_counts)), base::names(codon_counts))
  
  for (aa in base::unique(genetic_code)) {
    
    aa_codons <- base::names(genetic_code)[genetic_code == aa]
    n_i <- base::length(aa_codons)
    
    aa_counts <- base::sapply(aa_codons, function(codon) {
      ifelse(codon %in% base::names(codon_counts), codon_counts[codon], 0)
    })
    total <- base::sum(aa_counts) + pseudo_count * n_i
    rscu_values[aa_codons] <- (n_i * (aa_counts + pseudo_count)) / total
  }
  
  result_df <- base::data.frame(
    AA = genetic_code[base::names(codon_counts)],
    codon = base::names(codon_counts),
    eff = base::as.vector(codon_counts),
    RSCU = rscu_values,
    stringsAsFactors = FALSE
  )
  
  return(result_df)
}



#' Retrieve a Codon Table for a Genetic Code
#'
#' Returns the codon-to-amino acid mapping for a specified genetic code (codon table).
#'
#' @param codon_table_id An integer specifying the codon table. See vignette or NCBI documentation
#'   for available codes (e.g., 1 = Standard, 2 = Vertebrate Mitochondrial, etc.).
#'
#' @return A list with one element: \code{genetic_code}, a named character vector mapping codons to amino acids.
#'
#' @examples
#' get_codon_table(1) # Standard genetic code
#' get_codon_table(2) # Vertebrate mitochondrial code
#'
#' @seealso \code{\link{get_RSCU}}, \code{\link{get_RSCU_other}}, \code{\link{get_RSCU_other2}}
#' @export

get_codon_table <- function(codon_table_id) {
  
  #Standard
  genetic_code <- c(
    "ttt" = "Phe", "ttc" = "Phe", "tta" = "Leu", "ttg" = "Leu",
    "ctt" = "Leu", "ctc" = "Leu", "cta" = "Leu", "ctg" = "Leu",
    "att" = "Ile", "atc" = "Ile", "ata" = "Ile", "atg" = "Met",
    "gtt" = "Val", "gtc" = "Val", "gta" = "Val", "gtg" = "Val",
    "tct" = "Ser", "tcc" = "Ser", "tca" = "Ser", "tcg" = "Ser",
    "cct" = "Pro", "ccc" = "Pro", "cca" = "Pro", "ccg" = "Pro",
    "act" = "Thr", "acc" = "Thr", "aca" = "Thr", "acg" = "Thr",
    "gct" = "Ala", "gcc" = "Ala", "gca" = "Ala", "gcg" = "Ala",
    "tat" = "Tyr", "tac" = "Tyr", "taa" = "Stp", "tag" = "Stp",
    "cat" = "His", "cac" = "His", "caa" = "Gln", "cag" = "Gln",
    "aat" = "Asn", "aac" = "Asn", "aaa" = "Lys", "aag" = "Lys",
    "gat" = "Asp", "gac" = "Asp", "gaa" = "Glu", "gag" = "Glu",
    "tgt" = "Cys", "tgc" = "Cys", "tga" = "Stp", "tgg" = "Trp",
    "cgt" = "Arg", "cgc" = "Arg", "cga" = "Arg", "cgg" = "Arg",
    "agt" = "Ser", "agc" = "Ser", "aga" = "Arg", "agg" = "Arg",
    "ggt" = "Gly", "ggc" = "Gly", "gga" = "Gly", "ggg" = "Gly"
  ) 
  
  if (codon_table_id == 2) {  # Vertebrate Mitochondrial
    genetic_code[c("ata", "aga", "agg", "tga")] <- c("Met", "Stp", "Stp", "Trp")
  } 
  else if (codon_table_id == 3) { # Yeast Mitochondrial
    genetic_code[c("ata", "ctt", "ctc", "cta", "ctg", "tga", "cga", "cgc")] <- 
      c("Met", "Thr", "Thr", "Thr", "Thr", "Trp", "Stp", "Stp")
  }
  else if (codon_table_id == 4) {  # Mold Mitochondrial; Protozoan Mitochondrial; Coelenterate Mitochondrial; Mycoplasma; Spiroplasma
    genetic_code[c("tga", "ata")] <- c("Trp", "Met")
  }
  else if (codon_table_id == 5) {  # Invertebrate Mitochondrial
    genetic_code[c("aga", "agg", "tga")] <- c("Ser", "Ser", "Trp")
  }
  else if (codon_table_id == 6) {  # Ciliate Nuclear; Dasycladacean Nuclear; Hexamita Nuclear
    genetic_code[c("taa", "tag")] <- "Gln"
  }
  else if (codon_table_id == 9) {  # Echinoderm Mitochondrial; Flatworm Mitochondrial
    genetic_code[c("aaa", "aga")] <- c("Asn", "Ser")
  }
  else if (codon_table_id == 10) {  # Euplotid nuclear
    genetic_code["tga"] <- "Cys"
  }
  else if (codon_table_id == 11) {  # Bacterial, Archaeal and Plant Plastid
    genetic_code["ctg"] <- "Leu"
  }
  else if (codon_table_id == 12) {  # Alternative Yeast Nuclear
    genetic_code["ctg"] <- "Ser"
  }
  else if (codon_table_id == 13) {  # Ascidian Mitochondrial
    genetic_code[c("aga", "agg")] <- "Gly"
  }
  else if (codon_table_id == 14) {  # Flatworm Mitochondrial
    genetic_code[c("taa", "tag")] <- "Tyr"
  }
  else if (codon_table_id == 15) {  # Blepharisma Macronuclear
    genetic_code[c("taa", "tag")] <- "Glu"
  }
  else if (codon_table_id == 16) {  # Chlorophycean Mitochondrial
    genetic_code["tag"] <- "Leu"
  }
  else if (codon_table_id == 21) {  # Trematode Mitochondrial
    genetic_code["tga"] <- "Trp"
  }
  else if (codon_table_id == 22) {  # Scenedesmus obliquus Mitochondrial
    genetic_code["tca"] <- "Stp"
  }
  else if (codon_table_id == 23) {  # Thraustochytrium Mitochondrial
    genetic_code["tta"] <- "Stp"
  }
  else if (codon_table_id == 24) {  # Pterobranchia
    genetic_code["aga"] <- "Ser"
  }
  else if (codon_table_id == 25) {  # Candidate Division SR1 and Gracilibacteria
    genetic_code["tca"] <- "Ala"
  }
  else if (codon_table_id == 26) {  # Pachysolen tannophilus Nuclear
    genetic_code["ctg"] <- "Ala"
  }
    else if (codon_table_id == 1) {  # Standard
    genetic_code <- genetic_code
  }
  
  return(list(genetic_code = genetic_code))
}
