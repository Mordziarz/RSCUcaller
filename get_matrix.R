#' Calculating RSCU for multiple sequences
#'
#' @param get_RSCU_out the table was created using the get_RSCU function
#'
#' @return A `matrix` object containing matrix_RSCU.
#' @export
#'

get_matrix <- function(get_RSCU_out=get_RSCU_out){
  
  if (base::missing(get_RSCU_out)) {
    stop("The get_RSCU_out predictions are required. Please provide a valid argument.",
         call. = FALSE)
  }
  
  base::message(base::paste0("Generating matrix... "))
  
  RSCU_table <- get_RSCU_out[,c("Species","RSCU","AA","codon")]
  RSCU_table$AAcodon <- base::paste0(RSCU_table$AA," | ",RSCU_table$codon)
  RSCU_table$codon <- NULL
  RSCU_table$AA <- NULL 
  matrix_RSCU <- base::matrix(ncol=64,nrow=length(unique(RSCU_table$Species)))
  base::colnames(matrix_RSCU) <- base::unique(RSCU_table$AAcodon)
  base::rownames(matrix_RSCU) <- base::unique(RSCU_table$Species)
  matrix_RSCU <- as.data.frame(matrix_RSCU)
  
  for (i in 1:base::nrow(matrix_RSCU)) {
    j <- 64*i-63
    j <- base::as.numeric(j)
    k <- 64 * i
    k <- base::as.numeric(k)
    lista_rscu <- RSCU_table$RSCU[j:k]
    matrix_RSCU[i,] <- base::ifelse(
      base::rownames(matrix_RSCU)[i] %in% RSCU_table$Species & base::colnames(matrix_RSCU) %in% RSCU_table$AAcodon,
      lista_rscu, matrix_RSCU[i,])
  }
  base::message(base::paste0("Success"))
  return(matrix_RSCU)
}
