#' Creating a PCA based on the output from get_matrix function
#'
#' @param get_matrix_out the table was created using the get_RSCU function
#' @param grouping_table the table with Species and group column names
#' @return A plot object.
#' @export
#'

PCA_RSCU <- function(get_matrix_out=get_matrix_out,grouping_table=grouping_table){
  
  if (base::missing(get_matrix_out)) {
    stop("The get_matrix_out predictions are required. Please provide a valid argument.",
         call. = FALSE)
  }

  pca <- stats::prcomp(matrix_0, scale. = TRUE)
  pca_data <- base::as.data.frame(pca$x)
  pca_data$names <- base::rownames(matrix_0)
  pca_data <- base::merge(pca_data,grouping_table,by.x="names",by.y="Species",all.x = T)
  
  PCA_plot <- ggplot2::ggplot(pca_data, aes(x = PC1, y = PC2, label = names,color = group)) +
    ggplot2::geom_point(size = 3) +
    ggplot2::xlab(base::paste0("PC1 (", base::round(base::summary(pca)$importance[2,1]*100, 1), "%)")) +
    ggplot2::ylab(base::paste0("PC2 (", base::round(base::summary(pca)$importance[2,2]*100, 1), "%)")) +
    ggplot2::ggtitle("PCA") +
    ggplot2::theme_bw()+
    ggplot2::theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  
  return(PCA_plot)
  
  base::message(base::paste0("Success"))
  
}