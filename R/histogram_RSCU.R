#' Visualize Relative Synonymous Codon Usage (RSCU) Distribution
#' 
#' Generates a comprehensive histogram displaying codon usage patterns across samples,
#' integrating amino acid and codon annotations with sample grouping information.
#'
#' @param get_RSCU_out A data frame generated by \code{\link{get_RSCU}} containing
#'   RSCU values for multiple samples
#' @param title Main title for the visualization (optional)
#'
#' @return A \code{patchwork} object combining three plot elements:
#' \itemize{
#'   \item Left panel: Amino acid labels with codon annotations
#'   \item Center panel: Main RSCU histogram
#'   \item Right panel: Sample group legend
#' }
#' 
#' @details
#' This visualization provides a multi-layered view of codon usage patterns:
#' \itemize{
#'   \item Uses a reverse-sorted layout for easy comparison between samples
#'   \item Integrates amino acid and codon information in left margin
#'   \item Color-codes samples for group differentiation
#'   \item Automatically handles plot arrangement and scaling
#' }
#' 
#' The plot is particularly useful for identifying:
#' \itemize{
#'   \item Overrepresented codons within amino acid groups
#'   \item Sample-specific codon usage biases
#'   \item Grouping patterns across experimental conditions
#' }
#'
#' @examples
#' # Using built-in data
#' data("prepared_fasta", package = "RSCUcaller")
#' rscu_data <- get_RSCU(prepared_fasta)
#' 
#' # Visualization
#' hist_plot <- histogram_RSCU(rscu_data, title = "Codon Usage Distribution")
#' # Display plot
#' print(hist_plot)
#'
#' @seealso
#' Related visualization functions:
#' \code{\link{histogram_RSCU_double}} for comparative visualizations,
#' \code{\link{heatmap_RSCU}} for multivariate patterns
#' 
#' @export
#' @importFrom ggplot2 aes geom_tile geom_bar theme_minimal coord_flip ggtitle theme element_text element_blank
#' @importFrom patchwork plot_layout
#' @importFrom forcats fct_rev

histogram_RSCU <- function(get_RSCU_out=get_RSCU_out,title=""){

  if (base::missing(get_RSCU_out)) {
    stop("The get_RSCU_out predictions are required. Please provide a valid argument.",
         call. = FALSE)
  }

  base::message(base::paste0("Generating plot... "))

  get_RSCU_out$Col2 <- 1
  get_RSCU_out$Col2 <- as.factor(get_RSCU_out$Col2)

  p0<-ggplot2::ggplot(data=get_RSCU_out, aes(x=index,
                                             y=RSCU,
                                             fill=Col))+
    ggplot2::geom_bar(stat='identity')+
    ggplot2::theme_minimal()+
    ggplot2::coord_flip()+
    ggplot2::ggtitle(title)+
    ggplot2::theme(plot.title = element_text(hjust = 0.5)) +
    ggplot2::theme(axis.text.y =element_blank(),
                   axis.ticks.y=element_blank(),
                   axis.ticks.x=element_blank(),
                   axis.title.y=element_blank(),
                   panel.grid.major=element_blank(),
                   panel.border = element_blank(),
                   panel.background = element_blank(),
                   legend.position="none")


  p1<-ggplot2::ggplot(data=get_RSCU_out)+geom_tile(aes(x=AA,
                                                       y=fct_rev(Col),
                                                       fill=Col),
                                                   col="white")+
    ggplot2::geom_text(aes(x=AA,
                           y=Col,
                           label=codon),
                       col="white") +
    ggplot2::theme_void()+
    ggplot2::theme(legend.position="none")+
    ggplot2::coord_flip()

  p2<-ggplot2::ggplot(data=get_RSCU_out)+
    ggplot2::geom_tile(aes(x=AA,
                           y=fct_rev(Col2)),
                       col="white")+
    ggplot2::geom_text(aes(x=AA,
                           y=Col2,label=AA),
                       col="white",
                       fontface = "bold",angle = 90) +
    ggplot2::theme_void()+
    ggplot2::theme(legend.position="none") +
    ggplot2::coord_flip()+
    ggplot2::ggtitle("")


  p3<-ggplot2::ggplot(data=get_RSCU_out)+
    ggplot2::geom_tile(aes(x=index,
                           y=fct_rev(Col[1]),
                           fill=Species),
                       col="white") +
    ggplot2::theme_void()+
    ggplot2::theme(legend.position="right")+
    ggplot2::coord_flip()+
    ggplot2::ggtitle("")

  p=p2+p1+p0+p3+patchwork::plot_layout(ncol = 5, widths = c(0.2,1,2,0.5))

  base::message(base::paste0("Success"))
  return(p)
}
