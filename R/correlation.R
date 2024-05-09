#' Creating a csv table, boxplots and scatter plots table based on the table with the get_RSCU function
#'
#' @param get_RSCU_out the table was created using the get_RSCU function
#' @param Species_x The name from the get_RSCU table that you want to use to calculate the correlation will be on the x-axis
#' @param Species_y The name from the get_RSCU table that you want to use to calculate the correlation will be on the y-axis
#' @param xlab x axis title
#' @param ylab y axis titles
#'
#' @return A list with statistics and plot.
#' @export
#'

correlation <- function(get_RSCU_out,Species_x,Species_y,xlab,ylab) {

  if (base::missing(get_RSCU_out)) {
    stop("The get_RSCU_out predictions are required. Please provide a valid argument.",
         call. = FALSE)
  }

  if (base::missing(Species_x)) {
    stop("The Species_x predictions are required. Please provide a valid argument.",
         call. = FALSE)
  }

  if (base::missing(Species_y)) {
    stop("The Species_y predictions are required. Please provide a valid argument.",
         call. = FALSE)
  }

  get_RSCU_out_x <- get_RSCU_out[get_RSCU_out$Species %in% Species_x,]
  get_RSCU_out_y <- get_RSCU_out[get_RSCU_out$Species %in% Species_y,]
  correlation_xy <- base::merge(get_RSCU_out_x,get_RSCU_out_y,by="codon")
  res <- stats::cor.test(correlation_xy$RSCU.x, correlation_xy$RSCU.y, method = 'pearson')

  p <- ggplot2::ggplot(data = correlation_xy,mapping=aes(x = RSCU.x,
                                                         y = RSCU.y,
                                                         add="jitter"))+
    ggplot2::geom_point(shape=21,
                        fill='steelblue3',
                        color='white',
                        size = 3,
                        alpha=0.5)+
    ggplot2::geom_smooth(method=lm,
                         color='red4')+
    smplot2::sm_statCorr(show_text = F,
                         col="red4")+
    ggplot2::theme_bw()+
    ggplot2::xlab(xlab)+
    ggplot2::ylab(ylab)+
    ggplot2::ggtitle(base::paste0("R=",base::round(res$estimate,2),", ","p-value=",base::round(res$p.value,3),", ","95% conf.int=",round(res$conf.int,2)))

  return(list(res=res,plot=p))

}
