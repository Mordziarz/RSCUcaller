#' Creating a double barplot based on the table with the get_RSCU function
#'
#' @param get_RSCU_out_left the table was created using the get_RSCU function, which will be drawn on the left side of the chart.
#' @param get_RSCU_out_right the table was created using the get_RSCU function, which will be drawn on the right side of the chart.
#' @param title_left left chart title 
#' @param title_right right chart title 
#'
#' @return A plot object.
#' @export
#'

histogram_RSCU_double <- function(get_RSCU_out_left=get_RSCU_out,get_RSCU_out_right=get_RSCU_out,title_left="",title_right=""){
  
  if (base::missing(get_RSCU_out_left)) {
    stop("The get_RSCU_out_left predictions are required. Please provide a valid argument.",
         call. = FALSE)
  }
  
  if (base::missing(get_RSCU_out_right)) {
    stop("The get_RSCU_out_right predictions are required. Please provide a valid argument.",
         call. = FALSE)
  }
  
  base::message(base::paste0("Generating plot... "))
  
  get_RSCU_out_right$Col2 <- 1
  get_RSCU_out_right$Col2 <- as.factor(get_RSCU_out_right$Col2)
  
  p0<-ggplot2::ggplot(data=get_RSCU_out_right, aes(x=index, y=RSCU, fill=Col))+
      ggplot2::geom_bar(stat='identity')+
      ggplot2::theme_minimal()+
      ggplot2::coord_flip()+
      ggplot2::ggtitle(title_right)+
      ggplot2::theme(axis.text.y=element_blank(),
                     axis.ticks.y=element_blank(),
                     axis.ticks.x=element_blank(),
                     axis.title.y=element_blank(),
                     legend.position="none",
                     panel.grid.major = element_blank(),
                     panel.border = element_blank(),
                     panel.background = element_blank(),
                     plot.title = element_text(hjust = 0.5))
  
  
  p1<-ggplot2::ggplot(data=get_RSCU_out_right)+
      ggplot2::geom_tile(aes(x=AA, y=fct_rev(Col),fill=Col),col="white")+
      ggplot2::geom_text(aes(x=AA, y=Col,label=codon),col="white")+
      ggplot2::theme_void()+
      ggplot2::theme(legend.position="none")+ 
      ggplot2::coord_flip()
  
  p3<-ggplot2::ggplot(data=get_RSCU_out_left, aes(x=index, y=RSCU, fill=Col))+
      ggplot2::geom_bar(stat='identity')+
      ggplot2::theme_minimal()+
      ggplot2::ylab("RSCU")+ 
      ggplot2::scale_y_reverse()+
      ggplot2::coord_flip()+
      ggplot2::ggtitle(title_left)+
      ggplot2::theme(axis.text.y=element_blank(),
                     axis.ticks.y=element_blank(),
                     axis.ticks.x=element_blank(),
                     axis.title.y=element_blank(),
                     panel.grid.major = element_blank(),
                     panel.border = element_blank(),
                     panel.background = ggplot2::element_blank(),
                     plot.title = ggplot2::element_text(hjust = 0.5),
                     legend.position="none")
  
  p4<-ggplot2::ggplot(data=get_RSCU_out_right)+
      ggplot2::geom_tile(aes(x=AA, y=forcats::fct_rev(Col2)),col="white")+
      ggplot2::geom_text(aes(x=AA, y=Col2,label=AA),col="white", fontface = "bold",angle = 90)+
      ggplot2::theme_void()+
      ggplot2::theme(legend.position="none")+ 
      ggplot2::coord_flip()+
      ggplot2::ggtitle("")
  
  p5<-ggplot2::ggplot(data=get_RSCU_out_right)+
      ggplot2::geom_tile(aes(x=index, y=forcats::fct_rev(Col[1]),fill=Species),col="white")+
      ggplot2::theme_void()+theme(legend.position="right")+ 
      ggplot2::coord_flip()+
      ggplot2::ggtitle("")
  
  p<-p3+p4+p1+p0+p5+patchwork::plot_layout(ncol = 5, widths = c(2,0.2,1,2,0.5))
  
  base::message(base::paste0("Success"))
  return(p)
}
