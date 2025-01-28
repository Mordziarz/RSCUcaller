#' Creating a heatmap and grouping dendrogram based on the table with the get_RSCU function
#'
#' @param get_RSCU_out the table was created using the get_RSCU function
#' @param heatmap_color Selecting colors for the heatmap from among "red_green","green_red","blue_green","green_blue","blue_red","red_blue"
#' @param select the "heatmap" or "dendrogram" output
#'
#' @return A plot object.
#' @export
#'

heatmap_RSCU <- function(heatmap_color="",get_RSCU_out=get_RSCU_out,select=""){

  if (base::missing(get_RSCU_out)) {
    stop("The get_RSCU_out predictions are required. Please provide a valid argument.",
         call. = FALSE)
  }

  if (base::missing(select)) {
    stop("The select predictions are required. Please provide a valid argument.",
         call. = FALSE)
  }

  base::message(base::paste0("Generating plot... "))

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

  if (select=="heatmap"){
    valid_colors <- c("red_green",
                      "green_red",
                      "blue_green",
                      "green_blue",
                      "blue_red",
                      "red_blue")

    if (!(heatmap_color %in% valid_colors)) {
      stop("Invalid heatmap_color argument. Please choose from red_green, green_red, blue_green, green_blue, blue_red, red_blue valid options.")
    }

    if (heatmap_color=="green_red") {
      colors <- circlize::colorRamp2(c(0,0.2,0.4,0.6,0.8,1,1.2,1.4,1.6,1.8,2),
                                     c("#006d2c","#2ca25f","#66c2a4","#99d8c9","#ccece6","white","#fdd49e","#fdbb84", "#fc8d59","#e34a33","#b30000"))
    }

    if (heatmap_color=="red_green") {
      colors <- circlize::colorRamp2(c(0,0.2,0.4,0.6,0.8,1,1.2,1.4,1.6,1.8,2),
                                     c("#b30000","#e34a33","#fc8d59","#fdbb84","#fdd49e","white","#ccece6","#99d8c9","#66c2a4","#2ca25f","#006d2c"))
    }

    if (heatmap_color=="green_blue") {
      colors <- circlize::colorRamp2(c(0,0.2,0.4,0.6,0.8,1,1.2,1.4,1.6,1.8,2),
                                     c("#006d2c","#2ca25f","#66c2a4","#99d8c9","#ccece6","white","#c6dbef","#9ecae1","#6baed6","#3182bd","#08519c"))
    }

    if (heatmap_color=="blue_green") {
      colors <- circlize::colorRamp2(c(0,0.2,0.4,0.6,0.8,1,1.2,1.4,1.6,1.8,2),
                                     c("#08519c","#3182bd","#6baed6","#9ecae1","#c6dbef","white","#ccece6","#99d8c9","#66c2a4", "#2ca25f","#006d2c"))
    }

    if (heatmap_color=="blue_red") {
      colors <- circlize::colorRamp2(c(0,0.2,0.4,0.6,0.8,1,1.2,1.4,1.6,1.8,2),
                                     c("#08519c","#3182bd","#6baed6","#9ecae1","#c6dbef","white","#fdd49e","#fdbb84","#fc8d59","#e34a33","#b30000"))
    }

    if (heatmap_color=="red_blue") {
      colors <- circlize::colorRamp2(c(0,0.2,0.4,0.6,0.8,1,1.2,1.4,1.6,1.8,2),
                                     c("#b30000","#e34a33","#fc8d59","#fdbb84","#fdd49e","white","#c6dbef","#9ecae1","#6baed6","#3182bd","#08519c"))
    }


    matrix_RSCU <- base::as.matrix(matrix_RSCU)
    heatmap_rscu <- ComplexHeatmap::Heatmap(matrix_RSCU,name = "RSCU"),
                                            col = colors,
                                            column_dend_height = unit(3, "cm"),
                                            row_dend_width = unit(3, "cm"),
                                            cluster_rows = T,
                                            cluster_columns = T)


    heatmap_row <- ComplexHeatmap::rowAnnotation(italic_text = ComplexHeatmap::anno_text(rownames(matrix_RSCU),
                                                                                         gp = grid::gpar(fontface="italic")))
    heatmap_rscu <- heatmap_rscu + heatmap_row
    heatmap_rscu <- draw(heatmap_rscu)
    base::message(base::paste0("Success"))
    return(heatmap_rscu)

  }

  if (select=="dendogram") {
    matrix_RSCU <- base::as.matrix(matrix_RSCU)
    heatmap_rscu <- ComplexHeatmap::Heatmap(matrix_RSCU,name = "RSCU",
                                            rect_gp = grid::gpar(col = "black",
                                                                 lwd = 0.5),
                                            column_dend_height = unit(3, "cm"),
                                            row_dend_width = unit(3, "cm"),
                                            cluster_rows = T,
                                            cluster_columns = T)
    heatmap_rscu <- draw(heatmap_rscu)
    p <- ComplexHeatmap::row_dend(heatmap_rscu)
    p1 <- ggtree::ggtree(p)
    phylogram::write.dendrogram(p, file = "dendogram_from_heatmap.newick", append = FALSE, edges = TRUE)
    base::message(base::paste0("Success"))
    return(p1)

  }
}
