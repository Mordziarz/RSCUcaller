kmeans_RSCU <- function(get_matrix_out=get_matrix_out){

  if (base::missing(get_matrix_out)) {
    stop("The get_matrix_out predictions are required. Please provide a valid argument.",
         call. = FALSE)
  }
  
sd_rscu <- base::apply(get_matrix_out, 2, sd)
index <- base::which(sd_rscu != 0)
matrix_0 <- dane[, index]

p <- factoextra::fviz_nbclust(matrix_0, kmeans, method = "silhouette") +
  ggplot2::labs(subtitle = "Silhouette method")

number <- p$data[p$data$y %in% base::max(p$data$y),]
number <- base::as.numeric(number$clusters)
kmeans_result <- stats::kmeans(matrix_0, centers = number, nstart = 25)

clusters <- factoextra::fviz_cluster(kmeans_result, 
                                     data = matrix_0, 
                                     geom = "point",
                                     ellipse.type = "convex", 
                                     ggtheme = theme_bw()
)

matrix_0$cluster <- kmeans_result$cluster
data_frame_cluster <- base::table(rownames(matrix_0), matrix_0$cluster)
distance_matrix <- stats::dist(matrix_0, method = "euclidean")
hc <- stats::hclust(distance_matrix, method = "ward.D2")
dend <- ape::as.phylo(hc)
tree <- ggtree(dend)
ape::write.tree(dend, file = "k_means.newick")

return(list(table = data_frame_cluster, number_of_clusters = p, clusters = clusters, tree = tree))

base::message(base::paste0("Success"))

}