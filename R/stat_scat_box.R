#' Creating a csv table, boxplots and scatter plots table based on the table with the get_RSCU function
#'
#' @param get_RSCU_out the table was created using the get_RSCU function
#' @param width width of charts
#' @param height height of charts
#' @param res resolution of charts
#'
#' @return A table with statistics, boxplots and scatter plots.
#' @export
#'

stat_scat_box <-function(get_RSCU_out,width,height,res) {

  if (base::missing(get_RSCU_out)) {
    stop("The get_RSCU_out predictions are required. Please provide a valid argument.",
         call. = FALSE)
  }

  if (!base::file.exists("scatter_plots") & !base::file.exists("boxplots")){
    base::dir.create("scatter_plots")
    base::dir.create("boxplots")
  }

  if (base::file.exists("scatter_plots") & base::file.exists("boxplots")) {

    base::message(base::paste0("Calculating statistics and generating plots... "))
    aminoacids <- c("Ala","Arg","Asn","Asp","Cys","Gln","Glu","Gly","His","Ile","Leu","Lys","Phe","Pro","Ser","Stp","Thr","Tyr","Val")
    statistical_table <- base::data.frame(row.names = 1, .y. = NA, group1 = NA, group2 = NA, n1 = NA, n2 = NA, statistic = NA, p = NA, p.adj = NA, p.adj.signif = NA)

    for (i in 1:base::length(aminoacids)) {
      table_1 <- get_RSCU_out[get_RSCU_out$AA %in% aminoacids[i],]
      table_1 <- table_1[,c("RSCU", "codon")]
      stat <- stats::kruskal.test(RSCU ~ codon , data = table_1)
      cat("The differences between codons in the amino acid",aminoacids[i],"had a pvalue of:",stat$p.value,"\n")
      post_hoc <- base::as.data.frame(rstatix::dunn_test(RSCU ~ codon, data = table_1, p.adjust.method = "bonferroni"))
      statistical_table <- base::rbind(statistical_table,post_hoc)
      statistical_table <- statistical_table[!base::is.na(statistical_table$group1),]
      post_hoc <- dplyr::arrange(.data = post_hoc, p.adj)
      post_hoc_x_y <- rstatix::add_xy_position(post_hoc, x=post_hoc[2])
      png(paste0("scatter_plots/",aminoacids[i],".png"), width=width, height=height, units = "in", res=res)
      p <- ggpubr::ggscatter(table_1,
                             x="codon",
                             y="RSCU",
                             point = "FALSE") +
        ggplot2::geom_jitter(aes(color=codon),width = 0.1)+
        ggpubr::stat_pvalue_manual(post_hoc_x_y, label="p.adj.signif", hide.ns=TRUE) +
        ggplot2::theme(legend.position = "none") +
        ggplot2::ggtitle(base::paste0(aminoacids[i],", Kruskal-Wallis, p=",stat$p.value))
      base::print(p)
      dev.off()
      png(paste0("boxplots/",aminoacids[i],".png"), width=width, height=height, units = "in", res = res)
      p <- ggpubr::ggboxplot(table_1,
                             x="codon",
                             y="RSCU",
                             fill = "white",
                             color = "codon",
                             add = "jitter",
                             palette = c("dodgerblue3", "maroon2",  "forestgreen", "darkorange1", "blueviolet", "firebrick2")) +
        ggpubr::stat_pvalue_manual(post_hoc_x_y, label="p.adj.signif", hide.ns=TRUE) +
        ggplot2::theme(legend.position = "none") +
        ggplot2::ggtitle(base::paste0(aminoacids[i],", Kruskal-Wallis, p=",stat$p.value))
      base::print(p)
      dev.off()
    }
  }
  else (
    stop("Try to make directory scatter_plots & boxplots manually and run function again",
         call. = FALSE)
  )
  base::message(base::paste0("Done! Check scatter_plots & boxplots dirs !!!"))
  write.csv2(statistical_table,"Post_hoc_table_aminoacids.csv",row.names = F)
  return(statistical_table)
}
