#' Creating a boxplot and csv table based on the table with the get_RSCU function and grouping_table
#'
#' @param get_RSCU_out the table was created using the get_RSCU function
#' @param grouping_table the table with Species and group column names
#' @param xlab x axis title
#' @param width width of charts
#' @param height height of charts
#' @param res resolution of charts
#'
#' @return A plot and csv objects.
#' @export
#'

boxplot_between_groups <-function(get_RSCU_out,grouping_table,width,height,xlab,res,p.adjust.method = "bonferroni") {

  if (base::missing(get_RSCU_out)) {
    stop("The get_RSCU_out predictions are required. Please provide a valid argument.",
         call. = FALSE)
  }

  if (base::all(colnames(grouping_table)!=c("Species","group"))) {
    stop("Check the colnames of the grouping_table",
         call. = FALSE)
  }

  if (!base::file.exists("selected_species") & !base::file.exists("selected_species_barplots")){
    base::dir.create("selected_species")
    base::dir.create("selected_species_barplots")
  }

  if (base::file.exists("selected_species") & base::file.exists("selected_species_barplots")) {

    base::message(base::paste0("Calculating statistics and generating plots... "))
    codons <- c("aaa", "aac", "aag", "aat", "aca", "acc", "acg", "act", "aga", "agc", "agg", "agt", "ata", "atc", "att", "caa", "cac", "cag", "cat", "cca", "ccc", "ccg", "cct",
                "cga", "cgc", "cgg", "cgt", "cta", "ctc", "ctg", "ctt", "gaa", "gac", "gag", "gat", "gca", "gcc", "gcg", "gct", "gga", "ggc", "ggg", "ggt", "gta", "gtc", "gtg", "gtt",
                "taa", "tac", "tag", "tat", "tca", "tcc", "tcg", "tct", "tga", "tgc", "tgt", "tta", "ttc", "ttg", "ttt")
    statistical_table <- base::data.frame(row.names = 1, .y. = NA, group1 = NA, group2 = NA, n1 = NA, n2 = NA, statistic = NA, p=NA, p.adj=NA, p.adj.signif=NA)
    get_RSCU_out <- merge(get_RSCU_out,grouping_table,by="Species")
    get_RSCU_out$index2 <- paste0(get_RSCU_out$group,"_",get_RSCU_out$codon)
    get_RSCU_out <- get_RSCU_out[order(get_RSCU_out$group),]
    
    for (i in 1:base::length(codons)) {
      table_1 <- get_RSCU_out[get_RSCU_out$codon %in% codons[i],]
      table_1 <- table_1[,c("RSCU", "index2","group")]
      stat <- stats::kruskal.test(RSCU ~ index2 , data = table_1)
      post_hoc <- base::as.data.frame(rstatix::dunn_test(RSCU ~ index2, data = table_1, p.adjust.method = p.adjust.method))
      statistical_table <- base::rbind(statistical_table,post_hoc)
      statistical_table <- statistical_table[!base::is.na(statistical_table$group1),]
      post_hoc <- dplyr::arrange(.data = post_hoc, p.adj)
      post_hoc_x_y <- rstatix::add_xy_position(post_hoc, x=post_hoc[2])
      png(paste0("selected_species/",codons[i],".png"), width=width, height=height, units = "in", res=res)
      p <- ggpubr::ggboxplot(table_1,
                             x="group",
                             y="RSCU",
                             fill = "white",
                             color = "group",
                             xlab = xlab,
                             palette = c("dodgerblue3", "maroon2",  "forestgreen", "darkorange1", "blueviolet", "firebrick2",
                                         "deepskyblue", "orchid2", "chartreuse3", "gold", "slateblue1", "tomato" , "blue", "magenta", "green3",
                                         "yellow", "purple3", "red" ,"darkslategray1", "lightpink1", "lightgreen", "khaki1", "plum3", "salmon")) +
        ggpubr::stat_pvalue_manual(post_hoc_x_y, label="p.adj.signif", hide.ns=TRUE) +
        ggplot2::geom_jitter(aes(color = group), width = 0.2, height = 0) +
        ggplot2::theme(legend.position = "none") +
        ggplot2::scale_x_discrete() +
        ggplot2::ggtitle(base::paste0(codons[i],", Kruskal-Wallis, p=",stat$p.value))
      base::print(p)
      dev.off()
      table_1 <- table_1 %>% dplyr::group_by(group) %>% dplyr::summarize(`RSCU` = mean(RSCU))
      get_RSCU_out_groups <- get_RSCU_out[!duplicated(get_RSCU_out$group),]
      idx <- base::match(table_1$group,get_RSCU_out_groups$group)
      table_1 <- table_1[idx,]
      png(paste0("selected_species_barplots/", codons[i], ".png"), width = width, height = height, units = "in", res = res)
      p <- ggpubr::ggbarplot(table_1, 
                             x = "group", 
                             y = "RSCU", 
                             fill = "group", 
                             color = "group",
                             xlab = xlab,
                             palette = c("dodgerblue3", "maroon2",  "forestgreen", "darkorange1", "blueviolet", "firebrick2",
                                         "deepskyblue", "orchid2", "chartreuse3", "gold", "slateblue1", "tomato" , "blue", "magenta", "green3",
                                         "yellow", "purple3", "red" ,"darkslategray1", "lightpink1", "lightgreen", "khaki1", "plum3", "salmon")) + 
        ggpubr::stat_pvalue_manual(post_hoc_x_y,label = "p.adj.signif", hide.ns = TRUE) + 
        ggplot2::theme(legend.position = "none") + 
        ggplot2::scale_x_discrete() + 
        ggplot2::ggtitle(base::paste0(codons[i],", Kruskal-Wallis, p=", stat$p.value)) + ggplot2::ylab("mean(RSCU)")
      base::print(p)
      dev.off()
    }
  }
  else (
    stop("Try to make directory selected_species & selected_species_barplots manually and run function again",
         call. = FALSE)
  )
  base::message(base::paste0("Done! Check selected_species & selected_species_barplots dir !!!"))
  write.csv2(statistical_table,"Post_hoc_table_selected_species.csv",row.names = F)
  return(statistical_table)
}
