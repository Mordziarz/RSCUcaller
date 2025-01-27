# Installation 

install.packages("devtools")
library(devtools)
devtools::install_github('Mordziarz/RSCUcaller')
library(RSCUcaller)

#Other libraries

library(ggplot2)
library(dplyr)
library(seqinr)
library(rstatix)
library(ggpubr)
library(patchwork)
library(forcats)
library(phylogram)
library(circlize)
library(ComplexHeatmap)
library(smplot2)
library(ggtree)
library(stats)
library(stringr)
library(cluster)
library(factoextra)
library(ape)

set.seed(123) 

setwd("path/to/your/workdir")

################################################
################## Mitogenome ##################
################################################

#################################
###### read sample table ########
#################################

samples_table <- read.csv2("sample_table_mt.csv",sep = ";")

###################################
###### preprocessing fasta ########
###################################

RSCUcaller::prepare_fasta(samples_table = samples_table,path = "mitogenome_sequence.txt",file_out = "prepered_fasta.fasta")

#################################
###### calculate RSCU ###########
#################################

get_RSCU_out <- RSCUcaller::get_RSCU(merged_sequences = "prepered_fasta.fasta")

#################################
###### PR2 plot #################
#################################

png("Ex_PR2.png", width=4, height=4, units = "in", res = 300)
RSCUcaller::neutrality_pr2(get_RSCU_out = get_RSCU_out,select = "PR2_plot")
dev.off()

#################################
###### Neutrality plot ##########
#################################

png("Ex_neutrality.png", width=4, height=4, units = "in", res = 300)
RSCUcaller::neutrality_pr2(get_RSCU_out = get_RSCU_out,select = "neutrality_plot")
dev.off()

######################
###### matrix ########
######################

matrix <- RSCUcaller::get_matrix(get_RSCU_out = get_RSCU_out)

#######################
###### heatmap ########
#######################

png("Ex_heatmap.png", width=15, height=10, units = "in", res = 300)
RSCUcaller::heatmap_RSCU(get_RSCU_out = get_RSCU_out, select = "heatmap", heatmap_color = "red_blue")
dev.off()

##############################
###### dendogram/tree ########
##############################

RSCUcaller::heatmap_RSCU(get_RSCU_out = get_RSCU_out, select = "dendogram")
p1 <- ape::read.tree("dendogram_from_heatmap.newick")
p_tree <- ggtree(p1) + geom_tiplab(fontface="italic") + xlim(NA,0.9)

png("Ex_tree.png", width=5, height=5, units = "in", res = 300)
p_tree
dev.off()

##############################
###### histogram #############
##############################

get_RSCU_out_histogram <- get_RSCU_out[get_RSCU_out$Species %in% c("Apopellia_endiviifolia_A1",
                                                                   "Apopellia_endiviifolia_C1",
                                                                   "Pellia_neesiana_1",
                                                                   "Pellia_epiphylla_S1"),] # chose organisms

get_RSCU_out_histogram <- get_RSCU_out_histogram[order(get_RSCU_out_histogram$Species),] #order

png("Ex_histogram.png", width=8, height=8, units = "in", res = 300)
RSCUcaller::histogram_RSCU(get_RSCU_out = get_RSCU_out_histogram, title = "Mitogenome") # plot
dev.off()

###################
### Correlation ###
###################

png("Ex_correlation.png", width=5, height=5, units = "in", res = 300)
RSCUcaller::correlation(get_RSCU_out = get_RSCU_out, Species_x = "Apopellia_endiviifolia_A3", Species_y = "Pellia_neesiana_3", 
                        xlab = "Apopellia endiviifolia A3 \n Mitogenome", ylab = "Pellia neesiana 3 \n Mitogenome")
dev.off()

################################
### statistc between codons ###
################################

RSCUcaller::stat_scat_box(get_RSCU_out = get_RSCU_out, width = 6, height = 6, res = 300)

################################
### statistc between groups ###
################################

grouping_table <- read.csv2("grouping_table.csv",sep = ";") # read grouping table

RSCUcaller::boxplot_between_groups(get_RSCU_out = get_RSCU_out, grouping_table = grouping_table, 
                                   width = 6, height = 6, xlab = "Groups", res = 300) # calculate

################################
########### PCA ################
################################

png("Ex_PCA.png", width=5, height=5, units = "in", res = 300)
RSCUcaller::PCA_RSCU(get_matrix_out = matrix,grouping_table = grouping_table)
dev.off()

###############################################
########### k-means clustering ################
###############################################

object <- RSCUcaller::kmeans_RSCU(get_matrix_out = matrix)

object$table

png("Ex_clusters.png", width=6, height=6, units = "in", res = 300)
object$number_of_clusters
dev.off()

png("Ex_clusters_plot.png", width=6, height=6, units = "in", res = 300)
object$clusters
dev.off()

png("Ex_tree_k_m.png", width=6, height=6, units = "in", res = 300)
object$tree
dev.off()

#############################################
################## Plastid ##################
#############################################

# Path to sequence

path1 <- "OL654070.fasta"
path2 <- "OQ280817.txt"
path3 <- "OQ280829.txt"
path4 <- "OQ280824.txt"

samples_table_cp <- data.frame(sequence_path = c(path1,
                                              path2,
                                              path3,
                                              path4),
                            sample_name = c("1_Apopellia_endiviifolia_A1",
                                            "2_Apopellia_endiviifolia_C1",
                                            "3_Pellia_neesiana_1",
                                            "4_Pellia_epiphylla_S1")) #create samples_table

RSCUcaller::prepare_fasta(samples_table = samples_table_cp,file_out = "prepered_fasta_plastid.fasta") # prepare fasta file

get_RSCU_out_cp <- RSCUcaller::get_RSCU(merged_sequences = "prepered_fasta_plastid.fasta") #calculate RSCU

######################################################################################
############# Compare mitogenome and plastid data on a double histogram ##############
######################################################################################

get_RSCU_out_cp <- get_RSCU_out_cp[order(get_RSCU_out_cp$Species),] #order
get_RSCU_out_histogram <- get_RSCU_out_histogram[order(get_RSCU_out_histogram$Species),] #order

png("Ex_double_histogram.png", width=12, height=8, units = "in", res = 300)
RSCUcaller::histogram_RSCU_double(get_RSCU_out_left = get_RSCU_out_cp,get_RSCU_out_right = get_RSCU_out_histogram,
                                  title_left = "Plastid",title_right = "Mitogenome")
dev.off()
