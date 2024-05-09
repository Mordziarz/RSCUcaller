# RSCUcaller

The RSCUcaller package is designed to analyze multiple DNA sequences for their RSCU (relative synonymous codon usage) values. The package provides functions to obtain RSCU information, plot RSCU values, and calculate statistical significance.

# Installation

```{r}
#install.packages("devtools")
#devtools::install_github('Mordziarz/RSCUcaller')
#library(RSCUcaller)
```

# Input data 

The package saves tables and generates graphics to folders. It is recommended to use the setwd() function to control the program's output.

```{r}
#setwd("path/to/your/directory")
```

RSCUcaller requires DNA sequences from NCBI. In this tutorial, we will analyze CDS DNA sequences for 6 liverwort species. To start the analysis, the fasta file must be prepared appropriately. The prepare_fasta function is used for this purpose, which requires a table with 2 columns: sequence_path with the path to the fasta file and sample_name with the fasta name (these will be visible on the visualizations, so make sure the name is correct). However, if you downloaded multiple species from NCBI in a single fasta file, you can use the samples_table with the ID and GENBANK_ACCESSION columns, but you must also provide the path argument and the path to the fasta file. In both cases, you must specify the name of the fasta file you prepared, which will be saved in your working directory (argument file_out="name_of_your_output.fasta").

```{r}
# path1 <- "/dane/Program_RSCU/sequence(13).txt"
# path2 <- "/dane/Program_RSCU/sequence(15).txt"
# samples_table <- data.frame(sequence_path = c(path1,path2),
#                             sample_name = c("1_Riccia_fluitans","2_Marchantia_polymorpha"))
# 
# prepare_fasta(samples_table = samples_table,file_out = "merged_sequences1.fasta")
```

```{r}
# samples_table <- data.frame(ID = c("Aneura maxima M-1","Aneura maxima M-2"),
#                             GENBANK_ACCESSION = c("OQ700951","OQ700950"))
# prepare_fasta(samples_table = samples_table,path ="/dane/Aneura/CDS_cpDNA_Aneura.txt",file_out = "merged_sequences_lukasz.fasta")
```

# Calculating RSCU from multiple sequences

The main function in the package is get_RSCU(), which only requires the previously prepared sequences to work.

```{r}
#test <- get_RSCU(merged_sequences ="merged_sequences_lukasz.fasta")
```

# RSCU matrix

The get_matrix() function allows the user to create a matrix that the user can use in any way they like. The matrix will not be needed for the next steps. Simply use the result of the get_RSCU() function.

```{r}
#mat1 <- get_matrix(test)
```

# Heatmap

The heatmap_RSCU() function allows you to draw a heatmap and a dendrogram. The heatmap is created from the result of the get_RSCU() function. To call the heatmap, enter "heatmap" in the select argument and choose heatmap_color from: "red_green", "green_red", "blue_green", "green_blue", "blue_red", "red_blue". The names correspond to the colors used to color the heatmap.

```{r}
#heatmap_RSCU(get_RSCU_out = test,select = "heatmap", heatmap_color = "red_blue")
```

To get a dendrogram, you need to provide the result of the get_RSCU() function and enter "dendrogram" in the select argument. You can edit it, for example, in the ggtree package.

```{r}
#heatmap_RSCU(get_RSCU_out = test,select = "dendogram")
```

# Histograms

Our package allows you to draw histograms from the result of the get_RSCU() function. You can create a single histogram

```{r}
#histogram_RSCU(get_RSCU_out = test,title = "LUKASZ")
```

oraz histogram podwójny

```{r}
#histogram_RSCU_double(get_RSCU_out_left = test,get_RSCU_out_right = test,title_left = "LUKASZ",title_right = "MATEUSZ")
```

# Korelacja

RSCUcaller umozliwia zrobienie korelacji pomiedzy dwoma gatunkami

```{r}
#p <- correlation(get_RSCU_out = test,Species_x = "Pelia",Species_y = "Apopelia")
```

# Statystyka pomiedzy grupami poszczególnych codon

```{r}
#Boxplot_between_groups()
```

# Statystyka pomiedzy aminokwasami

```{r}
#stat_scat_box()
```

# Citation
Please cite RSCUcaller as: Mazdziarz M et al., RSCUcaller - R package for RSCU analysis.

# Support
Any issues connected with the RSCUcaller should be addressed to Mateusz Mazdziarz (mazdziarz.mateusz (at) uwm.edu.pl).

