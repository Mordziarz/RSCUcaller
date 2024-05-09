# RSCUcaller

The RSCUcaller package is designed to analyze multiple DNA sequences for their RSCU (relative synonymous codon usage) values. The package provides functions to obtain RSCU information, plot RSCU values, and calculate statistical significance.

# Installation

```{r}
install.packages("devtools")
devtools::install_github('Mordziarz/RSCUcaller')
library(RSCUcaller)
```

# Input data 

The package saves tables and generates graphics to folders. It is recommended to use the setwd() function to control the program's output.

```{r}
setwd("path/to/your/directory")
```

RSCUcaller requires DNA sequences from NCBI. To start the analysis, the fasta file must be prepared appropriately. The prepare_fasta() function is used for this purpose, which requires a table with 2 columns: sequence_path with the path to the fasta file and sample_name with the fasta name (these will be visible on the visualizations, and needs a number 1_ for example 1_Riccia_fluitans). However, if you downloaded multiple species from NCBI in a single fasta file, you can use the samples_table with the ID and GENBANK_ACCESSION columns, but you must also provide the path argument to the fasta file. In both cases, you must specify the name of the fasta file you prepared, which will be saved in your working directory (file_out="name_of_your_output.fasta").

```{r}
path1 <- "/path/to/your/fasta"
samples_table <- data.frame(sequence_path = path1,
                            sample_name = "1_fasta_name")
prepare_fasta(samples_table = samples_table,file_out = "your_fasta.fasta")
```

```{r}
samples_table <- data.frame(ID = "1_fasta_name",
                            GENBANK_ACCESSION = "gene_bank_accession_id")
prepare_fasta(samples_table = samples_table, path = "/path/to/multiple/sequence/fasta", file_out = "your_fasta.fasta")
```

# Calculating RSCU from multiple sequences

The main function in the package is get_RSCU(), which only requires the previously prepared sequences to work.

```{r}
get_RSCU_out <- get_RSCU(merged_sequences = "your_fasta.fasta")
```

# RSCU matrix

The get_matrix() function allows the user to create a matrix that the user can use in any way they like. The matrix will not be needed for the next steps. Simply use the result of the get_RSCU() function.

```{r}
matrix <- get_matrix(get_RSCU_out)
```

# Heatmap

The heatmap_RSCU() function allows you to draw a heatmap and a dendrogram. The heatmap is created from the result of the get_RSCU() function. To call the heatmap, enter "heatmap" in the select argument and choose heatmap_color from: "red_green", "green_red", "blue_green", "green_blue", "blue_red", "red_blue". The names correspond to the colors used to color the heatmap.

```{r}
heatmap_RSCU(get_RSCU_out = get_RSCU_out, select = "heatmap", heatmap_color = "red_blue")
```

To get a dendrogram, you need to provide the result of the get_RSCU() function and enter "dendrogram" in the select argument. You can edit result, for example, in the ggtree package.

```{r}
heatmap_RSCU(get_RSCU_out = get_RSCU_out, select = "dendogram")
```

# Histograms

Our package allows you to draw histograms from the result of the get_RSCU() function. You can create a single histogram

```{r}
histogram_RSCU(get_RSCU_out = get_RSCU_out, title = "graph title")
```

A double histogram is plotted using the function histogram_RSCU_double(). The function takes two input parameters: get_RSCU_out_left and get_RSCU_out_right. It also takes two optional title parameters: title_left and title_right.

```{r}
histogram_RSCU_double(get_RSCU_out_left = get_RSCU_out, get_RSCU_out_right = get_RSCU_out, title_left = "left title", title_right = "right title")
```

# Correlation

The RSCUcaller package allows you to perform Pearson correlation between two species. Specify the names of the sample column in the get_RSCU_out data frame that you want to correlate.

```{r}
correlation(get_RSCU_out = get_RSCU_out, x_name = "x name", y_name = "y name", xlab = "title of x lab", ylab = "title of y lab")
```

# Statystyka pomiedzy grupami poszczególnych codon

```{r}
Boxplot_between_groups()
```

# Statystyka pomiedzy aminokwasami

```{r}
stat_scat_box()
```

# Citation
Please cite RSCUcaller as: Mazdziarz M et al., RSCUcaller - R package for RSCU analysis.

# Support
Any issues connected with the RSCUcaller should be addressed to Mateusz Mazdziarz (mazdziarz.mateusz (at) uwm.edu.pl).

