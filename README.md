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

DNA sequences from NCBI are required by RSCUcaller. To initiate the analysis, the FASTA file needs to be prepared appropriately. This is achieved using the prepare_fasta() function, which requires a table with two columns:
    sequence_path: This column specifies the path to the FASTA file.
    sample_name: This column defines the FASTA name, which will be displayed on visualizations and should start with a number followed by an underscore (e.g., 1_Riccia_fluitans).

```{r}
path1 <- "/path/to/your/fasta"
samples_table <- data.frame(sequence_path = path1,
                            sample_name = "1_fasta_name")
prepare_fasta(samples_table = samples_table,file_out = "your_fasta.fasta")
```

However, if you downloaded multiple species from NCBI within a single FASTA file, you can utilize the samples_table with the ID and GENBANK_ACCESSION columns. In this case, you must also provide the path to the FASTA file using the path argument.

```{r}
samples_table <- data.frame(ID = "1_fasta_name",
                            GENBANK_ACCESSION = "gene_bank_accession_id")
prepare_fasta(samples_table = samples_table, path = "/path/to/multiple/sequence/fasta", file_out = "your_fasta.fasta")
```
Regardless of the method, the name of the prepared FASTA file needs to be specified. This file will be saved in your working directory (using the argument file_out="name_of_your_output.fasta").

# Calculating RSCU from multiple sequences

The main function in the package is get_RSCU(), which only requires the previously prepared sequences to work.

```{r}
get_RSCU_out <- get_RSCU(merged_sequences = "your_fasta.fasta")
```

# RSCU matrix

The get_matrix() function allows the user to create a matrix that the user can use in any way they like. The matrix will not be needed for the next steps. Simply use the result of the get_RSCU() function.

```{r}
get_matrix(get_RSCU_out = get_RSCU_out)
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

# Statistics between groups of individual codons

The boxplot_between_groups() function enables statistical analysis of RSCU values at the codon level. It utilizes the Kruskal-Wallis test to assess significant differences in RSCU values among groups, followed by Dunn's post-hoc test for pairwise comparisons.

To utilize this function, you need to prepare a table called grouping_table with two columns:
name: This column should contain the names of the sequences, corresponding to the name column from the get_RSCU_out() data frame.
group: This column should specify the group to which each sequence belongs.

```{r}
Boxplot_between_groups(get_RSCU_out = get_RSCU_out, grouping_table = grouping_table, width = 6, height = 6, xlab = "title of x lab", res = 300)
```

# Statistics between amino acids

The stat_scat_box() function performs statistical analysis and generates visualizations for your data. It utilizes the Kruskal-Wallis test to assess significant differences among groups, followed by Dunn's post-hoc test for pairwise comparisons. The function also generates table (Post_hoc_table_aminoacids.csv) and two folders (boxplots and scatter_plots) within your working directory:

```{r}
stat_scat_box(get_RSCU_out = get_RSCU_out, width = 6, height = 6, res = 300)
```

# Citation
Please cite RSCUcaller as: Mazdziarz M et al., RSCUcaller - R package for RSCU analysis.

# Support
Any issues connected with the RSCUcaller should be addressed to Mateusz Mazdziarz (mazdziarz.mateusz (at) uwm.edu.pl).

