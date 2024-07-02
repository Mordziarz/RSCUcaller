# RSCUcaller

RSCU (relative synonymous codon usage) values of multiple DNA sequences can be analyzed using the RSCUcaller package. The package offers functions to retrieve RSCU information, plot the values, and calculate statistical significance.

# Installation

```r
install.packages("devtools")
devtools::install_github('Mordziarz/RSCUcaller')
library(RSCUcaller)
```

# Input data 

Tables and graphics are saved to folders by the RSCUcaller package. To specify the desired output directory, it is recommended to use the setwd() function.

```r
setwd("path/to/your/directory")
```

RSCUcaller requires DNA sequences from NCBI. To begin the analysis, the FASTA file must be prepared appropriately using the prepare_fasta() function. This function requires a table with two columns:
"sequence_path" - This column specifies the path to the FASTA file.
"sample_name" - This column defines the FASTA name, which will be displayed on visualizations and should start with a number followed by an underscore (e.g., 1_Riccia_fluitans).

```r
path1 <- "/path/to/your/fasta"
samples_table <- data.frame(sequence_path = path1,
                            sample_name = "1_fasta_name")
prepare_fasta(samples_table = samples_table,file_out = "your_fasta.fasta")
```
Alternatively, if multiple species from NCBI were downloaded within a single FASTA file, the samples_table with the ID and GENBANK_ACCESSION columns can be used. In this scenario, the path to the FASTA file must also be provided using the path argument.

```r
samples_table <- data.frame(ID = "1_fasta_name",
                            GENBANK_ACCESSION = "gene_bank_accession_id")
prepare_fasta(samples_table = samples_table, path = "/path/to/multiple/sequence/fasta", file_out = "your_fasta.fasta")
```
The name of the prepared FASTA file needs to be specified regardless of the chosen method. This file will be created and saved in your working directory using the argument file_out="name_of_your_output.fasta".

# Calculating RSCU from multiple sequences

The main function of the package, get_RSCU(), can be used to calculate RSCU values directly from a set of previously prepared sequences.

```r
get_RSCU_out <- get_RSCU(merged_sequences = "your_fasta.fasta")
```

# RSCU matrix

The get_matrix() function allows the user to create a matrix that the user can use in any way they like. The matrix will not be needed for the next steps. Simply use the result of the get_RSCU() function.

```r
get_matrix(get_RSCU_out = get_RSCU_out)
```

# Heatmap

A heatmap and a dendrogram can be drawn using the heatmap_RSCU() function. The heatmap is generated from the output of the get_RSCU() function. To create the heatmap, specify "heatmap" in the select argument and choose a color scheme from "red_green", "green_red", "blue_green", "green_blue", "blue_red", or "red_blue".

```r
heatmap_RSCU(get_RSCU_out = get_RSCU_out, select = "heatmap", heatmap_color = "red_blue")
```

![Heatmap](https://github.com/Mordziarz/RSCUcaller/tree/main/graphs/heatmap.png , width=6, height=6)

A dendrogram can be obtained by providing the output of the get_RSCU() function and specifying "dendrogram" in the select argument. The result can be further edited using packages like ggtree. The dendrogram will also be saved as a file in your working directory named dendrogram_from_heatmap.newick by the system.

```r
heatmap_RSCU(get_RSCU_out = get_RSCU_out, select = "dendogram")
```
![Dendogram](https://github.com/Mordziarz/RSCUcaller/tree/main/graphs/dendogram.png)

# Histograms

Histograms can be created from the output of the get_RSCU() function using our package.

```r
histogram_RSCU(get_RSCU_out = get_RSCU_out, title = "graph title")
```

![Histogram](https://github.com/Mordziarz/RSCUcaller/tree/main/graphs/histogram.png)

Double histograms are plotted using the histogram_RSCU_double() function. This function requires two input parameters: get_RSCU_out_left and get_RSCU_out_right. Additionally, two optional title parameters, title_left and title_right, can be provided.

```r
histogram_RSCU_double(get_RSCU_out_left = get_RSCU_out, get_RSCU_out_right = get_RSCU_out, title_left = "left title", title_right = "right title")
```

![Histogram double](https://github.com/Mordziarz/RSCUcaller/tree/main/graphs/histogram_double.png)

# Correlation

Pearson correlation between two species can be performed using the RSCUcaller package. The names of the Species columns in the get_RSCU_out data frame to be correlated should be specified.

```r
correlation(get_RSCU_out = get_RSCU_out, Species_x = "Species_x", Species_y = "Species_y", xlab = "title of x lab", ylab = "title of y lab")
```

![Correlation](https://github.com/Mordziarz/RSCUcaller/tree/main/graphs/correlation.png)

# Statistics between groups of individual codons

Statistical analysis of RSCU values at the codon level can be performed using the boxplot_between_groups() function. This function employs the Kruskal-Wallis test to assess significant differences in RSCU values among groups, followed by Dunn's post-hoc test for pairwise comparisons.

To utilize this function, a table named grouping_table needs to be prepared with the following two columns:
"Species"-This column should contain the names of the sequences, corresponding to the Species column from the get_RSCU_out() data frame.
"group"-This column should specify the group to which each sequence belongs.

The function will generate graphical outputs in the "selected_species" folder and a table named Post_hoc_table_selected_species.csv in the working directory.

```r
boxplot_between_groups(get_RSCU_out = get_RSCU_out, grouping_table = grouping_table, width = 6, height = 6, xlab = "title of x lab", res = 300)
```

![Boxplot](https://github.com/Mordziarz/RSCUcaller/tree/main/graphs/Ala.png)

# Statistics between amino acids

Statistical analysis and visualizations can be generated for your data using the stat_scat_box() function. This function employs the Kruskal-Wallis test to assess significant differences among groups, followed by Dunn's post-hoc test for pairwise comparisons. A table named Post_hoc_table_aminoacids.csv and two folders, boxplots and scatter_plots, will also be created in your working directory.

```r
stat_scat_box(get_RSCU_out = get_RSCU_out, width = 6, height = 6, res = 300)
```
![Boxplot](https://github.com/Mordziarz/RSCUcaller/tree/main/graphs/Ala.png)

# Citation

# Support
Any issues connected with the RSCUcaller should be addressed to Mateusz Mazdziarz (mazdziarz.mateusz (at) uwm.edu.pl).

