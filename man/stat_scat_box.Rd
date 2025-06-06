% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stat_scat_box.R
\name{stat_scat_box}
\alias{stat_scat_box}
\title{Statistical Analysis and Visualization of RSCU Values Across Amino Acids}
\usage{
stat_scat_box(get_RSCU_out, width, height, res, p.adjust.method = "bonferroni")
}
\arguments{
\item{get_RSCU_out}{A data frame generated by \code{\link{get_RSCU}} containing RSCU values.}

\item{width}{Width of the output plots in inches.}

\item{height}{Height of the output plots in inches.}

\item{res}{Resolution of the output plots in DPI (dots per inch).}

\item{p.adjust.method}{Method for p-value adjustment in post-hoc tests. One of
\code{"holm"}, \code{"hochberg"}, \code{"hommel"}, \code{"bonferroni"}, \code{"BH"},
\code{"BY"}, or \code{"fdr"}. Default is \code{"bonferroni"}.}
}
\value{
A data frame with post-hoc statistical results. Side effects: 
  - PNG files with boxplots, barplots, and scatter plots are saved in the directories \code{boxplots}, \code{barplots}, and \code{scatter_plots}.
  - A CSV file \code{Post_hoc_table_aminoacids.csv} is saved in the working directory.
}
\description{
Performs comprehensive statistical analysis and generates publication-quality visualizations 
of Relative Synonymous Codon Usage (RSCU) values across amino acids. The function conducts:
\itemize{
  \item Omnibus tests (Kruskal-Wallis, Welch ANOVA, standard ANOVA)
  \item Post-hoc comparisons (Dunn's test, Pairwise t-tests, Tukey HSD)
  \item Multiple testing correction
}
Results are visualized and saved as PNG files, and a summary table is exported as CSV.
}
\examples{
withr::with_tempdir({data("prepared_fasta", package = "RSCUcaller")
rscu_data <- get_RSCU(prepared_fasta)
stat_scat_box(
  get_RSCU_out = rscu_data,
  width = 6,
  height = 6,
  res = 300,
  p.adjust.method = "bonferroni"
)})

}
\seealso{
\code{\link{get_RSCU}}, \code{\link{boxplot_between_groups}}
}
