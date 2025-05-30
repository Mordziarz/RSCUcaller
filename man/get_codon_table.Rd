% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_RSCU_other.R
\name{get_codon_table}
\alias{get_codon_table}
\title{Retrieve a Codon Table for a Genetic Code}
\usage{
get_codon_table(codon_table_id)
}
\arguments{
\item{codon_table_id}{An integer specifying the codon table. See vignette or NCBI documentation
for available codes (e.g., 1 = Standard, 2 = Vertebrate Mitochondrial, etc.).}
}
\value{
A list with one element: \code{genetic_code}, a named character vector mapping codons to amino acids.
}
\description{
Returns the codon-to-amino acid mapping for a specified genetic code (codon table).
}
\examples{
get_codon_table(1) # Standard genetic code
get_codon_table(2) # Vertebrate mitochondrial code

}
\seealso{
\code{\link{get_RSCU}}, \code{\link{get_RSCU_other}}, \code{\link{get_RSCU_other2}}
}
