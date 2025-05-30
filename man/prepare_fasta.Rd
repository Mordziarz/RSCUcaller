% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/prepare_fasta.R
\name{prepare_fasta}
\alias{prepare_fasta}
\title{Prepare FASTA files for RSCU analysis}
\usage{
prepare_fasta(samples_table = samples_table, path = "", file_out = "")
}
\arguments{
\item{samples_table}{A data frame specifying input sequences. Two formats are supported:
- Mode 1: Contains columns "sequence_path" (paths to individual FASTA files) 
  and "sample_name" (unique identifiers starting with numbers)
- Mode 2: Contains columns "ID" (unique sample IDs) and "GENBANK_ACCESSION" 
  (NCBI accession numbers) when using `path` parameter}

\item{path}{Path to master FASTA file (required for Mode 2 when using GenBank accessions)}

\item{file_out}{Output file path with .fasta extension}
}
\value{
Invisibly returns NULL. Writes processed FASTA file to disk.
}
\description{
Processes input sequences from multiple sources into standardized FASTA format
required for downstream analysis. Handles two input modes: multiple FASTA files
or a single FASTA file with GenBank accessions.
}
\details{
This function handles two common input scenarios:
1. Multiple individual FASTA files: Requires "sequence_path" and "sample_name" columns
2. Single FASTA with multiple accessions: Requires "ID" and "GENBANK_ACCESSION" columns

For Mode 2, the function will extract sequences from the master FASTA file using
the provided GenBank accessions. Sample IDs must be unique and follow numeric prefix
format (e.g., "1_Human", "2_Chimp").
}
\examples{
path1 <- system.file("extdata/OL654070.fasta", package = "RSCUcaller")
path2 <- system.file("extdata/OQ280817.txt", package = "RSCUcaller")
path3 <- system.file("extdata/OQ280829.txt", package = "RSCUcaller")
path4 <- system.file("extdata/OQ280824.txt", package = "RSCUcaller")
samples_table_cp <- data.frame(sequence_path = c(path1,
                                             path2,
                                             path3,
                                             path4),
                           sample_name = c("1_Apopellia_endiviifolia_A1",
                                           "2_Apopellia_endiviifolia_C1",
                                           "3_Pellia_neesiana_1",
                                           "4_Pellia_epiphylla_S1")) 


# Mode 1: Multiple FASTA files
withr::with_tempdir({prepare_fasta(samples_table_cp, file_out = "processed.fasta")})

# Mode 2: Single FASTA with GenBank accessions
mit_fasta <- system.file("extdata/mitogenome_sequence.txt", package = "RSCUcaller")
samples_path <- system.file("extdata/samples_table.csv", package = "RSCUcaller")
samples_table <- read.csv2(samples_path,sep = ";")
withr::with_tempdir({prepare_fasta(samples_table, 
              path = mit_fasta,
              file_out = "processed.fasta")})

}
\references{
For FASTA format specifications:
- NCBI FASTA format: \url{https://blast.ncbi.nlm.nih.gov/doc/fasta.html}
- Seqinr documentation: \url{https://seqinr.r-forge.r-project.org}
}
