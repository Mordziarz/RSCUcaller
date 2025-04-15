test_that("prepare_fasta processing", {

  library(seqinr)
  library(RSCUcaller)

  # Load test data
  data("samples_table_mt", package = "RSCUcaller")
  url <- "https://raw.githubusercontent.com/Mordziarz/RSCUcaller/main/test_data/mitogenome_sequence.txt"
  
  # Run function
  prepare_fasta(samples_table = samples_table,path = url ,file_out = "post_fasta.fasta")
  
  # Check output
  expect_true(file.exists("post_fasta.fasta"))
  
  # 🔻 Clean up any files generated (csv, rda, fasta)
  unlink(list.files(pattern = "\\.csv$|\\.rda$|\\.fasta$", full.names = TRUE), force = TRUE)
})
