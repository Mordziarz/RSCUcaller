test_that("Test get_RSCU and get_matrix functions", {
  
  library(seqinr)
  library(RSCUcaller)
  library(dplyr)
  
  # Load test data
  url <- "https://raw.githubusercontent.com/Mordziarz/RSCUcaller/main/test_data/OQ280817.txt"
  
  # Run function
  get_RSCU_out <- get_RSCU(merged_sequences = url)
  get_matrix_out <- get_matrix(get_RSCU_out = get_RSCU_out)
  
  # Ensure output is a dataframe
  expect_true(inherits(get_RSCU_out, "data.frame"))
  expect_true(inherits(get_matrix_out, "data.frame"))
  
  # 🔻 Clean up any files generated (csv, rda)
  unlink(list.files(pattern = "\\.csv$|\\.rda$", full.names = TRUE), force = TRUE)
})
