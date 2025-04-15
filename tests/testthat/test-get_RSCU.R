context("Testing RSCU calculation functions with prepared_fasta")
  
test_that("get_RSCU", {
  
  # Load test data
  data("prepared_fasta", package = "RSCUcaller")
  
  # Basic functionality test
  expect_error(
    rscu_std <- get_RSCU(merged_sequences = prepared_fasta),
    NA
  )
  
  # Check output structure
  expect_s3_class(rscu_std, "data.frame")
  expect_true(all(c("RSCU", "AA") %in% colnames(rscu_std)))
  
  # Check for non-negative values
  expect_true(all(rscu_std[,5] >= 0))
})
