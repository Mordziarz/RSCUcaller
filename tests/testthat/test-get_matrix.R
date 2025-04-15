context("Testing matrix creation functions")

test_that("get_matrix creates matrix from get_RSCU results", {

  data("prepared_fasta", package = "RSCUcaller")
  
  rscu_results <- get_RSCU(merged_sequences = prepared_fasta)
  expect_error(matrix_result <- get_matrix(get_RSCU_out = rscu_results), NA)
  expect_is(matrix_result, "data.frame")
  expected_rows <- length(unique(rscu_results$codon))
})

test_that("get_matrix throws error for invalid input", {
  invalid_df <- data.frame(
    Invalid_Column = c("A", "B")
  )
  expect_error(get_matrix(get_RSCU_out = invalid_df))
  empty_df <- data.frame()
  expect_error(get_matrix(get_RSCU_out = empty_df))
})
