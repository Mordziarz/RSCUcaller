context("Testing statistical functions")

test_that("correlation calculates correct correlation", {
  
  data("prepared_fasta", package = "RSCUcaller")
  
  rscu_results <- get_RSCU(merged_sequences = prepared_fasta)
  expect_error(result_corr <- correlation(
    get_RSCU_out = rscu_results,
    Species_x = "Apopellia_endiviifolia_B1",
    Species_y = "Apopellia_endiviifolia_A1",
    xlab = "Test X",
    ylab = "Test Y"), NA)
  expect_is(result_corr, "list")
})

test_that("neutrality_pr2 creates PR2 and neutrality plots", {
  
  data("prepared_fasta", package = "RSCUcaller")
  data("grouping_table", package = "RSCUcaller")
  
  rscu_results <- get_RSCU(merged_sequences = prepared_fasta)

  expect_error(result_pr2 <- neutrality_pr2(
    get_RSCU_out = rscu_results,
    select = "PR2_plot",
    grouping_table = grouping_table), NA)
  expect_is(result_pr2, "list")
  expect_error(result_neutrality <- neutrality_pr2(
    get_RSCU_out = rscu_results,
    select = "neutrality_plot",
    grouping_table = grouping_table), NA)
  expect_is(result_neutrality, "list")
})

test_that("PCA_RSCU performs PCA analysis", {

  data("prepared_fasta", package = "RSCUcaller")
  data("grouping_table", package = "RSCUcaller")
  
  rscu_results <- get_RSCU(merged_sequences = prepared_fasta)
  matrix_result <- get_matrix(get_RSCU_out = rscu_results)
  expect_error(result_pca <- PCA_RSCU(
    get_matrix_out = matrix_result,
    grouping_table = grouping_table), NA)
  expect_is(result_pca, "gg")
})
