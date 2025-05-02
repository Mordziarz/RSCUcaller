context("Testing visualization functions")

test_that("heatmap_RSCU creates correct plots", {

  withr::with_tempdir({
  data("prepared_fasta", package = "RSCUcaller")
  
  rscu_results <- get_RSCU(merged_sequences = prepared_fasta)
  expect_error(result_heatmap <- heatmap_RSCU(get_RSCU_out = rscu_results, 
                                              select = "heatmap", 
                                              heatmap_color = "red_blue"), NA)
  expect_s4_class(result_heatmap, "HeatmapList")
  expect_error(result_dendogram <- heatmap_RSCU(get_RSCU_out = rscu_results, 
                                                select = "dendogram"), NA)
  expect_s3_class(result_dendogram, "gg")
})})

test_that("histogram_RSCU creates correct histograms", {
  
  data("prepared_fasta", package = "RSCUcaller")
  
  rscu_results <- get_RSCU(merged_sequences = prepared_fasta)
  expect_error(result_histogram <- histogram_RSCU(get_RSCU_out = rscu_results, 
                                                  title = "Test Histogram"), NA)
  expect_is(result_histogram, "gg")
})

test_that("histogram_RSCU_double creates double histograms", {
  
  data("prepared_fasta", package = "RSCUcaller")
  
  rscu_results <- get_RSCU(merged_sequences = prepared_fasta)
  expect_error(result_double_histogram <- histogram_RSCU_double(
    get_RSCU_out_left = rscu_results,
    get_RSCU_out_right = rscu_results,
    title_left = "Test Left",
    title_right = "Test Right"), NA)
  expect_is(result_double_histogram, "patchwork")
})
