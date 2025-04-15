context("Testing prepare_fasta function")

test_that("prepare_fasta works with table containing sequence paths", {
  temp_file <- tempfile(fileext = ".fasta")
  write(">TEST_SEQ1\nATGGCTAGCTAG", temp_file)
  
  samples_table <- data.frame(
    sequence_path = rep(temp_file, 2),
    sample_name = c("1_test_sample", "2_test_sample")
  )
  
  output_file <- tempfile(fileext = ".fasta")
  
  expect_error(prepare_fasta(samples_table = samples_table, 
                             file_out = output_file), NA)
  expect_true(file.exists(output_file))
  file_content <- readLines(output_file)
  expect_true(any(grepl(">1_test_sample", file_content)))
  expect_true(any(grepl(">2_test_sample", file_content)))
  unlink(temp_file)
  unlink(output_file)
})

test_that("prepare_fasta works with table containing ID and GENBANK_ACCESSION", {
  temp_file <- tempfile(fileext = ".fasta")
  write(">gb_acc_1\nATGGCTAGCTAG\n>gb_acc_2\nATGGCTAGCTAG", temp_file)
  
  samples_table <- data.frame(
    ID = c("1_test_sample", "2_test_sample"),
    GENBANK_ACCESSION = c("gb_acc_1", "gb_acc_2")
  )
  
  output_file <- tempfile(fileext = ".fasta")
  
  expect_error(prepare_fasta(samples_table = samples_table, 
                             path = temp_file, 
                             file_out = output_file), NA)
  expect_true(file.exists(output_file))
  file_content <- readLines(output_file)
  expect_true(any(grepl(">1_test_sample", file_content)))
  expect_true(any(grepl(">2_test_sample", file_content)))
  unlink(temp_file)
  unlink(output_file)
})
