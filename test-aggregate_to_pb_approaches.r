test_that("assess_compute_corr", {
  expected <- to_pseudobulk(hagai_toy, "sum")
  expect_s3_class(expected, data)
})
