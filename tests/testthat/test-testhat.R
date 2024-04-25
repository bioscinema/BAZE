test_that("Data loads correctly", {
  data(ps)
  expect_true(inherits(ps, "phyloseq"))
})
