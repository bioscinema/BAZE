test_that("Data loads correctly", {
  data(age.rel)
  expect_true(inherits(age.rel, "phyloseq"))
})
