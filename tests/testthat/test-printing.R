test_that("show method works", {
  d <- S4Vectors::DataFrame(mtcars)
  dg <- group_by(d, am)
  out <- capture.output(dg)
  expect_identical(out[1], "DataFrame with 32 rows and 11 columns")
  expect_identical(out[2], "Groups:  am ")
})
