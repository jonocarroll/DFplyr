test_that("ungrouped show method works", {
    d <- S4Vectors::DataFrame(mtcars)
    out <- capture.output(d)
    expect_identical(out[1], "DataFrame with 32 rows and 11 columns")
    expect_identical(substr(out[2], 1, 27),  "                        mpg")
})

test_that("grouped show method works", {
    d <- S4Vectors::DataFrame(mtcars)
    dg <- group_by(d, am)
    out <- capture.output(dg)
    expect_identical(out[1], "DataFrame with 32 rows and 11 columns")
    expect_identical(out[2], "Groups:  am ")
    expect_identical(substr(out[3], 1, 27),  "                        mpg")
})
