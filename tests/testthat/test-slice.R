test_that("i and j slicing works", {
    m <- as(mtcars, "DataFrame")

    s <- m[1:4, 1:3]
    expect_s4_class(s, "DataFrame")
    expect_identical(nrow(s), 4L)
    expect_identical(ncol(s), 3L)
    expect_identical(names(s), c("mpg", "cyl", "disp"))

    s <- m[1:4, c("disp", "wt")]
    expect_s4_class(s, "DataFrame")
    expect_identical(nrow(s), 4L)
    expect_identical(ncol(s), 2L)
    expect_identical(names(s), c("disp", "wt"))

    s <- m[c("Datsun 710", "Mazda RX4"), 1:3]
    expect_s4_class(s, "DataFrame")
    expect_identical(nrow(s), 2L)
    expect_identical(ncol(s), 3L)
    expect_identical(names(s), c("mpg", "cyl", "disp"))
})

test_that("i only slicing works", {
  m <- as(mtcars, "DataFrame")

  s <- m[1:4]
  expect_s4_class(s, "DataFrame")
  expect_identical(nrow(s), 32L)
  expect_identical(ncol(s), 4L)
  expect_identical(names(s), c("mpg", "cyl", "disp", "hp"))

  s <- m[c("disp", "wt")]
  expect_s4_class(s, "DataFrame")
  expect_identical(nrow(s), 32L)
  expect_identical(ncol(s), 2L)
  expect_identical(names(s), c("disp", "wt"))
})
