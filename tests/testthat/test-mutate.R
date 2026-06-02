test_that("mutate works with regular columns", {
    d <- S4Vectors::DataFrame(mtcars)

    m <- mutate(d, newvar = cyl + mpg)
    expect_s4_class(m, "DataFrame")
    expect_true("newvar" %in% names(m))
    expect_identical(m$newvar, mtcars$cyl + mtcars$mpg)

    m <- mutate(d, newvar2 = cyl^2)
    expect_s4_class(m, "DataFrame")
    expect_true("newvar2" %in% names(m))
    expect_identical(m$newvar2, mtcars$cyl^2)

    carletters <- c(LETTERS, LETTERS[1:6])
    m <- mutate(d, newvar3 = paste0(carletters, cyl))
    expect_s4_class(m, "DataFrame")
    expect_true("newvar3" %in% names(m))
    expect_identical(m$newvar3, paste0(carletters, mtcars$cyl))

    m <- mutate_at(d, vars(starts_with("c")), ~ .^2)
    expect_s4_class(m, "DataFrame")
    expect_identical(names(m), names(d))
    expect_identical(m$cyl, mtcars$cyl^2)
    expect_identical(m$carb, mtcars$carb^2)
    expect_identical(m$mpg, mtcars$mpg)
})

test_that("mutate works with S4 columns", {
    skip_if_not_installed("GenomicRanges")
    skip_if_not_installed("IRanges")
    skip_if_not_installed("GenomeInfoDb")

    D <- S4Vectors::DataFrame(mtcars)
    D$gr <- GenomicRanges::GRanges("chrX", IRanges::IRanges(1:32, width = 10))

    m <- mutate(D, chr = GenomeInfoDb::seqnames(gr))
    expect_s4_class(m, "DataFrame")
    expect_identical(nrow(m), 32L)
    expect_true("chr" %in% names(m))
    expect_s4_class(m$chr, "Rle")
    expect_identical(m$chr, IRanges::RleList(factor(rep("chrX", 32)))[[1]])
})

test_that("mutate adds columns sequentially", {
  d <- S4Vectors::DataFrame(mtcars)
  m <- mutate(d, newvar = cyl * 2, newervar = newvar + am)
  expect_s4_class(m, "DataFrame")
  expect_true("newvar" %in% names(m))
  expect_true("newervar" %in% names(m))
  expect_identical(m$newvar, mtcars$cyl * 2)
  expect_identical(m$newervar, mtcars$am + mtcars$cyl * 2)
})


test_that("sequential mutation supports groups", {
  d <- S4Vectors::DataFrame(mtcars) |>
    group_by(gear)
  m <- mutate(d,
              newvar = cyl * 2,
              avggear = mean(gear) + 1,
              newervar = vs + newvar + avggear)
  expect_s4_class(m, "DataFrame")
  expect_s4_class(m, "GroupedDataFrame")
  expect_true("newvar" %in% names(m))
  expect_true("avggear" %in% names(m))
  expect_true("newervar" %in% names(m))

  # realign rows after grouping
  mm <- m[rownames(mtcars), ]
  expect_identical(mm$newvar, mtcars$cyl * 2)
  expect_identical(mm$avggear, mtcars$gear + 1)
  expect_identical(mm$newervar,
                   mtcars$vs +
                     (mtcars$gear + 1) +
                     (mtcars$cyl * 2))
})
