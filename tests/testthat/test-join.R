test_that("join works with regular columns", {
  da <- starwars[, c("name", "mass", "species")][1:10, ]
  db <- starwars[, c("name", "homeworld")]

  res_left <- left_join(da, db)
  res_right <- right_join(da, db)
  res_inner <- inner_join(da, db[1:3, ])
  res_full <- full_join(da, db[8:12,])

  Da <- as(da, "DataFrame")
  Db <- as(db, "DataFrame")

  Res_inner <- inner_join(Da, Db[1:3, ])
  Res_left <- left_join(Da, Db)
  Res_right <- right_join(Da, Db)
  Res_full <- full_join(Da, Db[8:12,])

  expect_s4_class(Res_left, "DataFrame")
  expect_s4_class(Res_right, "DataFrame")
  expect_s4_class(Res_inner, "DataFrame")
  expect_s4_class(Res_full, "DataFrame")

  expect_identical(Res_left, as(res_left, "DataFrame"))
  expect_identical(arrange(Res_right, name), arrange(as(res_right, "DataFrame"), name))
  expect_identical(Res_inner, as(res_inner, "DataFrame"))
  expect_identical(arrange(Res_full, name), arrange(as(res_full, "DataFrame"), name))
})
