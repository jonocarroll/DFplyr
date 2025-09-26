test_that("join preserves groups when possible", {
  da <- starwars[, c("name", "eye_color", "height", "mass")][1:10, ]
  db <- starwars[, c("name", "eye_color", "homeworld")]

  Da <- as(da, "DataFrame")
  Db <- as(db, "DataFrame")

  da <- group_by(da, eye_color)
  Da <- group_by(Da, eye_color)

  # groups preserved
  res_left <- left_join(da, db)
  res_right <- right_join(da, db)
  res_inner <- inner_join(da, db[1:3, ])
  res_full <- full_join(da, db[1:3, ])

  Res_left <- left_join(Da, Db)
  Res_right <- right_join(Da, Db)
  Res_inner <- inner_join(Da, Db[1:3, ])
  Res_full <- full_join(Da, Db[1:3, ])

  expect_s4_class(Res_left, "DataFrame")
  expect_s4_class(Res_right, "DataFrame")
  expect_s4_class(Res_inner, "DataFrame")
  expect_s4_class(Res_full, "DataFrame")

  expect_identical(Res_left, as(res_left, "DataFrame") %>% group_by(eye_color))
  expect_identical(
    arrange(Res_right, name),
    arrange(as(res_right, "DataFrame"), name) %>% group_by(eye_color)
  )
  expect_identical(
    Res_inner,
    as(res_inner, "DataFrame") %>% group_by(eye_color)
  )
  expect_identical(
    arrange(Res_full, name),
    arrange(as(res_full, "DataFrame"), name)
  )

  # groups altered
  res_left <- left_join(da, db, by = "name")
  res_right <- right_join(da, db, by = "name")
  res_inner <- inner_join(da, db[1:3, ], by = "name")
  res_full <- full_join(da, db[1:3, ], by = "name")

  # these currently fail because they try to group by a non-existent column
  Res_left <- left_join(Da, Db, by = "name")
  Res_right <- right_join(Da, Db, by = "name")
  Res_inner <- inner_join(Da, Db[1:3, ], by = "name")
  Res_full <- full_join(Da, Db[1:3, ], by = "name")

  expect_s4_class(Res_left, "DataFrame")
  expect_s4_class(Res_right, "DataFrame")
  expect_s4_class(Res_inner, "DataFrame")
  expect_s4_class(Res_full, "DataFrame")
})
