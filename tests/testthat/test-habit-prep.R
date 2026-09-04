test_that("habit_prep has the documented shape", {
  data(habit_prep, package = "frmtmb", envir = environment())

  expect_s3_class(habit_prep, "data.frame")
  expect_identical(nrow(habit_prep), 28940L)
  expect_identical(
    names(habit_prep),
    c("participant", "group", "subject_code", "trial", "prep_time",
      "remapped", "stimulus", "target_key", "response_key", "response")
  )
  expect_identical(levels(habit_prep$group), c("minimal", "4day", "20day"))
  expect_identical(levels(habit_prep$response), c("correct", "habit", "other"))
  expect_type(habit_prep$remapped, "logical")

  # 58 subject-conditions: experiment 1 measured twice, experiment 2 once
  cells <- unique(habit_prep[c("participant", "group")])
  expect_identical(nrow(cells), 58L)
  expect_identical(nlevels(habit_prep$participant), 36L)
})

test_that("habit_prep withholds the response category of unreadable trials", {
  data(habit_prep, package = "frmtmb", envir = environment())

  # a recorded key outside 1:4 is a hardware read failure, so the category is
  # not trustworthy and must not be scored as an ordinary error
  bad <- habit_prep$response_key < 1L | habit_prep$response_key > 4L
  expect_identical(sum(bad), 8L)
  expect_true(all(is.na(habit_prep$response[bad])))
  expect_false(any(is.na(habit_prep$response[!bad])))

  expect_false(any(is.na(habit_prep$prep_time)))
  expect_true(all(habit_prep$target_key %in% 1:4))
})

test_that("habit_prep splits the two stimulus classes by trial index", {
  data(habit_prep, package = "frmtmb", envir = environment())

  # every trial belongs to exactly one class, so the two counts partition the
  # rows; matching on preparation-time value instead would overlap them
  expect_identical(sum(habit_prep$remapped) + sum(!habit_prep$remapped),
                   nrow(habit_prep))
  by_cell <- tapply(habit_prep$remapped,
                    interaction(habit_prep$participant, habit_prep$group,
                                drop = TRUE), sum)
  expect_true(all(by_cell <= 250L))
  expect_true(all(by_cell >= 240L))
})
