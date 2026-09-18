source("C:/Users/adf44/source/r/frmtmb-wt-adefects/dev/adefects-rev-prelude.R")
suppressPackageStartupMessages({ library(frmtmb); library(testthat) })

addr <- function(x) {
  strsplit(trimws(utils::capture.output(.Internal(inspect(x)))[1L]),
           "[[:space:]]+")[[1L]][1L]
}

cat("== the recorded 3,696,708 figure, inside a test frame\n")
test_that("the terms environment drags the calling frame in", {
  set.seed(20260917)
  A <- diag(6); A[A == 0] <- 0.3; dimnames(A) <- list(1:6, 1:6)
  dd <- data.frame(g = factor(rep(1:6, each = 5)), x = rnorm(30))
  dd$y <- dd$x + rnorm(6)[dd$g] + rnorm(30)
  fa <- frm(y ~ x + (1 | gr(g, cov = A)), dd, data2 = list(A = A))
  cat("  serialize(fa$data) inside test_that :",
      length(serialize(fa$data, NULL)), "\n")
  cat("  columns alone                       :",
      length(serialize(lapply(fa$data, identity), NULL)), "\n")
  no_data <- fa; no_data$data <- NULL
  cat("  delta                               :",
      length(serialize(fa, NULL)) - length(serialize(no_data, NULL)), "\n")

  cat("\n== the sharing control: is data.frame(fa$data) an EQUAL copy?\n")
  cp_df <- data.frame(fa$data)
  cat("  identical(fa$data, data.frame(fa$data)):",
      identical(fa$data, cp_df), "\n")
  cat("  attributes dropped:",
      paste(setdiff(names(attributes(fa$data)), names(attributes(cp_df))),
            collapse = ", "), "\n")
  cp_true <- unserialize(serialize(fa$data, NULL))
  cat("  a TRUE equal copy, identical():",
      identical(fa$data, cp_true), "\n")
  cat("  its address differs           :",
      !identical(addr(fa$data), addr(cp_true)), "\n")
  expect_true(TRUE)
})
