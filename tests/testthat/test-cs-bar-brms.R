# `(cs(1) | g)` with brms attached after frmtmb (lane surface,
# 2026-10-06; the release review of 0.68.0, m10).
#
# Seen to fail on 0.68.1 (rellib-r6): with brms's cs() on the search
# path, parse_linpred() protected the call as .frm_cs(1), the refusal
# did not recognize it, and model.frame() stopped with the internal
# "variable lengths differ (found for '.frm_cs(1)')"; without brms the
# refusal was reached (dev/surface-repros.R item 6, dev/relrev-cs1g.R).
#
# A child process, because the defect is a property of the search path,
# which test-generic-collision.R explains.

refusal_cs1 <- paste("A category-specific effect inside a group-level",
                     "term, (cs(1) | g)")

run_child_cs <- function(code) {
  f <- tempfile(fileext = ".R")
  on.exit(unlink(f), add = TRUE)
  head <- sprintf(".libPaths(%s)",
                  paste0(deparse(.libPaths()), collapse = ""))
  writeLines(c(head, code), f)
  out <- suppressWarnings(system2(
    file.path(R.home("bin"), "Rscript"),
    c("--vanilla", shQuote(f)), stdout = TRUE, stderr = TRUE))
  paste(out, collapse = "\n")
}

cs_child <- function(attach, family) {
  c(attach,
    "set.seed(1)",
    "d <- data.frame(x = rnorm(200),",
    "                g = factor(sample(letters[1:6], 200, TRUE)))",
    "d$y <- sample(1:4, 200, TRUE)",
    sprintf(paste0("r <- tryCatch({frm(y ~ x + (cs(1) | g), ",
                   "family = %s(), data = d); 'FIT'}, ",
                   "error = conditionMessage)"), family),
    "cat('RESULT:', r, '\\n')",
    "cat('CHILDOK\\n')")
}

attach_both <- c("suppressPackageStartupMessages(library(frmtmb))",
                 "suppressPackageStartupMessages(library(brms))")

test_that("(cs(1) | g) gets the refusal with brms attached after frmtmb", {
  skip_on_cran()
  skip_if_not_installed("brms")
  for (fam in c("sratio", "cumulative")) {
    out <- run_child_cs(cs_child(attach_both, fam))
    expect_match(out, "CHILDOK", fixed = TRUE)
    expect_match(out, refusal_cs1, fixed = TRUE)
    expect_false(grepl(".frm_cs", out, fixed = TRUE))
  }
})

test_that("(cs(1) | g) gets the same refusal without brms (control)", {
  skip_on_cran()
  out <- run_child_cs(cs_child(
    "suppressPackageStartupMessages(library(frmtmb))", "sratio"))
  expect_match(out, "CHILDOK", fixed = TRUE)
  expect_match(out, refusal_cs1, fixed = TRUE)
})

test_that("a family without cs() gets brms's words with brms attached", {
  skip_on_cran()
  skip_if_not_installed("brms")
  out <- run_child_cs(cs_child(attach_both, "gaussian"))
  expect_match(out, paste("Category specific effects are not supported",
                          "for this family"), fixed = TRUE)
})
