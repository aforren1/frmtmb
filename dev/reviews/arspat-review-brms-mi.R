# Review of lane arspat: mi() cases the lane did not run, against brms
# 2.23.0's program integrated over Ymi, with the lane's own machinery
# (dev/arspat-l1b.R minus its driver loop) on new cases.
ex <- parse("dev/arspat-l1b.R")
drv <- vapply(ex, function(e) {
  txt <- paste(deparse(e), collapse = " ")
  grepl("^which_cases|^if \\(!length\\(which_cases|^for \\(nm in which_cases",
        txt)
}, NA)
for (e in ex[!drv]) eval(e, envir = globalenv())

cases <- list(
  mi_arma_se_edges = list(
    bf = "y | mi() + se(s, sigma = TRUE) ~ x + arma(t, g, cov = TRUE)",
    mi = c(1, 2, 6, 7, 12, 13, 14), marg = "laplace"),
  mi_ma_het_consec = list(
    bf = "bf(y | mi() ~ x + ma(t, g, cov = TRUE), sigma ~ x2)",
    mi = c(3, 4, 5, 30, 36), marg = "laplace"),
  mi_ar_student_edges = list(
    bf = "y | mi() ~ x + ar(t, g, cov = TRUE)", mi = c(1, 12, 13, 24),
    family = "student()", marg = "quad", tdf = 4, ng = 30),
  mi_unstr_se_edges = list(
    bf = "y | mi() + se(s, sigma = TRUE) ~ x + unstr(t, g)",
    mi = c(1, 6, 7, 8, 18), marg = "laplace")
)
which_cases <- commandArgs(trailingOnly = TRUE)
if (!length(which_cases)) which_cases <- names(cases)
for (nm in which_cases) {
  cat("\n######## ", nm, "\n", sep = "")
  r <- tryCatch(run_case(nm), error = function(e) {
    cat("ERROR: ", conditionMessage(e), "\n")
  })
}
