# Reviewer, priority 2: do predict()'s ndraws, re_formula,
# allow_new_levels and newdata behave as brms's do?
#
#   Rscript dev/shapes-rev-predargs.R

.libPaths(c("C:/Users/adf44/source/r/shapes-lib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
TREE <- "C:/Users/adf44/source/r/frmtmb-wt-shapes"
source(file.path(TREE, "dev/shapes-rev-fixtures.R"))
dd <- rev_data()
fm <- frm(bf(ymix ~ x + (1 | g)) + gaussian(), data = dd)
fg <- frm(bf(y ~ x + f) + gaussian(), data = dd)

p <- function(lab, e) {
  r <- tryCatch(suppressWarnings(e), error = function(c)
    paste0("ERROR: ", substr(conditionMessage(c), 1, 110)))
  cat(sprintf("  %-52s %s\n", lab,
              if (is.character(r) && length(r) == 1) r else
                paste(class(r)[1], paste(dim(r) %||% length(r),
                                         collapse = "x"))))
  invisible(r)
}
`%||%` <- function(a, b) if (is.null(a)) b else a

cat("== ndraws ==\n")
for (k in c(1L, 2L, 7L, 1000L))
  p(sprintf("dim(predict(summary = FALSE, ndraws = %d))", k),
    dim(predict(fg, summary = FALSE, ndraws = k)))
p("ndraws = 0", predict(fg, ndraws = 0))
p("ndraws = -1", predict(fg, ndraws = -1))
p("ndraws = 2.5", predict(fg, ndraws = 2.5))
p("ndraws = 'ten'", predict(fg, ndraws = "ten"))

cat("\n== does ndraws actually change the spread? (a bigger ndraws\n")
cat("   must not change the ANSWER systematically, only its noise)\n")
set.seed(1); a <- predict(fg, ndraws = 4000)[, "Est.Error"]
set.seed(2); b <- predict(fg, ndraws = 4000)[, "Est.Error"]
set.seed(3); c50 <- predict(fg, ndraws = 50)[, "Est.Error"]
cat(sprintf("   ndraws=4000 twice: max rel gap %.4f\n",
            max(abs(a / b - 1))))
cat(sprintf("   ndraws=50 vs 4000: max rel gap %.4f\n",
            max(abs(c50 / a - 1))))

cat("\n== re_formula ==\n")
p("re_formula = NULL  (all REs, brms default)",
  predict(fm, ndraws = 100))
p("re_formula = NA    (population level)",
  predict(fm, ndraws = 100, re_formula = NA))
p("re_formula = ~0", predict(fm, ndraws = 100, re_formula = ~0))
p("re_formula = ~(1|g)", predict(fm, ndraws = 100, re_formula = ~ (1 | g)))
p("re_formula = ~(1|nosuch)",
  predict(fm, ndraws = 100, re_formula = ~ (1 | nosuch)))
set.seed(11); e1 <- predict(fm, ndraws = 2000)[, "Estimate"]
set.seed(11); e2 <- predict(fm, ndraws = 2000, re_formula = NA)[, "Estimate"]
cat(sprintf("   NULL vs NA Estimate: max abs gap %.5f (must be > 0:\n",
            max(abs(e1 - e2))))
cat("     a population-level prediction drops the group modes)\n")
set.seed(12); s1 <- predict(fm, ndraws = 2000)[, "Est.Error"]
set.seed(12); s2 <- predict(fm, ndraws = 2000, re_formula = NA)[, "Est.Error"]
cat(sprintf("   NULL vs NA Est.Error: mean %.5f vs %.5f\n",
            mean(s1), mean(s2)))

cat("\n== allow_new_levels ==\n")
nd <- data.frame(x = c(0, 1), g = factor(c("99", "98")))
p("newdata with an unseen level, default (must refuse)",
  predict(fm, newdata = nd, ndraws = 50))
p("allow_new_levels = TRUE", predict(fm, newdata = nd, ndraws = 50,
                                     allow_new_levels = TRUE))
r <- suppressWarnings(tryCatch(predict(fm, newdata = nd, ndraws = 4000,
                                       allow_new_levels = TRUE),
                               error = function(c) NULL))
if (!is.null(r)) {
  nd2 <- data.frame(x = c(0, 1), g = factor(c("1", "2"),
                                            levels = levels(dd$g)))
  r2 <- predict(fm, newdata = nd2, ndraws = 4000)
  cat("   unseen level Est.Error: ", round(r[, "Est.Error"], 4), "\n")
  cat("   seen   level Est.Error: ", round(r2[, "Est.Error"], 4), "\n")
  cat("   (brms draws a fresh effect for an unseen level, so ITS\n")
  cat("    interval is wider; here the level is predicted at the\n")
  cat("    population level, so check whether the spread grows)\n")
}

cat("\n== newdata ==\n")
p("newdata = head(dd, 3)", predict(fg, newdata = head(dd, 3), ndraws = 50))
p("newdata = dd[0, ]", predict(fg, newdata = dd[0, ], ndraws = 50))
p("newdata missing a predictor",
  predict(fg, newdata = data.frame(x = 1), ndraws = 50))
p("newdata = list(...)", predict(fg, newdata = list(x = 1, f = "a"),
                                 ndraws = 50))
p("newdata with an unseen factor level of f",
  predict(fg, newdata = data.frame(x = 0, f = factor("zz")), ndraws = 50))

cat("\n== transform, robust, probs ==\n")
p("transform = log", predict(fg, ndraws = 50, transform = log))
p("robust = TRUE", predict(fg, ndraws = 50, robust = TRUE))
p("probs = NULL", predict(fg, ndraws = 50, probs = NULL))
p("probs = c(0.1, 0.5, 0.9)", predict(fg, ndraws = 50,
                                      probs = c(0.1, 0.5, 0.9)))
cat("   colnames at probs = c(0.1,0.5,0.9): ",
    paste(colnames(predict(fg, ndraws = 50,
                           probs = c(0.1, 0.5, 0.9))), collapse = ", "),
    "\n")
p("probs = 2", predict(fg, ndraws = 50, probs = 2))
