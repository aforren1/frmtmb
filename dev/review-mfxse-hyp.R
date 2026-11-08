# Review of lane mfxse, attack 3: hypotheses() on a mixed fit, both builds.
#   Rscript dev/review-mfxse-hyp.R [lane|before]
# Data seed 7411 (the lane test's mre_data()).
a <- commandArgs(TRUE)
arm <- if (length(a)) a[1] else "lane"
lane <- "C:/Users/adf44/source/r/wt-mfxse-lib"
base <- "C:/Users/adf44/source/r/rellib-r18"
user <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
.libPaths(switch(arm, lane = c(lane, base, user), before = c(base, user),
                 c(arm, base, user)))
Sys.setenv(NOT_CRAN = "true")
suppressPackageStartupMessages({library(frmtmb); library(marginaleffects)})
cat("lib:", find.package("frmtmb"), format(packageVersion("frmtmb")), "arm:", arm, "\n")
set.seed(7411)
J <- 12L; m <- 20L; n <- J * m
d <- data.frame(x = rnorm(n), f = factor(sample(c("a", "b"), n, TRUE)),
                g = factor(rep(seq_len(J), each = m)))
u0 <- rnorm(J, 0, 0.8); u1 <- rnorm(J, 0, 0.5)
d$cnt <- rpois(n, exp(0.5 + (0.4 + u1[d$g]) * d$x + 0.3 * (d$f == "b") + u0[d$g]))
fit <- suppressWarnings(frm(cnt ~ x + f + (1 + x | g), family = poisson(), data = d))
show <- function(lab, e) {
  ws <- character(0)
  r <- withCallingHandlers(tryCatch(e, error = function(err) {
    paste("ERROR:", conditionMessage(err))
  }), warning = function(w) {
    ws <<- c(ws, substr(conditionMessage(w), 1, 90))
    invokeRestart("muffleWarning")
  })
  cat("\n---", lab, "\n")
  if (is.character(r)) cat(r, "\n") else {
    df <- as.data.frame(r)
    keep <- intersect(c("term", "hypothesis", "estimate", "std.error",
                        "statistic", "p.value", "df1", "df2"), names(df))
    cat("rows:", nrow(df), "\n")
    print(utils::head(df[, keep, drop = FALSE], 6), digits = 6)
  }
  if (length(ws)) cat("warnings:", paste(unique(ws), collapse = " || "), "\n")
  invisible(r)
}
cat("get_coef length:", length(get_coef(fit)), "\n")
show("hypotheses(fit)", hypotheses(fit))
show("hypotheses(fit, 'b2 = 0')", hypotheses(fit, "b2 = 0"))
show("hypotheses(fit, 'b2 = b3')", hypotheses(fit, "b2 = b3"))
show("hypotheses(fit, hypothesis = ~ reference)", hypotheses(fit, hypothesis = ~ reference))
show("hypotheses(fit, hypothesis = ~ pairwise)", hypotheses(fit, hypothesis = ~ pairwise))
show("hypotheses(fit, hypothesis = ~ I(b2 - b3))", hypotheses(fit, hypothesis = ~ I(b2 - b3)))
show("hypotheses(fit, joint = TRUE)", hypotheses(fit, joint = TRUE))
show("hypotheses(fit, joint = c(2, 3))", hypotheses(fit, joint = c(2, 3)))
show("hypotheses(fit, joint = 'x')", hypotheses(fit, joint = "x"))
show("hypotheses(avg_slopes(fit))", hypotheses(avg_slopes(fit)))
show("hypotheses(avg_slopes(fit), 'b1 = b2')", hypotheses(avg_slopes(fit), "b1 = b2"))
show("hypotheses(avg_slopes(fit), joint = TRUE)", hypotheses(avg_slopes(fit), joint = TRUE))
show("hypotheses(fit, vcov = ~ g)", hypotheses(fit, vcov = ~ g))
show("hypotheses(fit, function)", hypotheses(fit, hypothesis = function(x) {
  data.frame(term = "sum", estimate = sum(coef(x)[1:2]))
}))
cat("\n--- reference: fixef\n"); print(fixef(fit))
# a leak of the row check: get_vcov() called alone arms it; the next
# predictions(vcov = FALSE) at a new level then warns where the help says
# vcov = FALSE gives no warning
nd <- data.frame(x = 0, f = "a", g = factor("new"))
invisible(get_vcov(fit))
show("predictions(vcov = FALSE, new level) after get_vcov()",
     predictions(fit, newdata = nd, vcov = FALSE, allow_new_levels = TRUE))
show("again, nothing armed", predictions(fit, newdata = nd, vcov = FALSE,
                                        allow_new_levels = TRUE))
# inferences(method = "simulation") draws over get_coef() x get_vcov()
set.seed(1)
show("inferences(simulation) on predictions",
     inferences(predictions(fit, newdata = d[1:3, ]), method = "simulation", R = 500))
p0 <- predictions(fit, newdata = d[1:3, ])
cat("delta se:", format(p0$std.error, digits = 5), "\n")
