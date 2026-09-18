source("C:/Users/adf44/source/r/frmtmb-wt-adefects/dev/adefects-rev-prelude.R")
suppressPackageStartupMessages({ library(frmtmb); library(brms) })

cat("== brms itself: validate_newdata with the grouping column removed\n")
bf1 <- get("brmsfit_example1", envir = asNamespace("brms"))
nd <- bf1$data[1:3, ]
nd$visit <- NULL
r <- tryCatch(brms:::validate_newdata(nd, bf1, allow_new_levels = TRUE),
              error = function(e) e)
if (inherits(r, "condition"))
  cat("  allow_new_levels = TRUE -> ERR:",
      substr(conditionMessage(r), 1, 90), "\n") else
  cat("  allow_new_levels = TRUE -> filled visit:",
      paste(r$visit, collapse = ","), " ncol:", ncol(r), "\n")
r2 <- tryCatch(brms:::validate_newdata(nd, bf1, allow_new_levels = FALSE),
               error = function(e) e)
cat("  allow_new_levels = FALSE ->",
    if (inherits(r2, "condition"))
      paste("ERR:", substr(conditionMessage(r2), 1, 90)) else "ACCEPTED",
    "\n")

cat("\n== frmtmb: the three routes must agree on the prediction\n")
set.seed(20260917)
n <- 80
d <- data.frame(x = rnorm(n), g = factor(rep_len(paste0("g", 1:8), n)))
d$y <- rnorm(n, 1 + 0.5 * d$x + rep_len(rnorm(8, 0, 0.7), n), 1)
fit <- frm(frmtmb::bf(y ~ x + (1 | g)), family = gaussian(), data = d)

nd0 <- data.frame(x = c(-0.5, 0.5))                      # no g at all
nd1 <- data.frame(x = c(-0.5, 0.5), g = factor(c("zz", "zz")))  # new level
p <- function(lbl, ...) {
  r <- tryCatch(predict(fit, ..., se.fit = TRUE), error = function(e) e)
  if (inherits(r, "condition")) {
    cat(sprintf("  %-38s ERR[%s]: %s\n", lbl,
                paste(class(r)[1:2], collapse = ","),
                substr(conditionMessage(r), 1, 80)))
    return(invisible(NULL))
  }
  cat(sprintf("  %-38s fit=%.10f,%.10f se=%.10f,%.10f\n", lbl,
              r$fit[1], r$fit[2], r$se.fit[1], r$se.fit[2]))
  invisible(r)
}
a <- p("no g, allow_new_levels = TRUE", newdata = nd0,
       allow_new_levels = TRUE)
b <- p("new level zz, allow_new_levels = TRUE", newdata = nd1,
       allow_new_levels = TRUE)
cc <- p("no g, re_formula = NA", newdata = nd0, re_formula = NA)
p("no g, no allow_new_levels", newdata = nd0)
if (!is.null(a) && !is.null(b)) {
  cat("  identical(fit a, fit b):", identical(a$fit, b$fit),
      " max|diff| =", format(max(abs(a$fit - b$fit)), digits = 17), "\n")
  cat("  identical(se a, se b):", identical(a$se.fit, b$se.fit),
      " max|diff| =", format(max(abs(a$se.fit - b$se.fit)), digits = 17),
      "\n")
}
if (!is.null(a) && !is.null(cc)) {
  cat("  a vs re_formula = NA: fit max|diff| =",
      format(max(abs(a$fit - cc$fit)), digits = 17),
      " se ratio =", format(a$se.fit / cc$se.fit, digits = 10), "\n")
}

cat("\n== which columns get filled, and which do not\n")
d2 <- d; d2$h <- factor(rep_len(paste0("h", 1:4), n))
f2 <- frm(frmtmb::bf(y ~ x + (1 | g) + (1 | h)), family = gaussian(), data = d2)
p2 <- function(lbl, nd, ...) {
  r <- tryCatch(predict(f2, newdata = nd, ...), error = function(e) e)
  cat(sprintf("  %-34s %s\n", lbl,
              if (inherits(r, "condition"))
                paste0("ERR: ", substr(conditionMessage(r), 1, 95)) else
                paste("OK", paste(round(r, 6), collapse = ","))))
}
p2("neither g nor h, no flag", data.frame(x = c(0, 1)))
p2("neither g nor h, flag", data.frame(x = c(0, 1)),
   allow_new_levels = TRUE)
p2("g only missing, flag",
   data.frame(x = c(0, 1), h = factor(c("h1", "h2"))),
   allow_new_levels = TRUE)
p2("re_formula = NA, neither", data.frame(x = c(0, 1)),
   re_formula = NA)
p2("re_formula = ~(1|h), g missing",
   data.frame(x = c(0, 1), h = factor(c("h1", "h2"))),
   re_formula = ~ (1 | h))

cat("\n== the guard's absent case: a grouping variable in the CALLING env\n")
gg <- factor(rep_len(paste0("g", 1:8), 2))
nd3 <- data.frame(x = c(0, 1))
r3 <- tryCatch(predict(fit, newdata = nd3), error = function(e) e)
cat("  g absent from newdata but `g` is a column name in d only:",
    if (inherits(r3, "condition"))
      paste0("ERR: ", substr(conditionMessage(r3), 1, 80)) else
      paste("OK", paste(round(r3, 6), collapse = ",")), "\n")
g <- factor(c("g1", "g2"))
r4 <- tryCatch(predict(fit, newdata = nd3), error = function(e) e)
cat("  with a global `g` of the right length:",
    if (inherits(r4, "condition"))
      paste0("ERR: ", substr(conditionMessage(r4), 1, 80)) else
      paste("OK", paste(round(r4, 6), collapse = ",")), "\n")
