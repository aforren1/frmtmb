# Reviewer: frmtmb's mo() interactions and group-level mo() terms on the
# lane build: names against brms's, equivalence with the indicator
# spelling, predictions on new factor levels and new groups, mo_search()
# on several simplexes and its cost.
# Usage: Rscript dev/moparity-rev-mo-probe.R [lib]
args <- commandArgs(TRUE)
lib <- if (length(args)) args[[1]] else "C:/Users/adf44/source/r/wt-moparity-lib"
.libPaths(c(lib, "C:/Users/adf44/source/r/rellib-r7",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), "\n")
set.seed(1)
n <- 90
d <- data.frame(x = sample(0:3, n, TRUE), f = factor(rep(c("a", "b", "c"), 30)),
                f2 = factor(rep(c("a", "b"), 45)), z = rnorm(n), w = rnorm(n),
                g = factor(rep(1:9, 10)), xl = rnorm(n) > 0)
d$y <- rnorm(n); d$x2 <- sample(0:2, n, TRUE)
d$fa <- as.numeric(d$f == "a"); d$fb <- as.numeric(d$f == "b")
d$fc <- as.numeric(d$f == "c")
# a signal so the simplexes are identified
eff <- c(0, 0.5, 0.7, 1.4)
d$y <- 1 + c(a = 0.4, b = -0.6, c = 1.1)[as.character(d$f)] * eff[d$x + 1] +
  0.5 * d$z + rnorm(n, 0, 0.5) + rnorm(9, 0, 0.4)[d$g]
nm <- function(fo) {
  r <- tryCatch({
    fit <- suppressWarnings(frm(bf(fo), family = gaussian(), data = d))
    gp <- get_prior(bf(fo), family = gaussian(), data = d)
    list(fit = fit,
         b = colnames(fit$frame$linpreds[[1]]$X),
         simo = gp$coef[gp$class == "simo"])
  }, error = function(e) e)
  if (inherits(r, "error")) {
    cat("==", deparse1(fo), ": ERROR", conditionMessage(r), "\n")
    return(invisible(NULL))
  }
  cat("==", deparse1(fo), "\n   X cols:", paste(r$b, collapse = " "),
      "\n   simo:", paste(r$simo, collapse = " "),
      "\n   logLik:", format(as.numeric(logLik(r$fit)), digits = 12), "\n")
  invisible(r$fit)
}
f1 <- nm(y ~ mo(x) * f)
f1i <- nm(y ~ fb + fc + mo(x) + mo(x):fb + mo(x):fc)
f2 <- nm(y ~ mo(x):f)
f2i <- nm(y ~ mo(x):fa + mo(x):fb + mo(x):fc)
f3 <- nm(y ~ f + mo(x):f)
f3i <- nm(y ~ fb + fc + mo(x):fa + mo(x):fb + mo(x):fc)
nm(y ~ mo(x) * f2)
nm(y ~ f2 + mo(x):f2)
nm(y ~ z + mo(x):z)
nm(y ~ z * mo(x))
nm(y ~ mo(x) * z)
nm(y ~ mo(x) * z * w)
nm(y ~ mo(x) * xl)
nm(y ~ (mo(x) | g))
nm(y ~ mo(x2) + (mo(x) | g))
nm(y ~ mo(x):mo(x2))
cat("\nlogLik gaps (factor vs indicator spelling):",
    format(c(as.numeric(logLik(f1)) - as.numeric(logLik(f1i)),
             as.numeric(logLik(f2)) - as.numeric(logLik(f2i)),
             as.numeric(logLik(f3)) - as.numeric(logLik(f3i))), digits = 4),
    "\n")

# predictions on new data
cat("\n== predictions, mo(x) * f\n")
LP <- function(...) { r <- frm_linpred(...); if (is.list(r)) r else as.numeric(r) }
p_in <- LP(f1)
nd <- d[d$f == "c", c("x", "f", "z")]
p_c <- LP(f1, newdata = nd)
cat("one level in newdata (factor), max |diff| vs fitted rows:",
    max(abs(as.numeric(p_c) - as.numeric(p_in[d$f == "c"]))), "\n")
nd2 <- nd
nd2$f <- as.character(nd2$f)
p_c2 <- tryCatch(LP(f1, newdata = nd2), error = function(e) e)
cat("one level as character:",
    if (inherits(p_c2, "error")) conditionMessage(p_c2) else
      max(abs(as.numeric(p_c2) - as.numeric(p_in[d$f == "c"]))), "\n")
nd3 <- nd
nd3$f <- factor("c", levels = c("c", "b", "a"))
p_c3 <- tryCatch(LP(f1, newdata = nd3), error = function(e) e)
cat("reordered factor levels:",
    if (inherits(p_c3, "error")) conditionMessage(p_c3) else
      max(abs(as.numeric(p_c3) - as.numeric(p_in[d$f == "c"]))), "\n")
nd4 <- nd[1:3, ]
nd4$f <- factor(c("d", "c", "c"))
p_c4 <- tryCatch(LP(f1, newdata = nd4), error = function(e) e)
cat("new level 'd':",
    if (inherits(p_c4, "error")) conditionMessage(p_c4) else
      paste(format(p_c4), collapse = " "), "\n")
p_se <- tryCatch(LP(f1, newdata = nd[1:4, ], se.fit = TRUE),
                 error = function(e) e)
cat("se.fit on newdata:",
    if (inherits(p_se, "error")) conditionMessage(p_se) else
      paste(format(p_se$se.fit, digits = 4), collapse = " "), "\n")
p_se0 <- LP(f1, se.fit = TRUE)
cat("se.fit in sample rows of level c:",
    paste(format(p_se0$se.fit[d$f == "c"][1:4], digits = 4), collapse = " "),
    "\n")

cat("\n== group-level mo(): y ~ mo(x) * f2 + (mo(x) | g)\n")
fr <- suppressWarnings(frm(bf(y ~ mo(x) * f2 + (mo(x) | g)),
                           family = gaussian(), data = d))
print(summary(fr))
p_in <- LP(fr)
p_nd <- LP(fr, newdata = d)
cat("predict(newdata = data) vs fitted, max |diff|:",
    max(abs(as.numeric(p_nd) - as.numeric(p_in))), "\n")
p_se_in <- LP(fr, se.fit = TRUE)
p_se_nd <- LP(fr, newdata = d, se.fit = TRUE)
cat("se.fit in sample vs newdata = data, max |diff|:",
    max(abs(p_se_in$se.fit - p_se_nd$se.fit)), "\n")
ndn <- d[1:4, ]
ndn$g <- factor(c("new1", "new1", "1", "2"))
pn <- tryCatch(LP(fr, newdata = ndn, allow_new_levels = TRUE),
               error = function(e) e)
cat("new groups, allow_new_levels:",
    if (inherits(pn, "error")) conditionMessage(pn) else
      paste(format(as.numeric(pn), digits = 6), collapse = " "), "\n")
pn0 <- LP(fr, newdata = ndn, re_formula = NA)
cat("same rows, re_formula = NA:",
    paste(format(as.numeric(pn0), digits = 6), collapse = " "), "\n")
pns <- tryCatch(LP(fr, newdata = ndn, allow_new_levels = TRUE,
                        se.fit = TRUE), error = function(e) e)
cat("new groups se.fit:",
    if (inherits(pns, "error")) conditionMessage(pns) else
      paste(format(pns$se.fit, digits = 4), collapse = " "), "\n")
pnf <- tryCatch(LP(fr, newdata = ndn), error = function(e) e)
cat("new groups without allow_new_levels:",
    if (inherits(pnf, "error")) substr(conditionMessage(pnf), 1, 120) else
      "no error", "\n")
# hand check: eta for row 1 of g = 1
est <- fr$estimates
re <- ranef(fr)$g
cat("ranef colnames:", paste(colnames(re), collapse = " "), "\n")

cat("\n== mo_search on several simplexes, cost\n")
for (ms in c(TRUE, FALSE)) {
  t0 <- proc.time()[[3]]
  ff <- suppressWarnings(frm(bf(y ~ mo(x) * f + (1 | g)), family = gaussian(),
                             data = d,
                             control = frmtmb_control(mo_search = ms)))
  cat("mo_search", ms, ": logLik", format(as.numeric(logLik(ff)), digits = 12),
      "| evals", ff$opt$evaluations[[1]] %||% NA,
      "| mo_search slot", paste(format(ff$opt$mo_search %||% NA), collapse = " "),
      "| secs", round(proc.time()[[3]] - t0, 2), "\n")
}
# a case with a face per term and opposite slopes
set.seed(7)
d$y2 <- 1 + c(a = 1, b = -1, c = 0.02)[as.character(d$f)] *
  c(0, 0, 1, 1)[d$x + 1] + rnorm(n, 0, 0.4)
for (ms in c(TRUE, FALSE)) {
  ff <- suppressWarnings(frm(bf(y2 ~ mo(x) * f), family = gaussian(),
                             data = d,
                             control = frmtmb_control(mo_search = ms,
                                                      verbose = FALSE)))
  W <- lapply(frmtmb::mo_frame_terms(ff), function(t) {
    round(frmtmb::mo_simplex(ff$estimates[[t$zeta]]), 4)
  })
  cat("faces y2, mo_search", ms, ": logLik",
      format(as.numeric(logLik(ff)), digits = 12),
      "| mo_search", paste(format(ff$opt$mo_search %||% NA), collapse = " "),
      "| weights", paste(vapply(W, paste, "", collapse = ","), collapse = " / "),
      "\n")
}
