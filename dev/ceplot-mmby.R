# Lane ceplot: an mm(g1, g2, by = cbind(f1, f2)) term at new levels,
# under band = "boot" and on draws, against the Wald band of the same
# display. The Wald band carries each new member's variance from its
# own by-level's block (lp_extra_var()), so a boot or draws band that
# draws the new levels correctly has about the Wald width; one that
# draws nothing is much narrower, and one that draws two members
# sharing a block and a value independently is narrower too.
#   Rscript dev/ceplot-mmby.R [base|lane] > dev/ceplot-log/mmby-<arm>.txt
# Seeds: data 45 (test-ce-levels.R), bootstrap 5, draws 1, curves 1.
arm <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(arm)) arm <- "lane"
libs <- c("C:/Users/adf44/source/r/rellib-r4",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
if (arm == "lane") libs <- c("C:/Users/adf44/source/r/wt-ceplot-lib", libs)
.libPaths(libs)
suppressMessages({
  library(frmtmb)
  library(frmtmb.sample)
})
cat("ARM", arm, "frmtmb from", find.package("frmtmb"), "\n")
say <- function(...) cat(sprintf(...), "\n", sep = "")
f4 <- function(v) paste(sprintf("%.4f", v), collapse = " ")
width <- function(d) d$upper__ - d$lower__
hand <- function(fit, n, seed = 1) {
  tpl <- fit$frame$par_template
  est <- unlist(lapply(names(tpl), function(cp) fit$estimates[[cp]]))
  set.seed(seed)
  M <- matrix(rep(est, each = n) + rnorm(n * length(est), 0, 0.05), n,
              dimnames = list(NULL, frmtmb::brms_par_labels(fit)))
  structure(list(stanfit = NULL,
                 draws = cbind(frmtmb.sample:::draws_to_natural(M, fit),
                               lp__ = 0),
                 fit = fit), class = "frmtmb_draws")
}
set.seed(45)
d <- data.frame(x = rnorm(300), g1 = factor(sample(1:10, 300, TRUE)),
                g2 = factor(sample(1:10, 300, TRUE)))
fl <- rep(c("a", "b"), each = 5)
d$f1 <- factor(fl[d$g1])
d$f2 <- factor(fl[d$g2])
u <- rnorm(10, 0, 1)
d$y <- rnorm(300, 1 + 0.5 * d$x + 0.5 * (u[d$g1] + u[d$g2]), 0.5)
fit <- frm(bf(y ~ x + (1 | mm(g1, g2, by = cbind(f1, f2)))),
           family = gaussian(), data = d)
ds <- hand(fit, 200)
say("block sd at the estimate (by-levels a, b): %s", f4(exp(fit$estimates$theta)))
ce <- function(o, cond, ...) {
  tryCatch(suppressMessages(conditional_effects(
    o, "x", resolution = 3, re_formula = NULL, conditions = cond,
    ...))$x, error = function(e) e)
}
run <- function(label, cond) {
  w <- ce(fit, cond)
  b <- ce(fit, cond, band = "boot", boot = 60, seed = 5)
  r <- ce(ds, cond, seed = 1)
  say("%s | wald: est %s | width %s", label, f4(w$estimate__), f4(width(w)))
  for (nm in c("boot60", "draws")) {
    z <- if (nm == "boot60") b else r
    if (inherits(z, "error")) {
      say("%s | %s: ERROR %s", label, nm,
          substr(conditionMessage(z), 1, 150))
    } else {
      say("%s | %s: est %s | width %s | width / wald %s", label, nm,
          f4(z$estimate__), f4(width(z)), f4(width(z) / width(w)))
    }
  }
}
run("both unset, f1 = a, f2 = b", list(f1 = "a", f2 = "b"))
run("both unset, f1 = f2 = a", list(f1 = "a", f2 = "a"))
run("g1 = 99, g2 = 98, f1 = a, f2 = b",
    list(g1 = "99", g2 = "98", f1 = "a", f2 = "b"))
run("g1 = 2 (seen, by a), g2 unset, f2 = b",
    list(g1 = "2", f1 = "a", f2 = "b"))
run("both seen, g1 = 2, g2 = 7", list(g1 = "2", g2 = "7", f1 = "a",
                                       f2 = "b"))
say("done")
