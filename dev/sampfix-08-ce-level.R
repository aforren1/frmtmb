# Lane sampfix, script 08: conditional_effects() at an observed group
# level, on full draws and on the same draws with the r_ columns removed.
#   Rscript dev/sampfix-08-ce-level.R
LANE <- "C:/Users/adf44/source/r/wt-sampfix-lib"
.libPaths(c(LANE, "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
suppressMessages(library(frmtmb.sample))
set.seed(9)
dd <- data.frame(x = rnorm(60), g = factor(rep(1:6, 10)))
dd$y <- rnorm(60, 1 + 0.5 * dd$x + rnorm(6, 0, 0.5)[dd$g], 1)
fit <- frm(bf(y ~ x + (1 | g)), family = gaussian(), data = dd)
lab <- frmtmb::brms_par_labels(fit)
tpl <- fit$frame$par_template
est <- unlist(lapply(names(tpl), function(cp) fit$estimates[[cp]]))
set.seed(1)
M <- matrix(rep(est, each = 5) + rnorm(5 * length(est), 0, 0.05), 5,
            dimnames = list(NULL, lab))
M <- cbind(frmtmb.sample:::draws_to_natural(M, fit), lp__ = 0)
full <- structure(list(stanfit = NULL, draws = M, fit = fit),
                  class = "frmtmb_draws")
lap <- full
lap$draws <- M[, !startsWith(colnames(M), "r_")]
for (cond in list(data.frame(g = 1), data.frame(g = factor(1, levels = 1:6)),
                  data.frame(g = "1"))) {
  cat("conditions g of class ", class(cond$g), "\n", sep = "")
  a <- tryCatch(conditional_effects(full, effects = "x", resolution = 3,
                                    re_formula = NULL, conditions = cond),
                error = function(e) e)
  if (inherits(a, "error")) { cat("  full ERROR ", conditionMessage(a), "\n"); next }
  b <- conditional_effects(full, effects = "x", resolution = 3,
                           re_formula = NA, conditions = cond)
  cat("  full, NULL vs NA estimate: ",
      paste(format(a$x$estimate__ - b$x$estimate__, digits = 4),
            collapse = " "), "\n", sep = "")
  cat("  g column in the frame: ", paste(a$x$g, collapse = " "), "\n")
  r <- tryCatch(conditional_effects(lap, effects = "x", resolution = 3,
                                    re_formula = NULL, conditions = cond),
                error = function(e) e)
  cat("  laplace: ", if (inherits(r, "error")) conditionMessage(r) else
    paste(format(r$x$estimate__ - a$x$estimate__, digits = 4),
          collapse = " "), "\n", sep = "")
}
cat("\n-- is level 1's effect read? r_g[1,Intercept] shifted by +10\n")
f2 <- full
f2$draws[, "r_g[1,Intercept]"] <- f2$draws[, "r_g[1,Intercept]"] + 10
cond <- data.frame(g = factor(1, levels = 1:6))
ce <- function(o) conditional_effects(o, effects = "x", resolution = 3,
                                      re_formula = NULL, conditions = cond,
                                      seed = 4)$x$estimate__
cat("  full: ", paste(format(ce(full), digits = 5), collapse = " "), "\n")
cat("  full, shifted: ", paste(format(ce(f2), digits = 5), collapse = " "),
    "\n")
cat("  laplace: ", paste(format(ce(lap), digits = 5), collapse = " "), "\n")
cat("-- the fit's own curve at g = 1 and at re_formula = NA\n")
cf <- function(rf) conditional_effects(fit, effects = "x", resolution = 3,
                                       re_formula = rf, conditions = cond)$x$estimate__
cat("  fit NULL: ", paste(format(cf(NULL), digits = 5), collapse = " "), "\n")
cat("  fit NA:   ", paste(format(cf(NA), digits = 5), collapse = " "), "\n")
cat("  ranef(fit) g = 1: ", format(ranef(fit)$g[1, 1], digits = 5), "\n")
cat("\n-- the same at level 2, which is not the new-level placeholder\n")
cond <- data.frame(g = factor(2, levels = 1:6))
f3 <- full
f3$draws[, "r_g[2,Intercept]"] <- f3$draws[, "r_g[2,Intercept]"] + 10
cat("  full: ", paste(format(ce(full), digits = 5), collapse = " "), "\n")
cat("  full, r_g[2] shifted: ", paste(format(ce(f3), digits = 5),
                                     collapse = " "), "\n")
r <- tryCatch(ce(lap), error = function(e) e)
cat("  laplace: ", if (inherits(r, "error")) substr(conditionMessage(r), 1, 60)
    else paste(format(r, digits = 5), collapse = " "), "\n")
