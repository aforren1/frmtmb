# Reviewer: fit method, conditions rows mixing an observed level with
# NA, allow_new_levels = TRUE (data seed 21 as in postfit2-rev-adv-fit.R).
#   Rscript dev/postfit2-rev-mixed-fit.R <base|lane>
args <- commandArgs(TRUE)
arm <- if (length(args)) args[1] else "lane"
libs <- c("C:/Users/adf44/source/r/wt-postfit2-lib",
          "C:/Users/adf44/source/r/rellib-r3",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
if (arm == "base") libs <- libs[-1]
.libPaths(libs)
suppressMessages(library(frmtmb))
say <- function(...) cat(sprintf(...), "\n", sep = "")
f3 <- function(v) paste(sprintf("%.4f", v), collapse = " ")
say("ARM %s", arm)
set.seed(21)
d <- data.frame(x = rnorm(80), g = factor(rep(1:8, 10)))
d$y <- rnorm(80, 1 + 0.5 * d$x + rnorm(8, 0, 2)[d$g], 1)
fit <- frm(bf(y ~ x + (1 | g)), family = gaussian(), data = d)
cm <- data.frame(g = factor(c("2", NA), levels = 1:8))
rownames(cm) <- c("lev2", "unset")
for (band in c("wald", "boot")) {
  r <- suppressMessages(conditional_effects(
    fit, "x", resolution = 3, re_formula = NULL, conditions = cm,
    allow_new_levels = TRUE, band = band, boot = 40, seed = 5))$x
  for (cv in c("lev2", "unset")) {
    s <- r[r$cond__ == cv, ]
    say("%s cond %s: est %s | width %s", band, cv, f3(s$estimate__),
        f3(s$upper__ - s$lower__))
  }
}
