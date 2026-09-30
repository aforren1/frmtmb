# Reviewer, punch round 1: the mm() term under re_formula = NULL (P3 of
# postfit2-rev-p1-attack.R, data seed 43). Does "nothing set" read the
# reference level as an observed group, and is the boot band at an
# observed mm group centered on its estimate?
#   Rscript dev/postfit2-rev-p1-mm.R <base|lane>
args <- commandArgs(TRUE); arm <- if (length(args)) args[1] else "lane"
libs <- c("C:/Users/adf44/source/r/wt-postfit2-lib", "C:/Users/adf44/source/r/rellib-r3",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
if (arm == "base") libs <- libs[-1]
.libPaths(libs)
suppressMessages(library(frmtmb))
f3 <- function(v) paste(sprintf("%.4f", v), collapse = " ")
cat("ARM", arm, "\n")
set.seed(43)
d3 <- data.frame(x = rnorm(200), g1 = factor(sample(1:10, 200, TRUE)),
                 g2 = factor(sample(1:10, 200, TRUE)))
u <- rnorm(10, 0, 1)
d3$y <- rnorm(200, 1 + 0.5 * d3$x + 0.5 * (u[d3$g1] + u[d3$g2]), 0.5)
f <- frm(bf(y ~ x + (1 | mm(g1, g2))), family = gaussian(), data = d3)
ce <- function(...) suppressMessages(conditional_effects(f, "x", resolution = 3, re_formula = NULL, ...))
a <- ce()$x; b <- ce(conditions = list(g1 = "1", g2 = "1"))$x
cat("nothing set: g1 g2 in the frame:", as.character(a$g1[1]), as.character(a$g2[1]), "; est", f3(a$estimate__), "\n")
cat("g1 = g2 = 1 set: est", f3(b$estimate__), "; identical to nothing set:", identical(a$estimate__, b$estimate__), "\n")
w <- ce(conditions = list(g1 = "2", g2 = "2"))$x
bo <- ce(conditions = list(g1 = "2", g2 = "2"), band = "boot", boot = 60, seed = 5)
m <- colMeans(attr(bo, "boot")$t, na.rm = TRUE)
cat("g1 = g2 = 2: wald est", f3(w$estimate__), "; boot refit mean", f3(m), "; boot width / wald width",
    f3((bo$x$upper__ - bo$x$lower__) / (w$upper__ - w$lower__)), "\n")
cat("VarCorr sd:", f3(VarCorr(f)[[1]]$sd[, 1]), "\n")
