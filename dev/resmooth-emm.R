# Lane wt-resmooth. emmeans builds its reference grid from the PLOTTABLE
# variables, which exclude a factor smooth's grouping factor, and from
# the group-level terms' grouping variables, which a smooth is not among.
# After this change the smooth is rebuilt at every re_formula and needs
# that column, so the grid could be missing it. Constructed here rather
# than assumed either way.
#   RESMOOTH_LIB=base Rscript dev/resmooth-emm.R > dev/resmooth-emm-before.txt
#   Rscript dev/resmooth-emm.R > dev/resmooth-emm-after.txt
base <- identical(Sys.getenv("RESMOOTH_LIB"), "base")
lane <- c("C:/Users/adf44/source/r/wt-resmooth-lib2",
          "C:/Users/adf44/source/r/wt-resmooth-lib")
.libPaths(c(if (!base) lane,
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true")
suppressMessages(library(frmtmb))
suppressMessages(library(emmeans))
cat("frmtmb from:", find.package("frmtmb"), "\n")
set.seed(5)
n <- 300
# g CROSSED with f: g cycling 1..10 against an alternating f gives both
# period 2, so g's parity fixes f and emmeans reads the design as nested
d <- data.frame(x = runif(n), f = factor(rep(c("a", "b"), length.out = n)),
                g = factor(rep(1:10, each = 30L)))
d$y <- sin(2 * pi * d$x) + 0.7 * (d$f == "b") +
  rnorm(10, 0, 0.5)[d$g] * d$x + rnorm(n, 0, 0.3)
show <- function(lab, expr) {
  v <- tryCatch({
    r <- force(expr)
    paste("OK:", paste(sprintf("%.4f", as.data.frame(r)[[2L]]),
                       collapse = " "))
  }, error = function(e) paste("ERROR:", substr(conditionMessage(e), 1, 160)))
  cat(sprintf("%-34s %s\n", lab, v))
}
fs <- suppressWarnings(frm(bf(y ~ f + s(x, g, bs = "fs", k = 5)), data = d))
show("emmeans(fs fit, ~f)", emmeans(fs, ~ f))
show("emmeans(fs fit, ~f, re_formula=NULL)",
     emmeans(fs, ~ f, re_formula = NULL))
re <- suppressWarnings(frm(bf(y ~ f + s(x) + s(g, bs = "re")), data = d))
show("emmeans(re fit, ~f)", emmeans(re, ~ f))
gg <- suppressWarnings(frm(bf(y ~ f + s(x) + (1 | g)), data = d))
show("emmeans((1|g) fit, ~f)", emmeans(gg, ~ f))
t2 <- suppressWarnings(frm(bf(y ~ f + t2(x, g, bs = c("cr", "re"))),
                           data = d))
show("emmeans(t2 fit, ~f)", emmeans(t2, ~ f))
