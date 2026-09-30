# emmeans() on the release build: which of poly() and offset() makes
# "undefined columns selected"? Run on rellib-r5 (0.67.0 candidate) and,
# with the argument "base", on rellib-r4 (0.66.0).
#
#   Rscript dev/rel067-emmoffset.R [base]
arm <- if (length(commandArgs(TRUE))) "base" else "release"
lib <- if (arm == "base") "C:/Users/adf44/source/r/rellib-r4" else
  "C:/Users/adf44/source/r/rellib-r5"
.libPaths(c(lib, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("ARM", arm, as.character(packageVersion("frmtmb")), "\n")
set.seed(20260930)
n <- 120
d <- data.frame(x = rnorm(n), z = rnorm(n), time = runif(n, 1, 5),
                f = factor(sample(c("a", "b", "c"), n, TRUE)))
d$yc <- rpois(n, d$time * exp(0.3 + 0.4 * d$x))
try1 <- function(lab, fo, spec = "f") {
  r <- tryCatch({
    fit <- frm(fo, data = d, family = poisson())
    e <- summary(emmeans::emmeans(fit, spec))
    paste(format(e$emmean, digits = 6), collapse = " ")
  }, error = function(e) paste("ERROR:", conditionMessage(e)))
  cat(sprintf("%-28s %s\n", lab, r))
}
try1("x + f", bf(yc ~ x + f))
try1("x + f + offset", bf(yc ~ x + f + offset(log(time))))
try1("poly(z, 2) + f", bf(yc ~ poly(z, 2) + f))
try1("poly(z, 2) + f + offset", bf(yc ~ poly(z, 2) + f + offset(log(time))))
try1("log(abs(z) + 1) + f", bf(yc ~ log(abs(z) + 1) + f))
try1("scale(z) + f", bf(yc ~ scale(z) + f))
