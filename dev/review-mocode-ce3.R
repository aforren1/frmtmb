## Fix round 1 check, frmtmb side of dev/review-mocode-brmsfit2.R (same
## data, seed 5): conditional_effects() reference and moderator values
## with two mo() terms; (mo(k) | g) SEs.
##   Rscript dev/review-mocode-ce3.R lane|base
arm <- commandArgs(trailingOnly = TRUE)[1L]
base <- c("C:/Users/adf44/source/r/rellib-r13",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
.libPaths(if (identical(arm, "lane")) {
  c("C:/Users/adf44/source/r/wt-mocode-lib", base)
} else base)
suppressMessages(library(frmtmb))
cat("frmtmb", arm, "lib:", dirname(system.file(package = "frmtmb")), "\n")
q <- function(e) suppressWarnings(suppressMessages(e))
f7 <- function(x) paste(format(signif(as.numeric(x), 7)), collapse = " ")
set.seed(5)
n <- 300
d <- data.frame(k = sample(1:4, n, replace = TRUE))
d$z <- rnorm(n)
d$g <- factor(sample(1:8, n, replace = TRUE))
d$j <- sample(2:6, n, replace = TRUE, prob = c(1, 1, 3, 1, 1))
d$y <- 2 + 4.5 * c(0, cumsum(c(0.6, 0.3, 0.1)))[d$k] +
  c(0, 0.5, 0.8, 1, 1.1)[d$j - 1L] + rnorm(8, sd = 0.4)[d$g] +
  d$k * rnorm(8, sd = 0.2)[d$g] + rnorm(n)
f2 <- q(frm(y ~ mo(k) + mo(j) + z, data = d))
for (eff in c("z", "k", "j", "z:k", "z:j", "k:j")) {
  ce <- tryCatch(q(conditional_effects(f2, effects = eff))[[1L]],
                 error = function(e) e)
  if (inherits(ce, "error")) {
    cat(sprintf("ce(%s): ERROR %s\n", eff, conditionMessage(ce)))
    next
  }
  cat(sprintf("ce(%s): k = {%s} j = {%s} z = {%s} effect2__ = {%s}\n", eff,
              f7(unique(ce$k)), f7(unique(ce$j)),
              f7(utils::head(unique(ce$z), 2)),
              paste(unique(as.character(ce$effect2__)), collapse = " ")))
  cat(sprintf("        estimate__ (first 4): %s\n", f7(head(ce$estimate__, 4))))
}
fre <- q(frm(y ~ mo(k) + (mo(k) | g), data = d))
cat("(mo(k) | g) g = 1, k = 1..4 se:",
    f7(q(fitted(fre, newdata = data.frame(k = 1:4, g = "1")))[, 2]), "\n")
cat("(mo(k) | g) re_formula = NA se:",
    f7(q(fitted(fre, newdata = data.frame(k = 1:4), re_formula = NA))[, 2]),
    "\n")
cat("(mo(k) | g) new level se:",
    f7(q(fitted(fre, newdata = data.frame(k = 1:4, g = "new"),
                allow_new_levels = TRUE))[, 2]), "\n")
cat("(mo(k) | g) g = 1 est:",
    f7(q(fitted(fre, newdata = data.frame(k = 1:4, g = "1")))[, 1]), "\n")
