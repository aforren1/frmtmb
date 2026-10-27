## Reviewer probe: conditional_effects() grid and reference values of a
## mo() variable that only a distributional parameter reads, and the
## fractional new-data cases, on the build in the first argument's
## library (lane or rv), seed 5.
##   Rscript dev/review-mocode-ce2.R lane|rv
arm <- commandArgs(trailingOnly = TRUE)[1L]
lp <- c("C:/Users/adf44/source/r/wt-mocode-lib",
        "C:/Users/adf44/source/r/rellib-r13",
        "C:/Users/adf44/AppData/Local/R/win-library/4.6")
if (identical(arm, "rv")) lp <- c("C:/Users/adf44/source/r/wt-mocode-rv-lib",
                                  lp)
.libPaths(lp)
suppressMessages(library(frmtmb))
cat("frmtmb", arm, "lib:", dirname(system.file(package = "frmtmb")), "\n")
q <- function(e) suppressWarnings(suppressMessages(e))
show <- function(lab, expr) {
  r <- tryCatch(q(expr), error = function(e) e)
  cat(sprintf("%-34s %s\n", lab, if (inherits(r, "error"))
    paste("ERROR:", substr(conditionMessage(r), 1, 140)) else r))
}
set.seed(5)
n <- 300
d <- data.frame(k = sample(1:4, n, replace = TRUE))
d$z <- rnorm(n)
d$y <- 2 + 4.5 * c(0, cumsum(c(0.6, 0.3, 0.1)))[d$k] + rnorm(n)
fs <- q(frm(bf(y ~ z, sigma ~ mo(k)), data = d))
show("sigma ~ mo(k): ce(k) grid length", {
  ce <- conditional_effects(fs, effects = "k", dpar = "sigma")[[1L]]
  paste(length(ce$k), "points:", paste(format(head(unique(ce$k), 5)),
                                       collapse = " "))
})
fz <- q(frm(y ~ mo(k) + z, data = d))
show("ce(z) held k", {
  ce <- conditional_effects(fz, effects = "z")[[1L]]
  paste(format(unique(ce$k)), "estimate at first z", format(ce$estimate__[1]))
})
show("ce(z:k) k values", {
  ce <- conditional_effects(fz, effects = "z:k")[[1L]]
  paste(format(unique(ce$k)), collapse = " ")
})
show("ce(k) estimates", {
  ce <- conditional_effects(fz, effects = "k")[[1L]]
  paste(format(ce$k), "->", format(signif(ce$estimate__, 6)), collapse = "; ")
})
for (x in c(0.6, 2.5, 4.4)) {
  show(paste("fitted(k =", x, ")"),
       fitted(fz, newdata = data.frame(k = x, z = 0))[, 1])
}
show("emmeans(~ z)", {
  s <- summary(emmeans::emmeans(fz, ~ z))
  paste(format(s$emmean), "at k =", format(s$k %||% NA))
})
