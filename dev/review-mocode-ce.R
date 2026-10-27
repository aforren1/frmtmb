## Reviewer probe: conditional_effects() of a mo() effect when the rows
## that hold the data's minimum are not fitted (missing response;
## subset()), seed 5.
##   Rscript dev/review-mocode-ce.R lane|base
arm <- commandArgs(trailingOnly = TRUE)[1L]
base <- c("C:/Users/adf44/source/r/rellib-r13",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
.libPaths(if (identical(arm, "lane")) {
  c("C:/Users/adf44/source/r/wt-mocode-lib", base)
} else base)
suppressMessages(library(frmtmb))
cat("frmtmb", format(packageVersion("frmtmb")), arm, "\n")
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
d$y2 <- d$y + rnorm(n)
dm <- d
dm$y[dm$k == 1L] <- NA
fm <- q(frm(y ~ mo(k), data = dm))
show("missing y: ce(k) grid", {
  ce <- conditional_effects(fm, effects = "k")[[1L]]
  paste(format(ce$k), collapse = " ")
})
fz <- q(frm(y ~ z + mo(k), data = dm))
show("missing y: ce(z) held k", {
  ce <- conditional_effects(fz, effects = "z")[[1L]]
  paste(format(unique(ce$k)), collapse = " ")
})
show("missing y: ce(z:k) values of k", {
  ce <- conditional_effects(fz, effects = "z:k")[[1L]]
  paste(format(unique(ce$k)), collapse = " ")
})
ds <- d
ds$s1 <- ds$k > 1L
ds$s2 <- TRUE
fs <- q(frm(bf(y | subset(s1) ~ mo(k)) + bf(y2 | subset(s2) ~ mo(k)) +
              set_rescor(FALSE), data = ds))
show("subset: ce(k) resp y", {
  ce <- conditional_effects(fs, effects = "k", resp = "y")[[1L]]
  paste(format(ce$k), collapse = " ")
})
show("subset: ce(k) resp y2", {
  ce <- conditional_effects(fs, effects = "k", resp = "y2")[[1L]]
  paste(format(ce$k), collapse = " ")
})
show("subset: ce(k) all", {
  ce <- conditional_effects(fs, effects = "k")
  paste(names(ce), vapply(ce, function(x) paste(format(x$k),
                                               collapse = " "), ""),
        collapse = " | ")
})
show("missing y: predict(fm) in-sample n", nrow(fitted(fm)))
show("missing y: fitted(newdata = dm)", {
  paste(dim(fitted(fm, newdata = dm[!is.na(dm$y), ])), collapse = "x")
})
show("missing y: fitted(newdata = dm all rows)", {
  paste(dim(fitted(fm, newdata = dm)), collapse = "x")
})
