# Lane predraw2: brms 2.23.0's posterior_predict() on the three weakly
# identified mixtures of backlog 4882 (the sweep's seed 101 data), against
# the known truth, beside frmtmb's predict() on the build given (4000
# draws each; brms 4 chains of 2000, seed 11, default priors).
# Usage: Rscript dev/predraw2-brms.R <lib or "base"> [family to run]
args <- commandArgs(TRUE)
lib <- if (args[1] != "base") args[1] else NULL
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
.libPaths(c(lib, "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({
  library(brms)
  library(frmtmb)
})
cat("brms", as.character(packageVersion("brms")), " frmtmb",
    find.package("frmtmb"), "\n")
n <- 300L
source("dev/ppdraw-sweep-gen.R")
src <- readLines("dev/predraw2-mixcov.R")
case <- NULL
summ <- function(label, p, d, fam, case) {
  case <<- case
  fam <<- fam
  tr <- truth_row(d)
  if (tr$type == "cat") {
    sh <- t(apply(p, 2L, function(col) tabulate(col, 4L) / nrow(p)))
    cat(sprintf("%-28s mae of shares %.4f\n", label, mean(abs(sh - tr$P))))
    return(invisible())
  }
  q <- apply(p, 2L, stats::quantile, c(0.01, 0.05, 0.95, 0.99),
             names = FALSE)
  c90 <- mean(vapply(seq_len(n), function(i) mass_in(tr, i, q[2, i],
                                                     q[3, i]), 0))
  c98 <- mean(vapply(seq_len(n), function(i) mass_in(tr, i, q[1, i],
                                                     q[4, i]), 0))
  tm <- if (tr$type == "mixn") tr$mu %*% tr$w else tr$lam %*% tr$w
  cat(sprintf("%-28s mean %.4f cov90 %.3f cov98 %.3f w98 %.3f mae %.4f\n",
              label, mean(p), c90, c98, mean(q[4, ] - q[1, ]),
              mean(abs(colMeans(p) - tm))))
}
fam <- NULL
eval(parse(text = src[grep("^truth_row <-", src):
                        (grep("^rows <- list", src) - 1L)]))
for (cs in list(c("mixture_gaussian", "boundary"),
                c("mixture_poisson", "boundary"),
                c("mixture_cumulative", "interior"))) {
  if (length(args) > 1 && cs[1] != args[2]) next
  d <- gen(cs[1], cs[2])
  cat("\n==", cs[1], cs[2], "\n")
  bfam <- switch(cs[1],
    mixture_gaussian = brms::mixture(gaussian, gaussian),
    mixture_poisson = brms::mixture(poisson, poisson),
    mixture_cumulative = brms::mixture(brms::cumulative(),
                                       brms::cumulative()))
  b <- brm(brms::bf(y ~ x), family = bfam, data = d, chains = 4,
           iter = 2000, seed = 11, refresh = 0)
  print(summary(b)$fixed[, c("Estimate", "Rhat")])
  set.seed(1)
  summ("brms posterior_predict", posterior_predict(b), d, cs[1], cs[2])
  s <- spec(cs[1])
  f <- suppressWarnings(frm(s[[1]], family = s[[2]], data = d))
  set.seed(1)
  pf <- suppressMessages(suppressWarnings(predict(f, ndraws = 4000,
                                                  summary = FALSE)))
  summ("frmtmb predict", pf, d, cs[1], cs[2])
  set.seed(1)
  p0 <- suppressWarnings(predict(f, ndraws = 4000, summary = FALSE,
                                 propagate_error = FALSE))
  summ("frmtmb plug-in", p0, d, cs[1], cs[2])
}
