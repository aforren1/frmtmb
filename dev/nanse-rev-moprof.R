# Reviewer: item 9, independently. For a few seeds of brms_monotonic's
# data code, the profile log-likelihood of ls ~ mo(income) * age over
# the two simplexes by an exhaustive grid (step 0.02 on each 2-simplex,
# both simplexes; given the simplexes the model is OLS), refined with
# Nelder-Mead on a softmax-free (projected) parameterization, against
# frmtmb's logLik on the base and lane builds.
#   Rscript dev/nanse-rev-moprof.R base|lane seeds
args <- commandArgs(trailingOnly = TRUE)
arm <- args[1]
seeds <- eval(parse(text = args[2]))
libs <- switch(arm,
  base = "C:/Users/adf44/source/r/rellib-r5",
  lane = c("C:/Users/adf44/source/r/wt-nanse-lib",
           "C:/Users/adf44/source/r/rellib-r5"))
.libPaths(c(libs, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("arm", arm, "from", find.package("frmtmb"), "\n")
mk <- function(s) {
  set.seed(s)
  lev <- c("below_20", "20_to_40", "40_to_100", "greater_100")
  income <- factor(sample(lev, 100, TRUE), levels = lev, ordered = TRUE)
  ls <- c(30, 60, 70, 75)[income] + rnorm(100, sd = 7)
  d <- data.frame(income, ls)
  d$age <- rnorm(100, mean = 40, sd = 10)
  d
}
grid_simplex <- function(h) {
  g <- seq(0, 1, by = h)
  m <- as.matrix(expand.grid(a = g, b = g))
  m <- m[rowSums(m) <= 1 + 1e-12, , drop = FALSE]
  cbind(m, pmax(0, 1 - rowSums(m)))
}
ll_at <- function(w1, w2, d) {
  code <- as.integer(d$income) - 1L
  m1 <- 3 * c(0, cumsum(w1))[code + 1L]
  m2 <- 3 * c(0, cumsum(w2))[code + 1L]
  X <- cbind(1, d$age, m1, m2 * d$age)
  r <- stats::lm.fit(X, d$ls)$residuals
  n <- nrow(d)
  -(n / 2) * (log(2 * pi * sum(r^2) / n) + 1)
}
G <- grid_simplex(0.05)
for (s in seeds) {
  d <- mk(s)
  # coarse: every pair on the 0.05 grid
  best <- -Inf
  for (i in seq_len(nrow(G))) {
    for (j in seq_len(nrow(G))) {
      v <- ll_at(G[i, ], G[j, ], d)
      if (v > best) {
        best <- v
        bi <- i
        bj <- j
      }
    }
  }
  # refine: projected coordinates, weights = |u| / sum |u|, which reach
  # the boundary exactly
  f <- function(u) {
    a <- abs(u[1:3]); b <- abs(u[4:6])
    -ll_at(a / sum(a), b / sum(b), d)
  }
  o <- stats::optim(c(G[bi, ], G[bj, ]) + 1e-9, f,
                    control = list(maxit = 5000, reltol = 1e-14))
  w <- abs(o$par)
  w1 <- w[1:3] / sum(w[1:3]); w2 <- w[4:6] / sum(w[4:6])
  fit <- suppressWarnings(suppressMessages(frm(ls ~ mo(income) * age,
                                               data = d)))
  llf <- as.numeric(logLik(fit))
  cat(sprintf(paste0("seed %3d frmtmb %.6f code %d | grid max %.6f, ",
                     "refined %.6f | gap %.4f | w1 %s w2 %s\n"),
              s, llf, fit$opt$convergence, best, -o$value, -o$value - llf,
              paste(signif(w1, 3), collapse = ","),
              paste(signif(w2, 3), collapse = ",")))
}
