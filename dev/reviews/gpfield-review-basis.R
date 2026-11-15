# Review of lane gpfield, frmtmb_gp_basis: how often it fires on ordinary
# approximate gp() fits, and whether each firing is a true
# non-identification. 30 gaussian fits: f drawn from an exact squared
# exponential GP with length scale ls0 on x in [0, 10], noise sd 0.3,
# k in 5, 8, 10, 15, 20, 30. For each fit: warnings, the length scale over
# the basis spacing, and (lane arm) the profile log-likelihood at fixed
# length scales of 0.05 to 4 spacings (others re-optimized).
# Usage: Rscript gpfield-review-basis.R <lane|base>
source("C:/Users/adf44/source/r/frmtmb-wt-gpfield/dev/gpfield-lib.R")
arm <- commandArgs(TRUE)[1]
gf_libs(arm)
ns <- asNamespace("frmtmb")
ks <- c(5, 8, 10, 15, 20, 30)
sets <- list(c(n = 60, ls0 = 0.7), c(n = 100, ls0 = 1.5), c(n = 100, ls0 = 3),
             c(n = 200, ls0 = 1), c(n = 50, ls0 = 2))
row <- 0L
for (si in seq_along(sets)) for (k in ks) {
  row <- row + 1L
  n <- sets[[si]][["n"]]; ls0 <- sets[[si]][["ls0"]]
  set.seed(9600 + row)
  x <- sort(runif(n, 0, 10))
  K <- exp(-outer(x, x, "-")^2 / (2 * ls0^2)) + diag(1e-8, n)
  f <- as.vector(t(chol(K)) %*% rnorm(n))
  d <- data.frame(x = x, y = f + 0.3 * rnorm(n))
  q <- gf_fit(bf(y ~ gp(x, k = k)), data = d, family = gaussian())
  fit <- q$fit
  if (inherits(fit, "error")) {
    cat("BASIS", row, "ERROR", conditionMessage(fit), "\n"); next
  }
  bk <- Filter(function(b) b$covstruct == "hsgp", fit$frame$re_blocks)[[1]]
  th <- fit$estimates$theta[bk$theta_idx]
  sp <- bk$gp_spacing %||% NA
  ll <- as.numeric(logLik(fit))
  ws <- unique(sub(":.*", "", q$warn))
  basis <- any(grepl("spacing its basis", q$warn))
  se_lost <- any(grepl("Standard errors are not available", q$warn))
  prof <- ""
  if (arm == "lane") {
    obj <- fit$obj
    par <- fit$opt$par
    pl <- obj$env$parList(par)
    # the outer position of theta's length scale
    tpos <- which(names(par) == "theta")[bk$theta_idx[2]]
    pv <- c()
    for (m in c(0.05, 0.25, 0.5, 1, 2, 4)) {
      p0 <- par
      p0[tpos] <- log(m * sp)
      free <- setdiff(seq_along(par), tpos)
      fn <- function(v) { p <- p0; p[free] <- v; obj$fn(p) }
      gr <- function(v) { p <- p0; p[free] <- v; obj$gr(p)[free] }
      best <- Inf
      # two starts: the fit's own others, and sd at the fit's sd^2 ls
      # product kept
      for (st in 1:2) {
        v0 <- par[free]
        if (st == 2) {
          spos <- which(free == which(names(par) == "theta")[bk$theta_idx[1]])
          v0[spos] <- th[1] + 0.5 * (th[2] - log(m * sp))
        }
        o <- tryCatch(suppressWarnings(stats::nlminb(v0, fn, gr)),
                      error = function(e) NULL)
        if (!is.null(o) && is.finite(o$objective)) best <- min(best, o$objective)
      }
      pv <- c(pv, -best)
    }
    prof <- paste(sprintf("%.3f", pv - ll), collapse = " ")
  }
  cat(sprintf(paste0("BASIS %2d %s n %3d ls0 %.1f k %2d | ll %.6f conv %d | ",
                     "ls/spacing %.3g sd %.3g | basis %d se_lost %d | ",
                     "prof(ls = .05 .25 .5 1 2 4 spacings) - ll: %s | ",
                     "gp2 %s | warn [%s]\n"),
              row, arm, n, ls0, k, ll, fit$opt$convergence,
              exp(th[2]) / sp, exp(th[1]), basis, se_lost, prof,
              paste(format(fit$opt$gp_second %||% NA, digits = 3),
                    collapse = ","), paste(ws, collapse = ";")))
}
