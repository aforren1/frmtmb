# A Laplace fit collapsed at sd = 0: the importance fit stays there (the lane's
# new test pins it: "settles where Laplace put it"). Is sd = 0 the exact maximum,
# or a local stop? Exact profile of the sd over a grid, from base R densities.
# The prior-on-Laplace sd_up_escape also exists, and importance_fit skips it.
source("C:/Users/adf44/source/r/frmtmb-wt-impsamp/dev/review-impsamp-common.R")
set_libs("after")
suppressPackageStartupMessages(library(frmtmb))
source("C:/Users/adf44/source/r/frmtmb-wt-impsamp/dev/review-impsamp-designs.R")
gh <- gh_rule(81L)
prof <- function(de, s_grid, np) {
  dl <- split(de$data, de$data$g)
  nll <- function(par) {
    sd <- de$sdof(par)
    tot <- 0
    for (d in dl) {
      fn <- function(u) de$ll(par, d, u) + dnorm(u, 0, sd, log = TRUE)
      o <- optimize(function(u) -fn(u), c(-8 * sd - 15, 8 * sd + 15), tol = 1e-7)
      m <- o$minimum
      h <- 1e-3
      cur <- -(fn(m + h) - 2 * fn(m) + fn(m - h)) / h^2
      s <- if (is.finite(cur) && cur > 1e-8) 1 / sqrt(cur) else sd
      nodes <- m + sqrt(2) * s * gh$x
      a <- fn(nodes) + gh$x^2 + log(gh$w)
      tot <- tot - (max(a) + log(sum(exp(a - max(a)))) + log(sqrt(2) * s))
    }
    tot
  }
  list(nll = nll, tab = NULL)
}
scan_case <- function(lab, de, draws) {
  fi <- suppressWarnings(frm(de$formula, family = de$family, data = de$data,
                             importance = draws))
  fl <- suppressWarnings(frm(de$formula, family = de$family, data = de$data))
  p0 <- unname(fi$opt$par)
  k <- length(p0)
  nll <- prof(de)$nll
  cat("\n==", lab, "fit par", format(p0, digits = 4), "lap par",
      format(unname(fl$opt$par), digits = 4), "ll imp",
      format(as.numeric(logLik(fi)), digits = 8), "\n")
  best <- Inf
  for (lt in c(-6, -4, -3, -2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2)) {
    o <- nlminb(p0[-k], function(b) nll(c(b, lt)))
    cat(sprintf("  log sd %5.1f  profile nll %.5f  (fit's value at its own %.5f)\n",
                lt, o$objective, nll(p0)))
    best <- min(best, o$objective)
  }
  cat("  exact nll at the fit minus best on the grid:", format(nll(p0) - best, digits = 4), "\n")
  # does any start above zero lead away from sd = 0 under the importance fit?
  st <- list()
  fi2 <- tryCatch(suppressWarnings(frm(de$formula, family = de$family, data = de$data,
                                       importance = draws,
                                       start = list(theta = 0))),
                  error = function(e) e)
  if (!inherits(fi2, "error")) {
    cat("  fit started at log sd 0:", format(unname(fi2$opt$par), digits = 4), "ll",
        format(as.numeric(logLik(fi2)), digits = 8), "\n")
  }
}
for (s in 1:12) {
  de <- design_list$unbal_bern(s)
  fl <- suppressWarnings(frm(de$formula, family = de$family, data = de$data))
  if (fl$opt$par[[3]] < -4) scan_case(paste("unbal_bern", s), de, de$draws)
}
for (s in 1:12) {
  de <- design_list$t_outlier(s)
  fl <- suppressWarnings(frm(de$formula, family = de$family, data = de$data))
  if (fl$opt$par[[4]] < -4) scan_case(paste("t_outlier", s), de, de$draws)
}
