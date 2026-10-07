# Lane setier: test-quadrature-defects.R's nested quadrature fits: which
# loses a standard error, and how the objective moves when that sd is
# pushed toward zero (the edge probe of se_at_edge()).
.libPaths(c("C:/Users/adf44/source/r/wt-setier-lib",
            "C:/Users/adf44/source/r/rellib-r6",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
set.seed(4)
ng <- 20
nt <- 5
n <- ng * nt
d <- data.frame(g = factor(rep(seq_len(ng), each = nt)),
                gb = factor(rep(c("i", "ii"), length.out = n)),
                x = rnorm(n))
d$ga <- d$g
d$eta <- 0.5 + 0.4 * d$x + rnorm(ng, 0, 0.4)[d$g] +
  rnorm(ng * 2, 0, 0.3)[as.integer(d$ga) * 2 + as.integer(d$gb) - 2]
cases <- list(
  poisson = list(fam = poisson(), y = function() rpois(n, exp(d$eta))),
  gamma = list(fam = Gamma(link = "log"), y = function() {
    rgamma(n, shape = 2, scale = exp(d$eta) / 2)
  }),
  beta = list(fam = Beta(), y = function() {
    p <- plogis(d$eta)
    rbeta(n, p * 5, (1 - p) * 5)
  }))
for (nm in names(cases)) for (re in c("(1 | g)", "(1 | ga/gb)")) {
  set.seed(11)
  dd <- d
  dd$y <- cases[[nm]]$y()
  fo <- stats::as.formula(paste("y ~ 1 + x +", re))
  fit <- suppressMessages(suppressWarnings(
    frm(bf(fo) + cases[[nm]]$fam, data = dd, quadrature = TRUE)))
  lost <- ns$sdr_of(fit)$se_lost
  cat(nm, re, "par", format(fit$opt$par, digits = 4), "lost",
      paste(names(lost), lost), "\n")
  if (!length(lost)) next
  p <- fit$opt$par
  j <- match(names(lost)[1], ns$outer_par_names(fit))
  f0 <- fit$obj$fn(p)
  for (s in c(-20, -5, -4, -3, -2, -1)) {
    q <- p
    q[j] <- q[j] + s
    cat("  step", s, "dll", fit$obj$fn(q) - f0, "\n")
  }
}
