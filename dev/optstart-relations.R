# Lane optstart: each family's relation between its scale parameter and
# the least-squares residuals (start_ls_rules), measured. Data from the
# family's own simulator with a covariate that explains most of the
# response; the cold start (par_template()) of the scale against the
# true value, and of the slope, on 20 seeds. Run on the lane build and on
# rellib-r8 (the start before the lane).
#   Rscript dev/optstart-relations.R <lib or ""> [n]
a <- commandArgs(TRUE)
.libPaths(c(if (nzchar(a[1])) a[1], "C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
n <- if (length(a) > 1) as.integer(a[2]) else 2000L
cat("lib", dirname(find.package("frmtmb")), "n", n, "\n")
reg <- frmtmb:::family_registry
fams <- list(
  gaussian = list(fam = gaussian(), dp = list(sigma = 1.5), link = "id",
                  b = c(10, 20)),
  student = list(fam = student(), dp = list(sigma = 1.5, nu = 4),
                 link = "id", b = c(10, 20)),
  skew_normal = list(fam = skew_normal(), dp = list(sigma = 1.5, alpha = 4),
                     link = "id", b = c(10, 20)),
  exgaussian = list(fam = exgaussian(), dp = list(sigma = 1.5, beta = 1.5),
                    link = "id", b = c(10, 20)),
  asym_laplace = list(fam = asym_laplace(),
                      dp = list(sigma = 1, quantile = 0.5), link = "id",
                      b = c(10, 20)),
  zero_inflated_asym_laplace = list(
    fam = zero_inflated_asym_laplace(),
    dp = list(sigma = 1, quantile = 0.5, zi = 0.3), link = "id",
    b = c(10, 20)),
  huber = list(fam = huber(), dp = list(sigma = 1.5), link = "id",
               b = c(10, 20)),
  lognormal = list(fam = lognormal(), dp = list(sigma = 0.4), link = "id",
                   b = c(1, 0.8)),
  shifted_lognormal = list(fam = shifted_lognormal(),
                           dp = list(sigma = 0.4, ndt = 0.3), link = "id",
                           b = c(1, 0.8)),
  hurdle_lognormal = list(fam = hurdle_lognormal(),
                          dp = list(sigma = 0.4, hu = 0.3), link = "id",
                          b = c(1, 0.8)),
  Gamma = list(fam = Gamma(link = "log"), dp = list(shape = 3),
               link = "log", b = c(1, 0.8)),
  weibull = list(fam = weibull(), dp = list(shape = 2), link = "log",
                 b = c(1, 0.8)),
  beta = list(fam = Beta(), dp = list(phi = 20), link = "logit",
              b = c(0.3, 1)),
  negbinomial = list(fam = negbinomial(), dp = list(shape = 3),
                     link = "log", b = c(1, 0.8))
)
inv <- list(id = identity, log = exp, logit = stats::plogis)
out <- list()
for (nm in names(fams)) {
  sp <- fams[[nm]]
  sim <- reg[[nm]]()$sim
  for (s in 1:20) {
    set.seed(s)
    x <- rnorm(n)
    mu <- inv[[sp$link]](sp$b[1] + sp$b[2] * x)
    if (is.null(sim)) {
      y <- mu + rnorm(n, 0, sp$dp$sigma)
    } else {
      dpl <- c(list(mu = mu), lapply(sp$dp, rep, n))
      y <- sim(dpl, list(), n)
    }
    d <- data.frame(x = x, y = y)
    st <- tryCatch(unclass(suppressMessages(par_template(
      bf(y ~ x), data = d, family = sp$fam))),
      error = function(e) NULL)
    if (is.null(st)) next
    for (dp in intersect(names(sp$dp), c("sigma", "shape", "phi", "beta"))) {
      v <- st$betad[[paste0(dp, "_(Intercept)")]]
      out[[length(out) + 1L]] <- data.frame(
        family = nm, seed = s, par = dp, truth = sp$dp[[dp]],
        start = exp(v), ratio = exp(v) / sp$dp[[dp]],
        slope = st$beta[["x"]], slope_truth = sp$b[2])
    }
  }
}
X <- do.call(rbind, out)
for (nm in unique(paste(X$family, X$par))) {
  x <- X[paste(X$family, X$par) == nm, ]
  cat(sprintf(paste0("%-34s start / true: median %.3f [%.3f, %.3f] | ",
                     "slope start / true median %.3f\n"),
              nm, median(x$ratio), min(x$ratio), max(x$ratio),
              median(x$slope / x$slope_truth)))
}
