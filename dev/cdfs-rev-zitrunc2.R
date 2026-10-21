# Reviewer: trace the NaN of zero_inflated_poisson / _negbinomial with
# trunc(ub = 30) along the optimizer's own path.
.libPaths(c("C:/Users/adf44/source/r/wt-cdfs-rev-lib",
            "C:/Users/adf44/source/r/rellib-r10",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
options(frmtmb.notices = FALSE)
ns <- asNamespace("frmtmb")
set.seed(5)
n <- 60
x <- rnorm(n)
d <- data.frame(x = x, yp = rpois(n, exp(1 + 0.3 * x)))
bad <- NULL
for (fm in list(zero_inflated_poisson(), zero_inflated_negbinomial())) {
  ob <- frm(yp | trunc(ub = 30) ~ x, family = fm, data = d,
            dry_run = "objective")$obj
  cat("par names:", names(ob$par), "\n")
  gr <- function(p) {
    g <- ob$gr(p)
    if (any(!is.finite(g)) && is.null(bad)) {
      bad <<- p
      cat("non-finite gradient at par", signif(p, 8), " fn", ob$fn(p), "\n")
    }
    g
  }
  o <- tryCatch(nlminb(ob$par, ob$fn, gr), error = function(e) conditionMessage(e))
  cat(ns$as_frmtmb_family(fm)$family, "nlminb:", if (is.list(o)) o$message else o, "\n")
  if (!is.null(bad)) {
    f <- ns$as_frmtmb_family(fm)
    eta <- bad[1] + bad[2] * d$x
    mu <- exp(eta)
    lz <- bad[length(bad)]
    zi <- plogis(lz)
    dp <- list(mu = mu, zi = zi)
    if (length(bad) == 4) dp$shape <- exp(bad[3])
    lpdf <- f$lpdf(d$yp, dp, list())
    F30 <- f$lcdf(rep(30, n), dp, list())
    cat("mu range", range(mu), " zi", zi, "\n")
    cat("lpdf range", range(lpdf), " any nonfinite", any(!is.finite(lpdf)), "\n")
    cat("F(30) range", range(F30), " rows with F(30) = 0:", sum(F30 == 0),
        " nonfinite:", sum(!is.finite(F30)), "\n")
    i <- which.min(F30)
    cat("worst row mu", mu[i], " y", d$yp[i], " F30", F30[i], "\n")
    cat("  base poisson lF at q = 30, mu", mu[i], ":",
        ns$poisson_cdf_log(30, mu[i])$lF, " stats", ppois(30, mu[i], log.p = TRUE), "\n")
    cat("  zi_cdf_log lF:", ns$zi_poisson_cdf_pair(30, list(mu = mu[i], zi = zi), list(),
                                                    ns$dpar_link("logit", "zi", "zero_inflated_poisson",
                                                                 ns$dpar_links_unit))$lF, "\n")
    bad <- NULL
  }
}
