# test-families-extreme.R's GEV reference moves from evd::dgev to
# extraDistr::dgev (review finding 10). On the test's own grid: the two
# references against each other, frmtmb's taped density against the
# new one under the test's own tolerance, and two wrong densities
# (sigma's log-Jacobian dropped; the sign of xi flipped) against the
# same tolerance, which must fail it.
#   Rscript dev/smallfix-gev.R <lib;lib>
# Log: dev/smallfix-out/gev.txt
a <- commandArgs(trailingOnly = TRUE)
.libPaths(c(strsplit(a[1], ";", fixed = TRUE)[[1]],
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("frmtmb from", dirname(find.package("frmtmb")), "| evd",
    format(packageVersion("evd")), "| extraDistr",
    format(packageVersion("extraDistr")), "\n")
eps <- .Machine$double.eps
tape_lpdf <- function(fam, y, eta, what = "lpdf") {
  dp <- names(eta)
  len <- vapply(eta, length, 1L)
  f <- function(p) {
    pars <- list()
    at <- 0L
    for (nm in dp) {
      e <- p[at + seq_len(len[[nm]])]
      at <- at + len[[nm]]
      pars[[nm]] <- fam$links[[nm]]$linkinv(e)
      pars[[paste0(".eta_", nm)]] <- e
    }
    fam[[what]](y, pars, list())
  }
  tp <- RTMB::MakeTape(f, unlist(eta, use.names = FALSE))
  tp(unlist(eta, use.names = FALSE))
}
g <- expand.grid(y = c(-2, 0.4, 3, 12), mu = c(-1, 0.5), sigma = c(0.5, 2),
                 xi = c(-0.3, 0.05, 0.6))
fam <- gen_extreme_value()
xi <- expm1(log1p(g$xi))
sg <- exp(log(g$sigma))
inside <- 1 + xi * (g$y - g$mu) / sg > 0
r_evd <- mapply(function(y, m, s, k) {
  evd::dgev(y, loc = m, scale = s, shape = k, log = TRUE)
}, g$y, g$mu, sg, xi)
r_ed <- mapply(function(y, m, s, k) {
  extraDistr::dgev(y, mu = m, sigma = s, xi = k, log = TRUE)
}, g$y, g$mu, sg, xi)
cat("rows", nrow(g), "inside", sum(inside), "\n")
cat("max |evd - extraDistr| / max(1, |evd|), inside:",
    format(max(abs(r_evd - r_ed)[inside] / pmax(1, abs(r_evd[inside]))),
           digits = 4), "\n")
cat("outside: evd", unique(r_evd[!inside]), "extraDistr",
    unique(r_ed[!inside]), "\n")
got <- tape_lpdf(fam, g$y, list(mu = g$mu, sigma = log(g$sigma),
                                xi = log1p(g$xi)))
zx <- abs(xi * (g$y - g$mu) / sg)
tol <- 64 * eps * pmax(1, abs(r_ed)) / pmin(1, zx)
pass <- function(x) all(abs(x - r_ed)[inside] <= tol[inside])
cat("frmtmb against extraDistr, the test's assertion:", pass(got),
    "| worst |got - ref| / tol:",
    format(max((abs(got - r_ed) / tol)[inside]), digits = 4), "\n")
w1 <- got + log(sg)
w2 <- mapply(function(y, m, s, k) {
  extraDistr::dgev(y, mu = m, sigma = s, xi = -k, log = TRUE)
}, g$y, g$mu, sg, xi)
cat("wrong: sigma's log-Jacobian dropped:", pass(w1), "\n")
cat("wrong: xi with its sign flipped (finite rows):",
    all(abs(w2 - r_ed)[inside & is.finite(w2)] <=
          tol[inside & is.finite(w2)]), "\n")
