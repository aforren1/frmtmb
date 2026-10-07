# Reviewer of lane optima, claim 4: ord_touch_nan() and the cs()
# threshold path leave a plain cumulative fit bitwise where it was.
# One process (pitfall 21): the lane's two ordinal functions against
# base's, swapped into the namespace.
#   Rscript dev/optima-rev-touch.R
.libPaths(c("C:/Users/adf44/source/r/wt-optima-lib",
            "C:/Users/adf44/source/r/rellib-r6",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
be <- new.env(parent = ns)
sys.source("dev/optima-rev-out/base-families.R", envir = be)
fns <- c("ord_cumulative_logpmf", "ord_cumulative_logpmf_cs")
lane <- mget(fns, envir = ns)
swap <- function(src) {
  for (f in fns) {
    unlockBinding(f, ns)
    g <- src[[f]]
    environment(g) <- ns
    assign(f, g, envir = ns)
    lockBinding(f, ns)
  }
}
mk <- function(seed, n = 300) {
  set.seed(seed)
  d <- data.frame(x = rnorm(n), z = rnorm(n))
  lat <- 1.2 * d$x + rlogis(n)
  d$y <- 1L + (lat > -1 + 0.2 * d$z) + (lat > 0) + (lat > 1 - 0.2 * d$z)
  d
}
cases <- list(
  cum_cs_logit = function(d) frm(bf(y ~ x + cs(z)), family = cumulative(),
                                 data = d),
  cum_cs_probit = function(d) frm(bf(y ~ x + cs(z)),
                                  family = cumulative("probit"), data = d),
  cum_cs_cloglog = function(d) frm(bf(y ~ x + cs(z)),
                                   family = cumulative("cloglog"), data = d),
  cum_cs_cauchit = function(d) frm(bf(y ~ x + cs(z)),
                                   family = cumulative("cauchit"), data = d),
  cum_plain = function(d) frm(bf(y ~ x + z), family = cumulative(), data = d),
  mix_nocs = function(d) frm(bf(y ~ x), family = mixture(cumulative(),
                                                         sratio()), data = d)
)
for (s in 1:5) {
  d <- mk(s)
  for (nm in names(cases)) {
    swap(lane)
    a <- tryCatch(suppressWarnings(cases[[nm]](d)), error = function(e) e)
    swap(be)
    b <- tryCatch(suppressWarnings(cases[[nm]](d)), error = function(e) e)
    swap(lane)
    if (inherits(a, "error") || inherits(b, "error")) {
      cat(nm, s, "ERROR\n")
      next
    }
    cat(sprintf("%-15s seed %d logLik identical %s par identical %s code %d\n",
                nm, s, identical(logLik(a), logLik(b)),
                identical(a$opt$par, b$opt$par), a$opt$convergence))
  }
}
