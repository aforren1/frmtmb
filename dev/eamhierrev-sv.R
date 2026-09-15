# REVIEW of lane eamhier, ranked item 1: is the `sv` collapse at seed
# 20260935 a flat likelihood or a failed optimization?
#
# The lane claims the former. This script decides it by CONSTRUCTION:
# a profile over `sv` held fixed on a grid, plus every optimizer
# frm_allfit() offers, plus a restart placed at the truth. A failed
# optimization has a higher likelihood somewhere on that grid; a flat
# one does not.
#
# It also re-fits seed 20260910 first, so that the toolchain is checked
# against a number this lane recorded before anything is concluded from
# a number it did not.
#
# Run: Rscript --vanilla dev/eamhierrev-sv.R

.libPaths(c("C:/Users/adf44/source/r/eamhierrev-lib",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({
  library(frmtmb)
  library(frmtmb.eam)
})

tr <- list(mu0 = 0.4, mu_cond = 0.9, bs = 1.4, ndt = 0.25, sv = 0.4)
mk <- function(seed, n = 12000L) {
  set.seed(seed)
  cond <- rep(0:1, each = n / 2L)
  d <- ddm_simulate(n, mu = tr$mu0 + tr$mu_cond * cond, bs = tr$bs,
                    ndt = tr$ndt, bias = 0.5, sv = tr$sv)
  d$cond <- factor(cond, labels = c("a", "b"))
  d
}
frm_sv <- function(d, ...) {
  frm(bf(rt | dec(upper) ~ cond, bs ~ 1, ndt ~ 1, bias = 0.5),
      family = wiener(variability = "sv"), data = d, ...)
}

## --------------------------------------------- the instrument, first
cat("== instrument: S1 seed 20260910 against its own record ==\n")
t0 <- Sys.time()
# already reproduced to every digit on the first pass; kept behind a
# switch so that a rerun of the expensive half does not pay for it
f10 <- if (nzchar(Sys.getenv("EAMHIERREV_SKIP_INSTRUMENT"))) NULL else
  frm_sv(mk(20260910L), se = TRUE)
if (is.null(f10)) cat("  skipped (reproduced -6657.979468 already)\n")
if (!is.null(f10)) {
cat(sprintf("  %.1f s\n", as.numeric(difftime(Sys.time(), t0,
                                              units = "secs"))))
b <- unlist(fixef(f10))
cat(sprintf("  logLik  %.6f   record -6657.979468\n",
            as.numeric(logLik(f10))))
cat(sprintf("  log sv  %.8f   record -0.83769782\n",
            unname(b[["sv.(Intercept)"]])))
cat(sprintf("  mu_cond %.8f   record 0.9382021016\n",
            unname(b[["mu.condb"]])))
}

## ---------------------------------------------------- the collapse
cat("\n== seed 20260935, sv free ==\n")
d <- mk(20260935L)
fit <- frm_sv(d, se = TRUE)
b <- unlist(fixef(fit))
ll_free <- as.numeric(logLik(fit))
cat(sprintf("  logLik %.6f  conv %d\n", ll_free, fit$opt$convergence))
print(b)
ci <- suppressWarnings(stats::confint(fit))
print(ci)
cat(sprintf("  sv on the natural scale: %.6e   record 4.392e-04\n",
            exp(unname(b[["sv.(Intercept)"]]))))
# the standard errors this fit reports, read off its own Wald
# interval, because that is the only route a user has
se_all <- (ci[, 2L] - ci[, 1L]) / (2 * stats::qnorm(0.975))
print(round(se_all, 6))
jsv <- grep("sv", names(se_all))[1L]
cat(sprintf("  se(log sv) %.2f ; largest OTHER se %.5f ; ratio %.0f\n",
            se_all[jsv], max(se_all[-jsv]),
            se_all[jsv] / max(se_all[-jsv])))
dg <- frmtmb::diagnose(fit, quiet = TRUE)
cat("  diagnose fields that are non-NULL: ",
    paste(names(dg)[!vapply(dg, is.null, logical(1))], collapse = ", "),
    "\n", sep = "")
cat(sprintf("  unbounded_dpar: %s | flat: %d | bad_se: %d | singular: %d\n",
            if (is.null(dg$unbounded_dpar)) "NULL" else "FIRED",
            length(dg$flat), length(dg$bad_se), length(dg$singular)))
cat(sprintf(paste0("  the check's own pair: abs(est) = %.4f (needs ",
                   "> 10), se = %.2f (needs > abs(est)): %s\n"),
            abs(unname(b[["sv.(Intercept)"]])),
            se_all[jsv],
            if (abs(unname(b[["sv.(Intercept)"]])) > 10) "BOTH" else
              "SECOND ONLY"))

## ------------------------------------------- no sv term at all
cat("\n== the same data with no sv term ==\n")
f0 <- frm(bf(rt | dec(upper) ~ cond, bs ~ 1, ndt ~ 1, bias = 0.5),
          family = wiener(), data = d, se = FALSE)
ll_0 <- as.numeric(logLik(f0))
cat(sprintf("  logLik %.8f | free - none = %.3e\n", ll_0,
            ll_free - ll_0))

## ------------------------------------------- the profile, the point
cat("\n== profile: sv HELD at a grid, everything else free ==\n")
grid <- c(1e-6, 1e-4, 1e-3, 0.01, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5)
prof <- vapply(grid, function(v) {
  f <- try(frm(bf(rt | dec(upper) ~ cond, bs ~ 1, ndt ~ 1,
                  bias = 0.5, sv = v),
               family = wiener(variability = "sv"), data = d,
               se = FALSE), silent = TRUE)
  if (inherits(f, "try-error")) return(NA_real_)
  as.numeric(logLik(f))
}, numeric(1))
for (i in seq_along(grid)) {
  cat(sprintf("  sv = %-8g logLik %14.6f   free - this = %10.6f\n",
              grid[i], prof[i], ll_free - prof[i]))
}
cat(sprintf("  best grid point: sv = %g at %.6f; free fit %.6f\n",
            grid[which.max(prof)], max(prof, na.rm = TRUE), ll_free))
above <- any(prof > ll_free + 1e-6, na.rm = TRUE)
cat(sprintf("  VERDICT: any grid point ABOVE the free fit? %s\n",
            if (above) "YES (failed optimization)" else
              "NO (flat/boundary)"))

## ------------------------------------------- restarts
cat("\n== restarts ==\n")
st <- try(frm_sv(d, se = FALSE,
                 start = list(beta = c("sv_(Intercept)" = log(0.4)))),
          silent = TRUE)
if (inherits(st, "try-error")) {
  cat("  start = list(beta = ...) refused: ",
      gsub("[\r\n]+", " ", as.character(st)), "\n", sep = "")
} else {
  cat(sprintf("  restart at sv = 0.4: logLik %.6f, sv = %.6e\n",
              as.numeric(logLik(st)),
              exp(unname(unlist(fixef(st))[["sv.(Intercept)"]]))))
}
af <- try(frmtmb::frm_allfit(fit), silent = TRUE)
if (inherits(af, "try-error")) {
  cat("  frm_allfit refused: ", gsub("[\r\n]+", " ",
                                     as.character(af)), "\n", sep = "")
} else {
  for (nm in names(af$fits)) {
    ff <- af$fits[[nm]]
    if (is.null(ff)) { cat("  ", nm, ": NULL\n", sep = ""); next }
    cat(sprintf("  %-14s logLik %14.6f  sv %.6e  conv %d\n", nm,
                as.numeric(logLik(ff)),
                exp(unname(unlist(fixef(ff))[["sv.(Intercept)"]])),
                ff$opt$convergence))
  }
}
cat("\ndone\n")
