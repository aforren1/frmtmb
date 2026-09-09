# rev-ndt: the idempotency guard, and whether the bound really is a
# property of the FITTED data on every family.
#
# ddm_ndt_install() guards a second wrap with `ndt_bound`, and the
# findings say the floor table is captured at finalize so that a refit
# or a prediction on rows that drop a group's fastest trial scores
# against the bound the fit was made with. The guard is checked with the
# guarded thing PRESENT (a second finalize) and ABSENT (the first one),
# and the claim is checked on wiener() and rdm() separately, because
# their family_finalize slots are not the same shape: wiener() rebuilds
# the family from its config and rdm() passes the family through.
#
# Seed 4242, worktree build.

.libPaths(c(Sys.getenv("REV_LIB",
                       "C:/Users/adf44/source/r/rev-ndt-lib2"),
            "C:/Users/adf44/source/r/reflib-r2",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(frmtmb); library(frmtmb.eam)})

set.seed(4242)
d <- ddm_simulate(400, mu = 1.0, bs = 1.4, ndt = 0.30, bias = 0.5)
k <- which.min(d$rt)
lo <- min(d$rt)
cat("min(rt) all rows", sprintf("%.6f", lo),
    " without the fastest", sprintf("%.6f", min(d$rt[-k])), "\n\n")

# ---- 1. a second family_finalize must not wrap twice ----------------
dbl <- function(fam, y, aterms, label) {
  f1 <- fam[["family_finalize"]](fam, y, aterms)
  f2 <- f1[["family_finalize"]](f1, y, aterms)
  dp <- list(mu = 1, bs = 1.4, ndt = 0.5, bias = 0.5)
  at <- c(aterms, list(ndt_floor = NULL))
  at$ndt_floor <- NULL
  v1 <- as.numeric(f1[["lpdf"]](0.6, dp, aterms))
  v2 <- as.numeric(f2[["lpdf"]](0.6, dp, aterms))
  cat(sprintf("%-8s lpdf after one finalize %.12f, after two %.12f  %s\n",
              label, sum(v1), sum(v2),
              if (identical(v1, v2)) "IDENTICAL" else "DIFFERS"))
  cat(sprintf("%-8s bound after one %.8f, after two %.8f\n", label,
              f1[["ndt_bound"]][["ub"]], f2[["ndt_bound"]][["ub"]]))
}
dbl(wiener(), d$rt, list(dec = as.numeric(d$upper)), "wiener")
set.seed(4242)
dr <- rdm_simulate(300, v = c(2.5, 1.5), A = 0.5, k = 0.5, ndt = 0.2)
dbl(rdm(2), dr$rt, list(vint1 = dr$choice), "rdm")

# ---- 2. a second finalize on DIFFERENT rows ------------------------
# the leave-one-out shape the findings name: does the bound follow the
# new rows, or stay the fitted one?
cat("\n")
moved <- function(fam, y, aterms, y2, aterms2, label) {
  f1 <- fam[["family_finalize"]](fam, y, aterms)
  f2 <- f1[["family_finalize"]](f1, y2, aterms2)
  cat(sprintf("%-8s bound on all rows %.8f, refinalized without the ",
              label, f1[["ndt_bound"]][["ub"]]))
  cat(sprintf("fastest %.8f  %s\n", f2[["ndt_bound"]][["ub"]],
              if (isTRUE(all.equal(f1[["ndt_bound"]][["ub"]],
                                   f2[["ndt_bound"]][["ub"]])))
                "KEPT" else "RE-DERIVED"))
}
moved(wiener(), d$rt, list(dec = as.numeric(d$upper)),
      d$rt[-k], list(dec = as.numeric(d$upper)[-k]), "wiener")
kr <- which.min(dr$rt)
moved(rdm(2), dr$rt, list(vint1 = dr$choice),
      dr$rt[-kr], list(vint1 = dr$choice[-kr]), "rdm")

# ---- 3. what a user sees: predict() on rows missing the fastest -----
cat("\n")
fit <- frm(bf(rt | dec(upper) ~ 1, bs ~ 1, ndt ~ 1, bias = 0.5),
           family = wiener(), data = d)
a <- suppressWarnings(predict(fit, newdata = d[1:5, ], dpar = "ndt",
                              type = "response"))
b <- suppressWarnings(predict(fit, newdata = d[-k, ][1:5, ],
                              dpar = "ndt", type = "response"))
cat("predict(ndt) on 5 rows, full data  ", sprintf("%.8f", a[1L]), "\n")
cat("predict(ndt) on 5 rows, fastest out", sprintf("%.8f", b[1L]), "\n")
fa <- suppressWarnings(predict(fit, newdata = d[1:5, ],
                               type = "response"))
fb <- suppressWarnings(predict(fit, newdata = d[-k, ][1:5, ],
                               type = "response"))
cat("fitted mean, full data             ", sprintf("%.8f", fa[1L]), "\n")
cat("fitted mean, fastest out           ", sprintf("%.8f", fb[1L]), "\n")
cat("ndt_time, full data                ",
    sprintf("%.8f", ndt_time(fit, d[1:5, ])[1L]), "\n")
cat("ndt_time, fastest out              ",
    sprintf("%.8f", ndt_time(fit, d[-k, ][1:5, ])[1L]), "\n")
