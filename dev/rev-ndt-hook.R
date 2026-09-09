# rev-ndt: was route 4 really closed?
#
# The lane rejected reporting the TIME through `post$dpar_response`
# because that hook is called as `value(dpars, dpar)` and receives no
# aterms, so it cannot reach a per-row bound. That reading is right for
# a per-row bound. It is NOT right for the SCALAR bound, which is every
# model without `ndt_group()`, that is the default and the whole
# backward-compatible majority: the scalar is known at finalize and can
# be captured in the closure.
#
# This builds that hook on top of the shipped family, without touching
# the package, and checks that predict(dpar = "ndt", type = "response")
# then reports seconds again, with a standard error, and that ndt_time()
# still agrees.
#
# Seed 4242, worktree build.

.libPaths(c("C:/Users/adf44/source/r/rev-ndt-lib",
            "C:/Users/adf44/source/r/reflib-r2",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(frmtmb); library(frmtmb.eam)})

set.seed(4242)
d <- ddm_simulate(400, mu = 1.0, bs = 1.4, ndt = 0.30, bias = 0.5)
d$s <- factor(rep(c("a", "b"), length.out = nrow(d)))
form <- bf(rt | dec(upper) ~ 1, bs ~ 1, ndt ~ 1, bias = 0.5)

# the shipped behavior
f0 <- frm(form, family = wiener(), data = d, se = TRUE)
p0 <- suppressWarnings(predict(f0, dpar = "ndt", type = "response",
                               se.fit = TRUE))
cat("shipped: predict(ndt, response) =", sprintf("%.6f", p0$fit[1L]),
    " se", sprintf("%.6f", p0$se.fit[1L]), "\n")
cat("shipped: ndt_time()             =",
    sprintf("%.6f", ndt_time(f0)[1L]), "\n")

# the same family, given the scalar reporting hook a family_finalize
# could install whenever `floors` is NULL
wiener_reporting <- function(...) {
  fam <- wiener(...)
  inner <- fam[["family_finalize"]]
  fam[["family_finalize"]] <- function(fam, y, aterms) {
    out <- inner(fam, y, aterms)
    bd <- out[["ndt_bound"]]
    if (is.null(bd[["floors"]])) {
      ub <- bd[["ub"]]
      out[["post"]][["dpar_response"]] <- list(
        dpars = "ndt",
        value = function(dpars, dpar) dpars[["ndt"]] * ub,
        deriv = function(dpars, dpar) {
          f <- dpars[["ndt"]] / ub
          ub * f * (1 - f)
        })
    }
    out
  }
  fam
}

f1 <- frm(form, family = wiener_reporting(), data = d, se = TRUE)
p1 <- suppressWarnings(predict(f1, dpar = "ndt", type = "response",
                               se.fit = TRUE))
cat("\nwith the hook: predict(ndt, response) =",
    sprintf("%.6f", p1$fit[1L]), " se", sprintf("%.6f", p1$se.fit[1L]),
    "\n")
cat("with the hook: ndt_time()             =",
    sprintf("%.6f", ndt_time(f1)[1L]), "\n")
cat("min(rt) =", sprintf("%.6f", min(d$rt)), "\n")
cat("log-likelihoods equal:",
    identical(as.numeric(logLik(f0)), as.numeric(logLik(f1))), "\n")

# the hook must NOT be installed where the bound is per group, because
# there the value is per row and the hook cannot see the row's group
f2 <- frm(bf(rt | dec(upper) + ndt_group(s) ~ 1, bs ~ 1, ndt ~ 1,
             bias = 0.5), family = wiener_reporting(), data = d,
          se = TRUE)
p2 <- suppressWarnings(predict(f2, dpar = "ndt", type = "response"))
cat("\nndt_group() model, hook not installed: predict =",
    sprintf("%.6f", p2[1L]), " (a fraction, as shipped)\n")
cat("ndt_group() model, ndt_time() =",
    sprintf("%.6f", ndt_time(f2)[1L]), "\n")

# and newdata still works through the hook, because the bound is scalar
nd <- d[1:3, , drop = FALSE]
cat("\nhook on newdata:",
    paste(sprintf("%.6f", suppressWarnings(
      predict(f1, newdata = nd, dpar = "ndt", type = "response"))),
      collapse = " "), "\n")
