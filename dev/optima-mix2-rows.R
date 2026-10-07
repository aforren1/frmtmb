# Lane optima, item 2: mixture(cumulative("probit"), acat()) on the
# review's data stops with a NaN gradient on the lane build. At the
# point where nlminb first gets one: which rows of component 1 carry
# it, and what their thresholds and linear predictor are.
.libPaths(c("C:/Users/adf44/source/r/wt-optima-lib",
            "C:/Users/adf44/source/r/rellib-r6",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
src <- readLines("dev/optima-mix3.R")
eval(parse(text = src[grep("^set.seed", src):grep("^d <- data.frame", src)]))
u <- frm(bf(y ~ x), family = mixture(cumulative("probit"), acat()),
         data = d, dry_run = "objective")
obj <- u$obj
last <- NULL
o <- tryCatch(stats::nlminb(obj$par, function(p) {
  v <- obj$fn(p); if (is.nan(v)) Inf else v
}, function(p) {
  g <- obj$gr(p); last <<- p; g
}, control = list(eval.max = 1000, iter.max = 1000)),
error = function(e) conditionMessage(e))
p0 <- last
pl <- obj$env$parList(p0)
tau <- frmtmb:::ord_tau_from_raw(pl$tau_raw1, TRUE)
cat("tau1:", format(tau, digits = 17), "\n")
u2 <- u
u2$estimates <- pl
class(u2) <- setdiff(class(u2), "frmtmb_unfitted")
resp <- names(u2$spec$responses)[1]
famo <- u2$spec$responses[[resp]]$family
dp <- frmtmb:::eval_dpars(u2)[[resp]]
av <- u2$frame[["aterm_values"]][[resp]]
ex <- frmtmb:::fit_extras(u2, resp)
y <- u2$frame[["y"]][[resp]]
f1 <- function(v) {
  dpv <- dp
  dpv[["mu1"]] <- v[1] * d$x
  exv <- ex
  exv$tau_raw1 <- v[2:4]
  famo[["mix"]][["comp_lpdf"]](y, dpv, av, 1L, exv)
}
v0 <- c(pl$beta[1], pl$tau_raw1)
bad <- which(vapply(seq_along(y), function(i) {
  any(!is.finite(RTMB::MakeTape(function(v) f1(v)[i], v0)$jacobian(v0)))
}, NA))
cat("rows with a non-finite gradient:", length(bad), "\n")
for (i in head(bad, 4)) {
  ti <- RTMB::MakeTape(function(v) f1(v)[i], v0)
  cat("row", i, "y", y[i], "eta", format(dp[["mu1"]][i], digits = 17),
      "value", format(ti(v0), digits = 10), "jac",
      format(ti$jacobian(v0)), "\n")
}
xK <- tau[3] - dp[["mu1"]][1]
lk <- frmtmb:::get_link("probit")
pieces <- list(
  pnorm_lo = function(x) RTMB::pnorm(-x, log.p = TRUE),
  pnorm_hi = function(x) RTMB::pnorm(x, log.p = TRUE),
  q = function(x) lk$logit_eta(x),
  l1m = function(x) frmtmb:::log1m_inv_logit(lk$logit_eta(x)),
  lil = function(x) frmtmb:::log_inv_logit(lk$logit_eta(x)),
  masked = function(x) 0 * frmtmb:::log1m_inv_logit(lk$logit_eta(x)))
for (nm in names(pieces)) {
  tp <- RTMB::MakeTape(pieces[[nm]], xK)
  cat(sprintf("%-9s at %.17g: value %s, derivative %s\n", nm, xK,
              format(tp(xK)), format(tp$jacobian(xK))))
}
for (gf in c(FALSE, TRUE)) for (yy in 1:4) {
  g <- function(v) {
    tau_v <- frmtmb:::ord_tau_from_raw_ad(v[2:4])
    frmtmb:::ord_cumulative_logpmf(yy, v[1] * d$x[1], tau_v, lk, 1,
                                   gap_floor = gf)
  }
  tp <- RTMB::MakeTape(g, v0)
  cat("gap_floor", gf, "y", yy, ": value", format(tp(v0), digits = 8),
      "jac", format(tp$jacobian(v0), digits = 4), "\n")
}
for (gf in c(FALSE, TRUE)) {
  g <- function(v) {
    tau_v <- frmtmb:::ord_tau_from_raw_ad(v[2:4])
    sum(frmtmb:::ord_cumulative_logpmf(y, v[1] * d$x, tau_v, lk, 1,
                                       gap_floor = gf))
  }
  tp <- RTMB::MakeTape(g, v0)
  cat("vectorized, gap_floor", gf, ": value", format(tp(v0), digits = 8),
      "jac", format(tp$jacobian(v0), digits = 4), "\n")
  for (k in 1:4) {
    gk <- function(v) {
      tau_v <- frmtmb:::ord_tau_from_raw_ad(v[2:4])
      sum(frmtmb:::ord_cumulative_logpmf(y[y == k], v[1] * d$x[y == k],
                                         tau_v, lk, 1, gap_floor = gf))
    }
    tk <- RTMB::MakeTape(gk, v0)
    cat("   rows y ==", k, ": jac", format(tk$jacobian(v0), digits = 4), "\n")
  }
}
