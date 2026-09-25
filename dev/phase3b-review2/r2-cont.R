# Reviewer 2, item 4: decision (a), contaminant_range =, and item 1's
# neighbour: contaminant + cens(). No optimization (dry_run).
source("dev/phase3b-review2/r2-prelude.R")
r2_lib("lane")
suppressPackageStartupMessages({library(frmtmb); library(frmtmb.eam); library(RWiener)})
options(width = 160)
try_obj <- function(expr) tryCatch(expr, error = function(e) e)
rep_err <- function(lab, x) {
  if (inherits(x, "error")) {
    cat(sprintf("%-44s REFUSED [%s] %s\n", lab, paste(class(x), collapse = ","),
                substr(conditionMessage(x), 1, 150)))
  } else {
    fam <- stats::family(x)
    f0 <- x$obj$fn(x$obj$par); g0 <- x$obj$gr(x$obj$par)
    cat(sprintf("%-44s built; range %s; fn %.6g; grad finite %s\n", lab,
                paste(format(fam$contaminant_range, digits = 6), collapse = ".."),
                f0, all(is.finite(g0))))
  }
  invisible(x)
}
set.seed(21)
d <- ddm_simulate(300, mu = 0.8, bs = 1.4, ndt = 0.3, bias = 0.5)
hit <- runif(300) < 0.08
d$rt[hit] <- runif(sum(hit), 0.1, 2.5); d$upper[hit] <- rbinom(sum(hit), 1, 0.5)
d <- d[d$rt < 2.5, ]
fb <- function(extra = "") stats::as.formula(sprintf(
  "rt | dec(upper)%s ~ 1", extra))
mk <- function(extra = "", ...) bf(fb(extra), bs ~ 1, ndt ~ 1, lambda ~ 1, bias = 0.5)

rep_err("no range, no trunc", try_obj(frm(mk(), family = wiener(contaminant = TRUE),
                                          data = d, dry_run = "objective")))
rep_err("no range, trunc(lb = 0.05) only", try_obj(frm(mk(" + trunc(lb = 0.05)"),
  family = wiener(contaminant = TRUE), data = d, dry_run = "objective")))
rep_err("no range, trunc(ub = 2.5)", try_obj(frm(mk(" + trunc(ub = 2.5)"),
  family = wiener(contaminant = TRUE), data = d, dry_run = "objective")))
rep_err("no range, trunc(ub = 2.5), max_ndt 0.5", try_obj(frm(mk(" + trunc(ub = 2.5)"),
  family = wiener(contaminant = TRUE, max_ndt = 0.5), data = d, dry_run = "objective")))
d$ub <- ifelse(seq_len(nrow(d)) %% 2 == 0, 2.5, 2.6)
rep_err("no range, trunc(ub = ub) varying", try_obj(frm(mk(" + trunc(ub = ub)"),
  family = wiener(contaminant = TRUE), data = d, dry_run = "objective")))
d$code <- as.integer(d$rt > 1.5); d$rt[d$code == 1] <- 1.5
rep_err("no range, cens()", try_obj(frm(mk(" + cens(code)"),
  family = wiener(contaminant = TRUE), data = d, dry_run = "objective")))
rep_err("range c(0, 2.5), cens()", oc <- try_obj(frm(mk(" + cens(code)"),
  family = wiener(contaminant = TRUE, contaminant_range = c(0, 2.5)), data = d,
  dry_run = "objective")))
rep_err("range c(0, 2.5), cens(), trunc(ub = 2.5)", try_obj(frm(
  mk(" + cens(code) + trunc(ub = 2.5)"),
  family = wiener(contaminant = TRUE, contaminant_range = c(0, 2.5)), data = d,
  dry_run = "objective")))
rep_err("no range, cens(), trunc(ub = 2.5)", try_obj(frm(
  mk(" + cens(code) + trunc(ub = 2.5)"),
  family = wiener(contaminant = TRUE), data = d, dry_run = "objective")))

cat("\n-- right-censored + contaminant against a hand-written mixture (RWiener)\n")
if (!inherits(oc, "error")) {
  fam <- stats::family(oc); ub <- fam$ndt_bound$ub
  cat("parameters:", names(oc$obj$par), " ndt ub", ub, "\n")
  hand <- function(p) {
    v <- p[1]; a <- exp(p[2]); t0 <- ub * plogis(p[3]); lam <- plogis(p[4])
    lo <- 0; hi <- 2.5
    resp <- ifelse(d$upper == 1, "upper", "lower")
    ll <- vapply(seq_len(nrow(d)), function(i) {
      y <- d$rt[i]
      if (d$code[i] == 1L) {
        log((1 - lam) * (1 - pwiener(y, a, t0, 0.5, v, resp = "both")) +
              lam * (1 - (y - lo) / (hi - lo)))
      } else {
        f <- if (y > t0) dwiener(y, a, t0, 0.5, v, resp = resp[i]) else 0
        log((1 - lam) * f + lam / (2 * (hi - lo)))
      }
    }, 0)
    -sum(ll)
  }
  for (p in list(oc$obj$par, oc$obj$par + c(0.2, 0.1, 0.3, 0.5))) {
    cat(sprintf("lane %.10f  hand %.10f  rel %.2e\n", oc$obj$fn(p), hand(p),
                abs(oc$obj$fn(p) - hand(p)) / abs(hand(p))))
  }
}

cat("\n-- a range that does not cover the data\n")
d2 <- d; d2$code <- NULL
cat("rows above 2.0:", sum(d2$rt > 2.0), " rows below 0.2:", sum(d2$rt < 0.2),
    " min rt", min(d2$rt), "\n")
o1 <- rep_err("range c(0.25, 2.0), default max_ndt", try_obj(frm(mk(),
  family = wiener(contaminant = TRUE, contaminant_range = c(0.25, 2.0)), data = d2,
  dry_run = "objective")))
o2 <- rep_err("range c(0.25, 2.0), max_ndt 0.5", try_obj(frm(mk(),
  family = wiener(contaminant = TRUE, contaminant_range = c(0.25, 2.0), max_ndt = 0.5),
  data = d2, dry_run = "objective")))
o3 <- rep_err("range c(0, 2.0), max_ndt 0.5", try_obj(frm(mk(),
  family = wiener(contaminant = TRUE, contaminant_range = c(0, 2.0), max_ndt = 0.5),
  data = d2, dry_run = "objective")))
# per-row values through the family, at ndt 0.3 and lambda 0.05
for (nm in c("o1", "o2")) {
  x <- get(nm); if (inherits(x, "error")) next
  fam <- stats::family(x)
  dp <- list(mu = 0.8, bs = 1.4, ndt = 0.3, bias = 0.5, lambda = 0.05,
             .eta_lambda = qlogis(0.05))
  lp <- fam$lpdf(d2$rt, dp, list(dec = d2$upper))
  cat(sprintf("%s: rows with log density < -1e3: %d (rt %s); min %.4g; any NaN %s\n", nm,
              sum(lp < -1e3), paste(format(d2$rt[lp < -1e3], digits = 3), collapse = " "),
              min(lp), anyNA(lp)))
}
