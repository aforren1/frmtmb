# Review of lane latperf: the sparse joint covariance against two references,
# in one process (dev/rtmb-pitfalls.md item 21), on designs the lane did not
# use.
#   ref   an independent dense V built here from joint_precision() alone:
#         solve(Q), or the block formula with the SE check's Vf where the
#         check kept one (no frmtmb reader involved)
#   dense the lane's own dense route: jcs_joint_cov() mocked to NULL, the
#         cache without Vjoint
# Also: the selected inverse against solve(Q_rr) at every stored entry (so a
# permutation bug shows), and how many rows fell to the direct solve.
# Usage: Rscript dev/review-latperf-algebra.R [design ...]  (seed 20261106)
.libPaths(c("C:/Users/adf44/source/r/wt-latperf-lib",
            "C:/Users/adf44/source/r/rellib-r18",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat(sprintf("lib: %s frmtmb %s\n", dirname(dirname(system.file(package =
  "frmtmb"))), format(packageVersion("frmtmb"))))
ns <- asNamespace("frmtmb")

mk <- function(seed = 20261106, G = 30, K = 14, sdid = 0.4) {
  set.seed(seed)
  d <- expand.grid(t = seq_len(K), g = factor(seq_len(G)))
  d <- d[order(d$g, d$t), ]
  n <- nrow(d)
  d$x <- stats::rnorm(n)
  d$z <- stats::runif(n)
  d$id <- factor(sample.int(8, n, TRUE))
  d$a <- factor(sample.int(7, n, TRUE))
  d$b <- factor(sample.int(9, n, TRUE))
  d$expo <- stats::runif(n, 0.5, 2)
  d$w <- stats::runif(n, 0.5, 1.5)
  d$sdx <- 0.3
  xt <- stats::rnorm(n)
  d$xo <- xt + stats::rnorm(n, 0, 0.3)
  e <- as.vector(replicate(G, stats::arima.sim(list(ar = 0.6), K, sd = 0.5)))
  e2 <- as.vector(replicate(G, stats::arima.sim(list(ar = 0.4), K, sd = 0.5)))
  u <- stats::rnorm(8, 0, sdid)
  ua <- stats::rnorm(7, 0, 0.3)
  ub <- stats::rnorm(9, 0, 0.3)
  sl <- stats::rnorm(8, 0, 0.2)
  eta <- 0.5 + 0.4 * d$x + e + u[d$id]
  d$cnt <- stats::rpois(n, exp(eta + sin(3 * d$z)))
  d$cnt0 <- stats::rpois(n, exp(0.5 + 0.4 * d$x + e))
  d$cnto <- stats::rpois(n, d$expo * exp(eta))
  d$cntx <- stats::rpois(n, exp(0.5 + 0.4 * xt + e))
  d$cntc <- stats::rpois(n, exp(0.5 + 0.4 * d$x + e + ua[d$a] + ub[d$b]))
  d$cnt2 <- stats::rpois(n, exp(0.2 + 0.3 * d$x + e2))
  d$tr <- 6L
  d$k <- stats::rbinom(n, d$tr, stats::plogis(eta - 0.5 + sl[d$id] * d$x))
  d$nb <- stats::rnbinom(n, mu = exp(eta), size = 3)
  d$pos <- exp(0.3 + 0.2 * d$x + e + stats::rnorm(n, 0, exp(-1 + d$z)))
  d
}

designs <- list(
  smooth = list(cnt ~ x + (1 | id) + s(z) + ar(t, g, cov = TRUE), poisson()),
  offwt = list(cnto | weights(w) ~ x + offset(log(expo)) + (1 | id) +
                 ar(t, g, cov = TRUE), poisson()),
  me = list(cntx ~ me(xo, sdx) + ar(t, g, cov = TRUE), poisson()),
  nlshape = list(bf(nb ~ x + ar(t, g, cov = TRUE)) + nlf(shape ~ exp(la)) +
                   lf(la ~ 1 + z), negbinomial()),
  sigma = list(bf(pos ~ x + ar(t, g, cov = TRUE), sigma ~ z), lognormal()),
  mv = list(mvbf(bf(cnt0 ~ x + ar(t, g, cov = TRUE)),
                 bf(cnt2 ~ x + ar(t, g, cov = TRUE))), list(poisson(),
                                                           poisson())),
  arma22 = list(cnt0 ~ x + arma(t, g, p = 2, q = 2, cov = TRUE), poisson()),
  ar2 = list(cnt0 ~ x + (1 | id) + ar(t, g, p = 2), poisson()),
  unstr = list(cnt0 ~ x + (1 | id) + unstr(t, g), poisson()),
  crossed = list(cntc ~ x + (1 | a) + (1 | b) + ar(t, g, cov = TRUE),
                 poisson()),
  slope = list(k | trials(tr) ~ x + (1 + x | id) + ar(t, g, cov = TRUE),
               binomial()),
  # no id effect in the data: (1 | id) runs to its boundary, so the SE check
  # keeps a propagation covariance (the repair branch)
  boundary = list(cnt0 ~ x + (1 | id) + ar(t, g, cov = TRUE), poisson())
)
args <- commandArgs(TRUE)
if (length(args)) designs <- designs[args]

dense_copy <- function(fit) {
  f2 <- fit
  keep <- setdiff(ls(fit$cache), "Vjoint")
  f2$cache <- list2env(mget(keep, envir = fit$cache))
  f2
}
with_dense <- function(expr) {
  old <- get("jcs_joint_cov", envir = ns)
  unlockBinding("jcs_joint_cov", ns)
  assign("jcs_joint_cov", function(fit) NULL, envir = ns)
  on.exit({
    assign("jcs_joint_cov", old, envir = ns)
    lockBinding("jcs_joint_cov", ns)
  })
  expr
}
nsolve <- 0L
old_solve <- get("jcs_rowvar_solve", envir = ns)
unlockBinding("jcs_rowvar_solve", ns)
assign("jcs_rowvar_solve", function(S, Ar, ri) {
  nsolve <<- nsolve + nrow(Ar)
  old_solve(S, Ar, ri)
}, envir = ns)

rel <- function(a, b) {
  a <- as.numeric(a)
  b <- as.numeric(b)
  ok <- is.finite(a) & is.finite(b) & b != 0
  if (!any(ok)) return(sprintf("NA (0, mis %d)", sum(is.finite(a) != is.finite(b))))
  sprintf("%.2e (%d, mis %d)", max(abs(a[ok] / b[ok] - 1)), sum(ok),
          sum(is.finite(a) != is.finite(b)))
}
ee <- function(v) {
  if (is.null(v)) return(NA)
  if (is.list(v) && !is.data.frame(v) && !is.null(v$se.fit)) return(v$se.fit)
  if (length(dim(v)) == 3L) return(v[, "Est.Error", ])
  if (is.matrix(v) && "Est.Error" %in% colnames(v)) return(v[, "Est.Error"])
  if (is.data.frame(v) && "Est.Error" %in% names(v)) return(v$Est.Error)
  NA
}
quiet <- function(expr) {
  suppressWarnings(suppressMessages(tryCatch(expr, error = function(e) {
    structure(conditionMessage(e), class = "err")
  })))
}
show <- function(x) if (inherits(x, "err")) paste("ERROR", substr(x, 1, 150)) else x

d <- mk()
nd <- data.frame(t = 1:4, g = factor(rep("new", 4)), x = c(-1, 0, 0.5, 1),
                 z = c(0.2, 0.4, 0.6, 0.8),
                 id = factor(c("1", "2", "99", "98")),
                 a = factor(c("1", "2", "77", "3")),
                 b = factor(c("1", "55", "2", "3")), expo = 1, w = 1,
                 xo = c(0, 0.5, -0.5, 1), sdx = 0.3, tr = 6L)

for (nm in names(designs)) {
  ds <- designs[[nm]]
  cat("\n== design", nm, "\n")
  fit <- quiet(frm(ds[[1]], family = ds[[2]], data = d))
  if (inherits(fit, "err")) {
    cat("  fit:", show(fit), "\n")
    next
  }
  sdr <- frmtmb:::sdr_of(fit)
  jc <- quiet(frmtmb:::get_joint_cov(fit))
  if (inherits(jc, "err")) {
    cat("  get_joint_cov:", show(jc), "\n")
    next
  }
  sparse <- frmtmb:::jc_is_sparse(jc)
  cat(sprintf("  n %d code %s sparse %s cov_fixed_prop %s se_lost %d nb %d\n",
              nrow(d), fit$opt$convergence, sparse,
              !is.null(sdr$cov_fixed_prop), length(sdr$se_lost),
              length(fit$estimates[["b"]])))
  if (!sparse) next
  S <- jc$V
  # ---- the independent dense reference ------------------------------------
  Q <- as.matrix(frmtmb:::joint_precision(fit))
  r <- S$r
  f <- S$f
  Vp <- sdr$cov_fixed_prop
  Wrr <- solve(Q[r, r])
  if (!is.null(Vp) && nrow(Vp) == length(f)) {
    J <- -Wrr %*% Q[r, f]
    Vref <- matrix(0, nrow(Q), nrow(Q))
    Vref[r, r] <- Wrr + J %*% Vp %*% t(J)
    Vref[r, f] <- J %*% Vp
    Vref[f, r] <- t(Vref[r, f])
    Vref[f, f] <- Vp
    cat("  reference: block formula with the SE check's Vf\n")
  } else {
    Vref <- solve(Q)
    cat("  reference: solve(Q)\n")
  }
  # ---- the selected inverse, at every stored entry -------------------------
  keys <- frmtmb:::jcs_selinv_keys(S)
  ki <- floor(keys$key / keys$N) + 1
  kj <- keys$key - (ki - 1) * keys$N + 1
  wv <- Wrr[cbind(ki, kj)]
  cat(sprintf("  selinv: %d stored entries, max |sel - solve| / max|W| %.2e\n",
              length(keys$x), max(abs(keys$x - wv)) / max(abs(Wrr))))
  # ---- whole V --------------------------------------------------------------
  Vs <- quiet(frm_joint_cov(fit)$V)
  if (!inherits(Vs, "err")) {
    ok <- is.finite(Vs)
    cat(sprintf("  frm_joint_cov V: max |sparse - ref| / max|ref| %.2e, nonfinite %d\n",
                max(abs(Vs[ok] - Vref[ok])) / max(abs(Vref)), sum(!ok)))
  } else cat("  frm_joint_cov:", show(Vs), "\n")
  # ---- fitted rows: linpred SE against the reference -----------------------
  nsolve <- 0L
  rsp <- if (nm == "mv") "cnt0" else NULL
  lp <- quiet(frm_linpred(fit, se.fit = TRUE, resp = rsp))
  cat("  rows sent to direct solve by frm_linpred:", nsolve, "\n")
  lb <- quiet(frm_lp_basis(fit, resp = rsp))
  if (!inherits(lb, "err") && !inherits(lp, "err")) {
    A <- as.matrix(lb$A)
    pos <- lb$coef_pos
    se_ref <- sqrt(pmax(rowSums((A %*% Vref[pos, pos]) * A), 0) + lb$extra_var)
    cat("  linpred se vs ref:", rel(lp$se.fit, se_ref), "\n")
  } else {
    cat("  linpred/lp_basis:", if (inherits(lp, "err")) show(lp) else "ok",
        "|", if (inherits(lb, "err")) show(lb) else "ok", "\n")
  }
  # ---- every reader against the lane's own dense route ----------------------
  fd <- dense_copy(fit)
  rd <- list(
    fitted = list(quiet(fitted(fit)), with_dense(quiet(fitted(fd)))),
    resid = list(quiet(residuals(fit)), with_dense(quiet(residuals(fd)))),
    linpred = list(lp, with_dense(quiet(frm_linpred(fd, se.fit = TRUE, resp = rsp)))),
    newdata = list(quiet(fitted(fit, newdata = nd, allow_new_levels = TRUE)),
                   with_dense(quiet(fitted(fd, newdata = nd,
                                           allow_new_levels = TRUE)))),
    newlp = list(quiet(frm_linpred(fit, newdata = nd, se.fit = TRUE, resp = rsp,
                                   allow_new_levels = TRUE)),
                 with_dense(quiet(frm_linpred(fd, newdata = nd, se.fit = TRUE, resp = rsp,
                                              allow_new_levels = TRUE)))),
    curve = list(quiet(frm_curve(fit, newdata = nd, resp = rsp, allow_new_levels = TRUE)),
                 with_dense(quiet(frm_curve(fd, newdata = nd, resp = rsp,
                                            allow_new_levels = TRUE)))))
  for (k in names(rd)) {
    a <- rd[[k]][[1]]
    b <- rd[[k]][[2]]
    if (inherits(a, "err") || inherits(b, "err")) {
      cat(sprintf("  %-8s sparse: %s | dense: %s\n", k,
                  if (inherits(a, "err")) show(a) else "ok",
                  if (inherits(b, "err")) show(b) else "ok"))
      next
    }
    av <- if (k == "curve") a$.se else ee(a)
    bv <- if (k == "curve") b$.se else ee(b)
    cat(sprintf("  %-8s vs dense: %s\n", k, rel(av, bv)))
  }
}
cat("\nDONE\n")
