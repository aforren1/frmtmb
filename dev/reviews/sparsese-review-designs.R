# Review of lane sparsese: the sparse route against the dense one on
# designs the lane did not test, small enough that both run.
#   Rscript sparsese-review-designs.R <design> [seed]
# Each design is fitted twice in one process, its readers called under
# options(frmtmb.joint_root_max = 1e6) ("dense": sdr_unwhiten()'s factor)
# and = 20 ("sparse": joint_root_rows() and jcs_joint_cov()), the option
# in force while the fit forms its sdreport and while every reader runs.
# Prints, per reader, the largest relative difference (or an error).
a <- commandArgs(trailingOnly = TRUE)
design <- a[1]
seed <- if (length(a) >= 2L) as.integer(a[2]) else 1L
.libPaths(c("C:/Users/adf44/source/r/wt-sparsese-lib",
            "C:/Users/adf44/source/r/rellib-r22",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("lib", dirname(find.package("frmtmb")), "design", design, "seed", seed,
    "\n")
warns <- character()
quiet <- function(expr) {
  withCallingHandlers(expr, warning = function(w) {
    m <- conditionMessage(w)
    if (!grepl("a product such as a' V a", m, fixed = TRUE)) {
      warns <<- c(warns, substr(m, 1, 160))
    }
    invokeRestart("muffleWarning")
  }, message = function(m) invokeRestart("muffleMessage"))
}
set.seed(seed)
ng <- 30L
n <- 300L
d <- data.frame(x = stats::runif(n, 1, 2),
                g = factor(rep(seq_len(ng), length.out = n)),
                h = factor(sample(rep(seq_len(12), length.out = n))),
                f = factor(rep(c("a", "b"), each = 1, length.out = n)))
for (k in 1:7) d[[paste0("x", k)]] <- d$x^k
ug <- stats::rnorm(ng, 0, 0.5)
sg <- stats::rnorm(ng, 0, 0.4)
uh <- stats::rnorm(12, 0, 0.3)
eta <- sin(3 * d$x) + ug[d$g] + 0.3 * (d$f == "b")
deg5 <- "x1 + x2 + x3 + x4 + x5"
deg7 <- "x1 + x2 + x3 + x4 + x5 + x6 + x7"
args <- list(REML = TRUE)
fo <- switch(design,
  base = {
    d$y <- eta + stats::rnorm(n, 0, 0.3)
    stats::as.formula(paste("y ~ f +", deg7, "+ (1 | g)"))
  },
  slope = {
    d$y <- eta + sg[d$g] * d$x + stats::rnorm(n, 0, 0.3)
    stats::as.formula(paste("y ~", deg5, "+ (1 + x1 | g)"))
  },
  crossed = {
    d$y <- eta + uh[d$h] + stats::rnorm(n, 0, 0.3)
    stats::as.formula(paste("y ~", deg5, "+ (1 | g) + (1 | h)"))
  },
  poisson_reml = {
    d$y <- stats::rpois(n, exp(0.5 + 0.5 * eta))
    args$family <- poisson()
    stats::as.formula(paste("y ~", deg5, "+ (1 | g)"))
  },
  poisson_profile = {
    d$y <- stats::rpois(n, exp(0.5 + 0.5 * eta))
    args <- list(family = poisson(), REML = FALSE,
                 control = frmtmb_control(profile = TRUE))
    stats::as.formula(paste("y ~", deg5, "+ (1 | g)"))
  },
  gauss_profile = {
    d$y <- eta + stats::rnorm(n, 0, 0.3)
    args <- list(REML = FALSE, control = frmtmb_control(profile = TRUE))
    stats::as.formula(paste("y ~", deg7, "+ (1 | g)"))
  },
  cumulative = {
    z <- 0.8 * eta + stats::rlogis(n)
    d$y <- factor(cut(z, c(-Inf, 0, 0.8, Inf), labels = c("lo", "mid", "hi")),
                  ordered = TRUE)
    args$family <- cumulative()
    stats::as.formula(paste("y ~", deg5, "+ (1 | g)"))
  },
  dpar = {
    d$y <- eta + stats::rnorm(n, 0, exp(-1.2 + 0.5 * d$x))
    bf(stats::as.formula(paste("y ~", deg5, "+ (1 | g)")), sigma ~ x1)
  },
  stop("unknown design"))
lims <- c(dense = 1e6, sparse = 20)
on <- function(arm, expr) {
  op <- options(frmtmb.joint_root_max = lims[[arm]])
  on.exit(options(op))
  expr
}
fits <- list()
for (arm in names(lims)) {
  fits[[arm]] <- on(arm, {
    f <- tryCatch(quiet(do.call(frm, c(list(fo, data = d), args))),
                  error = function(e) e)
    if (!inherits(f, "error")) invisible(quiet(frmtmb:::sdr_of(f)))
    f
  })
}
if (inherits(fits$dense, "error")) {
  cat("fit error:", conditionMessage(fits$dense), "\n")
  quit(save = "no")
}
fd <- fits$dense
fs <- fits$sparse
sd_ <- frmtmb:::sdr_of(fs)
cat("whitened:", !is.null(fs$beta_whiten),
    " joint parameters:", nrow(sd_$jointPrecision),
    " lost:", paste(names(sd_$se_lost), sd_$se_lost),
    " same estimates:", identical(fd$opt$par, fs$opt$par), "\n")
cat("rows attr dense:", length(attr(frmtmb:::sdr_of(fd)$joint_root, "rows")),
    " sparse:", length(attr(sd_$joint_root, "rows")),
    " route dense:", on("dense", class(frmtmb:::get_joint_cov(fd)$V)[1]),
    " sparse:", on("sparse", class(frmtmb:::get_joint_cov(fs)$V)[1]), "\n")
rel <- function(x, y) {
  x <- as.numeric(as.matrix(x))
  y <- as.numeric(as.matrix(y))
  if (length(x) != length(y)) return(NA_real_)
  if (all(is.na(x) == is.na(y)) && any(is.na(x))) {
    k <- !is.na(x)
    x <- x[k]
    y <- y[k]
  }
  max(abs(x - y) / pmax(abs(x), abs(y), 1e-300))
}
nd <- d[c(1:6, 151:156), ]
ndnew <- nd
ndnew$g <- factor(c(rep("new1", 6), rep("new2", 6)))
ndnew$h <- factor(rep(c("newh", "1"), 6))
readers <- list(
  linpred_noRE = function(f) {
    frm_linpred(f, newdata = nd, re_formula = NA, se.fit = TRUE)$se.fit
  },
  linpred_RE = function(f) {
    frm_linpred(f, newdata = nd, se.fit = TRUE)$se.fit
  },
  linpred_newlev = function(f) {
    frm_linpred(f, newdata = ndnew, se.fit = TRUE,
                allow_new_levels = TRUE)$se.fit
  },
  fitted_RE = function(f) fitted(f, newdata = nd)[, "Est.Error"],
  fitted_all = function(f) as.numeric(fitted(f, newdata = nd)),
  fitted_newlev = function(f) {
    fitted(f, newdata = ndnew, allow_new_levels = TRUE)[, "Est.Error"]
  },
  fitted_insample = function(f) fitted(f)[, "Est.Error"],
  vcov = function(f) vcov(f),
  vcov_full = function(f) vcov(f, full = TRUE),
  frm_joint_cov = function(f) frm_joint_cov(f)$V,
  ranef_condVar = function(f) {
    r <- ranef(f, condVar = TRUE)
    unlist(lapply(r, function(x) {
      cv <- attr(x, "condVar") %||% attr(x, "postVar")
      c(as.numeric(x), if (!is.null(cv)) as.numeric(cv))
    }))
  },
  hypothesis = function(f) {
    hypothesis(f, c("x1 + x2 = 0", "x5 = 0"))$hypothesis$Est.Error
  },
  lp_basis_V = function(f) {
    lb <- frm_lp_basis(f, newdata = nd)
    cat("      lp_basis: V", !is.null(lb$V), "V_root", !is.null(lb$V_root),
        "V_sparse", !is.null(lb$V_sparse), "\n")
    sqrt(rowSums((lb$A %*% lb$V) * lb$A))
  },
  lp_basis_root = function(f) {
    lb <- frm_lp_basis(f, newdata = nd)
    if (is.null(lb$V_root)) return(rep(NA_real_, nrow(nd)))
    sqrt(rowSums((lb$A %*% lb$V_root)^2) +
           rowSums((lb$A %*% lb$V_rest) * lb$A))
  },
  frm_curve = function(f) {
    cv <- frm_curve(f, newdata = nd, re_formula = NA, seed = 1)
    cat("      frm_curve cov_rel_error:",
        format(attr(cv, "check")$cov_rel_error, digits = 3), "\n")
    cv[[".se"]]
  },
  mfx_predictions = function(f) {
    marginaleffects::predictions(f, newdata = nd)$std.error
  },
  mfx_avg_predictions = function(f) {
    marginaleffects::avg_predictions(f)$std.error
  },
  mfx_avg_slopes = function(f) {
    marginaleffects::avg_slopes(f, variables = "x")$std.error
  },
  emmeans = function(f) {
    if (design != "base") return(NA_real_)
    as.data.frame(emmeans::emmeans(f, ~ f))$SE
  },
  draws_epred_var = function(f) {
    set.seed(11)
    dr <- fitted(f, newdata = nd, summary = FALSE, ndraws = 4000)
    s <- apply(dr, 2, stats::sd)
    an <- fitted(f, newdata = nd)[, "Est.Error"]
    cat(sprintf("      draw sd / analytic SE - 1: %s\n",
                paste(format(range(s / an - 1), digits = 3), collapse = " ")))
    s
  },
  predict_newlev_gauss = function(f) {
    set.seed(12)
    predict(f, newdata = ndnew, allow_new_levels = TRUE,
            sample_new_levels = "gaussian", ndraws = 4000)[, "Est.Error"]
  },
  predict_newlev_uncert = function(f) {
    set.seed(13)
    predict(f, newdata = ndnew, allow_new_levels = TRUE,
            sample_new_levels = "uncertainty", ndraws = 4000)[, "Est.Error"]
  }
)
for (nm in names(readers)) {
  warns <- character()
  vd <- tryCatch(on("dense", quiet(readers[[nm]](fd))),
                 error = function(e) structure(conditionMessage(e),
                                               class = "err"))
  wd <- warns
  warns <- character()
  vs <- tryCatch(on("sparse", quiet(readers[[nm]](fs))),
                 error = function(e) structure(conditionMessage(e),
                                               class = "err"))
  ws <- warns
  if (inherits(vd, "err") || inherits(vs, "err")) {
    cat(sprintf("%-22s dense: %s | sparse: %s\n", nm,
                if (inherits(vd, "err")) substr(vd, 1, 120) else "ok",
                if (inherits(vs, "err")) substr(vs, 1, 120) else "ok"))
  } else {
    cat(sprintf("%-22s n=%5d rel %.3e  NaN dense %d sparse %d\n", nm,
                length(as.numeric(as.matrix(vd))), rel(vd, vs),
                sum(is.nan(as.numeric(as.matrix(vd)))),
                sum(is.nan(as.numeric(as.matrix(vs))))))
  }
  if (length(wd) || length(ws)) {
    cat("      W dense:", paste(unique(wd), collapse = " | "), "\n")
    cat("      W sparse:", paste(unique(ws), collapse = " | "), "\n")
  }
}
