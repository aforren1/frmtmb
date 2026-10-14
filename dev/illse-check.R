# Lane illse: covariance and prediction accuracy on gaussian designs,
# against the exact covariance at the fit's own point (256-bit,
# dev/illse-common.R) and against lm() / lmer() at their own maxima.
# Per fit: code, warning tags, lost set, Newton steps, the shortfall
# to the exact maximum (lm()'s logLik; lmer()'s ML logLik with (1 | g)),
# the error of every entry of vcov(full = TRUE) (cor: over
# sqrt(Vii Vjj); rel: relative, entries with |correlation| >= 1e-3;
# diag: the variances), the SE ratio to lm()/lmer() over the
# coefficients, and fitted(newdata, re_formula = NA)$Est.Error over the
# exact one at x = 1, 1.25, 1.5, 1.75, 2.
#   ILLSE_ARM=base|lane Rscript dev/illse-check.R [seeds] [poly|other]
source("C:/Users/adf44/source/r/frmtmb-wt-illse/dev/illse-common.R")
source("C:/Users/adf44/source/r/frmtmb-wt-illse/dev/illse-designs.R")
args <- commandArgs(TRUE)
seeds <- if (length(args)) eval(parse(text = args[1])) else 1:3
suppressMessages(library(lme4))
tag <- function(w) {
  if (!length(w)) return("-")
  k <- c(illcond = "ill-conditioned", short = "stopped short",
         lost = "Standard errors are not available",
         grad = "gradient", flat = "flat along", conv = "false conv",
         prod = "product such as")
  paste(names(k)[vapply(k, function(s) any(grepl(s, w)), NA)],
        collapse = ",")
}
xr <- c(1, 1.25, 1.5, 1.75, 2)
fam_arg <- if (length(args) >= 2) args[2] else "poly"
designs <- list()
if (fam_arg == "poly") {
  for (deg in 5:8) for (re in c(FALSE, TRUE)) {
    designs[[length(designs) + 1]] <- des_poly(deg, re)
  }
  designs[[length(designs) + 1]] <- des_b5(FALSE)
  designs[[length(designs) + 1]] <- des_b5(TRUE)
} else {
  # covariates far from zero, near-collinear pairs, mixed scales
  for (re in c(FALSE, TRUE)) {
    for (off in c(1e5, 1e6, 3e6)) {
      designs[[length(designs) + 1]] <- des_offset(off, re)
    }
    for (dl in c(1e-5, 1e-6, 1e-7)) {
      designs[[length(designs) + 1]] <- des_collin(dl, re)
    }
    designs[[length(designs) + 1]] <- des_scales(re)
  }
}
for (mk in designs) for (seed in seeds) {
  ds <- mk(seed)
  re <- grepl("(1|g)", ds$lab, fixed = TRUE)
  r <- cap(frm(ds$formula, data = ds$data))
  if (inherits(r$v, "illse_err")) {
    cat(sprintf("%-30s s%d ERROR %s\n", ds$lab, seed, r$v)); next
  }
  f <- r$v
  p <- f$opt$par
  X <- as.matrix(f$frame$linpreds[[1]]$X)
  y <- ds$data$y
  Ve <- if (re) lmm_exact_cov(X, y, ds$data$g, p) else {
    gauss_exact_cov(X, y, p[names(p) == "beta"], p[names(p) == "betad"])
  }
  V <- suppressWarnings(vcov(f, full = TRUE))
  ee <- if (all(is.finite(V))) mp_cov_err(V, Ve) else
    c(cor = NaN, rel = NaN, diag = NaN)
  # reference fits at their own maxima
  if (re) {
    rf <- tryCatch(suppressWarnings(lmer(ds$formula, data = ds$data,
                                         REML = FALSE)),
                   error = function(e) NULL)
    ll_ref <- if (is.null(rf)) NA_real_ else as.numeric(logLik(rf))
    se_ref <- if (is.null(rf)) NA_real_ else sqrt(diag(as.matrix(vcov(rf))))
  } else {
    rf <- lm(ds$formula, data = ds$data)
    ll_ref <- as.numeric(logLik(rf))
    n <- nrow(ds$data)
    se_ref <- sqrt(diag(vcov(rf)) * (n - length(coef(rf))) / n)
  }
  short <- ll_ref - as.numeric(logLik(f))
  se_f <- sqrt(diag(V))[names(p) == "beta"]
  ser <- range(se_f / se_ref)
  # prediction rows at the population level
  nd <- if (fam_arg == "poly") {
    v <- ds$data[rep(1, length(xr)), , drop = FALSE]
    v$x <- xr
    v
  } else ds$data[seq_along(xr), , drop = FALSE]
  A <- model.matrix(delete.response(terms(lm(reformulas::nobars(ds$formula), data = ds$data,
                                             method = "model.frame"))),
                    nd)
  A <- A[, colnames(X), drop = FALSE]
  pb <- which(names(p) == "beta")
  Vbe <- Ve[pb, pb]
  Am <- mpfrArray(A, mp_prec, dim = dim(A))
  se_ex <- sqrt(asNumeric(rowSums((Am %*% Vbe) * Am)))
  fr <- cap(fitted(f, newdata = nd, re_formula = NA))
  pr <- if (inherits(fr$v, "illse_err")) rep(NaN, length(xr)) else
    fr$v[, "Est.Error"] / se_ex
  cat(sprintf(paste0("%-30s s%d code %d | warn %-18s | lost %-10s | ",
                     "steps %s | short %9.2e | cov cor %.1e rel %.1e ",
                     "diag %.1e | SE/ref %.4f-%.4f | fitted SE/exact %s",
                     " %s\n"),
              ds$lab, seed, f$opt$convergence, tag(r$w), lost_of(f),
              if (is.null(f$opt$polish)) 0 else f$opt$polish[["steps"]],
              short, ee[["cor"]], ee[["rel"]], ee[["diag"]],
              ser[1], ser[2],
              paste(sprintf("%.5f", pr), collapse = " "), tag(fr$w)))
}
