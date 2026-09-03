# dev/pdi-repro.R -- prior-dropping investigation, discriminating measurement.
#
# Does the prior term reach the SAMPLED objective? Build the augmented
# objective frm_sample() hands to tmbstan and evaluate its fn directly, so
# taping is separated from sampling. Run the identical script on the host and
# in the pkgcheck container and diff the numbers.
#
#   Rscript dev/pdi-repro.R
#
# Reads nothing, writes nothing; everything goes to stdout.

options(warn = 1, digits = 15)
suppressPackageStartupMessages(library(frmtmb))

cat("##### PDI repro\n")
cat("R          : ", R.version.string, " ", R.version$platform, "\n", sep = "")
for (p in c("frmtmb", "RTMB", "TMB", "Matrix", "rstan", "StanHeaders",
            "tmbstan", "numDeriv")) {
  v <- tryCatch(as.character(utils::packageVersion(p)), error = function(e) "-")
  cat(sprintf("%-12s: %s\n", p, v))
}
cat("BLAS       : ", sessionInfo()$BLAS %||% "?", "\n", sep = "")
cat("LAPACK     : ", sessionInfo()$LAPACK %||% "?", "\n", sep = "")

set.seed(1)
n <- 20L
d <- data.frame(x = rnorm(n), g = factor(rep(1:5, each = 4)))
d$y <- 1 + 0.8 * d$x + rep(rnorm(5, 0, 0.7), each = 4) + rnorm(n, 0, 0.5)

`%||%` <- function(a, b) if (is.null(a)) b else a

# ---------------------------------------------------------------- probe --
# obj_aug$fn(p) - obj_plain$fn(p) IS the prior term. Evaluated at two p that
# differ only in the priored coordinate it must CHANGE. Three outcomes:
#   delta == 0 everywhere      -> the prior never entered the objective
#   delta constant, non-zero   -> taped to a constant: invisible to NUTS,
#                                 because an additive constant in the log
#                                 density cancels in every Metropolis ratio
#                                 and has zero gradient. This is the shape
#                                 that makes a chain BYTE-IDENTICAL to flat.
#   delta varies as expected   -> the prior is live
probe <- function(label, fit, priors, comp, coefname) {
  cat("\n### ", label, "\n", sep = "")
  ri <- frmtmb:::resolve_prior_input(fit, priors)
  cat("entries: ", length(ri$entries), "\n", sep = "")
  for (e in ri$entries) {
    cat("  comp=", e$comp, " idx=", paste(e$idx, collapse = ","),
        " kind=", e$dist$kind, " scale=", e$scale, "\n", sep = "")
  }
  if (!length(ri$entries)) return(invisible(NULL))

  obj0 <- fit$obj
  obj1 <- frmtmb:::prior_augmented_obj(fit, ri$entries)
  cat("obj1 identical to fit$obj? ", identical(obj1, obj0), "\n", sep = "")

  p <- obj1$par
  cat("par names: ", paste(names(p), collapse = " "), "\n", sep = "")
  k <- which(names(p) == comp)
  # the priored coordinate inside the flat outer vector
  ent <- ri$entries[[1]]
  kk <- k[ent$idx[1]]
  cat("probing flat position ", kk, " (", names(p)[kk], ")\n", sep = "")

  vals <- c(-1, 0, 0.5, 2)
  cat(sprintf("%10s %22s %22s %22s\n", "value", "fn_plain", "fn_aug",
              "aug - plain"))
  deltas <- numeric(length(vals))
  for (i in seq_along(vals)) {
    q <- p
    q[kk] <- vals[i]
    f0 <- obj0$fn(q)
    f1 <- obj1$fn(q)
    deltas[i] <- f1 - f0
    cat(sprintf("%10.4f %22.12f %22.12f %22.12f\n", vals[i], f0, f1,
                f1 - f0))
  }
  cat("range(delta) = ", diff(range(deltas)), "\n", sep = "")

  # the same negative log prior, evaluated in plain R OUTSIDE any tape
  nlp <- frmtmb:::neg_log_prior_fn(ri$entries)
  est <- fit$estimates
  cat("untaped R-side nlp:\n")
  for (i in seq_along(vals)) {
    e2 <- est
    e2[[ent$comp]][ent$idx[1]] <- vals[i]
    cat(sprintf("%10.4f %22.12f\n", vals[i], nlp(e2)))
  }

  # gradients: a constant-taped prior leaves the gradient untouched
  q <- p
  q[kk] <- 0.5
  g0 <- obj0$gr(q)
  g1 <- obj1$gr(q)
  cat("gr_plain: ", paste(sprintf("%.8f", g0), collapse = " "), "\n", sep = "")
  cat("gr_aug  : ", paste(sprintf("%.8f", g1), collapse = " "), "\n", sep = "")
  cat("gr diff : ", paste(sprintf("%.8f", g1 - g0), collapse = " "), "\n",
      sep = "")
  invisible(NULL)
}

# (A) fixed effects only: no random block, the plainest possible tape
fitA <- frm(y ~ x, data = d, family = gaussian())
probe("A. gaussian y ~ x, prior_normal(0, 0.01) on x", fitA,
      list(x = prior_normal(0, 0.01)), "beta", "x")

# (B) a random intercept, prior on the sd: the shape test-reparam.R:627 uses
fitB <- frm(y ~ x + (1 | g), data = d, family = gaussian())
probe("B. y ~ x + (1|g), prior_normal(0, 0.01) on theta_1", fitB,
      list(theta_1 = prior_normal(0, 0.01)), "theta", "theta_1")

# (C) the FORMULA route's own resolution, which is what frm_sample() runs
cat("\n### C. sample_resolve_priors() default priors on B\n")
rpc <- frmtmb:::sample_resolve_priors(fitB, NULL)
cat("entries: ", length(rpc$ri$entries), "\n", sep = "")
for (e in rpc$ri$entries) {
  cat("  comp=", e$comp, " idx=", paste(e$idx, collapse = ","),
      " kind=", e$dist$kind, " scale=", e$scale,
      " loc=", e$dist$location %||% NA, " sc=", e$dist$scale %||% NA,
      " df=", e$dist$df %||% NA, "\n", sep = "")
}
if (length(rpc$ri$entries)) {
  objC <- frmtmb:::prior_augmented_obj(fitB, rpc$ri$entries)
  p <- objC$par
  cat("par names: ", paste(names(p), collapse = " "), "\n", sep = "")
  cat(sprintf("%10s %22s %22s %22s\n", "theta1", "fn_plain", "fn_aug",
              "aug - plain"))
  kk <- which(names(p) == "theta")[1]
  ds <- numeric(0)
  for (v in c(-1, 0, 0.5, 2)) {
    q <- p; q[kk] <- v
    f0 <- fitB$obj$fn(q); f1 <- objC$fn(q)
    ds <- c(ds, f1 - f0)
    cat(sprintf("%10.4f %22.12f %22.12f %22.12f\n", v, f0, f1, f1 - f0))
  }
  cat("range(delta) = ", diff(range(ds)), "\n", sep = "")
}

# ------------------------------------------------------ RTMB tape sanity --
# Does RTMB tape an added scalar term at all here? Minimal, frmtmb-free.
cat("\n### D. bare RTMB tape sanity\n")
f <- function(pars) {
  s <- sum(pars$a^2)
  s - sum(RTMB::dnorm(pars$a, 0, 0.01, log = TRUE))
}
o <- RTMB::MakeADFun(f, list(a = c(0.1, 0.2)), silent = TRUE)
g <- function(pars) sum(pars$a^2)
o2 <- RTMB::MakeADFun(g, list(a = c(0.1, 0.2)), silent = TRUE)
for (v in c(0, 0.5, 2)) {
  q <- c(v, 0.2)
  cat(sprintf("a1=%5.2f  plain=%18.10f  withprior=%18.10f  delta=%18.10f  analytic=%18.10f\n",
              v, o2$fn(q), o$fn(q), o$fn(q) - o2$fn(q),
              -sum(dnorm(q, 0, 0.01, log = TRUE))))
}
cat("gr with prior: ", paste(sprintf("%.6f", o$gr(c(0.5, 0.2))),
                             collapse = " "), "\n", sep = "")

cat("\n##### done\n")
