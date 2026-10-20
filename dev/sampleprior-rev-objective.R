# Reviewer of lane sampleprior: the log prior frm_sample() ACTUALLY adds
# in the objective it hands tmbstan, against brms 2.23.0's Stan lprior
# (read off brms::stancode(), dev/sampleprior-rev-log/dump*.txt) written
# out by hand here, carried to frmtmb's sampling coordinates by a
# NUMERICAL Jacobian of the natural parameters on those coordinates.
# Compared as differences between 6 points (point 1 the reference), so
# that normalizing constants (brms's lccdf terms) cancel.
#   Rscript dev/sampleprior-rev-objective.R [case ...]
source("C:/Users/adf44/source/r/frmtmb-wt-sampleprior/dev/sampleprior-rev-common.R")
source("C:/Users/adf44/source/r/frmtmb-wt-sampleprior/dev/sampleprior-rev-cases.R")
ns <- asNamespace("frmtmb.sample")
core <- asNamespace("frmtmb")
PATCH_MV <- identical(Sys.getenv("REV_PATCH_MV"), "1")
if (PATCH_MV) {
  # the reviewer's candidate fix for the shared nu: exempt it from the
  # "no resp" refusal (core's resp_missing_refusal())
  orig <- core$resp_missing_refusal
  utils::assignInNamespace("resp_missing_refusal", function(spec, frame, s) {
    if (!is.null(spec[["rescor_nu"]]) &&
        (identical(s$class, "nu") || identical(s$dpar, "nu"))) {
      return(NULL)
    }
    orig(spec, frame, s)
  }, "frmtmb")
}

lt <- function(v, df, mu, s) stats::dt((v - mu) / s, df, log = TRUE) - log(s)
linvg <- function(v, a, b) a * log(b) - lgamma(a) - (a + 1) * log(v) - b / v
mx <- mean(x)
# log |det d nat / d p| by central differences, nat a function of p
ljac <- function(fnat, p, h = 1e-5) {
  k <- length(fnat(p))
  J <- matrix(0, k, length(p))
  for (j in seq_along(p)) {
    a <- p; a[j] <- a[j] + h
    b <- p; b[j] <- b[j] - h
    J[, j] <- (fnat(a) - fnat(b)) / (2 * h)
  }
  as.numeric(determinant(J, logarithm = TRUE)$modulus)
}
# a scalar natural dpar: its log density plus the log Jacobian of
# `inv` (the reviewer's own inverse link) at `e`
nat1 <- function(e, inv, logd) logd(inv(e)) + ljac(function(z) inv(z), e)
plg <- stats::plogis
minY <- function(d, v) min(d[[v]])

B <- list(
  gauss_wide = function(p, d) lt(p[1] + p[2] * mx, 3, 33.4, 9.9) +
    nat1(p[3], exp, function(s) lt(s, 3, 0, 9.9)),
  sigma_x = function(p, d) lt(p[1] + p[2] * mx, 3, 6, 8.8) +
    lt(p[3] + p[4] * mx, 3, 0, 4.2),
  lognormal_wide = function(p, d) lt(p[1] + p[2] * mx, 3, 1.1, 4.7) +
    nat1(p[3], exp, function(s) lt(s, 3, 0, 4.7)),
  negbin = function(p, d) lt(p[1] + p[2] * mx, 3, 1.1, 2.5) +
    nat1(p[3], exp, function(s) linvg(s, 0.4, 0.3)),
  beta_phi_x = function(p, d) lt(p[1] + p[2] * mx, 3, 0, 2.5) +
    lt(p[3] + p[4] * mx, 3, 0, 2.5),
  student = function(p, d) lt(p[1] + p[2] * mx, 3, 4.1, 2.5) +
    nat1(p[3], exp, function(s) lt(s, 3, 0, 2.5)) +
    nat1(p[4], function(e) 1 + exp(e),
         function(v) stats::dgamma(v, 2, 0.1, log = TRUE)),
  zip_zi_x = function(p, d) lt(p[1] + p[2] * mx, 3, 0, 2.5) +
    stats::dlogis(p[3] + p[4] * mx, 0, 1, log = TRUE),
  hurdle_ln = function(p, d) lt(p[1] + p[2] * mx, 3, 0.2, 2.5) +
    nat1(p[3], exp, function(s) lt(s, 3, 0, 2.5)) +
    nat1(p[4], plg, function(h) stats::dbeta(h, 1, 1, log = TRUE)),
  hurdle_ln_hu1 = function(p, d) lt(p[1] + p[2] * mx, 3, 0.2, 2.5) +
    nat1(p[3], exp, function(s) lt(s, 3, 0, 2.5)) +
    stats::dlogis(p[4], 0, 1, log = TRUE),
  gev_xi1 = function(p, d) lt(p[1] + p[2] * mx, 3, 3.9, 2.5) +
    nat1(p[3], exp, function(s) lt(s, 3, 0, 2.5)) +
    stats::dnorm(p[4], 0, 4, log = TRUE),
  cum_thres_gr = function(p, d) {
    out <- 0
    for (gi in 0:1) {
      r <- p[1 + gi * 3 + 1:3]
      out <- out + sum(lt(cumsum(c(r[1], exp(r[2:3]))), 3, 0, 2.5)) +
        ljac(function(z) cumsum(c(z[1], exp(z[2:3]))), r)
    }
    out
  },
  cum = function(p, d) {
    r <- p[2:4]
    tau <- function(z) cumsum(c(z[1], exp(z[2:3])))
    sum(lt(tau(r) - mx * p[1], 3, 0, 2.5)) + ljac(tau, r)
  },
  mixture = function(p, d) lt(p[1] + p[2] * mx, 3, 5.3, 2.5) +
    lt(p[3] + p[4] * mx, 3, 5.3, 2.5) +
    nat1(p[5], exp, function(s) lt(s, 3, 0, 2.5)) +
    nat1(p[6], exp, function(s) lt(s, 3, 0, 2.5)) +
    # dirichlet(1) is flat on the simplex; the coordinate is the log
    # ratio of the second weight to the first
    ljac(function(z) plg(z), p[7]),
  poisson_sqrt = function(p, d) lt(p[1] + p[2] * mx, 3, 2.4, 2.5),
  beta_phi_id = function(p, d) lt(p[1] + p[2] * mx, 3, 0, 2.5) +
    stats::dgamma(p[3], 0.01, 0.01, log = TRUE),
  sigma_id_x = function(p, d) lt(p[1] + p[2] * mx, 3, 33.4, 9.9) +
    lt(p[3] + p[4] * mx, 3, 0, 9.9),
  shifted_ln = function(p, d) lt(p[1] + p[2] * mx, 3, -0.1, 2.5) +
    nat1(p[3], exp, function(s) lt(s, 3, 0, 2.5)) +
    nat1(p[4], exp, function(v) {
      if (v > 0 && v < minY(d, "y")) -log(minY(d, "y")) else -Inf
    }),
  wiener = function(p, d) {
    m <- minY(d, "rt")
    lt(p[1] + p[2] * mx, 3, 0.7, 2.5) +
      nat1(p[3], exp, function(v) stats::dgamma(v, 1, 1, log = TRUE)) +
      nat1(p[4], function(e) m * plg(e), function(v) -log(m)) +
      nat1(p[5], plg, function(v) stats::dbeta(v, 1, 1, log = TRUE))
  }
)

# cases whose covariance parameters go through core's own naming map
# (hyp_env_vals()): the natural quantities brms priors, as a function
# of the internal vector
B$re_sigma <- function(p, d, fit) {
  vo <- core$hyp_vals_only(fit)
  th <- p[names(p) == "theta"]
  nat <- function(t) {
    v <- vo$vals
    v[vo$comp == "theta"] <- t
    ev <- core$hyp_env_vals(fit, v, vo$comp)
    unlist(ev[c("sd_g__Intercept", "sd_g__x", "cor_g__Intercept__x",
                "sd_g__sigma_Intercept")])
  }
  q <- nat(th)
  lt(p[1] + p[2] * mx, 3, 33.4, 9.9) + lt(p[3], 3, 0, 2.5) +
    sum(lt(q[c(1, 2, 4)], 3, 0, 9.9)) + log(0.5) + ljac(nat, th)
}
B$mv_rescor <- function(p, d, fit) {
  vo <- core$hyp_vals_only(fit)
  bd <- p[names(p) == "betad"]
  tr <- p[names(p) == "thetar"]
  nat <- function(t) {
    v <- vo$vals
    v[vo$comp == "thetar"] <- t
    unlist(core$hyp_env_vals(fit, v, vo$comp)["rescor__y1__y2"])
  }
  lt(p[1] + p[2] * mx, 3, 4, 2.5) + lt(p[3] + p[4] * mx, 3, -0.9, 2.5) +
    nat1(bd[1], exp, function(s) lt(s, 3, 0, 2.5)) +
    nat1(bd[3], exp, function(s) lt(s, 3, 0, 2.5)) +
    nat1(bd[2], function(e) 1 + exp(e),
         function(v) stats::dgamma(v, 2, 0.1, log = TRUE)) +
    log(0.5) + ljac(nat, tr)
}

which <- commandArgs(TRUE)
if (!length(which)) which <- names(B)
set.seed(7)
res <- list()
for (nm in which) {
  cs <- CASES[[nm]]
  f <- eval(cs$frm)
  out <- tryCatch({
    # centered on both tapes: the non-centered retape reads z where the
    # flat tape reads b
    od <- cap_obj(f[[1]], f[[2]], f[[3]], reparameterize = FALSE)
    of <- cap_obj(f[[1]], f[[2]], f[[3]], prior = "flat",
                  reparameterize = FALSE)
    if (!identical(names(od$env$par), names(of$env$par))) {
      stop("default and flat objectives have different parameters")
    }
    fit <- suppressMessages(suppressWarnings(ns$sample_assemble(
      f[[1]], f[[2]], f[[3]], data2 = list(), start = NULL,
      fit_control = frmtmb_control(), na.action = stats::na.omit,
      REML = FALSE)))
    p0 <- od$env$par
    # points inside every support: around a neutral start, the ndt of
    # the shifted lognormal kept below min(y)
    P <- lapply(1:6, function(k) {
      q <- p0 + stats::rnorm(length(p0), 0, if (nm == "gev_xi1") 0.1 else 0.4)
      if (nm == "shifted_ln") q[4] <- log(minY(f[[2]], "y")) - 0.1 * k
      if (nm %in% c("beta_phi_id")) q[3] <- abs(p0[3]) * (0.5 + k / 4)
      if (nm %in% c("sigma_id_x")) {
        q[3] <- p0[3] - mx * q[4] + k   # a positive sigma at the mean
      }
      q
    })
    fl <- vapply(P, function(q) logprior_at(od, of, q), 1)
    bl <- vapply(P, function(q) {
      if (length(formals(B[[nm]])) == 3L) B[[nm]](q, f[[2]], fit) else
        B[[nm]](q, f[[2]])
    }, 1)
    dfl <- fl - fl[1]
    dbl <- bl - bl[1]
    sprintf("max |d frm - d brms| %.3g | d frm %s | d brms %s",
            max(abs(dfl - dbl)), paste(sprintf("%.4f", dfl[-1]), collapse = " "),
            paste(sprintf("%.4f", dbl[-1]), collapse = " "))
  }, error = function(e) paste("ERROR:", conditionMessage(e)))
  cat(sprintf("%-15s %s\n", nm, out))
}
