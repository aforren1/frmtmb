# Review of lane gpfield, item 1906: kernels against brms 2.23.0, cases the
# lane did not run: by = factor groups (exact and approximate), exact
# isotropic 2-D exponential, anisotropic matern52 exact (product kernel),
# approximate isotropic 2-D matern32, matern12 alias. The marginal at
# frm()'s estimates from brms's own R kernels (cov_gp(), spd_gp()) on
# standata()'s inputs. Also gp_field_dtheta() for Matern against a central
# difference of the field at fixed z.
source("C:/Users/adf44/source/r/frmtmb-wt-gpfield/dev/gpfield-lib.R")
gf_libs("lane")
suppressPackageStartupMessages(library(brms))
cat("brms", as.character(packageVersion("brms")), "\n")
ns <- asNamespace("frmtmb"); bns <- asNamespace("brms")
mk <- function(seed, n) {
  set.seed(seed)
  x1 <- round(runif(n, 0, 10), 2); x2 <- round(runif(n, 0, 5), 2)
  g <- factor(sample(c("a", "b", "c"), n, TRUE))
  y <- sin(x1) * (1 + (g == "b")) + 0.3 * x2 + 0.2 * rnorm(n)
  data.frame(x1, x2, g, y)
}
cases <- list(
  by_exact_m32 = 'y ~ g + gp(x1, by = g, cov = "matern32")',
  by_exact_exp = 'y ~ g + gp(x1, by = g, cov = "exponential")',
  by_hsgp_m52 = 'y ~ g + gp(x1, by = g, k = 10, cov = "matern52")',
  iso2_exact_expo = 'y ~ gp(x1, x2, cov = "exponential")',
  ani2_exact_m52 = 'y ~ gp(x1, x2, iso = FALSE, cov = "matern52")',
  iso2_hsgp_m32 = 'y ~ gp(x1, x2, k = 6, cov = "matern32")',
  ani2_hsgp_expo = 'y ~ gp(x1, x2, k = 6, iso = FALSE, cov = "exponential")',
  m12_alias = 'y ~ gp(x1, cov = "matern12")',
  ani2_c_hsgp_m52 = 'y ~ gp(x1, x2, k = 7, c = c(1.5, 2), iso = FALSE, cov = "matern52")')
i <- 0L
for (nm in names(cases)) {
  i <- i + 1L
  d <- mk(9400 + i, 45L)
  fml <- as.formula(cases[[nm]])
  q <- gf_fit(frmtmb::bf(fml), data = d, family = gaussian())
  f <- q$fit
  if (inherits(f, "error")) {
    cat(sprintf("KERNR %s ERROR %s\n", nm, conditionMessage(f))); next
  }
  sdat <- brms::standata(brms::bf(fml), data = d, family = gaussian())
  cv <- sub(".*cov = \"([a-z0-9_]+)\".*", "\\1", cases[[nm]])
  if (cv == "matern12") cv <- "exponential"
  bks <- Filter(function(b) b$covstruct %in% c("gp", "hsgp"),
                f$frame$re_blocks)
  sigma <- exp(f$estimates$betad[1])
  X <- model.matrix(if (grepl("~ g +", cases[[nm]], fixed = TRUE)) ~ g else ~ 1, d)
  mu <- as.vector(X %*% f$estimates$beta[seq_len(ncol(X))])
  V <- diag(sigma^2, nrow(d))
  # brms names per by level: Xgp_1_<k>, Igp_1_<k>, slambda_1_<k>
  sfx <- if (any(grepl("^Igp_1_", names(sdat)))) {
    sub("^Igp_1_", "", grep("^Igp_1_", names(sdat), value = TRUE))
  } else ""
  if (length(sfx) != length(bks)) {
    cat("KERNR", nm, "level count mismatch", length(sfx), length(bks), "\n")
  }
  # frmtmb's blocks per level in brms's level order
  for (k in seq_along(bks)) {
    bk <- bks[[k]]
    th <- f$estimates$theta[bk$theta_idx]
    s <- if (nzchar(sfx[1])) paste0("_", sfx[k]) else ""
    Xg <- sdat[[paste0("Xgp_1", s)]]
    I <- if (nzchar(s)) sdat[[paste0("Igp_1", s)]] else seq_len(nrow(d))
    J <- sdat[[paste0("Jgp_1", s)]]
    if (bk$covstruct == "gp") {
      Kb <- bns$cov_gp(Xg, sdgp = exp(th[1]),
                       lscale = exp(th[-1]) / bk$gp_lscale_div, cov = cv)
      if (is.null(J)) J <- seq_len(nrow(Kb))
      Kr <- Kb[J, J]
    } else {
      lam <- sdat[[paste0("slambda_1", s)]]
      spd <- bns$spd_gp(lam, sdgp = exp(th[1]),
                        lscale = matrix(exp(th[-1]), nrow = 1), cov = cv)
      Phi <- Xg
      if (!is.null(J)) Phi <- Phi[J, , drop = FALSE]
      Kr <- Phi %*% (as.numeric(spd) * t(Phi))
    }
    V[I, I] <- V[I, I] + Kr
  }
  ll_b <- mvtnorm::dmvnorm(d$y, mu, V, log = TRUE)
  ll_f <- as.numeric(logLik(f))
  cat(sprintf("KERNR %-16s %-11s blocks %d frm ll %.8f brms-kernel ll %.8f diff %.2e conv %d warn [%s]\n",
              nm, cv, length(bks), ll_f, ll_b, ll_f - ll_b,
              f$opt$convergence, paste(unique(sub(":.*", "", q$warn)),
                                       collapse = ";")))
  # derivative of an exact field at fixed z: analytic against central diff
  bk <- bks[[1]]
  if (bk$covstruct == "gp") {
    th <- f$estimates$theta[bk$theta_idx]
    z <- f$estimates$b[bk$b_idx]
    fld <- function(t) {
      L <- ns$gp_factor_ordered(t, bk)
      as.vector(exp(t[1]) * (L %*% z))
    }
    an <- ns$gp_field_dtheta(z, th, bk)
    an <- exp(th[1]) * (ns$gp_factor_ordered(th, bk) %*% an)
    err <- 0
    for (j in seq_along(th)[-1]) {
      h <- 1e-5
      tp <- th; tp[j] <- tp[j] + h; tm <- th; tm[j] <- tm[j] - h
      fd <- (fld(tp) - fld(tm)) / (2 * h)
      err <- max(err, max(abs(fd - an[, j])) / max(abs(fd)))
    }
    cat(sprintf("DTHETA %-16s rel err of d field / d log ls %.2e\n", nm, err))
  }
}
