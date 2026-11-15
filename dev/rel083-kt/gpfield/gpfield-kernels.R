# Lane gpfield, item 1906: brms's four gp() kernels against brms 2.23.0's
# own code. For each kernel and form, frm()'s gaussian fit, then the model's
# marginal log-likelihood at frm()'s estimates computed by hand from
# brms's covariance (brms:::cov_gp()) or its spectral density
# (brms:::spd_gp() at brms::standata()'s slambda and basis Xgp), with
# brms's scaled inputs. The two must read one number. Also the kernel
# matrix itself against brms's at one theta. Seed set.seed(9100 + case).
source("C:/Users/adf44/source/r/frmtmb-wt-rel083/dev/rel083-kt/gpfield/gpfield-lib.R")
gf_libs(commandArgs(TRUE)[1])
suppressPackageStartupMessages(library(brms))
cat("brms", as.character(packageVersion("brms")), "\n")
ns <- asNamespace("frmtmb")
bns <- asNamespace("brms")
kernels <- c("exp_quad", "matern52", "matern32", "exponential")
mk <- function(case, n) {
  set.seed(9100 + case)
  x1 <- round(runif(n, 0, 10), 2); x2 <- round(runif(n, 0, 5), 2)
  y <- sin(x1) + 0.3 * x2 + 0.2 * rnorm(n)
  data.frame(x1, x2, y)
}
forms <- list(
  exact1 = function(cv) y ~ gp(x1, cov = cv),
  exact2iso = function(cv) y ~ gp(x1, x2, cov = cv),
  exact2ani = function(cv) y ~ gp(x1, x2, cov = cv, iso = FALSE),
  hsgp1 = function(cv) y ~ gp(x1, k = 15, cov = cv),
  hsgp2ani = function(cv) y ~ gp(x1, x2, k = 6, cov = cv, iso = FALSE))
case <- 0L
for (fn in names(forms)) for (cv in kernels) {
  case <- case + 1L
  d <- mk(case, 40L)
  fml <- forms[[fn]](cv)
  q <- gf_fit(frmtmb::bf(fml), data = d, family = gaussian())
  f <- q$fit
  if (inherits(f, "error")) {
    cat(sprintf("KERN %s %s ERROR %s\n", fn, cv, conditionMessage(f)))
    next
  }
  bk <- Filter(function(b) b$covstruct %in% c("gp", "hsgp"),
               f$frame$re_blocks)[[1]]
  th <- f$estimates$theta[bk$theta_idx]
  sdgp <- exp(th[1])
  sigma <- exp(f$estimates$betad[1])
  b0 <- f$estimates$beta[1]
  sdat <- brms::standata(brms::bf(fml), data = d, family = gaussian())
  # brms's lscale is on its scaled inputs
  if (bk$covstruct == "gp") {
    ls_b <- exp(th[-1]) / bk$gp_lscale_div
    Xg <- sdat$Xgp_1
    Kb <- bns$cov_gp(Xg, sdgp = sdgp, lscale = ls_b, cov = cv)
    V <- Kb[sdat$Jgp_1, sdat$Jgp_1] + diag(sigma^2, nrow(d))
    # the kernel matrix alone, against frmtmb's at the same theta
    pos <- unname(as.matrix(Xg))
    Kf <- ns$gp_corr(th, bk)
    diag(Kf) <- diag(Kf) - ns$gp_nugget
    # frmtmb sorts its positions; match brms's rows by coordinates
    key_b <- apply(round(pos * bk$gp_lscale_div, 8), 1, paste, collapse = "_")
    fpos <- f$frame$gp_info %||% NULL
    gi <- NULL
    for (lp in f$frame$linpreds) for (g in lp$gps %||% list()) gi <- g
    key_f <- apply(round(gi$positions, 8), 1, paste, collapse = "_")
    o <- match(key_b, key_f)
    kerr <- max(abs(sdgp^2 * Kf[o, o] - Kb)) / sdgp^2
  } else {
    ls_b <- exp(th[-1])
    lam <- sdat$slambda_1
    spd <- bns$spd_gp(lam, sdgp = sdgp,
                      lscale = matrix(ls_b, nrow = 1), cov = cv)
    Phi <- sdat$Xgp_1
    if (!is.null(sdat$Jgp_1)) Phi <- Phi[sdat$Jgp_1, , drop = FALSE]
    V <- Phi %*% diag(as.numeric(spd), ncol(Phi)) %*% t(Phi) +
      diag(sigma^2, nrow(d))
    sf <- as.numeric(ns$hsgp_sds(th, bk$aux_omega, bk$gp_iso, cv))^2
    kerr <- max(abs(sort(sf) - sort(as.numeric(spd)))) / max(as.numeric(spd))
  }
  ll_b <- mvtnorm::dmvnorm(d$y, rep(b0, nrow(d)), V, log = TRUE)
  ll_f <- as.numeric(logLik(f))
  cat(sprintf(paste0("KERN %-9s %-11s frm ll %.8f brms-kernel ll %.8f ",
                     "diff %.2e | kernel rel diff %.2e | conv %d nwarn %d %s\n"),
              fn, cv, ll_f, ll_b, ll_f - ll_b, kerr, f$opt$convergence,
              length(q$warn), paste(unique(sub(":.*", "", q$warn)),
                                    collapse = ";")))
  flush.console()
}
cat("DONE\n")
