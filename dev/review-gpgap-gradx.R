# Short check of fix round 1, item 1: the tape's gradient and Stan's disagree at
# the target (review-gpgap-grad.R). Which is right? The density in 300-bit arithmetic,
# central differences, against each. Same points; data sorted, seed 6601, n 15.
#
#   bash dev/review-gpgap-run.sh dev/review-gpgap-gradx.R gradx.txt
.libPaths(c("C:/Users/adf44/source/r/wt-gpgap-lib",
            "C:/Users/adf44/source/r/rellib-r16",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(frmtmb); library(testthat); library(Rmpfr)})
H <- new.env(parent = asNamespace("frmtmb"))
for (h in c("helper-brms.R", "helper-stan.R")) {
  sys.source(file.path("tests/testthat", h), envir = H)
}
set.seed(6601)
n <- 15
x <- round(runif(n, 0, 6), 3)
e0 <- rnorm(n)
d <- data.frame(x = x, y = sin(2 * x) + 1e-3 * e0)
d <- d[order(d$x), ]
bb <- brms::bf(y ~ gp(x))
pr <- H$brms_flat_prior(bb, data = d, family = gaussian())
code <- brms::make_stancode(bb, data = d, family = gaussian(), prior = pr)
sdat <- suppressMessages(brms::make_standata(bb, data = d, family = gaussian(), prior = pr))
rel <- gsub("cov[n, n] += 1e-12;", "cov[n, n] += 1e-12 * square(sdgp);", code, fixed = TRUE)
sf <- suppressMessages(rstan::sampling(H$brms_stan_model(rel), data = sdat, chains = 0))
fit <- suppressWarnings(frm(bf(y ~ gp(x)), family = gaussian(), data = d))
fr <- fit$frame
fr$ncp_blocks <- 1L
fr$re_blocks[[1]]$gp_nugget <- 1e-12
fn <- frmtmb::build_objective(fr)
tpl <- fr$par_template
tpl$b[] <- 0
tape <- RTMB::MakeADFun(fn, tpl, map = fr$map, silent = TRUE)
dmax <- sdat$dmax_1
prec <- 300
mp <- function(v) mpfr(v, prec)
xm <- mp(d$x)
ym <- mp(d$y)
# log density in exact arithmetic: u = (b0, log sd, log rho, z, log sigma)
lpx <- function(u) {
  b0 <- u[1]; sd <- exp(u[2]); rho <- exp(u[3]); z <- u[4:(3 + n)]; sg <- exp(u[4 + n])
  L <- vector("list", n)
  Lm <- mpfrArray(0, prec, dim = c(n, n))
  for (j in seq_len(n)) {
    for (i in j:n) {
      kij <- sd^2 * exp(-(xm[i] - xm[j])^2 / (2 * rho^2)) +
        (if (i == j) mp(1e-12) * sd^2 else mp(0))
      s <- kij
      if (j > 1) for (q in seq_len(j - 1)) s <- s - Lm[i, q] * Lm[j, q]
      Lm[i, j] <- if (i == j) sqrt(s) else s / Lm[j, j]
    }
  }
  f <- mp(rep(0, n))
  for (i in seq_len(n)) {
    t0 <- mp(0)
    for (q in seq_len(i)) t0 <- t0 + Lm[i, q] * z[q]
    f[i] <- t0
  }
  r <- ym - b0 - f
  -sum(z^2) / 2 - n * u[4 + n] - sum(r^2) / (2 * sg^2) - n * log(2 * mp(pi)) / 2 -
    n * log(2 * mp(pi)) / 2
}
set.seed(44)
for (i in 1:8) {
  sg <- c(1.03e-3, 7e-4, 1.5e-3, 1.03e-3)[(i - 1) %% 4 + 1]
  sdg <- c(2.17, 1.5, 3.2, 2.9)[(i - 1) %% 4 + 1]
  ls <- c(0.2539, 0.22, 0.28, 0.27)[(i - 1) %% 4 + 1]
  z <- if (i <= 4) rnorm(n) else {
    K <- sdg^2 * exp(-outer(d$x, d$x, "-")^2 / (2 * (ls * dmax)^2)) + diag(1e-12 * sdg^2, n)
    V <- K + diag(sg^2, n)
    forwardsolve(t(chol(K)), as.numeric(K %*% solve(V, d$y - mean(d$y))))
  }
  if (i <= 4) next
  b0 <- mean(d$y)
  p <- list(Intercept = b0, sdgp_1 = array(sdg, 1L), lscale_1 = matrix(ls, 1, 1),
            zgp_1 = array(z, n), sigma = sg)
  u <- rstan::unconstrain_pars(sf, p)
  lp <- rstan::log_prob(sf, u, adjust_transform = FALSE, gradient = TRUE)
  gs <- attr(lp, "gradient")
  par <- tape$par
  nm <- names(par)
  par[nm == "beta"] <- b0
  par[nm == "betad"] <- log(sg)
  par[nm == "theta"] <- c(log(sdg), log(ls * dmax))
  par[nm == "b"] <- z
  gt <- -as.numeric(tape$gr(par))
  gt_st <- c(gt[nm == "beta"], gt[nm == "theta"], gt[nm == "b"], gt[nm == "betad"])
  u0 <- mp(c(b0, log(sdg), log(ls * dmax), z, log(sg)))
  h <- mp(1e-40)
  gx <- vapply(seq_along(u0), function(k) {
    e <- u0; e[k] <- e[k] + h
    m <- u0; m[k] <- m[k] - h
    as.numeric((lpx(e) - lpx(m)) / (2 * h))
  }, 0)
  lab <- c("b0", "logsd", "logrho", paste0("z", 1:n), "logsig")
  cat(sprintf("GX pt %d sigma %.2e sdgp %.2f | max|exact| %.3e | max |tape - exact| %.3e | max |stan - exact| %.3e\n",
              i, sg, sdg, max(abs(gx)), max(abs(gt_st - gx)), max(abs(gs - gx))))
  if (i == 5) {
    print(signif(data.frame(exact = gx, tape = gt_st, stan = as.numeric(gs),
                            row.names = lab), 6))
  }
}
