# Short check of fix round 1, item 1: same data sorted, so both factors run over
# one order and z means one thing. Is the sampler's tape the same function as
# brms's program with its line at the sampler's relative 1e-12, value and gradient,
# at points near the target? A gap here would be in this change; none means the
# remaining difference from brms is the order of the factor. Seed 6601, n 15.
#
#   bash dev/review-gpgap-run.sh dev/review-gpgap-grad.R grad.txt
.libPaths(c("C:/Users/adf44/source/r/wt-gpgap-lib",
            "C:/Users/adf44/source/r/rellib-r16",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(frmtmb); library(testthat)})
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
cat("tape order:", paste(unique(names(tape$par)), collapse = " "), "\n")
dmax <- sdat$dmax_1
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
  b0 <- mean(d$y)
  p <- list(Intercept = b0, sdgp_1 = array(sdg, 1L), lscale_1 = matrix(ls, 1, 1),
            zgp_1 = array(z, n), sigma = sg)
  u <- rstan::unconstrain_pars(sf, p)
  lp <- rstan::log_prob(sf, u, adjust_transform = FALSE, gradient = TRUE)
  gs <- attr(lp, "gradient")
  # the tape's vector, by name
  par <- tape$par
  nm <- names(par)
  par[nm == "beta"] <- b0
  par[nm == "betad"] <- log(sg)
  par[nm == "theta"] <- c(log(sdg), log(ls * dmax))
  par[nm == "b"] <- z
  ft <- tape$fn(par)
  gt <- as.numeric(tape$gr(par))
  # Stan's order: Intercept, sdgp, lscale, z, sigma; tape's by name
  gt_st <- c(gt[nm == "beta"], gt[nm == "theta"], gt[nm == "b"], gt[nm == "betad"])
  # both drop the same constants, up to one number that is the same at every point
  cat(sprintf(paste0("GRAD pt %d %s sigma %.2e sdgp %.2f | value: tape %.6f stan %.6f, diff %.3e | ",
                       "max |grad diff| %.3e, max |grad| %.3e, relative %.2e\n"),
              i, if (i <= 4) "typical z" else "at the field", sg, sdg, -ft, as.numeric(lp),
              -ft - as.numeric(lp),
              max(abs(-gt_st - gs)), max(abs(gs)), max(abs(-gt_st - gs)) / max(abs(gs))))
}
ee <- 1000L
t0 <- proc.time()[["elapsed"]]
for (k in seq_len(ee)) tape$gr(par)
t1 <- proc.time()[["elapsed"]]
for (k in seq_len(ee)) rstan::grad_log_prob(sf, u, adjust_transform = FALSE)
t2 <- proc.time()[["elapsed"]]
cat(sprintf("GRAD cost per gradient, n %d: tape %.1f us, Stan through R %.1f us\n", n,
            1e6 * (t1 - t0) / ee, 1e6 * (t2 - t1) / ee))
