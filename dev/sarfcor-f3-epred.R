# Fix round 1, F3: frmtmb's fitted() on same-row newdata against brms
# 2.23.0's posterior_epred() on the same newdata, with brms's draws held
# at frmtmb's estimates (Stan's Fixed_param). Rscript <this>
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win",
           FRMTMB_STAN_CACHE =
             "C:/Users/adf44/source/r/frmtmb-wt-sarfcor/dev/stan-cache",
           FRMTMB_BRMS_FIT_TESTS = "true", NOT_CRAN = "true")
.libPaths(c("C:/Users/adf44/source/r/wt-sarfcor-lib",
            "C:/Users/adf44/source/r/rellib-r12",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb.spatial))
wt <- "C:/Users/adf44/source/r/frmtmb-wt-sarfcor"
env <- new.env(parent = asNamespace("frmtmb"))
sys.source(file.path(wt, "extensions/frmtmb.spatial/tests/testthat",
                     "helper-brms.R"), envir = env)
epred_at <- function(bform, data, data2, pars, newdata) {
  fam <- gaussian()
  prior <- env$brms_flat_prior(bform, data = data, family = fam,
                               data2 = data2)
  code <- brms::make_stancode(bform, data = data, family = fam,
                              prior = prior, data2 = data2)
  sdat <- env$brms_standata(bform, data = data, family = fam,
                            prior = prior, data2 = data2)
  sc <- suppressMessages(brms::brm(bform, data = data, family = fam,
                                   prior = prior, data2 = data2,
                                   empty = TRUE))
  sf <- suppressMessages(rstan::sampling(
    env$brms_stan_model(code), data = sdat, chains = 1, iter = 2,
    warmup = 0, algorithm = "Fixed_param", init = list(pars),
    refresh = 0, seed = 1))
  sc$fit <- sf
  bfit <- brms::rename_pars(sc)
  if (is.null(newdata)) return(bfit)
  brms::posterior_epred(bfit, newdata = newdata)
}
set.seed(2)
k <- 7
g <- expand.grid(r = seq_len(k), c = seq_len(k))
n <- nrow(g)
A <- outer(seq_len(n), seq_len(n), function(i, j) {
  as.numeric(abs(g$r[i] - g$r[j]) + abs(g$c[i] - g$c[j]) == 1)
})
dg <- 1 / sqrt(rowSums(A))
Ws <- A * outer(dg, dg)
d <- data.frame(x = rnorm(n))
d$y <- as.vector(solve(diag(n) - 0.5 * Ws, 1 + 0.5 * d$x + rnorm(n, 0, 0.7)))
V <- 0.6^abs(outer(1:n, 1:n, "-"))
nd <- d
nd$x <- nd$x + 1
for (type in c("lag", "error")) {
  f <- frm(as.formula(paste0("y ~ x + sar(Ws, type = '", type, "')")),
           data = d, data2 = list(Ws = Ws))
  b <- unname(coef(f))
  ac <- f$frame$autocor[[1]]
  rho <- unname(frmtmb:::autocor_natural(f$estimates$thetaac, ac))
  pars <- list(b = array(b[2], 1), Intercept = b[1] + b[2] * mean(d$x),
               sigma = sigma(f))
  pars[[paste0(type, "sar")]] <- rho
  bfit <- epred_at(brms::bf(as.formula(paste0(
    "y ~ x + sar(Ws, type = '", type, "')"))), d, list(Ws = Ws), pars,
    NULL)
  be <- tryCatch(brms::posterior_epred(bfit, newdata = nd),
                 error = function(e) NULL)
  if (is.null(be)) {
    # brms 2.23.0 on an error SAR: every method that reaches
    # posterior_epred_gaussian() fails, with or without newdata
    for (call in c("brms::posterior_epred(bfit)", "stats::fitted(bfit)",
                   "brms::posterior_predict(bfit)", "brms::log_lik(bfit)")) {
      r <- tryCatch({eval(parse(text = call), list(bfit = bfit)); "ok"},
                    error = function(e) conditionMessage(e))
      cat(sprintf("sar %-5s %-30s %s\n", type, call, r))
    }
    next
  }
  fe <- fitted(f, newdata = nd)[, "Estimate"]
  cat(sprintf("sar %-5s draws %d  max |brms - frmtmb| / max |frmtmb| = %.3g  (moved from fitted rows %.3g)\n",
              type, nrow(be), max(abs(be[1, ] - fe)) / max(abs(fe)),
              max(abs(fe - fitted(f)[, "Estimate"]))))
}
d$yf <- 1 + 0.5 * d$x + as.vector(t(chol(V)) %*% rnorm(n)) * 0.8
f <- frm(yf ~ x + fcor(V), data = d, data2 = list(V = V))
b <- unname(coef(f))
pars <- list(b = array(b[2], 1), Intercept = b[1] + b[2] * mean(d$x),
             sigma = sigma(f))
be <- epred_at(brms::bf(yf ~ x + fcor(V)), d, list(V = V), pars, nd)
fe <- fitted(f, newdata = nd)[, "Estimate"]
cat(sprintf("fcor       draws %d  max |brms - frmtmb| / max |frmtmb| = %.3g\n",
            nrow(be), max(abs(be[1, ] - fe)) / max(abs(fe))))
# brms on rows that differ: the same refusal frmtmb gives
r <- tryCatch(epred_at(brms::bf(yf ~ x + fcor(V)), d, list(V = V), pars,
                       nd[1:5, ]), error = function(e) conditionMessage(e))
cat("brms on 5 rows:", if (is.character(r)) r else "no error", "\n")
