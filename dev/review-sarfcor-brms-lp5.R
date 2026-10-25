# Reviewer: brms 2.23.0's own Stan log density for sar() lag/error and
# fcor() under student() at nu = 5 (and gaussian), at hand-set
# parameters away from any estimate, against frmtmb's objective at the
# same point and a dense mvtnorm reference. Uses the lane's helper for
# the cached Stan models (flat priors). Rscript dev/review-sarfcor-brms-lp5.R
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win",
           FRMTMB_STAN_CACHE =
             "C:/Users/adf44/source/r/frmtmb-wt-sarfcor/dev/stan-cache")
.libPaths(c("C:/Users/adf44/source/r/wt-sarfcor-lib",
            "C:/Users/adf44/source/r/rellib-r12",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(frmtmb.spatial); library(mvtnorm)})
wt <- "C:/Users/adf44/source/r/frmtmb-wt-sarfcor"
env <- new.env(parent = asNamespace("frmtmb"))
sys.source(file.path(wt, "extensions/frmtmb.spatial/tests/testthat",
                     "helper-brms.R"), envir = env)
cat("lib:", find.package("frmtmb"), "\n")
set.seed(2)
k <- 7
gr <- expand.grid(r = seq_len(k), c = seq_len(k))
n <- nrow(gr)
A <- outer(seq_len(n), seq_len(n), function(i, j) {
  as.numeric(abs(gr$r[i] - gr$r[j]) + abs(gr$c[i] - gr$c[j]) == 1)
})
W <- A / rowSums(A)
ds <- data.frame(x = rnorm(n))
IA <- diag(n) - 0.5 * W
ds$ylag <- as.vector(solve(IA, 1 + 0.5 * ds$x + rnorm(n, 0, 0.7)))
ds$yerr <- 1 + 0.5 * ds$x + as.vector(solve(IA, rnorm(n, 0, 0.7)))
dg <- 1 / sqrt(rowSums(A))
Ws <- A * outer(dg, dg)
V <- 0.6^abs(outer(1:n, 1:n, "-"))
ds$yfc <- 1 + 0.5 * ds$x + as.vector(t(chol(V)) %*% rnorm(n)) * 0.8

b0 <- 0.8; b1 <- 0.4; sg <- 0.9; rho <- 0.35
X <- cbind(1, ds$x)
mu <- as.vector(X %*% c(b0, b1))
brms_lp <- function(form, fam, resp_rho = NULL, nu = NULL, data2) {
  prior <- env$brms_flat_prior(form, data = ds, family = fam, data2 = data2)
  code <- brms::make_stancode(form, data = ds, family = fam, prior = prior,
                              data2 = data2)
  sdat <- env$brms_standata(form, data = ds, family = fam, prior = prior,
                            data2 = data2)
  mod <- env$brms_stan_model(code)
  sf <- suppressMessages(rstan::sampling(mod, data = sdat, chains = 0))
  pars <- list(b = array(b1, 1), Intercept = b0 + b1 * mean(ds$x),
               sigma = sg)
  if (!is.null(nu)) pars$nu <- nu
  if (!is.null(resp_rho)) pars[[resp_rho]] <- rho
  up <- rstan::unconstrain_pars(sf, pars)
  rstan::log_prob(sf, up, adjust_transform = FALSE, gradient = FALSE)
}
frm_lp <- function(fit, nu = NULL) {
  p <- fit$obj$env$last.par.best
  p[names(p) == "beta"] <- c(b0, b1)
  bd <- which(names(p) == "betad")
  p[bd[1]] <- log(sg)
  if (!is.null(nu)) p[bd[2]] <- log(nu - 1)
  ac <- fit$frame[["autocor"]][[1L]]
  if (length(ac$theta_idx)) {
    aux <- ac$aux
    p[names(p) == "thetaac"] <- qlogis((rho - aux$lo) / (aux$hi - aux$lo))
  }
  -fit$obj$fn(p)
}
row <- function(lab, b, f, dense, const) {
  cat(sprintf("%-22s brms %.10f frmtmb %.10f dense %.10f | brms - frmtmb - const %.3g | frmtmb - dense %.3g\n",
              lab, b, f, dense, b - f - const, f - dense))
}
c2 <- 0.5 * n * log(2 * pi); cp <- 0.5 * n * log(pi)
for (type in c("lag", "error")) {
  y <- if (type == "lag") ds$ylag else ds$yerr
  yn <- if (type == "lag") "ylag" else "yerr"
  Aw <- diag(n) - rho * Ws
  m <- if (type == "lag") solve(Aw, mu) else mu
  S <- sg^2 * solve(crossprod(Aw))
  form <- brms::bf(as.formula(paste(yn, "~ x + sar(Ws, type = '", type, "')", sep = "")))
  fform <- bf(as.formula(paste(yn, "~ x + sar(Ws, type = '", type, "')", sep = "")))
  pn <- paste0(type, "sar")
  fg <- frm(fform, data = ds, data2 = list(Ws = Ws))
  row(paste("gauss", type), brms_lp(form, gaussian(), pn, NULL, list(Ws = Ws)),
      frm_lp(fg), dmvnorm(y, m, S, log = TRUE), c2)
  fs <- frm(fform, family = student(), data = ds, data2 = list(Ws = Ws))
  for (nu in c(3, 5, 30)) {
    row(paste("student", type, "nu", nu),
        brms_lp(form, brms::student(), pn, nu, list(Ws = Ws)),
        frm_lp(fs, nu), dmvt(y, m, S, df = nu, log = TRUE), cp)
  }
}
fg <- frm(yfc ~ x + fcor(V), data = ds, data2 = list(V = V))
row("gauss fcor", brms_lp(brms::bf(yfc ~ x + fcor(V)), gaussian(), NULL, NULL,
                          list(V = V)),
    frm_lp(fg), dmvnorm(ds$yfc, mu, sg^2 * V, log = TRUE), 0)
fs <- frm(yfc ~ x + fcor(V), family = student(), data = ds,
          data2 = list(V = V))
for (nu in c(3, 5, 30)) {
  row(paste("student fcor nu", nu),
      brms_lp(brms::bf(yfc ~ x + fcor(V)), brms::student(), NULL, nu,
              list(V = V)),
      frm_lp(fs, nu), dmvt(ds$yfc, mu, sg^2 * V, df = nu, log = TRUE), 0)
}
