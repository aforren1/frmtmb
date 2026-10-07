# Reviewer of lane optima, task 4: sample the same objective frm_sample()
# samples (bare likelihood + frm_sample's default Intercept/sigma priors),
# plus brms's simo ~ dirichlet(1) (a constant on the simplex) PLUS the
# log-Jacobian of zeta -> (w_1, ..., w_{D-1}) of the arm's own map. If the
# coordinates are the only difference, both arms then target brms's
# posterior.
#   Rscript dev/optima-rev-sample4-dirichlet.R base|lane
arm <- commandArgs(TRUE)[1]
libs <- switch(arm,
  base = "C:/Users/adf44/source/r/rellib-r6",
  lane = c("C:/Users/adf44/source/r/wt-optima-lib",
           "C:/Users/adf44/source/r/rellib-r6"))
.libPaths(c(libs, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressMessages({library(frmtmb); library(frmtmb.sample)})
cat("arm", arm, "frmtmb from", find.package("frmtmb"), "\n")
out <- "C:/Users/adf44/source/r/frmtmb-wt-optima/dev/optima-rev-out"
ns <- asNamespace("frmtmb"); nss <- asNamespace("frmtmb.sample")
D <- 3L

if (arm == "lane") {
  ms <- ns$mo_simplex
  B <- ns$mo_sphere_basis(D); P <- -1 / sqrt(D)
  # log |det d(w_1..w_{D-1}) / dzeta| of the chart:
  #   (D-1) log 2 + sum_j log|u_j| + (D-1) log(1 + sum(u)/sqrt(D)),
  # and 1 + sum(u)/sqrt(D) = 2 / (1 + r2) on the sphere
  logJ <- function(z) {
    r2 <- sum(z * z)
    acc <- (D - 1) * log(2) + (D - 1) * (log(2) - log(1 + r2))
    for (j in seq_len(D)) {
      uj <- (2 * sum(B[j, ] * z) + (r2 - 1) * P) / (r2 + 1)
      acc <- acc + 0.5 * log(uj * uj)
    }
    acc
  }
} else {
  ms <- function(z) { x <- exp(c(0, z)); x / sum(x) }
  # softmax: sum_j log w_j
  logJ <- function(z) {
    s <- 1
    for (k in seq_along(z)) s <- s + exp(z[k])
    sum(z) - D * log(s)
  }
}
# numerical check of logJ against a finite-difference Jacobian
fdJ <- function(z, h = 1e-6) {
  J <- sapply(seq_along(z), function(k) {
    e <- replace(numeric(length(z)), k, h)
    (ms(z + e)[-D] - ms(z - e)[-D]) / (2 * h)
  })
  log(abs(det(J)))
}
set.seed(5)
for (i in 1:5) {
  z <- rnorm(D - 1, sd = 1.5)
  cat(sprintf("logJ check z = (%s): analytic %.8f  finite-diff %.8f\n",
              paste(sprintf("%.3f", z), collapse = ", "), logJ(z), fdJ(z)))
}

dl <- readRDS(file.path(out, "data.rds"))
cases <- list(weak = list(f = bf(y ~ mo(x)), d = dl$weak),
              strong = list(f = bf(ls ~ mo(income)), d = dl$strong))
res <- list()
for (nm in names(cases)) {
  cs <- cases[[nm]]
  fit <- frm(cs$f, data = cs$d, family = gaussian())
  pl <- nss$default_priors_for(fit)
  ent <- ns$resolve_prior_input(fit, pl)$entries
  nll <- ns$build_objective(fit$frame)
  nlp <- ns$neg_log_prior_fn(ent)
  zi <- grep("^zeta", names(fit$estimates), value = TRUE)
  stopifnot(length(fit$estimates[[zi]]) == D - 1)
  obj <- RTMB::MakeADFun(function(pars) {
    nll(pars) + nlp(pars) - logJ(pars[[zi]])   # dirichlet(1) is constant
  }, fit$estimates, map = fit$frame[["map"]], silent = TRUE)
  # start every chain at the ML point but with the simplex moved to an
  # interior point: the lane's ML weak simplex is ON a face, where the
  # Jacobian is 0 and the log density -Inf
  set.seed(1)
  inits <- lapply(1:2, function(ch) {
    p <- fit$estimates
    p[[zi]] <- rnorm(D - 1, sd = 0.25)
    p
  })
  s <- withCallingHandlers(
    tmbstan::tmbstan(obj, chains = 2, iter = 2000, warmup = 1000, seed = 1,
                     cores = 1, refresh = 0, init = inits),
    warning = function(w) { cat("WARN:", conditionMessage(w), "\n"); invokeRestart("muffleWarning") })
  a <- rstan::extract(s, permuted = FALSE)
  sm <- rstan::summary(s)$summary
  cat("\n==", nm, "+ dirichlet(1) + logJ ==\n")
  print(round(sm[, c("mean", "sd", "n_eff", "Rhat")], 4))
  sp <- rstan::get_sampler_params(s, inc_warmup = FALSE)
  div <- sapply(sp, function(x) sum(x[, "divergent__"]))
  cat("divergent per chain:", div, "\n")
  res[[nm]] <- list(arr = a, summary = sm, div = div)
}
saveRDS(res, file.path(out, paste0("draws-dir-", arm, ".rds")))
