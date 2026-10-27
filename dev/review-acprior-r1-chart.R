# Reviewer, fix round 1, item 3: the cosy chart.
#  (a) the stored draws are the sampler's v carried back: cosy from the
#      stanfit's raw column through L + (U - L) plogis(v) equals the
#      draws' cosy column;
#  (b) cosy + mo(): simplex draws on the simplex, and cosy, the simplex
#      and bsp against brms 2.23.0 on the same data (4 x 2000, seed 3);
#  (c) a user prior and user bounds on cosy; lb = -1/(d - 1) lifts the
#      default bound with no warning;
#  (d) log_lik() and loo() on charted draws.
# Usage: Rscript dev/review-acprior-r1-chart.R frm|brms|summary
# Log: dev/review-acprior-log/r1-chart-<arm>.txt
arm <- commandArgs(TRUE)[1]
wt <- "C:/Users/adf44/source/r/frmtmb-wt-acprior"
lg <- file.path(wt, "dev", "review-acprior-log")
ulib <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
.libPaths(c(if (identical(arm, "frm")) "C:/Users/adf44/source/r/wt-acprior-lib",
            "C:/Users/adf44/source/r/rellib-r13", ulib))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
Sys.setenv(FRMTMB_STAN_CACHE = file.path(wt, "dev", "stan-cache"))
set.seed(31)
ng <- 40; k <- 4
d <- data.frame(t = rep(seq_len(k), ng), g = rep(seq_len(ng), each = k),
                xo = sample(1:4, ng * k, replace = TRUE), x = rnorm(ng * k))
zeta <- c(0, 0.6, 0.9, 1)
d$y <- 1 + 0.4 * d$x + 1.5 * zeta[d$xo] + rnorm(ng, 0, sqrt(0.02))[d$g] +
  rnorm(ng * k, 0, sqrt(0.98))
d$xo <- factor(d$xo, ordered = TRUE)
fo <- y ~ x + mo(xo) + cosy(time = t, gr = g)
CH <- 4; IT <- 2000; SEED <- 3
catch <- function(expr) {
  w <- character(0)
  v <- withCallingHandlers(suppressMessages(expr), warning = function(x) {
    w <<- c(w, conditionMessage(x))
    invokeRestart("muffleWarning")
  })
  list(value = v, warnings = w)
}
if (identical(arm, "frm")) {
  suppressPackageStartupMessages({
    library(frmtmb)
    library(frmtmb.sample)
  })
  cat("frmtmb from", dirname(find.package("frmtmb")), "\n")
  fit <- frm(fo, data = d)
  cat("ML cosy", format(unlist(frmtmb:::autocor_natural(
    fit$estimates$thetaac, fit$frame$autocor[[1]])), digits = 4), "\n")
  r <- catch(frm_sample(fit, chains = CH, iter = IT, seed = SEED,
                        refresh = 0))
  ds <- r$value
  cat("warnings:", if (length(r$warnings)) r$warnings else "none", "\n")
  sp <- rstan::get_sampler_params(ds$stanfit, inc_warmup = FALSE)
  cat("divergent:", sum(unlist(lapply(sp, function(m) m[, "divergent__"]))),
      "\n")
  mm <- as.matrix(ds)
  cat("columns:", paste(colnames(mm), collapse = " "), "\n")
  # (a) the raw chart coordinate against the stored cosy
  raw <- rstan::extract(ds$stanfit, permuted = FALSE)
  rawm <- do.call(rbind, lapply(seq_len(dim(raw)[2]), function(ch) raw[, ch, ]))
  ch <- frmtmb::autocor_bound_chart(fit, list(thetaac_1 = qlogis(
    (0 + 1 / 3) / (1 + 1 / 3))), list())
  cat("chart terms:", length(ch), "L", ch[[1]]$L, "U", ch[[1]]$U, "\n")
  jcol <- which(apply(rawm, 2, function(v) {
    isTRUE(all.equal(ch[[1]]$L + (ch[[1]]$U - ch[[1]]$L) * plogis(v),
                     unname(mm[, "cosy"]), tolerance = 1e-10))
  }))
  cat("raw columns whose chart image equals the cosy draws:", jcol, "\n")
  cat(sprintf("cosy min %.3g mean %.4f\n", min(mm[, "cosy"]), mean(mm[, "cosy"])))
  simo <- mm[, grep("^simo", colnames(mm)), drop = FALSE]
  cat(sprintf("simplex: %d columns, max |sum - 1| %.2e, min %.3g\n",
              ncol(simo), max(abs(rowSums(simo) - 1)), min(simo)))
  saveRDS(mm, file.path(lg, "r1-chart-frm.rds"))
  # (d)
  ll <- tryCatch(log_lik(ds), error = function(e) conditionMessage(e))
  if (is.matrix(ll)) {
    cat(sprintf("log_lik: %d x %d, all finite %s\n", nrow(ll), ncol(ll),
                all(is.finite(ll))))
    lo <- tryCatch(suppressWarnings(loo::loo(ll)),
                   error = function(e) conditionMessage(e))
    if (inherits(lo, "loo")) {
      cat(sprintf("loo elpd %.2f (se %.2f)\n",
                  lo$estimates["elpd_loo", 1], lo$estimates["elpd_loo", 2]))
    } else cat("loo:", lo, "\n")
  } else cat("log_lik:", ll, "\n")
  # (c) user prior and user bounds, cosy alone
  fc <- frm(y ~ x + cosy(time = t, gr = g), data = d)
  for (nm in c("normal", "box", "lift")) {
    pr <- switch(nm,
      normal = set_prior("normal(0, 0.1)", class = "cosy"),
      box = set_prior("", class = "cosy", lb = 0.02, ub = 0.3),
      lift = set_prior("", class = "cosy", lb = -1/3))
    r <- catch(frm_sample(fc, prior = pr, chains = 2, iter = 1000,
                          seed = SEED, refresh = 0))
    cz <- as.matrix(r$value)[, "cosy"]
    sp <- rstan::get_sampler_params(r$value$stanfit, inc_warmup = FALSE)
    cat(sprintf("%-6s cosy min %.4g max %.4g mean %.4f; divergent %d; warnings: %s\n",
                nm, min(cz), max(cz), mean(cz),
                sum(unlist(lapply(sp, function(m) m[, "divergent__"]))),
                if (length(r$warnings)) paste(r$warnings, collapse = " | ")
                else "none"))
  }
} else if (identical(arm, "brms")) {
  suppressPackageStartupMessages(library(brms))
  mod <- rstan::stan_model(model_code = stancode(fo, data = d))
  sf <- rstan::sampling(mod, data = standata(fo, data = d), chains = CH,
                        iter = IT, seed = SEED, refresh = 0)
  cat("divergent:", rstan::get_num_divergent(sf), "\n")
  mm <- as.matrix(sf)
  mm <- mm[, grep("^(cosy|simo|bsp|sigma|b\\[|b_Intercept)", colnames(mm)),
           drop = FALSE]
  cat("columns:", paste(colnames(mm), collapse = " "), "\n")
  saveRDS(mm, file.path(lg, "r1-chart-brms.rds"))
} else {
  .libPaths(ulib)
  F <- readRDS(file.path(lg, "r1-chart-frm.rds"))
  B <- readRDS(file.path(lg, "r1-chart-brms.rds"))
  fs <- grep("^simo", colnames(F), value = TRUE)
  bs <- grep("^simo", colnames(B), value = TRUE)
  pairs <- rbind(c("cosy", "cosy"), cbind(fs, bs),
                 c(grep("^bsp", colnames(F), value = TRUE)[1],
                   grep("^bsp", colnames(B), value = TRUE)[1]))
  for (i in seq_len(nrow(pairs))) {
    f <- F[, pairs[i, 1]]; b <- B[, pairs[i, 2]]
    mf <- posterior::mcse_mean(matrix(f, ncol = CH))
    mb <- posterior::mcse_mean(matrix(b, ncol = CH))
    cat(sprintf("%-14s %-14s frm %.4f brms %.4f z %6.2f sd ratio %.3f\n",
                pairs[i, 1], pairs[i, 2], mean(f), mean(b),
                (mean(f) - mean(b)) / sqrt(mf^2 + mb^2), sd(f) / sd(b)))
  }
}
