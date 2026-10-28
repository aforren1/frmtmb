# Reviewer rerun of frm_sample() against brms 2.23.0's own draws, with
# a seed the lane did not use (11) and two models it did not sample:
#   stuar  student() + ar(cov = TRUE), ragged groups, default priors
#   armasx gaussian arma(cov = TRUE) + sigma ~ x2, default priors
# Arms: after (lane build), brms (raw stanfit through the Stan cache),
# summary (z = difference of means over the pooled MCSE).
# Usage: Rscript dev/review-sigmainn-agree.R after|before|brms|summary
arm <- commandArgs(TRUE)[1]
wt <- "C:/Users/adf44/source/r/frmtmb-wt-sigmainn"
out <- file.path(wt, "dev", "review-sigmainn-out")
.libPaths(c(if (identical(arm, "after")) "C:/Users/adf44/source/r/wt-sigmainn-lib",
            if (!identical(arm, "brms")) "C:/Users/adf44/source/r/rellib-r14",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win",
           FRMTMB_STAN_CACHE = file.path(wt, "dev", "stan-cache"))
CH <- 4; IT <- 3000; SEED <- as.integer(Sys.getenv("REV_SEED", "11"))
TAG <- Sys.getenv("REV_TAG", "")
set.seed(20261030)
mk <- function(ng, nt, ar, ma, sx) {
  do.call(rbind, lapply(seq_len(ng), function(i) {
    k <- sample(5:nt, 1)
    e <- as.numeric(arima.sim(list(ar = ar, ma = ma), k, sd = 1))
    data.frame(g = i, t = seq_len(k), e = e)
  })) -> d
  d$x <- rnorm(nrow(d)); d$x2 <- rnorm(nrow(d))
  d$y <- 1 + 0.5 * d$x + exp(-0.7 + sx * d$x2) * d$e
  d
}
dat <- list(stuar = mk(16, 9, 0.6, 0, 0), armasx = mk(16, 9, 0.6, 0.3, 0.3))
dat$stuar$y <- dat$stuar$y + 0.3 * rt(nrow(dat$stuar), 4)
fo <- list(stuar = y ~ x + ar(time = t, gr = g, cov = TRUE),
           armasx = y ~ x + arma(time = t, gr = g, cov = TRUE))
dp <- list(armasx = sigma ~ x2)
pars_of <- function(m) {
  if (m == "stuar") c("ar[1]", "sigma", "nu", "b_Intercept", "b_x") else
    c("ar[1]", "ma[1]", "b_sigma_Intercept", "b_sigma_x2", "b_Intercept", "b_x")
}
if (arm %in% c("before", "after")) {
  suppressPackageStartupMessages({
    library(frmtmb)
    library(frmtmb.sample)
  })
  cat("frmtmb from", dirname(find.package("frmtmb")), "\n")
  res <- list()
  for (m in names(fo)) {
    form <- if (!is.null(dp[[m]])) bf(fo[[m]], dp[[m]]) else bf(fo[[m]])
    fam <- if (m == "stuar") student() else gaussian()
    ds <- frm_sample(form, family = fam, data = dat[[m]], chains = CH,
                     iter = IT, seed = SEED, refresh = 0)
    mm <- as.matrix(ds)
    res[[m]] <- mm[, pars_of(m), drop = FALSE]
    cat(m, "divergent:", sum(rstan::get_divergent_iterations(ds$stanfit)), "\n")
  }
  saveRDS(res, file.path(out, paste0("agree-", arm, TAG, ".rds")))
} else if (identical(arm, "brms")) {
  suppressPackageStartupMessages(library(brms))
  res <- list()
  for (m in names(fo)) {
    form <- if (!is.null(dp[[m]])) brms::bf(fo[[m]], dp[[m]]) else brms::bf(fo[[m]])
    fam <- if (m == "stuar") brms::student() else gaussian()
    code <- stancode(form, data = dat[[m]], family = fam)
    sdat <- standata(form, data = dat[[m]], family = fam)
    f <- tempfile(fileext = ".stan")
    writeLines(c(code, paste0("// rstan ", packageVersion("rstan"))), f)
    path <- file.path(Sys.getenv("FRMTMB_STAN_CACHE"),
                      paste0(unname(tools::md5sum(f)), ".rds"))
    mod <- if (file.exists(path)) readRDS(path) else {
      mm <- rstan::stan_model(model_code = code, save_dso = TRUE)
      saveRDS(mm, path)
      mm
    }
    sf <- rstan::sampling(mod, data = sdat, chains = CH, iter = IT,
                          seed = SEED, refresh = 0)
    cat(m, "divergent", rstan::get_num_divergent(sf), "\n")
    pn <- pars_of(m)
    raw <- sub("^b_x$", "b[1]", sub("^b_sigma_x2$", "b_sigma[1]", pn))
    mm <- as.matrix(sf, pars = raw)
    colnames(mm) <- pn
    res[[m]] <- mm
  }
  saveRDS(res, file.path(out, paste0("agree-brms", TAG, ".rds")))
} else {
  suppressPackageStartupMessages(library(posterior))
  B <- readRDS(file.path(out, paste0("agree-brms", TAG, ".rds")))
  for (a in c("after", "before")) {
    fn <- file.path(out, paste0("agree-", a, TAG, ".rds"))
    if (!file.exists(fn)) next
    A <- readRDS(fn)
    cat("####", a, "\n")
    for (m in names(B)) {
      zz <- vapply(colnames(B[[m]]), function(p) {
        x <- A[[m]][, p]; y <- B[[m]][, p]
        mx <- mcse_mean(matrix(x, ncol = CH)); my <- mcse_mean(matrix(y, ncol = CH))
        c(z = (mean(x) - mean(y)) / sqrt(mx^2 + my^2), r = sd(x) / sd(y),
          ma = mean(x), mb = mean(y))
      }, numeric(4))
      cat(m, paste(sprintf("%s z %.2f (sd ratio %.3f; %.4f vs %.4f)", colnames(zz),
                           zz["z", ], zz["r", ], zz["ma", ], zz["mb", ]),
                   collapse = "; "), "\n")
    }
  }
}
