# frm_sample() against brms 2.23.0's own draws (no weights, nothing
# converted) for sigma, ar/ma and b under a cov = TRUE ar(), ma() and
# arma(), under the default priors and under user priors, and under
# sigma ~ x2 with a user prior on the sigma intercept and slope.
#
# Arms:
#   before  rellib-r14 (frmtmb 0.76.0 / frmtmb.sample 0.23.0), its draws
#           as reported (sigma the marginal sd there)
#   after   the lane library
#   brms    brms's raw stanfit from stancode()/standata(), compiled
#           through the Stan cache
#   summary z-scores as dev/acprior-innov.R reports them: the difference
#           of the means over the pooled Monte Carlo standard error, and
#           the sd ratio
# 4 chains of 4000, seed 3, on every side. Data (seed 20261027): ar, ma,
# arma as in dev/acprior-innov.R; sigx as in dev/sigmainn-route.R.
# Usage: Rscript dev/sigmainn-agree.R before|after|brms|summary
# Log: dev/sigmainn-out/agree-<arm>.txt, draws agree-<arm>.rds
arm <- commandArgs(TRUE)[1]
wt <- "C:/Users/adf44/source/r/frmtmb-wt-sigmainn"
out <- file.path(wt, "dev", "sigmainn-out")
.libPaths(c(if (identical(arm, "after")) "C:/Users/adf44/source/r/wt-sigmainn-lib",
            if (!identical(arm, "brms")) "C:/Users/adf44/source/r/rellib-r14",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
Sys.setenv(FRMTMB_STAN_CACHE = file.path(wt, "dev", "stan-cache"))
src <- readLines(file.path(wt, "dev", "acprior-agree.R"))
i0 <- grep("^make_data <- function", src)
i1 <- i0 + which(src[(i0 + 1):length(src)] == "}")[1]
eval(parse(text = src[i0:i1]))
sim_arma <- function(ng, k, ar, ma, seed) {
  set.seed(seed)
  d <- data.frame(t = rep(seq_len(k), ng), g = rep(seq_len(ng), each = k),
                  x = rnorm(ng * k))
  e <- unlist(lapply(seq_len(ng), function(i) {
    as.numeric(stats::arima.sim(list(ar = if (ar != 0) ar,
                                      ma = if (ma != 0) ma), n = k))
  }))
  d$y <- 1 + 0.5 * d$x + e
  d
}
mk_sigx <- function() {
  set.seed(20261027)
  ng <- 20; k <- 6
  d <- data.frame(t = rep(seq_len(k), ng), g = rep(seq_len(ng), each = k),
                  x = rnorm(ng * k), x2 = rnorm(ng * k))
  si <- exp(-0.8 + 0.3 * d$x2)
  e <- numeric(nrow(d))
  for (i in seq_len(nrow(d))) {
    e[i] <- if (d$t[i] == 1) rnorm(1, 0, 1 / sqrt(1 - 0.49)) else
      0.7 * e[i - 1] + rnorm(1)
  }
  d$y <- 1 + 0.5 * d$x + si * e
  d
}
dat <- list(ar = make_data()$ar,
            ma = sim_arma(10, 6, 0, 0.6, 20261027),
            arma = sim_arma(12, 6, 0.7, 0.3, 20261027))
dat$ar_user <- dat$ar
dat$sigx <- mk_sigx()
fo <- list(ar = y ~ x + ar(time = t, gr = g, cov = TRUE),
           ma = y ~ x + ma(time = t, gr = g, cov = TRUE),
           arma = y ~ x + arma(time = t, gr = g, cov = TRUE),
           ar_user = y ~ x + ar(time = t, gr = g, cov = TRUE),
           sigx = y ~ x + ar(time = t, gr = g, cov = TRUE))
dp <- list(sigx = sigma ~ x2)
# user priors written once, for brms and for frmtmb (the same strings)
up <- list(
  ar_user = list(list("normal(0, 0.5)", "sigma", "", NA, 1),
                 list("normal(0.5, 0.3)", "ar", "", NA, NA)),
  sigx = list(list("normal(-1, 0.5)", "Intercept", "sigma", NA, NA),
              list("normal(0, 0.5)", "b", "sigma", NA, NA)))
CH <- 4; IT <- 4000; SEED <- 3
pars_of <- function(m) {
  c(if (m %in% c("ar", "arma", "ar_user", "sigx")) "ar[1]",
    if (m %in% c("ma", "arma")) "ma[1]",
    if (m == "sigx") c("b_sigma_Intercept", "b_sigma_x2") else "sigma",
    "b_Intercept", "b_x")
}

if (arm %in% c("before", "after")) {
  suppressPackageStartupMessages({
    library(frmtmb)
    library(frmtmb.sample)
  })
  cat("frmtmb", format(packageVersion("frmtmb")), "from",
      dirname(find.package("frmtmb")), "; frmtmb.sample from",
      dirname(find.package("frmtmb.sample")), "\n")
  res <- list()
  for (m in names(fo)) {
    pr <- if (!is.null(up[[m]])) do.call(c, lapply(up[[m]], function(u) {
      frmtmb::set_prior(u[[1]], class = u[[2]], dpar = u[[3]], lb = u[[4]],
                        ub = u[[5]])
    }))
    form <- if (!is.null(dp[[m]])) bf(fo[[m]], dp[[m]]) else bf(fo[[m]])
    msgs <- character(0)
    ds <- tryCatch(withCallingHandlers(
      frm_sample(form, family = gaussian(), data = dat[[m]], prior = pr,
                 chains = CH, iter = IT, seed = SEED, refresh = 0),
      message = function(c) {
        msgs <<- c(msgs, conditionMessage(c))
        invokeRestart("muffleMessage")
      },
      warning = function(w) {
        msgs <<- c(msgs, paste("WARNING:", conditionMessage(w)))
        invokeRestart("muffleWarning")
      }), error = function(e) {
        msgs <<- c(msgs, paste("ERROR:", conditionMessage(e)))
        NULL
      })
    cat("\n--", m, ":\n", paste(msgs, collapse = "\n"), "\n")
    if (is.null(ds)) next
    mm <- as.matrix(ds)
    res[[m]] <- mm[, pars_of(m), drop = FALSE]
    cat("divergent:", sum(rstan::get_divergent_iterations(ds$stanfit)),
        "\n")
  }
  saveRDS(res, file.path(out, paste0("agree-", arm, ".rds")))
} else if (identical(arm, "brms")) {
  suppressPackageStartupMessages(library(brms))
  cat("brms", format(packageVersion("brms")), "\n")
  res <- list()
  for (m in names(fo)) {
    pr <- if (!is.null(up[[m]])) do.call(c, lapply(up[[m]], function(u) {
      brms::set_prior(u[[1]], class = u[[2]], dpar = u[[3]], lb = u[[4]],
                      ub = u[[5]])
    })) else brms::empty_prior()
    form <- if (!is.null(dp[[m]])) brms::bf(fo[[m]], dp[[m]]) else
      brms::bf(fo[[m]])
    code <- stancode(form, data = dat[[m]], prior = pr)
    sdat <- standata(form, data = dat[[m]], prior = pr)
    cat("\n--", m, "\n")
    print(grep("sigma|lprior", strsplit(code, "\n")[[1]], value = TRUE))
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
    # the raw stanfit names the slopes b[1] and b_sigma[1]
    raw <- sub("^b_x$", "b[1]", sub("^b_sigma_x2$", "b_sigma[1]", pn))
    mm <- as.matrix(sf, pars = raw)
    colnames(mm) <- pn
    res[[m]] <- mm
  }
  saveRDS(res, file.path(out, "agree-brms.rds"))
} else {
  .libPaths("C:/Users/adf44/AppData/Local/R/win-library/4.6")
  B <- readRDS(file.path(out, "agree-brms.rds"))
  for (tg in c("before", "after")) {
    p <- file.path(out, paste0("agree-", tg, ".rds"))
    if (!file.exists(p)) next
    Fm <- readRDS(p)
    cat("\n#### frm_sample", tg, "against brms's own draws\n")
    for (m in names(B)) {
      cat("\n=====", m, "\n")
      if (is.null(Fm[[m]])) {
        cat("(no draws: refused, see agree-", tg, ".txt)\n", sep = "")
        next
      }
      for (pn in colnames(B[[m]])) {
        f <- Fm[[m]][, pn]
        b <- B[[m]][, pn]
        mf <- posterior::mcse_mean(matrix(f, ncol = CH))
        mb <- posterior::mcse_mean(matrix(b, ncol = CH))
        z <- (mean(f) - mean(b)) / sqrt(mf^2 + mb^2)
        ef <- posterior::ess_bulk(matrix(f, ncol = CH))
        cat(sprintf("%-18s frm mean %.4f sd %.4f | brms mean %.4f sd %.4f | z %7.2f | sd ratio %.3f | frm ess_bulk %.0f\n",
                    pn, mean(f), sd(f), mean(b), sd(b), z, sd(f) / sd(b),
                    ef))
      }
    }
  }
}
