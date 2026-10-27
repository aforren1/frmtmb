# Fix round 1, item 1: frm_sample() under its defaults against brms
# 2.23.0's OWN posterior (no weights), for a first-order cov = TRUE ar,
# ma and arma, after frmtmb.sample's sigma default moved to brms's
# innovation sd. Every parameter: ar or ma, sigma (brms's innovation
# sd; frmtmb's draws are converted, sigma_i = sigma * c), b_Intercept,
# b_x. 4 chains of 4000, seed 3, on both sides.
# Data (seed 20261027): ar, the AR(1) 0.9 data of dev/acprior-agree.R;
# ma, an MA(1) with theta 0.6 in 10 groups of 6; arma, ARMA(1, 1) with
# ar 0.7, ma 0.3 in 12 groups of 6.
# Usage: Rscript dev/acprior-innov.R frm|brms|summary  [before]
# Log: dev/acprior-out/innov.txt (summary)
arm <- commandArgs(TRUE)[1]
before <- identical(commandArgs(TRUE)[2], "before")
wt <- "C:/Users/adf44/source/r/frmtmb-wt-acprior"
out <- file.path(wt, "dev", "acprior-out")
.libPaths(c(if (!before) "C:/Users/adf44/source/r/wt-acprior-lib",
            "C:/Users/adf44/source/r/rellib-r13",
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
dat <- list(ar = make_data()$ar,
            ma = sim_arma(10, 6, 0, 0.6, 20261027),
            arma = sim_arma(12, 6, 0.7, 0.3, 20261027))
fo <- list(ar = y ~ x + ar(time = t, gr = g, cov = TRUE),
           ma = y ~ x + ma(time = t, gr = g, cov = TRUE),
           arma = y ~ x + arma(time = t, gr = g, cov = TRUE))
CH <- 4; IT <- 4000; SEED <- 3
cfac <- function(m, ar, ma) {
  ar <- if (is.null(ar)) 0 else ar
  ma <- if (is.null(ma)) 0 else ma
  sqrt((1 - ar^2) / (1 + ma^2 + 2 * ar * ma))
}
tag <- if (before) "before" else "after"
if (identical(arm, "frm")) {
  suppressPackageStartupMessages({
    library(frmtmb)
    library(frmtmb.sample)
  })
  cat("frmtmb from", dirname(find.package("frmtmb")), "\n")
  res <- list()
  for (m in names(fo)) {
    fit <- frm(fo[[m]], data = dat[[m]])
    msgs <- character(0)
    ds <- withCallingHandlers(
      frm_sample(fit, chains = CH, iter = IT, seed = SEED, refresh = 0),
      message = function(c) {
        msgs <<- c(msgs, conditionMessage(c))
        invokeRestart("muffleMessage")
      },
      warning = function(w) {
        msgs <<- c(msgs, paste("WARNING:", conditionMessage(w)))
        invokeRestart("muffleWarning")
      })
    cat("\n--", m, ":\n", paste(msgs, collapse = "\n"), "\n")
    mm <- as.matrix(ds)
    ar <- if ("ar[1]" %in% colnames(mm)) mm[, "ar[1]"]
    ma <- if ("ma[1]" %in% colnames(mm)) mm[, "ma[1]"]
    sig <- mm[, "sigma"] * cfac(m, ar, ma)
    res[[m]] <- cbind(mm[, intersect(c("ar[1]", "ma[1]"), colnames(mm)),
                         drop = FALSE], sigma = sig,
                      b_Intercept = mm[, "b_Intercept"], b_x = mm[, "b_x"])
  }
  saveRDS(res, file.path(out, paste0("innov-frm-", tag, ".rds")))
} else if (identical(arm, "brms")) {
  suppressPackageStartupMessages(library(brms))
  res <- list()
  for (m in names(fo)) {
    code <- stancode(fo[[m]], data = dat[[m]])
    sdat <- standata(fo[[m]], data = dat[[m]])
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
    pars <- c(intersect(c("ar[1]", "ma[1]"),
                        dimnames(as.array(sf))[[3]]), "sigma",
              "b_Intercept", "b_x")
    # the raw stanfit names the slope b[1]; brms renames it b_x
    mm <- as.matrix(sf, pars = sub("^b_x$", "b[1]", pars))
    colnames(mm) <- pars
    res[[m]] <- mm
  }
  saveRDS(res, file.path(out, "innov-brms.rds"))
} else {
  .libPaths("C:/Users/adf44/AppData/Local/R/win-library/4.6")
  R <- readRDS(file.path(out, "innov-brms.rds"))
  for (tg in c("before", "after")) {
    p <- file.path(out, paste0("innov-frm-", tg, ".rds"))
    if (!file.exists(p)) next
    Fm <- readRDS(p)
    cat("\n#### frm_sample", tg, "against brms's own draws\n")
    for (m in names(R)) {
      cat("\n=====", m, "\n")
      for (pn in colnames(R[[m]])) {
        f <- Fm[[m]][, pn]
        b <- R[[m]][, pn]
        mf <- posterior::mcse_mean(matrix(f, ncol = CH))
        mb <- posterior::mcse_mean(matrix(b, ncol = CH))
        z <- (mean(f) - mean(b)) / sqrt(mf^2 + mb^2)
        cat(sprintf("%-12s frm mean %.4f sd %.4f | brms mean %.4f sd %.4f | z %6.2f | sd ratio %.3f\n",
                    pn, mean(f), sd(f), mean(b), sd(b), z, sd(f) / sd(b)))
      }
    }
  }
}
