# A prior of the user's own on an autocorrelation class, under
# frm_sample() (before and after the lane) and under brms 2.23.0, on the
# data of dev/acprior-agree.R (seed 20261027):
#   ar_normal: ar(cov = TRUE), set_prior("normal(0.5, 0.2)", class = "ar")
#   ar_lb:     ar(cov = TRUE), set_prior("", class = "ar", lb = 0.5),
#              a bounds-only prior: brms is flat on [0.5, 1]
#   cosy_norm: cosy(), set_prior("normal(0, 0.1)", class = "cosy"): brms
#              keeps its declared [0, 1] under a user density
# brms's ar draws are importance-weighted to frmtmb's marginal-sigma
# model as in dev/acprior-agree.R.
# Usage: Rscript dev/acprior-agree-up.R before|after|brms|summary [unweighted]
# Log: dev/acprior-out/agree-up.txt (summary), agree-up-<arm>.txt
arm <- commandArgs(TRUE)[1]
wt <- "C:/Users/adf44/source/r/frmtmb-wt-acprior"
out <- file.path(wt, "dev", "acprior-out")
ulib <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
lib <- c(if (identical(arm, "after")) "C:/Users/adf44/source/r/wt-acprior-lib",
         "C:/Users/adf44/source/r/rellib-r13", ulib)
.libPaths(lib)
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
Sys.setenv(FRMTMB_STAN_CACHE = file.path(wt, "dev", "stan-cache"))
# the data of dev/acprior-agree.R, regenerated the same way
src <- readLines(file.path(wt, "dev", "acprior-agree.R"))
i0 <- grep("^make_data <- function", src)
i1 <- i0 + which(src[(i0 + 1):length(src)] == "}")[1]
eval(parse(text = src[i0:i1]))
dat <- make_data()
CH <- 4; IT <- 4000; SEED <- 3
cases <- list(
  ar_normal = list(m = "ar", f = y ~ x + ar(time = t, gr = g, cov = TRUE),
                   p = "normal(0.5, 0.2)", cl = "ar", lb = NA),
  ar_lb = list(m = "ar", f = y ~ x + ar(time = t, gr = g, cov = TRUE),
               p = "", cl = "ar", lb = 0.5),
  cosy_norm = list(m = "cosy", f = y ~ x + cosy(time = t, gr = g),
                   p = "normal(0, 0.1)", cl = "cosy", lb = NA))

if (arm %in% c("before", "after")) {
  suppressPackageStartupMessages({
    library(frmtmb)
    library(frmtmb.sample)
  })
  cat("frmtmb from", dirname(find.package("frmtmb")), "; frmtmb.sample from",
      dirname(find.package("frmtmb.sample")), "\n")
  res <- list()
  for (nm in names(cases)) {
    cs <- cases[[nm]]
    pr <- set_prior(cs$p, class = cs$cl, lb = cs$lb)
    fit <- frm(cs$f, data = dat[[cs$m]])
    msgs <- character(0)
    ds <- withCallingHandlers(
      frm_sample(fit, prior = pr, chains = CH, iter = IT, seed = SEED,
                 refresh = 0),
      message = function(c) {
        msgs <<- c(msgs, conditionMessage(c))
        invokeRestart("muffleMessage")
      },
      warning = function(w) {
        msgs <<- c(msgs, paste("WARNING:", conditionMessage(w)))
        invokeRestart("muffleWarning")
      })
    cat("\n--", nm, "messages:\n", paste(msgs, collapse = "\n"), "\n")
    print(prior_summary(ds))
    col <- if (cs$m == "ar") "ar[1]" else "cosy"
    res[[nm]] <- as.matrix(ds)[, col]
  }
  saveRDS(res, file.path(out, paste0("agree-up-", arm, ".rds")))
} else if (identical(arm, "brms")) {
  suppressPackageStartupMessages(library(brms))
  res <- list()
  for (nm in names(cases)) {
    cs <- cases[[nm]]
    pr <- set_prior(cs$p, class = cs$cl, lb = cs$lb)
    code <- stancode(cs$f, data = dat[[cs$m]], prior = pr)
    cat("\n--", nm, "\n")
    cat(grep(paste0("lprior.*", cs$cl, "|> ", cs$cl, ";"),
             strsplit(code, "\n")[[1]], value = TRUE), sep = "\n")
    sdat <- standata(cs$f, data = dat[[cs$m]], prior = pr)
    f <- tempfile(fileext = ".stan")
    writeLines(c(code, paste0("// rstan ", packageVersion("rstan"))), f)
    key <- unname(tools::md5sum(f))
    path <- file.path(Sys.getenv("FRMTMB_STAN_CACHE"), paste0(key, ".rds"))
    mod <- if (file.exists(path)) readRDS(path) else {
      mm <- rstan::stan_model(model_code = code, save_dso = TRUE)
      saveRDS(mm, path)
      mm
    }
    sf <- rstan::sampling(mod, data = sdat, chains = CH, iter = IT,
                          seed = SEED, refresh = 0)
    a <- as.matrix(sf, pars = c(if (cs$m == "ar") "ar[1]" else "cosy",
                                "sigma"))
    res[[nm]] <- list(par = a[, 1], sigma = a[, 2])
  }
  saveRDS(res, file.path(out, "agree-up-brms.rds"))
} else {
  .libPaths(ulib)
  rd <- function(a) readRDS(file.path(out, paste0("agree-up-", a, ".rds")))
  B <- rd("before"); A <- rd("after"); R <- rd("brms")
  pt <- function(s) dt(s / 2.5, 3) / 2.5
  ess <- function(x) posterior::ess_bulk(matrix(x, ncol = CH))
  stat <- function(x, w = rep(1, length(x)), essv) {
    mu <- sum(w * x) / sum(w)
    sdv <- sqrt(sum(w * (x - mu)^2) / sum(w))
    o <- order(x); cw <- cumsum(w[o]) / sum(w)
    q <- vapply(c(0.05, 0.5, 0.95), function(p) x[o][which(cw >= p)[1]], 0)
    c(mean = mu, sd = sdv, q05 = q[1], q50 = q[2], q95 = q[3], min = min(x),
      ess = essv, mcse = sdv / sqrt(essv))
  }
  for (nm in names(cases)) {
    cs <- cases[[nm]]
    w <- rep(1, length(R[[nm]]$par))
    if (cs$m == "ar" && !identical(commandArgs(TRUE)[2], "unweighted")) {
      s <- sqrt(1 - R[[nm]]$par^2)
      w <- pt(R[[nm]]$sigma / s) / (pt(R[[nm]]$sigma) * s)
    }
    ew <- 1 / sum((w / sum(w))^2)
    tab <- rbind(before = stat(B[[nm]], essv = ess(B[[nm]])),
                 after = stat(A[[nm]], essv = ess(A[[nm]])),
                 brms_ref = stat(R[[nm]]$par, w,
                                 essv = min(ew, ess(R[[nm]]$par))))
    cat("\n=====", nm, "\n")
    print(round(tab, 4))
    z <- function(a) (tab[a, "mean"] - tab["brms_ref", "mean"]) /
      sqrt(tab[a, "mcse"]^2 + tab["brms_ref", "mcse"]^2)
    cat(sprintf("z(after - brms) %.2f; z(before - brms) %.2f\n", z("after"),
                z("before")))
  }
}
