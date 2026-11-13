# Lane ordstart, item 3: what brms 2.23.0 finds on a few of the 40 fits.
# For one family and seed of dev/ordstart-starts.R's data:
# 1. frmtmb's candidates from the default and 10 jittered starts (the
#    generator of dev/optima-ordmix-starts.R): the default, the best
#    proper point and the best degenerate one (the warning's rule);
# 2. brms's log posterior under its DEFAULT priors (no Jacobian, the
#    function rstan::optimizing() maximizes) and the log-likelihood (the
#    same program with flat priors, less lgamma(K)) at each candidate;
# 3. rstan::optimizing() from 10 random starts, on the default-prior
#    program (brms's posterior mode) and on the flat one (maximum
#    likelihood by L-BFGS);
# 4. rstan::sampling(): 4 chains, brms's random inits, 1000 + 1000
#    iterations; per chain the median log-likelihood of 100 thinned
#    draws and the median largest |b|.
#   Rscript dev/ordstart-brms.R family seed out.txt
args <- commandArgs(trailingOnly = TRUE)
famname <- args[1]
s <- as.integer(args[2])
.libPaths(c("C:/Users/adf44/source/r/rellib-r20",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressPackageStartupMessages({
  library(frmtmb)
  library(brms)
  library(rstan)
})
cat("family", famname, "seed", s, "frmtmb", format(packageVersion("frmtmb")),
    find.package("frmtmb"), "brms", format(packageVersion("brms")),
    "rstan", format(packageVersion("rstan")), "\n")
src <- readLines("dev/optima-csmix.R")
eval(parse(text = src[grep("^mk <- ", src):(grep("^one <- ", src) - 1L)]))
tsrc <- readLines("tests/testthat/test-ordinal-mixture-brms.R")
eval(parse(text = tsrc[grep("^omb_par_names <- function", tsrc):
                         (grep("^omb_check <- function", tsrc) - 1L)]))
environment(omb_translate) <- asNamespace("frmtmb")
cache <- normalizePath("dev/stan-cache", winslash = "/")
stan_model_cached <- function(code) {
  f <- tempfile(fileext = ".stan")
  writeLines(c(code, paste0("// rstan ", packageVersion("rstan"))), f)
  key <- unname(tools::md5sum(f))
  path <- file.path(cache, paste0(key, ".rds"))
  if (file.exists(path)) {
    cat("program", key, "from the cache\n")
    return(readRDS(path))
  }
  cat("program", key, "compiling\n")
  mod <- rstan::stan_model(model_code = code, save_dso = TRUE)
  tmp <- paste0(path, ".part", Sys.getpid())
  saveRDS(mod, tmp)
  file.rename(tmp, path)
  mod
}
F <- asNamespace("frmtmb")
famf <- switch(famname,
  cum_cum = F$mixture(F$cumulative(), F$cumulative()),
  cum_sratio = F$mixture(F$cumulative(), F$sratio()))
famb <- suppressMessages(switch(famname,
  cum_cum = brms::mixture(brms::cumulative(), brms::cumulative()),
  cum_sratio = brms::mixture(brms::cumulative(), brms::sratio())))
d <- mk(s, mix = TRUE)
fitq <- function(start = NULL) {
  w <- character()
  f <- withCallingHandlers(
    frmtmb::frm(frmtmb::bf(y ~ x + z), family = famf, data = d,
                start = start),
    warning = function(x) {
      w <<- c(w, conditionMessage(x))
      invokeRestart("muffleWarning")
    })
  list(fit = f, ll = as.numeric(logLik(f)),
       degen = any(grepl("degenerate boundary", w)))
}
runs <- list(fitq())
tpl <- frmtmb::frm(frmtmb::bf(y ~ x + z), family = famf, data = d,
                   dry_run = "objective")$estimates
set.seed(1000 + s)
for (i in 1:10) {
  st <- lapply(tpl, function(v) v + stats::rnorm(length(v), 0, 0.5))
  dk <- grepl("^disc", names(tpl$betad))
  st$betad[dk] <- tpl$betad[dk]
  runs[[i + 1L]] <- tryCatch(fitq(st), error = function(e) NULL)
}
runs <- Filter(Negate(is.null), runs)
ll <- vapply(runs, `[[`, 0, "ll")
dg <- vapply(runs, `[[`, NA, "degen")
cand <- list(default = runs[[1]])
if (any(!dg)) cand$best_proper <- runs[[which(!dg)[which.max(ll[!dg])]]]
if (any(dg)) cand$best_degenerate <- runs[[which(dg)[which.max(ll[dg])]]]

fb <- brms::bf(y ~ x + z)
quiet <- function(e) suppressWarnings(suppressMessages(e))
code_def <- quiet(as.character(brms::stancode(fb, data = d, family = famb)))
sdat_def <- quiet(brms::standata(fb, data = d, family = famb))
pr <- quiet(brms::get_prior(fb, data = d, family = famb))
cat("brms default priors:\n")
print(pr[, c("prior", "class", "coef", "dpar")])
pr$prior <- ""
code_flat <- quiet(as.character(brms::stancode(fb, data = d, family = famb,
                                               prior = pr)))
sdat_flat <- quiet(brms::standata(fb, data = d, family = famb, prior = pr))
mod_def <- stan_model_cached(code_def)
mod_flat <- stan_model_cached(code_flat)
sf_def <- quiet(rstan::sampling(mod_def, data = sdat_def, chains = 0))
sf_flat <- quiet(rstan::sampling(mod_flat, data = sdat_flat, chains = 0))
lconst <- lgamma(2)
pars <- omb_par_names(code_flat)
# log posterior (default priors, no Jacobian) and log-likelihood at a
# constrained parameter list
at <- function(pl) {
  pl <- pl[pars]
  c(lpost = rstan::log_prob(sf_def, rstan::unconstrain_pars(sf_def, pl),
                            adjust_transform = FALSE),
    ll = rstan::log_prob(sf_flat, rstan::unconstrain_pars(sf_flat, pl),
                         adjust_transform = FALSE) - lconst)
}
maxb <- function(pl) {
  max(abs(unlist(pl[grepl("^(b|Intercept)", names(pl))])))
}
cat("\nfrmtmb candidates (ll = frmtmb logLik; brms_ll and lpost from",
    "brms's programs at the translated point)\n")
for (nm in names(cand)) {
  f <- cand[[nm]]$fit
  pl <- omb_translate(f, f$estimates, code_flat, sdat_flat)
  v <- at(pl)
  cat(sprintf("CAND %-16s ll %.4f brms_ll %.4f lpost %.4f maxb %.2f\n", nm,
              cand[[nm]]$ll, v[["ll"]], v[["lpost"]], maxb(pl)))
}
cat("\nrstan::optimizing from 10 random inits\n")
for (which in c("default", "flat")) {
  mod <- if (which == "default") mod_def else mod_flat
  sdat <- if (which == "default") sdat_def else sdat_flat
  for (i in 1:10) {
    o <- tryCatch(quiet(rstan::optimizing(mod, data = sdat, init = "random",
                                          seed = i, as_vector = FALSE,
                                          iter = 10000)),
                  error = function(e) NULL)
    if (is.null(o)) {
      cat(sprintf("OPT %-7s %2d error\n", which, i))
      next
    }
    v <- tryCatch(at(o$par), error = function(e) c(lpost = NA, ll = NA))
    cat(sprintf("OPT %-7s %2d code %d ll %.4f lpost %.4f maxb %.2f\n", which,
                i, o$return_code, v[["ll"]], v[["lpost"]], maxb(o$par)))
  }
}
cat("\nrstan::sampling, default priors, 4 chains, random inits\n")
t0 <- proc.time()[["elapsed"]]
sf <- quiet(rstan::sampling(mod_def, data = sdat_def, chains = 4,
                            iter = 2000, warmup = 1000, seed = 1,
                            refresh = 0))
cat("sampling seconds", proc.time()[["elapsed"]] - t0, "\n")
skel <- sf@inits[[1]][pars]
A <- rstan::extract(sf, pars = pars, permuted = FALSE)
for (ch in seq_len(dim(A)[2])) {
  keep <- seq(10, dim(A)[1], by = 10)
  lls <- vapply(keep, function(it) {
    pl <- utils::relist(A[it, ch, ], skel)
    pl <- lapply(seq_along(skel), function(j) {
      v <- pl[[j]]
      if (is.null(dim(skel[[j]]))) array(as.numeric(v), length(v)) else
        array(as.numeric(v), dim(skel[[j]]))
    })
    names(pl) <- pars
    at(pl)[["ll"]]
  }, 0)
  mb <- apply(A[keep, ch, grepl("^(b|Intercept)", dimnames(A)[[3]]),
                drop = FALSE], 1, function(v) max(abs(v)))
  th <- A[keep, ch, grepl("^theta", dimnames(A)[[3]]), drop = FALSE]
  cat(sprintf(paste("CHAIN %d ll median %.4f [q10 %.4f q90 %.4f] max ll",
                    "%.4f; max|b| median %.2f max %.2f; theta1 median %.3f\n"),
              ch, stats::median(lls), stats::quantile(lls, 0.1),
              stats::quantile(lls, 0.9), max(lls), stats::median(mb),
              max(mb), stats::median(th[, 1, 1])))
}
sm <- rstan::summary(sf)$summary
cat("Rhat max", max(sm[, "Rhat"], na.rm = TRUE), "; divergent",
    sum(vapply(rstan::get_sampler_params(sf, inc_warmup = FALSE),
               function(x) sum(x[, "divergent__"]), 0)), "\n")
print(round(sm[grepl("^(b_|Intercept|theta)", rownames(sm)),
               c("mean", "sd", "2.5%", "97.5%", "Rhat")], 3))
