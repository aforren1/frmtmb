# Review of lane ordstart: where does brms 2.23.0's posterior mode (its
# default priors, rstan::optimizing from 10 random inits) lie relative to
# the lane's answer and 0.81.0's answer on a seed the lane moved?
#   Rscript ordstart-review-brms.R family seed
args <- commandArgs(trailingOnly = TRUE)
famname <- args[1]
s <- as.integer(args[2])
.libPaths(c("C:/Users/adf44/source/r/wt-ordstart-lib",
            "C:/Users/adf44/source/r/rellib-r20",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressPackageStartupMessages({
  library(frmtmb)
  library(brms)
  library(rstan)
})
cat("family", famname, "seed", s, "frmtmb", find.package("frmtmb"),
    "brms", format(packageVersion("brms")), "\n")
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
famf <- switch(famname,
  cum_cum = frmtmb::mixture(frmtmb::cumulative(), frmtmb::cumulative()),
  cum_sratio = frmtmb::mixture(frmtmb::cumulative(), frmtmb::sratio()))
famb <- suppressMessages(switch(famname,
  cum_cum = brms::mixture(brms::cumulative(), brms::cumulative()),
  cum_sratio = brms::mixture(brms::cumulative(), brms::sratio())))
d <- mk(s, mix = TRUE)
q <- function(e) suppressWarnings(suppressMessages(e))
f_lane <- q(frmtmb::frm(frmtmb::bf(y ~ x + z), family = famf, data = d))
f_base <- q(frmtmb::frm(frmtmb::bf(y ~ x + z), family = famf, data = d,
                        control = frmtmb::frmtmb_control(starts = 0)))
fb <- brms::bf(y ~ x + z)
code_def <- q(as.character(brms::stancode(fb, data = d, family = famb)))
sdat_def <- q(brms::standata(fb, data = d, family = famb))
pr <- q(brms::get_prior(fb, data = d, family = famb))
pr$prior <- ""
code_flat <- q(as.character(brms::stancode(fb, data = d, family = famb,
                                           prior = pr)))
sdat_flat <- q(brms::standata(fb, data = d, family = famb, prior = pr))
mod_def <- stan_model_cached(code_def)
mod_flat <- stan_model_cached(code_flat)
sf_def <- q(rstan::sampling(mod_def, data = sdat_def, chains = 0))
sf_flat <- q(rstan::sampling(mod_flat, data = sdat_flat, chains = 0))
pars <- omb_par_names(code_flat)
at <- function(pl) {
  pl <- pl[pars]
  c(lpost = rstan::log_prob(sf_def, rstan::unconstrain_pars(sf_def, pl),
                            adjust_transform = FALSE),
    ll = rstan::log_prob(sf_flat, rstan::unconstrain_pars(sf_flat, pl),
                         adjust_transform = FALSE) - lgamma(2))
}
# the parameter vector, label-invariant for two components of one family:
# the smaller of the distances under the identity and under the swap
vec <- function(pl) {
  list(c1 = c(pl$b_mu1, pl$Intercept_mu1, pl$theta1),
       c2 = c(pl$b_mu2, pl$Intercept_mu2, pl$theta2))
}
dist <- function(a, b) {
  va <- vec(a)
  vb <- vec(b)
  d1 <- sqrt(sum((va$c1 - vb$c1)^2) + sum((va$c2 - vb$c2)^2))
  if (famname != "cum_cum") return(d1)
  d2 <- sqrt(sum((va$c1 - vb$c2)^2) + sum((va$c2 - vb$c1)^2))
  min(d1, d2)
}
maxb <- function(pl) {
  max(abs(unlist(pl[grepl("^(b|Intercept)", names(pl))])))
}
pl_lane <- omb_translate(f_lane, f_lane$estimates, code_flat, sdat_flat)
pl_base <- omb_translate(f_base, f_base$estimates, code_flat, sdat_flat)
for (nm in c("lane", "base")) {
  pl <- if (nm == "lane") pl_lane else pl_base
  f <- if (nm == "lane") f_lane else f_base
  # a collapsed threshold pair is not a point of brms's ordered space
  v <- tryCatch(at(pl), error = function(e) c(lpost = NA, ll = NA))
  cat(sprintf("FRMTMB %-4s ll %.4f brms_ll %.4f lpost %.4f maxb %.2f\n", nm,
              as.numeric(logLik(f)), v[["ll"]], v[["lpost"]], maxb(pl)))
}
cat(sprintf("distance lane-base %.3f\n", dist(pl_lane, pl_base)))
for (i in 1:10) {
  o <- tryCatch(q(rstan::optimizing(mod_def, data = sdat_def,
                                    init = "random", seed = i,
                                    as_vector = FALSE, iter = 10000)),
                error = function(e) NULL)
  if (is.null(o)) {
    cat(sprintf("MODE %2d error\n", i))
    next
  }
  v <- at(o$par)
  cat(sprintf(paste("MODE %2d code %d ll %.4f lpost %.4f maxb %.2f",
                    "dist_lane %.3f dist_base %.3f\n"),
              i, o$return_code, v[["ll"]], v[["lpost"]], maxb(o$par),
              dist(o$par, pl_lane), dist(o$par, pl_base)))
}
cat("DONE\n")
