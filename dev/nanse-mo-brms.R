# Defect 8 against brms 2.23.0: `ls ~ mo(income) * age` on given seeds
# of brms_monotonic's data code. brms with its default priors (flat on
# b, Dirichlet(1) on each simplex) sampled 4 x 2000; frmtmb's standard
# errors from the build under test. Stan programs come from
# dev/stan-cache, keyed as tests/testthat/helper-brms.R keys them.
#   Rscript dev/nanse-mo-brms.R [lib] [seeds]
args <- commandArgs(trailingOnly = TRUE)
lib <- if (length(args)) args[1] else "C:/Users/adf44/source/r/wt-nanse-lib"
seeds <- if (length(args) > 1) eval(parse(text = args[2])) else c(1, 7, 8, 11)
.libPaths(unique(c(lib, "C:/Users/adf44/source/r/rellib-r5",
                   "C:/Users/adf44/AppData/Local/R/win-library/4.6")))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressMessages({
  library(frmtmb)
  library(brms)
  library(rstan)
})
cat("frmtmb", as.character(packageVersion("frmtmb")), "from",
    find.package("frmtmb"), "; brms", as.character(packageVersion("brms")),
    "\n")
cache <- normalizePath("dev/stan-cache", winslash = "/")
stan_cached <- function(code) {
  f <- tempfile(fileext = ".stan")
  writeLines(c(code, paste0("// rstan ", packageVersion("rstan"))), f)
  key <- unname(tools::md5sum(f))
  path <- file.path(cache, paste0(key, ".rds"))
  if (file.exists(path)) return(readRDS(path))
  m <- rstan::stan_model(model_code = code, save_dso = TRUE)
  saveRDS(m, path)
  m
}
mk <- function(s) {
  set.seed(s)
  lev <- c("below_20", "20_to_40", "40_to_100", "greater_100")
  income <- factor(sample(lev, 100, TRUE), levels = lev, ordered = TRUE)
  ls <- c(30, 60, 70, 75)[income] + rnorm(100, sd = 7)
  d <- data.frame(income, ls)
  d$age <- rnorm(100, mean = 40, sd = 10)
  d
}
for (s in seeds) {
  d <- mk(s)
  f <- suppressWarnings(frm(ls ~ mo(income) * age, data = d))
  fe <- fixef(f)
  bf0 <- brm(ls ~ mo(income) * age, data = d, empty = TRUE)
  mod <- stan_cached(stancode(bf0))
  sf <- rstan::sampling(mod, data = standata(bf0), chains = 4,
                        iter = 2000, seed = s, refresh = 0)
  bf0$fit <- sf
  bf0 <- brms::rename_pars(bf0)
  bfe <- fixef(bf0)
  sig <- posterior::as_draws_df(bf0)$sigma
  rn <- rownames(fe)
  cat(sprintf("\n== seed %d (frmtmb code %d, se_lost: %s)\n", s,
              f$opt$convergence,
              paste(names(frmtmb:::sdr_of(f)$se_lost), collapse = ",")))
  cat(sprintf("  %-14s %10s %10s %10s %10s %7s\n", "coef", "frm est",
              "brms mean", "frm SE", "brms sd", "ratio"))
  for (r in rn) {
    j <- match(r, rownames(bfe))
    cat(sprintf("  %-14s %10.4g %10.4g %10.4g %10.4g %7.3f\n", r,
                fe[r, "Estimate"], bfe[j, "Estimate"], fe[r, "Est.Error"],
                bfe[j, "Est.Error"], fe[r, "Est.Error"] / bfe[j, "Est.Error"]))
  }
  ss <- summary(f)$spec_pars
  cat(sprintf("  %-14s %10.4g %10.4g %10.4g %10.4g %7.3f\n", "sigma",
              ss["sigma", "Estimate"], mean(sig), ss["sigma", "Est.Error"],
              stats::sd(sig), ss["sigma", "Est.Error"] / stats::sd(sig)))
  sm <- summary(sf, pars = c("simo_1", "simo_2"))$summary
  cat("  brms simplexes (mean, sd):\n")
  print(round(sm[, c("mean", "sd")], 3))
  z <- f$estimates[grepl("^zeta", names(f$estimates))]
  cat("  frmtmb simplexes:\n")
  for (zn in names(z)) {
    x <- exp(c(0, z[[zn]]))
    cat("   ", zn, signif(x / sum(x), 3), "\n")
  }
  cat("  max Rhat", round(max(sm[, "Rhat"]), 3), "\n")
}
