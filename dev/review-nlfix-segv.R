# The hurdle_cumulative + disc segfault. Which layer: brms's program, its data,
# or rstan? One step per argument, one process per step. Seed 4952.
args <- commandArgs(trailingOnly = TRUE)
step <- args[1]
.libPaths(c("C:/Users/adf44/source/r/rellib-r16", "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win", FRMTMB_STAN_CACHE = normalizePath("dev/stan-cache", winslash = "/"))
suppressPackageStartupMessages({ library(frmtmb); library(rstan) })
E <- new.env(parent = asNamespace("frmtmb"))
E$`%||%` <- function(a, b) if (is.null(a)) b else a
for (f in c("dev/equiv-lib.R", "dev/equiv-specs.R")) sys.source(f, envir = E)
d <- E$eq_dat("ordinal", seed = 4952)
set.seed(4952)
d$y <- ordered(ifelse(runif(nrow(d)) < 0.2, 0L, as.integer(d$y)), levels = 0:4)
variants <- list(
  disc_hu = brms::bf(y ~ x, disc ~ 0 + z, hu ~ z),
  disc_only = brms::bf(y ~ x, disc ~ 0 + z),
  hu_only = brms::bf(y ~ x, hu ~ z),
  plain = brms::bf(y ~ x),
  disc_const = brms::bf(y ~ x, disc = 1.5))
v <- args[2]
fb <- variants[[v]]
fam <- brms::hurdle_cumulative()
prior <- E$brms_flat_prior(fb, data = d, family = fam)
code <- brms::make_stancode(fb, data = d, family = fam, prior = prior)
sdat <- brms::standata(fb, data = d, family = fam, prior = prior)
cat("variant", v, "step", step, "\n")
cat("Y range", range(sdat$Y), " nthres", sdat$nthres, " names", paste(names(sdat), collapse = ","), "\n")
if (step == "code") {
  ln <- strsplit(as.character(code), "\n")[[1]]
  i <- grep("hurdle_cumulative", ln)
  a <- grep("real hurdle_cumulative", ln)
  cat(ln, sep = "\n")
  quit("no")
}
m <- E$brms_stan_model(code)
cat("compiled\n"); flush.console()
if (step == "logprob") {
  sf <- suppressMessages(rstan::sampling(m, data = sdat, chains = 0))
  cat("zero-chain object\n"); flush.console()
  n <- rstan::get_num_upars(sf)
  set.seed(1)
  u <- rep(0, n)
  for (k in seq_len(5)) {
    u <- rnorm(n, 0, 0.5)
    cat("log_prob", k, rstan::log_prob(sf, u, adjust_transform = FALSE), "\n"); flush.console()
  }
  cat("grad", head(rstan::grad_log_prob(sf, u), 3), "\n")
}
if (step == "sample") {
  s <- suppressWarnings(rstan::sampling(m, data = sdat, chains = 1, iter = 20, refresh = 0, seed = 1))
  cat("sampled", dim(as.matrix(s)), "\n")
}
if (step == "opt") {
  s <- rstan::optimizing(m, data = sdat, seed = 1, iter = 50, as_vector = FALSE)
  cat("optimizing value", s$value, s$return_code, "\n")
}
if (step == "brm") {
  f <- brms::brm(fb, data = d, family = fam, chains = 1, iter = 20, refresh = 0, backend = "rstan", seed = 1, prior = prior)
  print(dim(as.matrix(f)))
}
