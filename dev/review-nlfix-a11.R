# Is the surviving ac_refit mutant equivalent? (Set REV_MUT.) And the hurdle off-by-one
# against brms's own brm(), with and without the top category.
a <- commandArgs(trailingOnly = TRUE)
step <- a[1]
if (step == "acrefit") {
  Sys.setenv(REV_MUT = a[2])
  .libPaths(c("C:/Users/adf44/source/r/rev-nlfix-lib", "C:/Users/adf44/source/r/rellib-r16", "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
  suppressPackageStartupMessages(library(frmtmb))
  set.seed(4401); n <- 200
  d <- data.frame(x = rnorm(n), g = factor(rep(1:10, each = 20)), t = rep(1:20, 10))
  d$y <- 1 + 0.5 * d$x + ave(rnorm(n, 0, 0.5), d$g, FUN = function(z) as.numeric(stats::filter(z, 0.5, "recursive")))
  fn <- suppressWarnings(frm(bf(y ~ a + b * x, a + b ~ 1, nl = TRUE, autocor = ~ ar(time = t, gr = g, cov = TRUE)), data = d, start = list(beta = c("a_(Intercept)" = 1, "b_(Intercept)" = 0.5))))
  u <- suppressWarnings(update(fn, data = d[d$g != 1, ]))
  cat("mutant", a[2], "original ll", as.numeric(logLik(fn)), " update(data =) ll", as.numeric(logLik(u)), " autocor kept:", length(u$frame[["autocor"]]) > 0, "\n")
}
if (step == "hurdle") {
  .libPaths(c("C:/Users/adf44/source/r/rellib-r16", "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
  Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win", FRMTMB_STAN_CACHE = normalizePath("dev/stan-cache", winslash = "/"))
  suppressPackageStartupMessages({ library(frmtmb); library(rstan) })
  E <- new.env(parent = asNamespace("frmtmb"))
  E$`%||%` <- function(a, b) if (is.null(a)) b else a
  for (f in c("dev/equiv-lib.R", "dev/equiv-specs.R")) sys.source(f, envir = E)
  d <- E$eq_dat("ordinal", seed = 4952)
  set.seed(4952)
  d$y <- ordered(ifelse(runif(nrow(d)) < 0.2, 0L, as.integer(d$y)), levels = 0:4)
  fb <- brms::bf(y ~ x, disc ~ 0 + z, hu ~ z)
  fam <- brms::hurdle_cumulative()
  prior <- E$brms_flat_prior(fb, data = d, family = fam)
  code <- brms::make_stancode(fb, data = d, family = fam, prior = prior)
  m <- E$brms_stan_model(code)
  sdat <- brms::make_standata(fb, data = d, family = fam, prior = prior)
  cat("top category rows:", sum(sdat$Y == sdat$nthres + 1L), " nthres", sdat$nthres, " max Y", max(sdat$Y), " code's top test: y == nthres + 2\n")
  lp <- function(sd, lab) {
    sf <- suppressMessages(rstan::sampling(m, data = sd, chains = 0))
    set.seed(1); u <- rnorm(rstan::get_num_upars(sf), 0, 0.3)
    r <- tryCatch(rstan::log_prob(sf, u, adjust_transform = FALSE), error = function(e) paste("ERROR:", substr(conditionMessage(e), 1, 170)))
    cat(lab, ":", r, "\n")
  }
  lp(sdat, "all rows")
  keep <- sdat$Y != sdat$nthres + 1L
  sd2 <- sdat; sd2$N <- sum(keep); sd2$Y <- sdat$Y[keep]; sd2$X <- sdat$X[keep, , drop = FALSE]; sd2$X_disc <- sdat$X_disc[keep, , drop = FALSE]; sd2$X_hu <- sdat$X_hu[keep, , drop = FALSE]
  lp(sd2, "without the top category")
  f <- tryCatch(suppressWarnings(suppressMessages(brms::brm(fb, data = d, family = fam, chains = 1, iter = 30, refresh = 0, backend = "rstan", seed = 1, prior = prior))), error = function(e) e)
  cat("brm() on all rows:", if (inherits(f, "error")) paste("ERROR", substr(conditionMessage(f), 1, 200)) else paste("draws", paste(dim(as.matrix(f)), collapse = "x")), "\n")
}
