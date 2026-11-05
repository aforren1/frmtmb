# Measure every core autocor term x core family class against brms 2.23.0.
#
# One row per (term, family): what brms's generator does (the error in its
# own words, if any), whether the program compiles through the suite's
# content-addressed cache, whether it evaluates (log_prob finite at one
# random point), and what frm() does on the build under test.
#
# Usage: Rscript dev/latres-grid.R <worker> <n_workers> <lib> <outdir>
#   <lib> is the frmtmb library under test, ";"-separated, tried first.
# Seed: set.seed(20261029) for the data. Compiling every family is not
# needed: a family's class is represented by the `compile = TRUE` rows,
# and the others run the generator and frm() only.
args <- commandArgs(TRUE)
worker <- as.integer(args[[1]])
workers <- as.integer(args[[2]])
lib <- strsplit(args[[3]], ";", fixed = TRUE)[[1]]
outdir <- args[[4]]
.libPaths(c(lib, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
root <- "C:/Users/adf44/source/r/frmtmb-wt-latres"
Sys.setenv(FRMTMB_STAN_CACHE = file.path(root, "dev/stan-cache"))
suppressPackageStartupMessages({
  library(frmtmb)
  library(rstan)
})
frm_ns <- function(txt) eval(parse(text = txt), envir = asNamespace("frmtmb"))
brm_ns <- function(txt) eval(parse(text = txt), envir = asNamespace("brms"))
sys.source(file.path(root, "tests/testthat/helper-brms.R"),
           envir = environment())

terms <- c(
  ar1t = "ar(t, g, p = 1, cov = TRUE)", ar1f = "ar(t, g, p = 1)",
  ar2t = "ar(t, g, p = 2, cov = TRUE)", ar2f = "ar(t, g, p = 2)",
  ma1t = "ma(t, g, q = 1, cov = TRUE)", ma1f = "ma(t, g, q = 1)",
  ma2t = "ma(t, g, q = 2, cov = TRUE)", ma2f = "ma(t, g, q = 2)",
  arma1t = "arma(t, g, p = 1, q = 1, cov = TRUE)",
  arma1f = "arma(t, g, p = 1, q = 1)",
  arma2t = "arma(t, g, p = 2, q = 2, cov = TRUE)",
  arma2f = "arma(t, g, p = 2, q = 2)",
  cosy = "cosy(t, g)", unstr = "unstr(t, g)", fcor = "fcor(M)")

# class, family text, whether the program is compiled, and the response
fams <- list(
  list("natural", "gaussian()", TRUE, "real"),
  list("natural", "student()", TRUE, "real"),
  list("count", "poisson()", TRUE, "count"),
  list("count", "negbinomial()", TRUE, "count"),
  list("count", "geometric()", FALSE, "count"),
  list("count", "zero_inflated_poisson()", TRUE, "count"),
  list("count", "hurdle_poisson()", TRUE, "count"),
  list("binary", "bernoulli()", TRUE, "bin"),
  list("binomial", "binomial()", TRUE, "trials"),
  list("binomial", "beta_binomial()", FALSE, "trials"),
  list("positive", "Gamma(link = 'log')", TRUE, "pos"),
  list("positive", "Gamma()", TRUE, "pos"),
  list("positive", "lognormal()", TRUE, "pos"),
  list("positive", "weibull()", FALSE, "pos"),
  list("positive", "exponential()", FALSE, "pos"),
  list("positive", "inverse.gaussian()", FALSE, "pos"),
  list("real_other", "skew_normal()", FALSE, "real"),
  list("real_other", "exgaussian()", FALSE, "pos"),
  list("real_other", "asym_laplace()", FALSE, "real"),
  list("bounded", "Beta()", TRUE, "unit"),
  list("bounded", "zero_one_inflated_beta()", FALSE, "unit"),
  list("ordinal", "cumulative()", TRUE, "ord"),
  list("ordinal", "sratio()", FALSE, "ord"),
  list("ordinal", "cratio()", FALSE, "ord"),
  list("ordinal", "acat()", FALSE, "ord"),
  list("categorical", "categorical()", TRUE, "cat"),
  list("multi_mu", "multinomial()", TRUE, "mult"),
  list("multi_mu", "dirichlet()", FALSE, "dir"),
  list("mixture", "mixture(poisson(), poisson())", TRUE, "count"),
  list("mixture", "mixture(gaussian(), gaussian())", FALSE, "real")
)

set.seed(20261029)
d <- expand.grid(t = 1:5, g = factor(1:8))
n <- nrow(d)
d$x <- rnorm(n)
eta <- 0.3 + 0.4 * d$x + rep(rnorm(8, 0, 0.4), each = 5)
d$real <- eta + rnorm(n)
d$count <- rpois(n, exp(eta))
d$bin <- rbinom(n, 1, plogis(eta))
d$trials <- rbinom(n, 10, plogis(eta))
d$tr <- 10L
d$pos <- rgamma(n, 2, 2 / exp(eta))
d$unit <- rbeta(n, 2 * plogis(eta), 2 * (1 - plogis(eta)))
d$ord <- factor(findInterval(eta + rlogis(n), c(-0.5, 0.5, 1.5)) + 1L,
                ordered = TRUE)
d$cat <- factor(c("a", "b", "c")[findInterval(eta + rlogis(n), c(0, 1)) + 1L])
m <- t(rmultinom(n, 10, c(0.2, 0.3, 0.5)))
colnames(m) <- c("c1", "c2", "c3")
d$mult <- m
dm <- matrix(rgamma(3 * n, 2), n, 3)
dm <- dm / rowSums(dm)
colnames(dm) <- c("d1", "d2", "d3")
d$dir <- dm
M <- diag(n)
for (i in seq_len(n - 1L)) M[i, i + 1L] <- M[i + 1L, i] <- 0.3

lhs <- function(resp) switch(resp,
  trials = "trials | trials(tr)", mult = "mult | trials(10)", resp)

grid <- expand.grid(term = names(terms), fam = seq_along(fams),
                    stringsAsFactors = FALSE)
grid <- grid[(seq_len(nrow(grid)) - 1L) %% workers == worker - 1L, ]

clip <- function(x) {
  x <- gsub("[\t\r\n]+", " ", paste(x, collapse = " "))
  if (nchar(x) > 400) paste0(substr(x, 1, 400), "...") else x
}

out <- file.path(outdir, sprintf("grid-%02d.tsv", worker))
hdr <- c("term", "class", "family", "brms", "brms_error", "compiled",
         "evaluates", "latent", "frm", "frm_error")
cat(paste(hdr, collapse = "\t"), "\n", sep = "", file = out)

for (k in seq_len(nrow(grid))) {
  tm <- grid$term[k]
  fm <- fams[[grid$fam[k]]]
  f_txt <- paste0(lhs(fm[[4]]), " ~ x + ", terms[[tm]])
  fam_b <- brm_ns(fm[[2]])
  bf_b <- brms::bf(as.formula(f_txt))
  rec <- list(brms = NA, brms_error = "", compiled = NA, evaluates = NA,
              latent = NA, frm = NA, frm_error = "")
  dd <- d
  code <- tryCatch({
    pr <- brms_flat_prior(bf_b, data = dd, family = fam_b,
                          data2 = list(M = M))
    brms::make_stancode(bf_b, data = dd, family = fam_b, prior = pr,
                        data2 = list(M = M))
  }, error = function(e) e)
  if (inherits(code, "error")) {
    rec$brms <- "refuses"
    rec$brms_error <- clip(conditionMessage(code))
  } else {
    rec$brms <- "builds"
    rec$latent <- grepl("zerr", code)
    if (isTRUE(fm[[3]])) {
      mod <- tryCatch(brms_stan_model(code), error = function(e) e)
      if (inherits(mod, "error")) {
        rec$compiled <- FALSE
        rec$brms_error <- clip(conditionMessage(mod))
      } else {
        rec$compiled <- TRUE
        sdat <- brms::make_standata(bf_b, data = dd, family = fam_b,
                                    prior = pr, data2 = list(M = M))
        lp <- tryCatch({
          sf <- suppressMessages(rstan::sampling(mod, data = sdat,
                                                 chains = 0))
          set.seed(1)
          u <- stats::rnorm(rstan::get_num_upars(sf), 0, 0.3)
          rstan::log_prob(sf, u, adjust_transform = FALSE)
        }, error = function(e) e)
        rec$evaluates <- if (inherits(lp, "error")) {
          clip(conditionMessage(lp))
        } else {
          is.finite(lp)
        }
      }
    }
  }
  fit <- tryCatch({
    fam_f <- frm_ns(fm[[2]])
    suppressWarnings(suppressMessages(frmtmb::frm(
      frmtmb::bf(as.formula(f_txt)), family = fam_f, data = dd,
      data2 = list(M = M))))
  }, error = function(e) e)
  if (inherits(fit, "error")) {
    rec$frm <- "refuses"
    rec$frm_error <- clip(conditionMessage(fit))
  } else {
    rec$frm <- sprintf("fits logLik %.6f", as.numeric(logLik(fit)))
  }
  cat(paste(c(tm, fm[[1]], fm[[2]], unlist(rec)), collapse = "\t"), "\n",
      sep = "", file = out, append = TRUE)
}
