# Lane sampfront: one reproducer per item of the brief.
#   Rscript dev/sampfront-repro.R [base] > dev/sampfront-out/repro-<arm>.txt
# "base" drops the lane library (rellib-r20, frmtmb 0.81.0). Data seed 11
# (items 1-5, 7) and 77 (item 6); sampler seed 1; lap_pair() seed 1.
args <- commandArgs(TRUE)
base <- length(args) && identical(args[1], "base")
.libPaths(c(if (!base) "C:/Users/adf44/source/r/wt-sampfront-lib",
            "C:/Users/adf44/source/r/rellib-r20",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win",
           FRMTMB_STAN_CACHE =
             "C:/Users/adf44/source/r/frmtmb-wt-sampfront/dev/stan-cache")
suppressPackageStartupMessages({library(frmtmb); library(frmtmb.sample)})
cat("frmtmb", format(packageVersion("frmtmb")), find.package("frmtmb"),
    "\nfrmtmb.sample", format(packageVersion("frmtmb.sample")),
    find.package("frmtmb.sample"), "\n")
short <- function(m) substr(gsub("\\s+", " ", m), 1L, 330L)
go <- function(lab, expr) {
  cat("\n-- ", lab, "\n", sep = "")
  ws <- character()
  r <- withCallingHandlers(
    tryCatch(expr, error = function(e) {
      m <- conditionMessage(e)
      # the reason a subset is refused ends the message
      why <- if (grepl("This newdata", m, fixed = TRUE)) {
        paste("\n  ... |",
              sub(".*(This newdata)", "\\1", gsub("\\s+", " ", m)))
      }
      structure(paste("ERROR:", short(m), why), class = "err")
    }),
    warning = function(w) {
      ws <<- c(ws, short(conditionMessage(w)))
      invokeRestart("muffleWarning")
    },
    message = function(m) invokeRestart("muffleMessage"))
  for (w in unique(ws)) {
    if (!grepl("R-hat|Effective Samples|divergent|pairs\\(\\)", w)) {
      cat("  WARNING:", w, "\n")
    }
  }
  if (inherits(r, "err")) cat(" ", r, "\n") else print(r)
  invisible(r)
}
smp <- function(...) {
  suppressWarnings(suppressMessages(
    frm_sample(..., chains = 1, iter = 200, refresh = 0, seed = 1)))
}

set.seed(11)
d <- data.frame(g = factor(rep(1:6, each = 8)), x = rnorm(48),
                t = rep(1:8, 6))
d$y <- 1 + 0.5 * d$x + rnorm(6)[d$g] * 0.7 + rnorm(48)
d$z <- -1 + 0.3 * d$x + rnorm(48)

cat("\n==== 1. 4717: newdata that is a subset of the fitted rows\n")
ds <- smp(bf(y ~ x + (1 | g)), data = d)
ll <- log_lik(ds)
go("log_lik(newdata = d[1:10, ]) vs fitted columns", {
  ls <- log_lik(ds, newdata = d[1:10, ])
  c(dim = dim(ls), identical = identical(unclass(ls)[, 1:10],
                                         unclass(ll)[, 1:10]))
})
go("log_lik(newdata = d[c(5, 3, 20), ]) identical to columns 5, 3, 20", {
  ls <- log_lik(ds, newdata = d[c(5, 3, 20), ])
  identical(unname(unclass(ls)[, 1:3]), unname(unclass(ll)[, c(5, 3, 20)]))
})
go("loo(newdata = d[1:10, ]) pointwise vs loo() rows 1:10", {
  l0 <- loo(ds)
  l1 <- loo(ds, newdata = d[1:10, ])
  c(n = nrow(l1$pointwise), maxdiff = max(abs(
    l1$pointwise[, "elpd_loo"] - l0$pointwise[1:10, "elpd_loo"])))
})
go("waic(newdata = d[1:10, ])", {
  w1 <- waic(ds, newdata = d[1:10, ])
  nrow(w1$pointwise)
})
go("psis(ds, d[1:10, ])", dim(psis(ds, d[1:10, ])$log_weights))
dmv <- d[1:10, ]
dmv$x[2] <- dmv$x[2] + 1
go("a moved x in a subset row is refused", log_lik(ds, newdata = dmv))
dnr <- d[11:20, ]
rownames(dnr) <- NULL
go("row names that do not match the values are refused",
   log_lik(ds, newdata = dnr))
dsa <- smp(bf(y ~ x + ar(time = t, gr = g)), data = d)
go("ar(): a subset is refused, naming the term",
   log_lik(dsa, newdata = d[1:10, ]))

cat("\n==== 2. 4584: log_lik(combine = FALSE) without rescor\n")
dmv2 <- smp(mvbf(bf(y ~ x), bf(z ~ x), rescor = FALSE), data = d)
go("log_lik(combine = FALSE)", {
  a <- log_lik(dmv2, combine = FALSE)
  list(dim = dim(a), dimnames = dimnames(a),
       sum_equals = max(abs(apply(a, c(1, 2), sum) - log_lik(dmv2))))
})

cat("\n==== 3. 3215: frm_sample(formula, gen_extreme_value())\n")
set.seed(4)
rgev <- function(n, mu, s, xi) mu + s * ((-log(runif(n)))^(-xi) - 1) / xi
dv <- data.frame(y = rgev(100, 0, 1, 0.1))
go("formula route", {
  dg <- smp(y ~ 1, family = gen_extreme_value(), data = dv)
  round(colMeans(as.matrix(dg))[c("b_Intercept", "sigma", "xi")], 3)
})

cat("\n==== 4. 3594: posterior_linpred() with an unknown argument\n")
go("posterior_linpred(ds, junk = 1)", dim(posterior_linpred(ds, junk = 1)))

cat("\n==== 5. 4370: prior_samples()\n")
go("prior_samples exists", exists("prior_samples"))
go("prior_samples(ds) vs prior_draws(ds)", {
  dp <- smp(bf(y ~ x), data = d, prior = set_prior("normal(0, 2)", "b"),
            sample_prior = "yes")
  identical(prior_samples(dp), prior_draws(dp))
})

cat("\n==== 7a. 1547: as.array(chains =)\n")
go("as.array(ds, chains = 1)", dim(as.array(ds, chains = 1)))
cat("\n==== 7b. 1655: drop_unused_levels\n")
go("frm_sample(drop_unused_levels = FALSE)", {
  d3 <- d
  d3$g <- factor(d3$g, levels = 1:7)
  x <- smp(bf(y ~ x + (1 | g)), data = d3, drop_unused_levels = FALSE)
  c(levels_in_data = nlevels(x$data$g),
    r_g = sum(grepl("^r_g\\[", variables(x))))
})
cat("\n==== 7c. 1559: posterior_average(), pp_average(), model_weights()\n")
ds0 <- smp(bf(y ~ 1 + (1 | g)), data = d)
go("posterior_average missing variable",
   posterior_average(ds, ds0, variable = c("b_x", "nope"),
                     weights = c(0.5, 0.5)))
go("posterior_average missing variable, missing = 0",
   posterior_average(ds, ds0, variable = "nope", missing = 0,
                     weights = c(0.5, 0.5)))
go("pp_average exported",
   "pp_average" %in% getNamespaceExports("frmtmb.sample"))
go("model_weights exported",
   "model_weights" %in% getNamespaceExports("frmtmb.sample"))
cat("\n==== 7d. 1533: REML with laplace = TRUE on a prior route\n")
go("frm_sample(REML fit, laplace = TRUE) conditions", {
  fr <- frm(bf(y ~ x + (1 | g)), data = d, REML = TRUE)
  ms <- character()
  withCallingHandlers(
    frm_sample(fr, laplace = TRUE, chains = 1, iter = 100, refresh = 0,
               seed = 1),
    message = function(m) {
      ms <<- c(ms, short(conditionMessage(m)))
      invokeRestart("muffleMessage")
    },
    warning = function(w) invokeRestart("muffleWarning"))
  grep("REML", ms, value = TRUE)
})

cat("\n==== 6. 1456, 1530: the laplace probe and watch\n")
lap_pair <- function(fit, n = 6L, seed = 1L, tweak = NULL) {
  tpl <- fit$frame[["par_template"]]
  est <- unlist(lapply(names(tpl), function(cp) fit$estimates[[cp]]))
  lab <- frmtmb::brms_par_labels(fit)
  set.seed(seed)
  M <- matrix(rep(est, each = n) + stats::rnorm(n * length(est), 0, 0.05),
              n, dimnames = list(NULL, lab))
  if (!is.null(tweak)) M <- tweak(M)
  M <- cbind(frmtmb.sample:::draws_to_natural(M, fit), lp__ = 0)
  inner <- setdiff(lab, frmtmb::brms_par_labels(fit, include_random = FALSE))
  full <- structure(list(stanfit = NULL, draws = M, fit = fit),
                    class = "frmtmb_draws")
  lap <- full
  lap$draws <- M[, setdiff(colnames(M), inner), drop = FALSE]
  list(full = full, lap = lap)
}
set.seed(77)
G <- 8; n <- 10
dd <- data.frame(g = factor(rep(seq_len(G), each = n)))
dd$x <- rnorm(nrow(dd))
u <- rnorm(G, 0, 0.6)[dd$g]
dd$y <- 0.5 + 0.4 * dd$x + u + rnorm(nrow(dd), 0, 0.7)
dd$yn <- 0.5 * dd$x + exp(0.3 + u)^0.8 + rnorm(nrow(dd), 0, 0.3)
q <- function(...) suppressWarnings(suppressMessages(frm(...)))
fnl <- q(bf(yn ~ log(c1) * x + exp(a)^k, c1 ~ 1, a ~ 1 + (1 | g), k ~ 1,
            nl = TRUE), family = gaussian(), data = dd,
         prior = set_prior("normal(1.6, 1)", nlpar = "c1"))
tw <- function(M) {
  M[1, "b_c1_Intercept"] <- -1
  M[1, "b_k_Intercept"] <- 0
  M
}
p <- lap_pair(fnl, tweak = tw)
go("1456 posterior_epred(laplace), non-finite cells", {
  a <- posterior_epred(p$lap)
  c(nonfinite = sum(!is.finite(a)), cells = length(a))
})
go("1456 full draws non-finite cells", {
  b <- posterior_epred(p$full)
  c(nonfinite = sum(!is.finite(b)), cells = length(b))
})
fx <- q(bf(y ~ x + (0 + x | g)), family = gaussian(), data = dd)
px <- lap_pair(fx)
nd0 <- data.frame(x = c(0, 0), g = factor(c(1, 2), levels = 1:G))
go("1530 posterior_epred(laplace, newdata x = 0) vs full", {
  a <- posterior_epred(px$lap, newdata = nd0)
  b <- posterior_epred(px$full, newdata = nd0)
  identical(a, b)
})
go("1530 control: x = 1 still refused",
   posterior_epred(px$lap, newdata = transform(nd0, x = 1)))
cat("\nDONE\n")
