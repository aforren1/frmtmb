# Attack 7: the grid's claim, sampled on the cells the lane's own grid only built
# ("B", no compile): make_stancode() AND rstan::stanc() on brms's side, a real fit
# on frm()'s side, not a dry run. A stanc failure does not count as a build.
source("C:/Users/adf44/source/r/frmtmb-wt-latres/dev/review-latres-common.R")
suppressPackageStartupMessages({library(brms); library(rstan); library(frmtmb.survival)})
set.seed(20261102)
d <- expand.grid(t = 1:5, g = factor(1:10))
d <- d[order(d$g, d$t), ]
n <- nrow(d)
d$x <- rnorm(n)
e <- as.vector(replicate(10, arima.sim(list(ar = .5), 5, sd = .4)))
d$cnt <- rpois(n, exp(1 + .3 * d$x + e))
d$bi <- rbinom(n, 8, plogis(.3 * d$x + e))
d$tr <- 8L
d$pos <- rgamma(n, 4, 4 / exp(.5 + .3 * d$x + e))
d$un <- pmin(pmax(rbeta(n, 4 * plogis(.3 * d$x + e), 4 * (1 - plogis(.3 * d$x + e))), 1e-3), 1 - 1e-3)
d$un01 <- d$un
d$un01[sample(n, 8)] <- 0
d$un01[sample(n, 8)] <- 1
d$gs <- .3 * d$x + e + rnorm(n, 0, .5)
d$ang <- pmin(pmax(d$gs / 3, -3), 3)
d$ord <- factor(findInterval(.3 * d$x + e + rlogis(n), c(-1, 0, 1)) + 1L, ordered = TRUE)
fam_cells <- list(
  list("geometric()", "cnt"), list("beta_binomial()", "bi | trials(tr)"), list("weibull()", "pos"),
  list("exponential()", "pos"), list("inverse.gaussian()", "pos"), list("skew_normal()", "gs"),
  list("exgaussian()", "gs"), list("zero_one_inflated_beta()", "un01"), list("sratio()", "ord"),
  list("cratio()", "ord"), list("acat()", "ord"), list("asym_laplace()", "gs"),
  list("von_mises()", "ang"), list("frechet()", "pos"), list("gen_extreme_value()", "gs"),
  list("shifted_lognormal()", "pos"), list("zero_inflated_beta()", "un"), list("xbeta()", "un"),
  list("zero_inflated_negbinomial()", "cnt"), list("hurdle_negbinomial()", "cnt"),
  list("hurdle_lognormal()", "pos"), list("zero_inflated_binomial()", "bi | trials(tr)"),
  list("cumulative('probit')", "ord"), list("negbinomial('sqrt')", "cnt"), list("poisson('identity')", "cnt"),
  list("Gamma('identity')", "pos"), list("cox()", "pos"))
terms <- c(ar1t = "ar(t, g, cov = TRUE)", ar2f = "ar(t, g, p = 2)", ma1t = "ma(t, g, cov = TRUE)",
           arma1t = "arma(t, g, cov = TRUE)", cosy = "cosy(t, g)", unstr = "unstr(t, g)")
rows <- list()
for (fc in fam_cells) for (tn in names(terms)) {
  f <- as.formula(paste(fc[[2]], "~ x +", terms[[tn]]))
  fb <- tryCatch(eval(parse(text = fc[[1]]), asNamespace("brms")), error = function(e) NULL)
  ff <- tryCatch(eval(parse(text = fc[[1]]), asNamespace("frmtmb")), error = function(e) NULL)
  b <- "?"
  bm <- ""
  if (!is.null(fb)) {
    r <- msg_of(brms::stancode(f, d, family = fb))
    if (!is.null(r$err)) {
      b <- "R"
      bm <- r$err
    } else {
      chk <- msg_of(rstan::stanc(model_code = r$val, verbose = FALSE))
      if (!is.null(chk$err)) {
        b <- "SF"
        bm <- chk$err
      } else b <- "B"
    }
  } else {
    b <- "nofam"
  }
  m <- msg_of(frm(f, data = d, family = ff))
  fr <- if (!is.null(m$err)) "r" else if (m$val$opt$convergence == 0) "f0" else "f1"
  rows[[length(rows) + 1]] <- data.frame(fam = fc[[1]], term = tn, brms = b, frm = fr,
                                         ll = if (is.null(m$err)) as.numeric(logLik(m$val)) else NA,
                                         brms_msg = substr(gsub("\\s+", " ", bm), 1, 90), frm_msg = substr(gsub("\\s+", " ", m$err %||% ""), 1, 120), stringsAsFactors = FALSE)
}
o <- do.call(rbind, rows)
print(table(o$brms, o$frm))
print(xtabs(~ fam + paste0(brms, "/", frm), o))
bad <- o[(o$brms == "B" & o$frm == "r") | (o$brms != "B" & o$frm != "r") | abs(o$ll) > 1e6 & !is.na(o$ll), ]
print(bad[, c("fam", "term", "brms", "frm", "ll", "brms_msg", "frm_msg")], row.names = FALSE, right = FALSE)
write.table(o, "dev/review-latres-out/grid.tsv", sep = "\t", row.names = FALSE)
