# Lane sampfix, script 12 (R3): on every watched path, is a read of a
# value laplace draws do not hold NA (or NaN), and never +-Inf?
#
#   Rscript dev/sampfix-12-fillpaths.R      (data seed 77, draws seed 1)
#
# The probe is switched off and the watch replaced by a recorder, so each
# call runs every draw with the integrated values at NA, as it would
# without the refusal. For each watched value v the recorder counts cells
# that are Inf, and cells that are non-finite in the laplace result but
# finite in the full draws' result at the same draw (the reads). A read
# that surfaces as Inf would show as inf_read > 0.
.libPaths(c("C:/Users/adf44/source/r/wt-sampfix-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(frmtmb); library(frmtmb.sample); library(testthat)})
src <- readLines("C:/Users/adf44/source/r/frmtmb-wt-sampfix/dev/sampfix-rev-r2-01-watch.R")
eval(parse(text = src[seq_len(grep("^res <- function", src)[1L] - 1L)]))

rec <- new.env()
record <- function(x, what, at_na = NULL) {
  function(v) {
    v <- suppressWarnings(as.numeric(unlist(v, use.names = FALSE)))
    rec$vals[[length(rec$vals) + 1L]] <- v
    invisible(NULL)
  }
}
run <- function(obj, f) {
  rec$vals <- list()
  set.seed(3)
  r <- tryCatch(suppressWarnings(suppressMessages(
    with_mocked_bindings(f(obj),
      draws_laplace_probe = function(...) invisible(NULL),
      draws_laplace_watch = record, .package = "frmtmb.sample"))),
    error = function(e) e)
  list(r = r, vals = rec$vals)
}
tally <- function(lab, p, f) {
  a <- run(p$lap, f)
  b <- run(p$full, f)
  if (inherits(a$r, "error") || !length(a$vals)) {
    cat(sprintf("%-46s %s\n", lab, if (inherits(a$r, "error"))
      paste("ERROR", substr(conditionMessage(a$r), 1, 60)) else
        "no watched value"))
    return(invisible(NULL))
  }
  la <- unlist(a$vals); lb <- unlist(b$vals)
  read <- !is.finite(la) & is.finite(lb)
  cat(sprintf("%-46s cells %6d  reads %6d  as NA/NaN %6d  as Inf %d\n",
              lab, length(la), sum(read), sum(read & is.na(la)),
              sum(read & is.infinite(la))))
}

set.seed(77)
G <- 8; n <- 10
dd <- data.frame(g = factor(rep(seq_len(G), each = n)))
dd$x <- rnorm(nrow(dd))
u <- rnorm(G, 0, 0.6)[dd$g]
dd$y <- 0.5 + 0.4 * dd$x + u + rnorm(nrow(dd), 0, 0.7)
dd$yp <- rpois(nrow(dd), exp(0.3 + 0.4 * dd$x + u))
dd$yb <- rbinom(nrow(dd), 1, plogis(0.2 + 0.8 * dd$x + u))
dd$pos <- exp(0.2 + 0.1 * dd$x + 0.3 * u + rnorm(nrow(dd), 0, 0.2))
dd$o <- factor(cut(0.8 * dd$x + u + rlogis(nrow(dd)),
                   c(-Inf, -0.7, 0.6, Inf), labels = FALSE), ordered = TRUE)
dd$yn <- 0.5 * dd$x + exp(0.3 + u)^0.8 + rnorm(nrow(dd), 0, 0.3)
dd$ym <- c(rnorm(40, -2), rnorm(40, 3)) + u
dm <- dd
dm$x[c(3, 11, 19)] <- NA
q <- function(...) suppressWarnings(suppressMessages(frm(...)))
nd <- data.frame(x = c(-1, 0, 1), g = factor(c(1, 2, 3), levels = 1:G))

ep <- function(d) posterior_epred(d)
lp <- function(d) posterior_linpred(d)
pp <- function(d) posterior_predict(d)
epn <- function(d) posterior_epred(d, newdata = nd)
ppn <- function(d) posterior_predict(d, newdata = nd)

mods <- list(
  gaussian = q(bf(y ~ x + (1 | g)), family = gaussian(), data = dd),
  poisson = q(bf(yp ~ x + (1 | g)), family = poisson(), data = dd),
  bernoulli = q(bf(yb ~ x + (1 | g)), family = bernoulli(), data = dd),
  lognormal_mu_sigma = q(bf(pos ~ x + (1 | g), sigma ~ (1 | g)),
                         family = lognormal(), data = dd),
  cumulative = q(bf(o ~ x + (1 | g)), family = cumulative(), data = dd),
  smooth = q(bf(y ~ s(x, k = 5)), family = gaussian(), data = dd),
  nl_exp = q(bf(yn ~ c0 + exp(a)^k, c0 ~ 1 + x, a ~ 1 + (1 | g), k ~ 1,
                nl = TRUE), family = gaussian(), data = dd))
for (nm in names(mods)) {
  cat("== ", nm, "\n", sep = "")
  p <- lap_pair(mods[[nm]])
  tally("posterior_epred", p, ep)
  tally("posterior_linpred", p, lp)
  tally("posterior_predict (dpars)", p, pp)
  if (nm != "smooth") {
    tally("posterior_epred newdata", p, epn)
    tally("posterior_predict newdata (dpars)", p, ppn)
  }
}
cat("== smooth + (1 | g): conditional_effects at re_formula = NULL\n")
p <- lap_pair(mods$smooth <- q(bf(y ~ s(x, k = 5) + (1 | g)),
                               family = gaussian(), data = dd))
tally("conditional_effects re_formula = NULL", p, function(d)
  conditional_effects(d, effects = "x", resolution = 5, re_formula = NULL,
                      seed = 2))
tally("conditional_effects default", p, function(d)
  conditional_effects(d, effects = "x", resolution = 5))
cat("== mi() predictor\n")
p <- lap_pair(q(bf(y ~ mi(x)) + bf(x | mi() ~ 1) + set_rescor(FALSE),
                family = gaussian(), data = dm))
tally("posterior_epred resp y", p, function(d) posterior_epred(d, resp = "y"))
tally("posterior_predict resp y (dpars)", p,
      function(d) posterior_predict(d, resp = "y"))
cat("== grouped mixture\n")
p <- lap_pair(q(bf(ym ~ 1 + (1 | g)),
                family = frmtmb::mixture(gaussian(), gaussian()), data = dd))
tally("pp_mixture", p, function(d) pp_mixture(d, summary = FALSE))
tally("posterior_epred", p, ep)
cat("== hypothesis through the fit (sd_ of a smooth + (1 | g) model)\n")
p <- lap_pair(q(bf(y ~ x + (1 | g)), family = gaussian(), data = dd))
tally("hypothesis sd_g__Intercept", p, function(d)
  hypothesis(d, "sd_g__Intercept > 0", class = NULL))
cat("DONE\n")
