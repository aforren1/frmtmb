# Short check: the mm flag hole, names' order, the bootstrap on a fit with
# dropped rows, the warning question, chain rows by reader share. Seed 11500.
source("dev/review-refuse-lib.R")
suppressPackageStartupMessages(library(frmtmb.sample))
quiet <- function(e) suppressMessages(suppressWarnings(e))
set.seed(11500)
n <- 120; lev <- paste0("L", 1:8)
dm <- data.frame(g1 = sample(lev, n, TRUE), g2 = sample(lev, n, TRUE), x = rnorm(n), w = rnorm(n))
re <- rnorm(8); res <- 0.7 * re + rnorm(8, 0, 0.2); names(re) <- names(res) <- lev
dm$y <- 1 + (re[dm$g1] + re[dm$g2]) / 2 + rnorm(n, 0, exp(-0.2 + (res[dm$g1] + res[dm$g2]) / 2))
for (txt in c("bf(y ~ x + (1 | mm(g1, g2, id = 'p')), sigma ~ (1 | mm(g1, g2, id = 'p', cor = FALSE)))",
              "bf(y ~ x + (1 | p | mm(g1, g2)), sigma ~ (1 | p | mm(g1, g2, cor = FALSE)))",
              "bf(y ~ x + (1 | mm(g1, g2, id = 'p', cor = FALSE)), sigma ~ (1 | mm(g1, g2, id = 'p')))")) {
  a <- tr(frm(eval(parse(text = txt), envir = asNamespace("frmtmb")), data = dm))
  co <- tr(brms::make_stancode(eval(parse(text = txt), envir = asNamespace("brms")), data = dm))
  pr("%-80s frm %s, blocks %s, logLik %s; brms %s, has cor_1 %s", txt, msg(a, 90),
     if (is_err(a)) "-" else paste(vapply(a$frame$re_blocks, function(k) paste0(k$covstruct, k$dim), ""), collapse = ","),
     if (is_err(a)) "-" else format(as.numeric(logLik(a)), digits = 7), msg(co, 50), if (is_err(co)) "-" else grepl("cor_1", co))
}
a <- tr(frm(bf(y ~ x + (1 | p | mm(g1, g2)), sigma ~ (1 | p | mm(g1, g2))), data = dm))
pr("the correlated fit: %s, logLik %s", msg(a, 80), if (is_err(a)) "-" else format(as.numeric(logLik(a)), digits = 7))

# the order of the draw names
d <- data.frame(x = rnorm(60), y = rgamma(60, 2, 1)); d$y[c(5, 25, 45)] <- NA
s <- quiet(frm_sample(bf(y | mi() ~ x), family = stats::Gamma(link = "log"), data = d, chains = 1, iter = 60, seed = 1, refresh = 0))
pr("frm_sample order: %s", paste(names(posterior::as_draws_df(s)), collapse = " "))

# the bootstrap on dropped rows: same sd as on the complete-case fit?
set.seed(11100)
n <- 120; x <- rnorm(n)
d <- data.frame(x = x, y = rgamma(n, 2, 2 / exp(0.2 + 0.4 * x))); M <- seq(4, n, by = 12); d$y[M] <- NA
Gam <- stats::Gamma(link = "log")
a <- frm(bf(y | mi() ~ x), family = Gam, data = d)
b <- frm(bf(y ~ x), family = Gam, data = d[-M, ])
ba <- quiet(frm_bootstrap(a, nsim = 300, seed = 5)); bb <- quiet(frm_bootstrap(b, nsim = 300, seed = 5))
ta <- matrix(ba$t, nrow = 300); tb <- matrix(bb$t, nrow = 300)
pr("bootstrap sd: dropped rows %s; complete case %s", paste(round(apply(ta, 2, sd), 4), collapse = " "), paste(round(apply(tb, 2, sd), 4), collapse = " "))
# a refit on simulated data: are the unread rows out of it? the refit's nobs
sim <- simulate(a, nsim = 1, seed = 1)
pr("simulate fills the unread rows: %d", sum(!is.na(sim[[1]][M])))
pf <- tr(get("frm_bootstrap", asNamespace("frmtmb"))(a, FUN = function(f) attr(logLik(f), "nobs"), nsim = 3, seed = 2))
pr("bootstrap refit's nobs attr: %s %s", msg(pf, 120), if (is_err(pf)) "" else paste(unlist(pf$t), collapse = " "))

# when does a read row act unread? the share: 60 of 200 read by a reader whose
# slope is r, exgaussian; frm against the observed rows on y's own parameters,
# which the observed rows nearly recover when the reader carries no signal
set.seed(11501)
for (r in c(0, 0.02, 0.05, 0.1, 0.2, 0.4, 0.7)) for (mm in c(20, 60)) {
  shs <- numeric(8)
  for (k in 1:8) {
    n <- 200; x <- rnorm(n)
    y <- 0.3 + 0.5 * x + rnorm(n, 0, 0.5) + rexp(n, 1 / 2) - 2
    d <- data.frame(x = x, y = y); d$z <- 0.5 + r * y + rnorm(n, 0, 1)
    rows <- sample(n, mm); d$y[rows] <- NA
    j <- tr(frm(bf(y | mi() ~ x) + exgaussian() + bf(z ~ mi(y)) + gaussian(), data = d))
    cc <- tr(frm(bf(y ~ x), family = exgaussian(), data = d[-rows, ]))
    if (is_err(j) || is_err(cc)) { shs[k] <- NA; next }
    H <- stats::optimHess(cc$opt$par, cc$obj$fn, cc$obj$gr); se <- sqrt(diag(solve(H)))
    ej <- j$estimates$betad[grep("^y_", names(j$estimates$betad))]
    shs[k] <- max(abs((ej - cc$estimates$betad) / se[3:4]))
  }
  pr("exgaussian, reader slope %.2f, %2d read of 200: y's dpar shift from the observed rows, SE: median %.2f max %.2f (NA %d)", r, mm, median(shs, na.rm = TRUE), max(shs, na.rm = TRUE), sum(is.na(shs)))
}
