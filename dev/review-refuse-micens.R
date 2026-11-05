# Item 7: brms's mi() x cens() program, closely. Seed 10004.
source("dev/review-refuse-lib.R")
set.seed(10004); n <- 180; x <- rnorm(n)
base <- data.frame(x = x, y = rnorm(n, 0.5 * x), ce = 0L)
base$y[seq(5, n, 9)] <- NA
rows <- which(is.na(base$y))
form <- brms::bf(y | mi() + cens(ce) ~ x)
pri <- E$brms_flat_prior(form, data = base, family = stats::gaussian())
code <- brms::make_stancode(form, data = base, prior = pri)
writeLines(strsplit(code, "\n")[[1]][38:50])
mod <- E$brms_stan_model(code)
pts <- list(b = array(0.5, 1L), Intercept = 0, sigma = 1)
for (cc in c(0L, 1L, -1L)) {
  dd <- base; dd$ce[rows] <- cc
  sdat <- E$brms_standata(form, data = dd, prior = pri)
  pr("coded %+d: rows in Jevent %d, Jrcens %d, Jlcens %d of %d; Nmi %d; Y there %s; Nevent %d Nrcens %d Nlcens %d", cc,
     sum(rows %in% sdat$Jevent), sum(rows %in% sdat$Jrcens), sum(rows %in% sdat$Jlcens), length(rows), sdat$Nmi, unique(sdat$Y[rows]), sdat$Nevent, sdat$Nrcens, sdat$Nlcens)
  sf <- suppressMessages(rstan::sampling(mod, data = sdat, chains = 0))
  p1 <- c(pts, list(Ymi = array(rep(0.2, sdat$Nmi))))
  p2 <- c(pts, list(Ymi = array(rep(-1, sdat$Nmi))))
  lp <- function(p) rstan::log_prob(sf, rstan::unconstrain_pars(sf, p), adjust_transform = FALSE)
  pr("   log_prob at two values of Ymi: %s, %s", format(lp(p1)), format(lp(p2)))
  s <- tr(suppressWarnings(rstan::sampling(mod, data = sdat, chains = 1, iter = 200, refresh = 0, seed = 1,
     init = function() list(b = array(0.5, 1L), Intercept = 0, sigma = 1, Ymi = array(rep(0.1, sdat$Nmi))))))
  pr("   sampling, explicit inits: %s; draws kept %s", msg(s, 150), if (is_err(s)) "-" else sum(s@sim$n_save))
  if (!is_err(s) && length(s@sim$samples)) {
    e <- rstan::extract(s)
    if (length(e$Ymi)) pr("   Ymi[1] draw sd %.2f", sd(e$Ymi[, 1]))
  }
}
# the no-NA control: brms's cens() is read; a response with no missing value
d0 <- base; d0$y[rows] <- 0.1
sd0 <- brms::make_standata(form, data = d0)
pr("no NA rows: Nmi %s, Nevent %s", sd0$Nmi, sd0$Nevent)
# what frm() says in the same two cases
a <- tr(frm(bf(y | mi() + cens(ce) ~ x), data = base)); pr("frm with NA rows: %s", msg(a, 140))
a <- tr(frm(bf(y | mi() + cens(ce) ~ x), data = d0)); pr("frm with no NA rows: %s", msg(a, 140))
