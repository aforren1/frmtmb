# Item 8: the changed compat rows against what both packages do, on designs of
# mine; mi() and se() over every family of frm_compat(), so a row that says
# works and does not is seen; a fresh, uncached compile for a sample of the
# sweep's "builds and compiles". Seed 10100.
source("dev/review-refuse-lib.R")
library(frmtmb.sample)
set.seed(10100)
n <- 100
d <- data.frame(x = rnorm(n), z = rnorm(n), ce = 0L, s = 0.3, tr = 12L, w = 1)
d$y <- rnorm(n); d$ypos <- rgamma(n, 2, 1); d$yu <- runif(n, 0.05, 0.95)
d$ycnt <- rpois(n, 3); d$yo <- findInterval(d$x + rnorm(n), c(-1, 0, 1)) + 1L
d$ybin <- rbinom(n, 1, 0.5)
d$yc <- runif(n, -3, 3)
d$g <- factor(rep(1:10, 10))
cp <- frm_compat()
print(names(cp))
fams <- unique(cp$feature_a[cp$kind_a == "family"])
pr("family count in the table: %d", length(fams))
rowof <- function(a, b) {
  r <- cp[(cp$feature_a == a & cp$feature_b == b) | (cp$feature_a == b & cp$feature_b == a), ]
  if (nrow(r)) r$status[1] else "none"
}
for (pair in list(c("cens()", "disc"), c("cens()", "equidistant"), c("vint()", "sum_to_zero"), c("trials()", "cs_pred()"),
                  c("s()", "prior"), c("t2()", "prior"), c("mm()", "|ID|"), c("mm()", "gr_cov"), c("|ID|", "mmc()"), c("gr_cov", "mmc()"),
                  c("mi()", "se()"), c("se()", "skew_normal"), c("cens()", "mixture"), c("trunc()", "mixture"), c("thres()", "mvbf"),
                  c("mi()", "asym_laplace"), c("mi()", "cens()"), c("mi()", "trunc()"), c("mixture", "mi()"), c("mi()", "von_mises"), c("mi()", "gaussian"), c("mi()", "poisson"))) {
  pr("row %-12s x %-14s %s", pair[1], pair[2], rowof(pair[1], pair[2]))
}

# frmtmb family object and response for each table name; brms's by name
pick_y <- function(f) {
  if (grepl("binomial|bernoulli", f)) return(if (grepl("bernoulli", f)) "ybin" else "ycnt")
  if (grepl("poisson|negbinomial|nbinom|geometric|discrete_weibull|com_poisson", f)) return("ycnt")
  if (grepl("beta", f)) return("yu")
  if (grepl("von_mises", f)) return("yc")
  if (grepl("cumulative|sratio|cratio|acat|categorical|multinomial", f)) return("yo")
  if (grepl("Gamma|gamma|weibull|exponential|frechet|inverse|lognormal|gaussian$", f) && !identical(f, "gaussian")) return("ypos")
  "y"
}
fam_frm <- function(f) tryCatch({
  if (f == "Gamma") return(stats::Gamma(link = "log"))
  if (f == "inverse.gaussian") return(stats::inverse.gaussian(link = "log"))
  if (f == "beta") return(frmtmb::Beta())
  get(f, asNamespace("frmtmb"))()
}, error = function(e) NULL)
fam_brm <- function(f) tryCatch({
  if (f == "Gamma") return(brms::brmsfamily("gamma", "log"))
  if (f == "beta") return(brms::Beta())
  brms::brmsfamily(gsub("[.]", "_", f))
}, error = function(e) tryCatch(get(f, asNamespace("brms"))(), error = function(e2) NULL))
bad <- 0L
for (f in fams) {
  fo <- fam_frm(f); fb <- fam_brm(f)
  yv <- pick_y(f)
  dd <- d
  dd$yy <- dd[[yv]]; dd$yy[3] <- NA
  st <- rowof("mi()", f)
  extra <- if (grepl("binomial", f)) " + trials(tr)" else ""
  ff <- as.formula(sprintf("yy | mi()%s ~ x", extra))
  a <- if (is.null(fo)) structure(list(msg = "no family object"), class = "tr_err") else tr(frm(bf(ff), family = fo, data = dd, dry_run = "frame"))
  b <- if (is.null(fb)) structure(list(msg = "no brms family"), class = "tr_err") else tr(brms::make_stancode(brms::bf(ff), family = fb, data = dd))
  a_ok <- !is_err(a); b_ok <- !is_err(b)
  flag <- ""
  if ((identical(st, "works") && !a_ok) || (identical(st, "refused") && a_ok)) flag <- "ROW!=FRM"
  if (b_ok != a_ok && !f %in% c("asym_laplace")) flag <- paste(flag, "FRM!=BRMS")
  if (nzchar(flag)) {
    bad <- bad + 1L
    pr("mi  %-24s row=%-10s frm=%-5s brms=%-5s %s | %s | %s", f, st, a_ok, b_ok, flag, msg(a, 70), msg(b, 70))
  }
  # se()
  st <- rowof("se()", f)
  a <- if (is.null(fo)) structure(list(msg = "no family object"), class = "tr_err") else tr(frm(bf(as.formula(paste0(yv, " | se(s) ~ x"))), family = fo, data = d, dry_run = "frame"))
  b <- if (is.null(fb)) structure(list(msg = "no brms family"), class = "tr_err") else tr(brms::make_stancode(brms::bf(as.formula(paste0(yv, " | se(s) ~ x"))), family = fb, data = d))
  a_ok <- !is_err(a); b_ok <- !is_err(b)
  flag <- ""
  if ((identical(st, "works") && !a_ok) || (identical(st, "refused") && a_ok)) flag <- "ROW!=FRM"
  if (b_ok != a_ok) flag <- paste(flag, "FRM!=BRMS")
  if (nzchar(flag)) pr("se  %-24s row=%-10s frm=%-5s brms=%-5s %s | %s | %s", f, st, a_ok, b_ok, flag, msg(a, 70), msg(b, 70))
  # cens() x this family inside a mixture
}
pr("mi() flags: %d", bad)

# cens x mixture: the components that have no CDF; the row says "when every
# component has a CDF"
for (fx in c("mixture(gaussian(), gaussian(), order = 'none')", "mixture(lognormal(), Gamma(link = 'log'), order = 'none')",
             "mixture(poisson(), poisson(), order = 'none')", "mixture(weibull(), exponential(), order = 'none')",
             "mixture(student(), skew_normal(), order = 'none')", "mixture(hurdle_gamma(), gamma(), order = 'none')",
             "mixture(von_mises(), von_mises(), order = 'none')", "mixture(Beta(), Beta(), order = 'none')",
             "mixture(gaussian(), wiener(), order = 'none')")) {
  yy <- if (grepl("Beta", fx)) "yu" else if (grepl("poisson", fx)) "ycnt" else if (grepl("von_mises", fx)) "yc" else if (grepl("lognormal|weibull|hurdle", fx)) "ypos" else "y"
  ftxt <- as.formula(sprintf("%s | cens(ce) ~ 1", yy))
  a <- tr(frm(bf(ftxt), family = eval(parse(text = fx), envir = asNamespace("frmtmb")), data = d))
  fxb <- gsub("gamma\\(\\)", "brms::hurdle_gamma()", fx)
  b <- tr(brms::make_stancode(brms::bf(ftxt), family = eval(parse(text = gsub("Gamma\\(link = 'log'\\)", "brmsfamily('gamma', 'log')", fx)), envir = asNamespace("brms")), data = d))
  pr("cens x %-55s frm %s | brms %s", fx, msg(a, 90), msg(b, 90))
}

# a fresh, uncached compile of a sample of the sweep's "builds and compiles"
sw <- list(
  list("cens_mixture", brms::bf(y | cens(ce) ~ 1, mu1 ~ x, mu2 ~ x), brms::mixture(gaussian(), gaussian(), order = "none")),
  list("trunc_mixture", brms::bf(ypos | trunc(lb = 0) ~ 1, mu1 ~ x, mu2 ~ x), brms::mixture(gaussian(), gaussian(), order = "none")),
  list("mi_trunc", brms::bf(y | mi() + trunc(lb = -50) ~ x), gaussian()),
  list("mi_se", brms::bf(y | mi() + se(s) ~ x), gaussian()),
  list("mixture_se", brms::bf(y | se(s) ~ 1, mu1 ~ x, mu2 ~ x), brms::mixture(gaussian(), gaussian(), order = "none")),
  list("rescor_weights", brms::bf(brms::mvbind(y, z) | weights(w) ~ x) + brms::set_rescor(TRUE), gaussian()))
Sys.unsetenv("FRMTMB_STAN_CACHE")
for (s in sw) {
  dd <- d; dd$y[5] <- if (grepl("^mi", s[[1]])) NA else dd$y[5]
  co <- tr(brms::make_stancode(s[[2]], family = s[[3]], data = dd))
  if (is_err(co)) { pr("sweep %-16s stancode: %s", s[[1]], msg(co, 150)); next }
  t0 <- Sys.time()
  m <- tr(rstan::stan_model(model_code = co, auto_write = FALSE))
  sd0 <- tr(brms::make_standata(s[[2]], family = s[[3]], data = dd))
  ss <- if (is_err(m)) m else tr(rstan::sampling(m, data = sd0, chains = 1, iter = 60, refresh = 0, seed = 1))
  pr("sweep %-16s fresh compile %s (%.0f s); a short sampling run %s", s[[1]], msg(m, 120), as.numeric(difftime(Sys.time(), t0, units = "secs")), msg(ss, 120))
}
