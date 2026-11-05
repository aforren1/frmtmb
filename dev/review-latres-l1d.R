# Short check 2, item 2: is the "no Cholesky factor" of the L1 cell in the fit or
# in the translator? Read Lcortime and the map through us_chol_L() directly, not
# through chol() of the rebuilt correlation, on the test's data and on the sweep's
# Gamma() cell.
source("C:/Users/adf44/source/r/frmtmb-wt-latres/dev/review-latres-common.R")
Sys.setenv(FRMTMB_STAN_CACHE = "C:/Users/adf44/source/r/frmtmb-wt-latres/dev/stan-cache", FRMTMB_BRMS_FIT_TESTS = "true")
suppressPackageStartupMessages(library(rstan))
hel <- new.env(parent = asNamespace("frmtmb"))
for (f in c("helper-brms.R", "helper-brms-methods.R", "helper-warnings.R")) sys.source(file.path("tests/testthat", f), envir = hel)
assignInNamespace("latres_refuse_margin", function(fit) invisible(NULL), ns = "frmtmb")
t0 <- parse(text = readLines("tests/testthat/test-latres-brms.R"))
for (e in t0) if (is.call(e) && identical(e[[1]], as.name("<-")) && as.character(e[[2]]) %in% c("lrb_data")) eval(e)
d <- lrb_data()
orig <- hel$brms_latres_chol
hel$brms_latres_chol <- function(a, th, k, J = NULL) {
  if (identical(a$kind, "unstr") && identical(as.integer(J), seq_len(a$d))) {
    return(as.matrix(fi("us_chol_L")(th, a$d)))
  }
  orig(a, th, k, J)
}
for (cc in list(list("Beta", brms::brmsfamily("beta"), Beta(), "unit"), list("Gamma()", brms::brmsfamily("Gamma"), Gamma(), "pos"))) {
  f <- stats::as.formula(paste(cc[[4]], "~ x + unstr(t, g)"))
  fit <- suppressWarnings(frm(frmtmb::bf(f), family = cc[[3]], data = d))
  r <- tryCatch(hel$brms_lp_check(brms::bf(f), cc[[2]], d, fit, joint = TRUE), error = function(e) conditionMessage(e))
  cat(cc[[1]], ": ", if (is.character(r)) substr(r, 1, 200) else paste("const", r$measured_const, " ours", r$ours), "\n")
}
cat("\n== frm_lp_basis at 8000 rows, memory of the matrix\n")
set.seed(1)
N <- 8000
dx <- data.frame(t = seq_len(N), g = factor(1), x = rnorm(N))
dx$y <- rpois(N, exp(.5 + .3 * dx$x + as.numeric(arima.sim(list(ar = .5), N, sd = .5))))
fb <- suppressWarnings(frm(y ~ x + ar(t, g, cov = TRUE), family = poisson(), data = dx))
t2 <- system.time(z <- frm_lp_basis(fb))[[3]]
cat("seconds", t2, " class", class(z)[1], "\n")
str(z, max.level = 1, give.attr = FALSE, list.len = 8)
cat("MB", as.numeric(object.size(z)) / 1e6, "\n")
