# Fix round 2, (1): what brms 2.23.0 does with newdata when the fit
# dropped rows for missing values. 36 rows, x missing in rows 5 and 20,
# fcor(V) and sar(Ws) (lag), brms's draws held at frmtmb's estimates
# (Fixed_param). Rscript dev/sarfcor-fr2-rows.R [frmtmb]
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win",
           FRMTMB_STAN_CACHE =
             "C:/Users/adf44/source/r/frmtmb-wt-sarfcor/dev/stan-cache",
           FRMTMB_BRMS_FIT_TESTS = "true", NOT_CRAN = "true")
.libPaths(c("C:/Users/adf44/source/r/wt-sarfcor-lib",
            "C:/Users/adf44/source/r/rellib-r12",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb.spatial))
wt <- "C:/Users/adf44/source/r/frmtmb-wt-sarfcor"
env <- new.env(parent = asNamespace("frmtmb"))
sys.source(file.path(wt, "extensions/frmtmb.spatial/tests/testthat",
                     "helper-brms.R"), envir = env)
fixed_fit <- function(bform, data, data2, pars) {
  fam <- gaussian()
  prior <- env$brms_flat_prior(bform, data = data, family = fam,
                               data2 = data2)
  code <- brms::make_stancode(bform, data = data, family = fam,
                              prior = prior, data2 = data2)
  sdat <- env$brms_standata(bform, data = data, family = fam,
                            prior = prior, data2 = data2)
  sc <- suppressMessages(brms::brm(bform, data = data, family = fam,
                                   prior = prior, data2 = data2,
                                   empty = TRUE))
  sc$fit <- suppressMessages(rstan::sampling(
    env$brms_stan_model(code), data = sdat, chains = 1, iter = 2,
    warmup = 0, algorithm = "Fixed_param", init = list(pars),
    refresh = 0, seed = 1))
  brms::rename_pars(sc)
}
set.seed(2)
k <- 6
g <- expand.grid(r = seq_len(k), c = seq_len(k))
n <- nrow(g)
A <- outer(seq_len(n), seq_len(n), function(i, j) {
  as.numeric(abs(g$r[i] - g$r[j]) + abs(g$c[i] - g$c[j]) == 1)
})
dg <- 1 / sqrt(rowSums(A))
Ws <- A * outer(dg, dg)
V <- 0.6^abs(outer(1:n, 1:n, "-"))
d <- data.frame(x = rnorm(n))
d$y <- as.vector(solve(diag(n) - 0.4 * Ws, 1 + 0.5 * d$x + rnorm(n, 0, 0.7)))
d$x[c(5, 20)] <- NA
keep <- setdiff(seq_len(n), c(5, 20))
cases <- list(
  full36 = transform(d, x = x + 1),
  kept34 = transform(d[keep, ], x = x + 1),
  filled36 = transform(d, x = ifelse(is.na(x), 0, x) + 1),
  othermissing36 = local({
    z <- transform(d, x = ifelse(is.na(x), 0, x) + 1)
    z$x[c(7, 30)] <- NA
    z
  }),
  ynull36 = transform(d, x = x + 1, y = NA))
for (term in c("fcor", "sar")) {
  form <- if (term == "fcor") y ~ x + fcor(V) else y ~ x + sar(Ws)
  d2 <- if (term == "fcor") list(V = V) else list(Ws = Ws)
  f <- suppressMessages(frm(form, data = d, data2 = d2))
  b <- unname(coef(f))
  pars <- list(b = array(b[2L], 1),
               Intercept = b[1L] + b[2L] * mean(d$x[keep]),
               sigma = sigma(f))
  if (term == "sar") {
    ac <- f$frame$autocor[[1L]]
    pars$lagsar <- unname(frmtmb:::autocor_natural(f$estimates$thetaac, ac))
  }
  bfit <- suppressWarnings(fixed_fit(brms::bf(form), d, d2, pars))
  for (nm in names(cases)) {
    nd <- cases[[nm]]
    rb <- tryCatch({
      e <- brms::posterior_epred(bfit, newdata = nd)[1L, ]
      sprintf("%d values, %d NA", length(e), sum(is.na(e)))
    }, error = function(e) paste("ERROR:", conditionMessage(e)))
    rf <- tryCatch({
      e <- fitted(f, newdata = nd)[, "Estimate"]
      out <- sprintf("%d values, %d NA", length(e), sum(is.na(e)))
      if (!startsWith(rb, "ERROR")) {
        be <- brms::posterior_epred(bfit, newdata = nd)[1L, ]
        if (length(be) == length(e)) {
          ok <- !is.na(e)
          out <- sprintf("%s; NA rows as brms: %s; max rel diff %s", out,
                         identical(is.na(be), is.na(unname(e))),
                         if (any(ok)) format(max(abs(be[ok] - e[ok])) /
                                               max(abs(e[ok])), digits = 3)
                         else "-")
        }
      }
      out
    }, error = function(e) paste("ERROR:", substr(conditionMessage(e), 1, 90)))
    cat(sprintf("%-5s %-15s brms: %-60s | frmtmb: %s\n", term, nm,
                substr(rb, 1, 60), rf))
  }
}
