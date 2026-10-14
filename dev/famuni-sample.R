# frm_sample() (frmtmb.sample) and as_tmbstan() on the three families:
# that they sample, what default priors they announce against brms's
# default_prior(), and the draws methods (log_lik, posterior_predict,
# posterior_epred, pp_check). Seed 20261017, n = 300, 2 chains x 1000.
# Output: dev/famuni-out/sample.txt
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
.libPaths(c("C:/Users/adf44/source/r/wt-famuni-lib",
            "C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
stopifnot(pkgbuild::has_build_tools())
suppressMessages({
  library(frmtmb)
  library(frmtmb.sample)
})
cat("frmtmb", format(packageVersion("frmtmb")), "from", find.package("frmtmb"),
    "; frmtmb.sample from", find.package("frmtmb.sample"), "\n")
set.seed(20261017)
n <- 300
d <- data.frame(x = rnorm(n), z = rnorm(n))
d$yf <- exp(1 + 0.3 * d$x) / gamma(1 - 1 / 4) * rexp(n)^(-1 / 4)
d$yg <- 1 + 0.5 * d$x + exp(-0.3 + 0.2 * d$z) * (rexp(n)^(-0.15) - 1) / 0.15
q <- plogis(1 + 0.4 * d$x)
d$yd <- pmax(ceiling((log1p(-runif(n)) / log(q))^(1 / 1.2) - 1), 0)
specs <- list(
  frechet = list(f = bf(yf ~ x), fam = frechet(), bfam = brms::frechet()),
  gev = list(f = bf(yg ~ x, sigma ~ z, xi ~ z), fam = gen_extreme_value(),
             bfam = brms::gen_extreme_value()),
  gev_scalar = list(f = bf(yg ~ x, sigma ~ z), fam = gen_extreme_value(),
                    bfam = brms::gen_extreme_value()),
  dw = list(f = bf(yd ~ x), fam = discrete_weibull(),
            bfam = brms::brmsfamily("discrete_weibull")))
for (nm in names(specs)) {
  s <- specs[[nm]]
  cat("\n=====", nm, "\n")
  # brms reads its own bf(); handed frmtmb's, it dropped the dpar formulas
  bfo <- do.call(brms::bf, unname(lapply(c(list(s$f$formula), s$f$pforms),
                                         stats::as.formula)))
  bp <- suppressWarnings(brms::default_prior(bfo, data = d, family = s$bfam))
  cat("brms default_prior (rows with a prior):\n")
  print(as.data.frame(bp)[nzchar(bp$prior), c("prior", "class", "coef",
                                               "dpar", "lb")],
        row.names = FALSE)
  cat("frmtmb get_prior(route = \"sample\"):\n")
  gp <- get_prior(s$f, data = d, family = s$fam, route = "sample")
  print(as.data.frame(gp)[gp$prior != "(flat)", c("prior", "class", "dpar")],
        row.names = FALSE)
  fit <- frm(s$f, family = s$fam, data = d)
  msgs <- character(0)
  dr <- withCallingHandlers(
    tryCatch(frm_sample(fit, chains = 2, iter = 1000, seed = 1),
             error = function(e) e),
    message = function(m) {
      msgs <<- c(msgs, conditionMessage(m))
      invokeRestart("muffleMessage")
    },
    warning = function(w) {
      cat("  WARN:", conditionMessage(w), "\n")
      invokeRestart("muffleWarning")
    })
  cat("frm_sample() messages:\n", paste(" ", msgs, collapse = "\n"), "\n")
  if (inherits(dr, "error")) {
    cat("  ERROR:", conditionMessage(dr), "\n")
    next
  }
  print(summary(dr))
  ll <- tryCatch(log_lik(dr), error = function(e) conditionMessage(e))
  cat("log_lik dim:", if (is.character(ll)) ll else dim(ll), "\n")
  pp <- tryCatch(posterior_predict(dr, ndraws = 20),
                 error = function(e) conditionMessage(e))
  cat("posterior_predict dim:", if (is.character(pp)) pp else dim(pp), "\n")
  pe <- tryCatch(posterior_epred(dr, ndraws = 20),
                 error = function(e) conditionMessage(e))
  cat("posterior_epred dim:", if (is.character(pe)) pe else dim(pe), "\n")
  pc <- tryCatch(class(pp_check(dr, ndraws = 10))[1],
                 error = function(e) conditionMessage(e))
  cat("pp_check:", pc, "\n")
  ts <- tryCatch({
    st <- as_tmbstan(fit, chains = 1, iter = 400, seed = 1, refresh = 0)
    class(st)[1]
  }, error = function(e) conditionMessage(e))
  cat("as_tmbstan:", ts, "\n")
}
