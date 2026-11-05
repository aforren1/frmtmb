# Which families brms 2.23.0 accepts for an mi() response and for an
# mi() predictor: every family of brms::brmsfamily's list, stancode plus
# a cached compile. Seed 169. Writes dev/refuse-brms-mi-out/mi-families-pred.tsv.
.libPaths(c("C:/Users/adf44/source/r/rellib-r16",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(FRMTMB_STAN_CACHE = normalizePath("dev/stan-cache", winslash = "/"),
           R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressPackageStartupMessages({library(frmtmb); library(testthat)})
E <- new.env(parent = asNamespace("frmtmb"))
for (f in c("dev/equiv-lib.R", "dev/equiv-specs.R", "dev/equiv-flow.R")) {
  sys.source(f, envir = E)
}
dir.create("dev/refuse-brms-mi-out", showWarnings = FALSE)
fams <- ls(asNamespace("brms"), pattern = "^[.]family_", all.names = TRUE)
fams <- sub("^[.]family_", "", fams)
fams <- fams[fams != "custom"]
set.seed(169)
n <- 120
base <- data.frame(x = rnorm(n), z = rnorm(n))
pos <- exp(rnorm(n))
draw <- function(f) {
  switch(f,
    beta = , xbeta = , zero_inflated_beta = , zero_one_inflated_beta = ,
    logistic_normal = plogis(rnorm(n)),
    von_mises = runif(n, -3, 3),
    bernoulli = , zero_inflated_binomial = , binomial = ,
    beta_binomial = rbinom(n, 1, .5),
    poisson = , negbinomial = , negbinomial2 = , geometric = , com_poisson = ,
    discrete_weibull = , zero_inflated_poisson = ,
    zero_inflated_negbinomial = , hurdle_poisson = , hurdle_negbinomial = ,
    zero_inflated_beta_binomial = rpois(n, 3),
    categorical = , multinomial = , cumulative = , sratio = , cratio = ,
    acat = , hurdle_cumulative = sample(1:3, n, TRUE),
    gaussian = , student = , skew_normal = , asym_laplace = ,
    gen_extreme_value = rnorm(n),
    pos)
}
out <- list()
one <- function(role, f) {
  d <- base
  v <- draw(f)
  d$v <- v
  d$vm <- v
  d$vm[sample(n, 10)] <- NA
  d$w <- rnorm(n)
  fam <- tryCatch(brms::brmsfamily(f), error = function(e) NULL)
  if (is.null(fam)) return(data.frame(role, family = f, built = NA,
                                      compiled = NA,
                                      note = "no brmsfamily"))
  lhs <- if (f %in% c("binomial", "beta_binomial", "zero_inflated_binomial", "zero_inflated_beta_binomial")) " | trials(10)" else ""
  bform <- tryCatch(switch(role,
    resp = brms::bf(vm | mi() ~ x, family = fam),
    pred = brms::bf(stats::as.formula(paste0("v", lhs, " ~ mi(w2) + x")), family = fam) +
      brms::bf(w2 | mi() ~ z, family = stats::gaussian()) +
      brms::set_rescor(FALSE),
    both = brms::bf(y ~ mi(vm) + x) + brms::bf(vm | mi() ~ z, family = fam) +
      brms::set_rescor(FALSE)), error = function(e) e)
  if (inherits(bform, "error")) return(data.frame(role, family = f,
    built = FALSE, compiled = NA, note = conditionMessage(bform)))
  d$w2 <- d$w
  d$w2[sample(n, 10)] <- NA
  d$y <- rnorm(n)
  code <- tryCatch(suppressWarnings(brms::stancode(bform, data = d)),
                   error = function(e) e)
  if (inherits(code, "error")) return(data.frame(role, family = f,
    built = FALSE, compiled = NA,
    note = substr(gsub("[\r\n\t]+", " ", conditionMessage(code)), 1, 200)))
  cc <- tryCatch({suppressWarnings(E$brms_stan_model(code)); TRUE},
                 error = function(e) e)
  if (inherits(cc, "error")) return(data.frame(role, family = f, built = TRUE,
    compiled = FALSE,
    note = substr(gsub("[\r\n\t]+", " ", conditionMessage(cc)), 1, 200)))
  data.frame(role, family = f, built = TRUE, compiled = TRUE, note = "")
}
for (role in c("pred")) {
  for (f in fams) {
    r <- one(role, f)
    out[[length(out) + 1L]] <- r
    cat(role, f, r$built, r$compiled, r$note, "\n")
    utils::write.table(do.call(rbind, out),
                       "dev/refuse-brms-mi-out/mi-families-pred.tsv",
                       sep = "\t", row.names = FALSE)
  }
}
