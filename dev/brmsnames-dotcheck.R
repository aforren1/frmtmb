## For the frmtmb_draws methods whose own arguments run out before
## brms's (dev/brmsnames-beyond.R), whether an extra POSITIONAL argument
## is refused or answered. One call per method, the extra argument a
## value brms would read in its next slot.
##   Rscript dev/brmsnames-dotcheck.R lane
arm <- commandArgs(trailingOnly = TRUE)[1L]
source("dev/brmsnames-libs.R")
brmsnames_libs(arm)
q <- function(e) suppressWarnings(suppressMessages(e))
q(library(frmtmb)); q(library(frmtmb.sample))
ds <- readRDS(sprintf("dev/stan-cache/brmsnames-draws-%s.rds", arm))$ds
calls <- list(
  "bridge_sampler(ds, TRUE)" = quote(bridge_sampler(ds, TRUE)),
  "expose_functions(ds, TRUE)" = quote(expose_functions(ds, TRUE)),
  "family(ds, 'y')" = quote(family(ds, "y")),
  "log_lik(ds, NULL, NULL, NULL, NULL, NULL, TRUE)" =
    quote(log_lik(ds, NULL, NULL, NULL, NULL, NULL, TRUE)),
  "loo_moment_match(ds, NULL)" = quote(loo_moment_match(ds, NULL)),
  "nobs(ds, 'y')" = quote(nobs(ds, "y")),
  "nsamples(ds, 1:5)" = quote(nsamples(ds, 1:5)),
  "nuts_params(ds, 'stepsize__')" = quote(nuts_params(ds, "stepsize__")),
  "plot(ds, '^b_')" = quote(plot(ds, "^b_")),
  "posterior_epred(ds, ..., sort = TRUE) positional" =
    quote(posterior_epred(ds, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
                          NULL, TRUE)),
  "posterior_samples(ds, NA, FALSE)" = quote(posterior_samples(ds, NA, FALSE)),
  "print(ds, 3)" = quote(print(ds, 3)),
  "prior_summary(ds, FALSE)" = quote(prior_summary(ds, FALSE)),
  "reloo(ds, NULL)" = quote(reloo(ds, NULL)),
  "stancode(ds, TRUE)" = quote(stancode(ds, TRUE)),
  "standata(ds, NULL)" = quote(standata(ds, NULL)),
  "update(ds, y ~ 1)" = quote(update(ds, y ~ 1))
)
for (nm in names(calls)) {
  r <- tryCatch({ q(utils::capture.output(v <- eval(calls[[nm]]))); "ANSWERED" },
                error = function(e) paste("refused:",
                                          substr(conditionMessage(e), 1, 60)))
  cat(sprintf("%-50s %s\n", nm, r))
}
cat("DONE\n")
