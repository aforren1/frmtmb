## Recheck, priority 5: sample port-ledger rows independently of
## dev/famlink-port-ledger.R. Each expression is evaluated directly and
## its value or error message printed, so whether a pass is genuine is
## read from the message, not from a classifier.
## Usage: Rscript dev/famlink-rev2-ledger-sample.R <lane|base>
ARM <- commandArgs(trailingOnly = TRUE)[1]
source("dev/famlink-rev-common.R")
set.seed(1)
dat <- data.frame(y = rnorm(10), x = rnorm(10), g = rep(1:5, 2))
rows <- list(
  `6 negbinomial(inverse) ~ "inverse"` = quote(negbinomial(inverse)),
  `11 weibull()$family == weibull` = quote(weibull()$family),
  `17 hurdle_negbinomial("inverse")` = quote(hurdle_negbinomial("inverse")),
  `21 zero_inflated_poisson(list(1)) ~ zip` = quote(zero_inflated_poisson(list(1))),
  `26 zero_inflated_binomial(y ~ x) ~ zib` = quote(zero_inflated_binomial(y ~ x)),
  `28 categorical(probit) ~ "probit"` = quote(categorical(probit)),
  `42 beta_binomial(link_phi = "logit")` = quote(beta_binomial(link_phi = "logit")),
  `70 family_names(mix)` = quote(brms:::family_names(mix)),
  `77 mixture(poisson, binomial, order = "x")` = quote(mixture(poisson, binomial, order = "x")),
  `88 frm(family = "ordinal")` = quote(frm(y ~ x, dat, family = "ordinal")),
  `89 standata cratio on factor(-1:1)` = quote(frm(y ~ 1, data = data.frame(y = factor(-1:1)), family = "cratio", dry_run = "frame"))
)
for (nm in names(rows)) {
  v <- tryCatch({ x <- eval(rows[[nm]]); paste("VALUE:", paste(format(if (inherits(x, "frmtmb_family")) "<family object>" else x), collapse = " ")) },
                error = function(e) paste("ERROR:", conditionMessage(e)))
  cat(sprintf("%-44s %s\n", nm, substr(gsub("\n", " ", v), 1, 120)))
}
