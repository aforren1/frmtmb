# Reviewer of lane sampleprior, fix round 1: a scalar-xi GEV from a
# formula, and from a fit, on this arm (REV_ARM=r9 for the base).
source("C:/Users/adf44/source/r/frmtmb-wt-sampleprior/dev/sampleprior-rev-common.R")
source("C:/Users/adf44/source/r/frmtmb-wt-sampleprior/dev/sampleprior-rev-cases.R")
fam <- frmtmb::brmsfamily("gen_extreme_value")
for (route in c("formula", "fit")) {
  r <- tryCatch({
    m <- if (route == "formula") y ~ x else suppressWarnings(frm(y ~ x, data = D$gev, family = fam))
    ds <- suppressMessages(suppressWarnings(frm_sample(m, data = if (route == "formula") D$gev,
      family = if (route == "formula") fam, chains = 1, iter = 200, refresh = 0, seed = 1)))
    sprintf("OK, xi mean %.3f", mean(ds$draws[, "xi"]))
  }, error = function(e) paste("ERROR", substr(conditionMessage(e), 1, 90)))
  cat(route, ":", r, "\n")
}
