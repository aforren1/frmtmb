# Lane wt-mvprior: frm_sample()'s default rows on the several-location
# and multivariate cases, per arm, beside brms's in brms-probe.txt.
source("C:/Users/adf44/source/r/frmtmb-wt-mvprior/dev/mvprior-prelude.R")
suppressMessages({library(frmtmb); library(frmtmb.sample)})
source(file.path(MVPRIOR_ROOT, "dev/mvprior-cases.R"))
cat("arm", mvprior_arm, "\n")
for (case in c("catre", "mixre", "mvcat")) {
  m <- mvprior_model(case, "frmtmb")
  tab <- tryCatch(as.data.frame(default_prior(m[[1]], data = m[[2]],
                                              family = m[[3]],
                                              route = "sample")),
                  error = function(e) conditionMessage(e))
  cat("==", case, "\n")
  if (is.data.frame(tab)) {
    tab <- tab[tab$prior != "(flat)" & tab$class != "theta", ]
    print(tab[, c("prior", "class", "group", "resp", "dpar")],
          row.names = FALSE)
  } else cat(tab, "\n")
  fit <- suppressMessages(frm(m[[1]], data = m[[2]], family = m[[3]],
                              dry_run = "objective"))
  cat(tryCatch({
    frmtmb.sample:::announce_default_priors(
      frmtmb.sample:::default_priors_for(fit), character(0))
    "announced"
  }, message = function(msg) conditionMessage(msg),
  error = function(e) paste("ERROR", conditionMessage(e))), "\n")
}
