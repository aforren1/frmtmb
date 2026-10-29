## Review rescorpp: the lane's frmtmb.sample on rellib-r15's core (0.77.0).
##   Rscript dev/review-rescorpp-floor.R rv
source("C:/Users/adf44/source/r/frmtmb-wt-rescorpp/dev/review-rescorpp-common.R")
cat("DESCRIPTION floor:", grep("frmtmb \\(", readLines(file.path(
  find.package("frmtmb.sample"), "DESCRIPTION")), value = TRUE), "\n")
source(file.path(wt, "dev/rescorpp-data.R"))
d <- rescorpp_data("g3")
fit <- quiet(frm(eval(rescorpp_frm_forms[["g3"]]), family = "gaussian", data = d))
ds <- quiet(frm_sample(fit, chains = 1, iter = 300, refresh = 0, seed = 2))
show("posterior_predict (joint)", dim(posterior_predict(ds)))
show("log_lik resp c(y1, y3)", dim(log_lik(ds, resp = c("y1", "y3"))))
show("log_lik resp y1", dim(log_lik(ds, resp = "y1")))
show("posterior_predict resp y1", dim(posterior_predict(ds, resp = "y1")))
