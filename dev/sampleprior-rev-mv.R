# Reviewer of lane sampleprior: frm_sample() on a Student-t mvbf with
# rescor, both spellings, really sampling (no capture).
source("C:/Users/adf44/source/r/frmtmb-wt-sampleprior/dev/sampleprior-rev-common.R")
source("C:/Users/adf44/source/r/frmtmb-wt-sampleprior/dev/sampleprior-rev-cases.R")
run <- function(lab, expr) {
  r <- tryCatch({
    ds <- suppressMessages(suppressWarnings(eval(expr)))
    paste("OK, draws", nrow(ds$draws), "x", ncol(ds$draws), ":",
          paste(grep("nu|rescor", colnames(ds$draws), value = TRUE),
                collapse = " "))
  }, error = function(e) paste("ERROR:", conditionMessage(e)))
  cat(lab, "->", r, "\n")
}
run("mvbind student rescor", quote(frm_sample(
  bf(mvbind(y1, y2) ~ x) + set_rescor(TRUE), data = D$mv,
  family = frmtmb::brmsfamily("student"), chains = 1, iter = 200,
  refresh = 0, seed = 1)))
run("mvbf student rescor", quote(frm_sample(
  mvbf(bf(y1 ~ x), bf(y2 ~ x)) + set_rescor(TRUE), data = D$mv,
  family = frmtmb::brmsfamily("student"), chains = 1, iter = 200,
  refresh = 0, seed = 1)))
run("fit route mvbind student rescor", quote(frm_sample(
  frm(bf(mvbind(y1, y2) ~ x) + set_rescor(TRUE), data = D$mv,
      family = frmtmb::brmsfamily("student")), chains = 1, iter = 200,
  refresh = 0, seed = 1)))
run("mvbind gaussian rescor", quote(frm_sample(
  bf(mvbind(y1, y2) ~ x) + set_rescor(TRUE), data = D$mv,
  family = gaussian(), chains = 1, iter = 200, refresh = 0, seed = 1)))
dp <- tryCatch(default_prior(bf(mvbind(y1, y2) ~ x) + set_rescor(TRUE),
                             data = D$mv,
                             family = frmtmb::brmsfamily("student"),
                             route = "sample"),
               error = function(e) conditionMessage(e))
print(dp)
