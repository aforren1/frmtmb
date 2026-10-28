# Lane reslot: the random-effect and rescor double prior slot, measured.
# Usage: Rscript dev/reslot-double.R before|after
# Seed 1 (re_data()), block (1 + x + z | g), 8 groups of 10.
source("C:/Users/adf44/source/r/frmtmb-wt-reslot/dev/reslot-common.R")
suppressPackageStartupMessages(library(frmtmb.sample))
dd <- re_data(1L)
fit <- frm(bf(y ~ x + z + (1 + x + z | g)), family = gaussian(), data = dd)
cat("theta names:", names(fit$estimates$theta), "\n")
pars <- fit$obj$env$parList(fit$obj$env$last.par.best)
nlp <- function(pl) {
  ri <- frmtmb:::resolve_prior_input(fit, pl)
  frmtmb:::neg_log_prior_fn(ri$entries)(pars)
}
lkj <- set_prior("lkj(2)", class = "cor")
th4 <- set_prior("normal(0, 1)", class = "theta", coef = "theta_4")
try_show("frm: cor lkj(2) alone", frmtmb:::resolve_prior_input(fit, lkj))
try_show("frm: cor lkj(2) then theta_4 normal(0, 1)",
         frmtmb:::resolve_prior_input(fit, lkj + th4))
try_show("frm: theta_4 then cor", frmtmb:::resolve_prior_input(fit, th4 + lkj))
cat("neg log prior at the ML estimate: lkj alone",
    sprintf("%.12f", nlp(lkj)), " lkj + theta_4",
    tryCatch(sprintf("%.12f", nlp(lkj + th4)), error = function(e) "ERROR"),
    "\n")
sdx <- set_prior("student_t(3, 0, 1)", class = "sd", group = "g", coef = "x")
th2 <- set_prior("normal(0, 1)", class = "theta", coef = "theta_2")
try_show("frm: sd g x then theta_2", frmtmb:::resolve_prior_input(fit,
                                                                   sdx + th2))
try_show("frm: class theta, no coef",
         frmtmb:::resolve_prior_input(fit, set_prior("normal(0, 2)",
                                                     class = "theta")))
try_show("frm: theta_5 lb/ub box on one correlation",
         frmtmb:::resolve_prior_input(fit, set_prior("", class = "theta",
                                                     coef = "theta_5",
                                                     lb = -1, ub = 1)))
# MAP fits: the order of the two specifications decides the model
for (pr in list(lkj + th4, th4 + lkj)) {
  m <- tryCatch(frm(bf(y ~ x + z + (1 + x + z | g)), family = gaussian(),
                    data = dd, prior = pr), error = function(e) e)
  if (inherits(m, "error")) {
    cat("MAP ERROR:", conditionMessage(m), "\n")
  } else {
    cat("MAP objective", sprintf("%.10f", m$opt$objective), " theta:",
        sprintf("%.6f", m$estimates$theta), "\n")
  }
}

cat("\n== frm_sample's resolved stack (no sampling)\n")
sr <- function(pr, f = fit) {
  r <- suppressMessages(frmtmb.sample:::sample_resolve_priors(f, pr))
  show_entries(r$ri)
  invisible(NULL)
}
try_show("sample: defaults only", invisible(sr(NULL)))
try_show("sample: theta_4 normal(0, 1)", invisible(sr(th4)))
try_show("sample: sd g x + theta_2", invisible(sr(sdx + th2)))

cat("\n== default_prior() tables\n")
print(default_prior(bf(y ~ x + z + (1 + x + z | g)), data = dd,
                    family = gaussian()))
print(default_prior(bf(y ~ x + z + (1 + x + z | g)), data = dd,
                    family = gaussian(), route = "sample"))

cat("\n== rescor\n")
mv <- frm(bf(y ~ x) + bf(y2 ~ x) + set_rescor(TRUE), data = dd)
cat("thetar names:", names(mv$estimates$thetar), "\n")
rl <- set_prior("lkj(2)", class = "rescor")
tr1 <- set_prior("normal(0, 1)", class = "theta", coef = "thetar_1")
try_show("frm: rescor lkj(2) then thetar_1",
         frmtmb:::resolve_prior_input(mv, rl + tr1))
try_show("frm: thetar_1 then rescor", frmtmb:::resolve_prior_input(mv,
                                                                   tr1 + rl))
try_show("sample: thetar_1 alone", sr(tr1, mv))
print(default_prior(bf(y ~ x) + bf(y2 ~ x) + set_rescor(TRUE), data = dd))
