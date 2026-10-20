# Reviewer of lane sampleprior, fix round 2: where the bound warning
# fires (brms 2.23.0 and frmtmb), and how many times per call.
source("C:/Users/adf44/source/r/frmtmb-wt-sampleprior/dev/sampleprior-rev-common.R")
source("C:/Users/adf44/source/r/frmtmb-wt-sampleprior/dev/sampleprior-rev-cases.R")
d <- D$wide
count <- function(expr) {
  n <- 0L
  r <- withCallingHandlers(tryCatch(expr, error = function(e)
    paste("ERROR", substr(conditionMessage(e), 1, 70))),
    warning = function(w) {
      if (grepl("bounded prior", conditionMessage(w))) n <<- n + 1L
      invokeRestart("muffleWarning")
    })
  list(n = n, err = if (is.character(r) && length(r) == 1L &&
                          startsWith(r, "ERROR")) r else "")
}
show <- function(lab, expr) {
  x <- count(expr)
  cat(sprintf("%-46s bound warnings %d %s\n", lab, x$n, x$err))
}
show("brms: bare set_prior(uniform on b)",
     brms::set_prior("uniform(2.2, 2.6)", class = "b", coef = "x"))
bp <- brms::set_prior("uniform(2.2, 2.6)", class = "b", coef = "x")
show("brms: validate_prior()", brms::validate_prior(bp, y ~ x, data = d))
show("brms: stancode()", brms::stancode(y ~ x, data = d, prior = bp))
show("frmtmb: bare set_prior()",
     set_prior("uniform(2.2, 2.6)", class = "b", coef = "x"))
pr <- set_prior("uniform(-10, 10)", class = "b", coef = "x")
show("frmtmb: frm()", frm(y ~ x, data = d, prior = pr))
show("frmtmb: validate_prior()", validate_prior(pr, y ~ x, data = d))
show("frmtmb: par_template()", par_template(y ~ x, data = d, prior = pr))
f0 <- frm(y ~ x, data = d)
show("frmtmb: update(fit, prior =)", update(f0, prior = pr))
show("frmtmb: frm_sample(formula) [captured]", cap_obj(y ~ x, d, gaussian(), prior = pr))
show("frmtmb: frm_sample(fit) [captured]", {
  capture_on(); on.exit(capture_off())
  tryCatch(frm_sample(f0, prior = pr, chains = 1, iter = 10, refresh = 0),
           error = function(e) NULL)
})
show("frmtmb: frm_simulate()", frm_simulate(y ~ x, data = d, prior = pr,
                                            nsim = 1, seed = 1))
fm <- suppressWarnings(frm(y ~ x, data = d, prior = pr))
show("frmtmb: update(MAP fit) with no new prior", update(fm, data = d))
