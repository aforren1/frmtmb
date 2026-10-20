# Reviewer of lane sampleprior, fix round 2: frm_sample() from a formula,
# warnings not suppressed, tmbstan captured.
source("C:/Users/adf44/source/r/frmtmb-wt-sampleprior/dev/sampleprior-rev-common.R")
source("C:/Users/adf44/source/r/frmtmb-wt-sampleprior/dev/sampleprior-rev-cases.R")
n <- 0L
capture_on()
invisible(withCallingHandlers(suppressMessages(tryCatch(frm_sample(y ~ x, data = D$wide,
  family = gaussian(), prior = set_prior("uniform(-10, 10)", class = "b", coef = "x"),
  chains = 1, iter = 10, refresh = 0), error = function(e) NULL)),
  warning = function(w) { if (grepl("bounded prior", conditionMessage(w))) n <<- n + 1L
    invokeRestart("muffleWarning") }))
capture_off()
cat("frm_sample(formula): bound warnings", n, "\n")
