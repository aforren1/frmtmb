# Lane optima: what the SE-check surface says on the flat fixture of
# tests/testthat/helper-mo-flat.R, to write the replacement assertions
# from a run rather than from expectation.
.libPaths(c("C:/Users/adf44/source/r/wt-optima-lib",
            "C:/Users/adf44/source/r/rellib-r6",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
source("tests/testthat/helper-mo-flat.R")
cap <- function(expr) {
  w <- character()
  v <- withCallingHandlers(expr, warning = function(x) {
    w <<- c(w, conditionMessage(x))
    invokeRestart("muffleWarning")
  })
  list(value = v, warnings = w)
}
for (fo in list(ls ~ mo(inc) + age, ls ~ mo(inc) * age)) {
  cat("=====", deparse(fo), "\n")
  r <- cap(frm(fo, data = mo_flat_data(1)))
  fit <- r$value
  cat("warnings:", length(r$warnings), "\n")
  cat(r$warnings, sep = "\n")
  print(frmtmb:::sdr_of(fit)$se_lost)
  print(fixef(fit))
  ci <- confint(fit)
  print(ci)
  print(diagnose(fit, quiet = TRUE)$flat)
  for (v in c("age", "inc", "inc:age")) {
    rr <- tryCatch(cap(conditional_effects(fit, v)),
                   error = function(e) list(value = NULL,
                                            warnings = conditionMessage(e)))
    se <- if (!is.null(rr$value)) rr$value[[1]]$se__ else NA
    cat("conditional_effects", v, ": se finite", sum(is.finite(se)), "of",
        length(se), "; inc values", if (!is.null(rr$value))
          paste(unique(rr$value[[1]]$inc), collapse = " "),
        "; warnings", length(rr$warnings), substr(rr$warnings, 1, 80), "\n")
  }
  s <- capture.output(print(summary(fit)))
  cat(grep("without a standard error|zeta", s, value = TRUE), sep = "\n")
  cat("vcov warnings:", length(cap(vcov(fit))$warnings), "\n")
}
r <- cap(frm(ls ~ mo(inc) + age, data = mo_flat_data(1),
             control = frmtmb_control(check_se = "ignore")))
cat("check_se = ignore warnings:", length(r$warnings), "\n")
r <- tryCatch(frm(ls ~ mo(inc) + age, data = mo_flat_data(1),
                  control = frmtmb_control(check_se = "stop")),
              error = function(e) conditionMessage(e))
cat("check_se = stop:", substr(r, 1, 80), "\n")
r <- cap(frm(ls ~ mo(inc) + age, data = mo_flat_data(1), se = TRUE))
cat("se = TRUE warnings:", length(r$warnings), "\n")
r <- cap(frm(ls ~ mo(inc) + age, data = mo_flat_data(1),
             control = frmtmb_control(optCtrl = list(iter.max = 3,
                                                     eval.max = 5),
                                      restarts = 0)))
cat("stopped short: warnings", length(r$warnings), ":",
    substr(r$warnings, 1, 60), "\n")
# the healthy absent case: the same data with every category observed
set.seed(7)
f7 <- cap(frm(ls ~ mo(inc), data = transform(mo_flat_data(1),
                                             inc = sample(0:2, 100, TRUE))))
cat("all categories observed: warnings", length(f7$warnings), "\n")
