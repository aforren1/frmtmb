# Does frm() say anything on p082's hsgp fit (lscale at 2e-05)?
source("dev/samplers-designs.R")
samplers_lib()
eq <- samplers_eq_env()
spec <- samplers_design("p082", eq)
w <- character(); m <- character()
fit <- withCallingHandlers(eq$eq_frm_fit(spec),
  warning = function(x) { w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning") },
  message = function(x) { m <<- c(m, conditionMessage(x)); invokeRestart("muffleMessage") })
cat("warnings:", length(w), "\n"); print(w); cat("messages:", length(m), "\n"); print(m)
cat("class:", class(fit), "\n")
print(summary(fit))
