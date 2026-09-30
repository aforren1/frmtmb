source("dev/defects-pre.R")
fit1 <- brms_fixture(1)
nd <- fit1$data[1:5, ]
print(str(nd))
options(warn = 2)
r <- tryCatch(fitted(fit1, newdata = nd), error = function(e) {
  print(conditionCall(e))
  print(sys.calls())
})
withCallingHandlers(fitted(fit1, newdata = nd), warning = function(w) {
  cs <- sys.calls()
  print(lapply(cs[max(1, length(cs) - 12):length(cs)],
               function(x) deparse1(x)[1]))
  invokeRestart("muffleWarning")
})
