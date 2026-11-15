source("dev/refits-lib.R")
set.seed(5)
dm <- data.frame(x = rnorm(200), N = 10L)
pr <- cbind(1, exp(0.3 + 0.5 * dm$x), exp(-0.2 * dm$x))
pr <- pr / rowSums(pr)
dm$Y <- t(sapply(1:200, function(i) rmultinom(1, 10, pr[i, ])))
colnames(dm$Y) <- c("a", "b", "c")
fm <- frm(Y | trials(N) ~ x, family = multinomial(), data = dm)
md <- marginaleffects:::get_modeldata(fm)
str(md)
r <- tryCatch(withCallingHandlers(marginaleffects::predictions(fm),
  error = function(e) {
    cs <- sys.calls()
    for (k in seq_along(cs)) cat(k, substr(deparse1(cs[[k]]), 1, 150), "\n")
  }), error = function(e) conditionMessage(e))
print(r)
# brms's own route: does it read the response matrix?
cat("\ninsight::get_data:\n")
str(insight::get_data(fm))
