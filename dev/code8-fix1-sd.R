# Fixture 1: why the smooth's sd at exp(-354) is "nonfinite" rather than
# at its edge.
source("dev/code8-common.R")
fit0 <- get("brmsfit_example1", envir = asNamespace("brms"))
d <- as.data.frame(fit0$data); attr(d, "terms") <- NULL
form <- bf(count ~ Trt * Age + mo(Exp) + s(Age) + volume + offset(Age) +
             (1 + Trt | visit) + arma(visit, patient), sigma ~ Trt)
f <- quiet(frm(form, family = student(), data = d))
cat("log_sd_theta_index:", fns$log_sd_theta_index(f), "\n")
cat("se_at_edge theta_4:", fns$se_at_edge(f, "theta_4", -1), "\n")
g <- fns$se_sd_gain_up(f, 4)
cat("se_sd_gain_up(4):", format(g), " moved:", format(attr(g, "moved")),
    " at:", format(attr(g, "at")), "\n")
cat("se_edge_tol:", fns$se_edge_tol, "\n")
cat("fitted sigma with conditions:\n")
r <- withCallingHandlers(fitted(f, dpar = "sigma"),
  warning = function(w) { cat("W:", conditionMessage(w), "\n");
    invokeRestart("muffleWarning") })
print(head(r, 3))
s <- withCallingHandlers(summary(f), warning = function(w) {
  cat("W:", conditionMessage(w), "\n"); invokeRestart("muffleWarning") })
print(s)
