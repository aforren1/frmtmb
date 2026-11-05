# band = "boot" (the refusal's remedy), the large-response nl start, both arms. Seed 9901.
source("dev/review-nlfix-common.R")
set.seed(9901)
n <- 200
d <- data.frame(x = rnorm(n), z = rnorm(n), g = factor(rep(1:10, each = 20)))
d$yo <- cut(0.8 * d$x + rlogis(n), c(-Inf, -1, 0, 1, Inf), labels = 1:4, ordered_result = TRUE)
d$y <- 1 + 0.5 * d$x + rnorm(n)
d$xs <- d$x * 1e-10
d$ybig <- 1e6 + 5e4 * d$x + rnorm(n, 0, 1e3)
ll <- function(r) if (iserr(r)) paste("ERR", substr(r$msg, 1, 200)) else paste0(format(as.numeric(logLik(r)), digits = 10), " conv=", r$opt$convergence, " beta=", paste(format(r$estimates$beta, digits = 4), collapse = ","), " w=", length(attr(r, "w")))
cat("--- large response\n")
cat("linear y ~ x    :", ll(tt(frm(ybig ~ x, data = d))), "\n")
cat("nl a + b * x    :", ll(tt(frm(bf(ybig ~ a + b * x, a + b ~ 1, nl = TRUE), data = d))), "\n")
cat("nl a + b * xs   :", ll(tt(frm(bf(ybig ~ a + b * xs, a + b ~ 1, nl = TRUE), data = d))), "\n")
cat("--- band = boot\n")
try_boot <- function(lab, fit, ...) {
  x <- tryCatch(withCallingHandlers(conditional_effects(fit, effects = "x", band = "boot", boot = 8, seed = 1, ...), message = function(m) invokeRestart("muffleMessage")), error = function(e) e)
  cat(sprintf("%-26s %s\n", lab, if (inherits(x, "error")) paste("ERR", conditionMessage(x), "| call:", deparse1(conditionCall(x))[1]) else paste("ok rows", nrow(x[[1]]))))
}
if (arm == "after") {
  fa <- frm(bf(yo ~ b1 * x, b1 ~ 1, nl = TRUE), family = cumulative(), data = d)
  try_boot("ordinal body", fa)
  try_boot("ordinal body categorical", fa, categorical = TRUE)
}
fl <- frm(yo ~ 0 + x, family = cumulative(), data = d)
try_boot("ordinal linear", fl)
gg <- suppressWarnings(frm(bf(y ~ b0 + b1 * x, b0 + b1 ~ 1, nl = TRUE), data = d))
try_boot("gaussian body", gg)
gl <- frm(y ~ x, data = d)
try_boot("gaussian linear", gl)
pb <- frm(bf(cnt ~ a + b * x, a + b ~ 1, nl = TRUE), family = poisson(), data = transform(d, cnt = rpois(n, 2)))
try_boot("poisson body", pb)
cat("--- the refusal's text, ordinal body, in full\n")
if (arm == "after") {
  x <- tryCatch(conditional_effects(fa, effects = "x"), error = function(e) conditionMessage(e)); cat(x, "\n")
  cat("--- method = predict, named too\n")
  x <- tryCatch(conditional_effects(fa, effects = "x", method = "predict"), error = function(e) conditionMessage(e)); cat(if (is.character(x)) x else "ok", "\n")
  x <- tryCatch(conditional_effects(gg, effects = "x", method = "predict"), error = function(e) conditionMessage(e)); cat(if (is.character(x)) x else "ok gaussian body", "\n")
}
