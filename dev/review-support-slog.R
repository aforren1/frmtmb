# Reviewer probe (lane support): shifted_lognormal under a non-binding
# trunc(lb <= 0), and the Frechet's ub = 1e6 alone. Seed 402 as in
# dev/review-support-trunc.R.
lib <- Sys.getenv("FRMTMB_PORT_LIB", "C:/Users/adf44/source/r/rellib-r19")
.libPaths(c(strsplit(lib, ";")[[1]],
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("lib:", dirname(system.file(package = "frmtmb")), "\n")
set.seed(402)
n <- 300
e <- rexp(n)
d <- data.frame(y = 2 * e + 0.05, lb = -1, lb0 = 0, ub = 1e6)
f0 <- tryCatch(frm(y ~ 1, family = shifted_lognormal(), data = d),
               error = function(e) conditionMessage(e))
print(if (is.character(f0)) f0 else logLik(f0))
for (lbn in c("lb", "lb0")) {
  f1 <- tryCatch(frm(as.formula(paste0("y | trunc(lb = ", lbn, ") ~ 1")),
                     family = shifted_lognormal(), data = d),
                 error = function(e) conditionMessage(e))
  print(if (is.character(f1)) f1 else c(logLik(f1), fixef(f1)[, 1]))
}
# the slot itself at a q below ndt, its value and its taped gradient
fam <- shifted_lognormal()
print(fam$cdf_log(c(-1, 0, 0.01), list(mu = 0, sigma = 1, ndt = 0.02), list()))
g <- RTMB::MakeTape(function(p) {
  sum(fam$cdf_log(c(-1, 0), list(mu = p[1], sigma = p[2], ndt = p[3]),
                  list())$lS)
}, c(0, 1, 0.02))
print(g$jacobian(c(0, 1, 0.02)))
fr0 <- frm(y ~ 1, family = frechet(), data = d)
fr1 <- frm(y | trunc(ub = ub) ~ 1, family = frechet(), data = d)
cat("frechet ub-only dlogLik", as.numeric(logLik(fr1) - logLik(fr0)),
    " nu", exp(fixef(fr0)[2, 1]), "\n")
print(fixef(fr0))
