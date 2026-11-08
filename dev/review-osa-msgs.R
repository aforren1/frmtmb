# Review lane osa: refusal wording for bernoulli()/nbinom1() with cens() and
# trunc(), what brms 2.23.0 says on the same calls, quadrature = TRUE
# fits, and glmmTMB's reading of weights in simulate().
#
#   Rscript dev/review-osa-msgs.R before|after
arm <- commandArgs(TRUE)[1]
.libPaths(c(if (arm == "after") "C:/Users/adf44/source/r/wt-osa-lib",
            "C:/Users/adf44/source/r/rellib-r18",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib: frmtmb", as.character(packageVersion("frmtmb")),
    dirname(find.package("frmtmb")), "\n")
cat("BEGIN review-osa-msgs", arm, "\n")
msg <- function(expr) {
  tryCatch({ force(expr); "no error" }, error = function(e) {
    gsub("\n", " ", conditionMessage(e))
  })
}
set.seed(1)
d <- data.frame(y = rep(0:1, 20), x = rnorm(40), cc = 0, n = rpois(40, 3))
cat("bernoulli cens():", msg(frm(y | cens(cc) ~ x, family = bernoulli(), data = d)), "\n")
cat("bernoulli trunc(ub = 1):", msg(frm(y | trunc(ub = 1) ~ x, family = bernoulli(), data = d)), "\n")
cat("nbinom1 cens():", msg(frm(n | cens(cc) ~ x, family = nbinom1(), data = d)), "\n")
cat("nbinom1 trunc(lb = 0):", msg(frm(n | trunc(lb = 0) ~ x, family = nbinom1(), data = d)), "\n")
cat("compat cens() x bernoulli:", frm_compat("cens()", "bernoulli")$status, "|",
    frm_compat("cens()", "bernoulli")$note, "\n")
cat("compat trunc() x nbinom1:", frm_compat("trunc()", "nbinom1")$status, "|",
    frm_compat("trunc()", "nbinom1")$note, "\n")
cat("compat residuals_osa x von_mises:", frm_compat("von_mises", "residuals_osa")$status, "\n")
cat("compat residuals_osa x tweedie:", frm_compat("tweedie", "residuals_osa")$status, "\n")
cat("compat residuals_osa x huber:",
    msg(print(frm_compat("huber", "residuals_osa")$status)), "\n")

# brms 2.23.0, no compile
suppressPackageStartupMessages(library(brms))
cat("brms", as.character(packageVersion("brms")), "\n")
cat("brms bernoulli cens():", msg(brms::stancode(y | cens(cc) ~ x, data = d, family = bernoulli())), "\n")
cat("brms bernoulli trunc(ub = 1):", msg(brms::stancode(y | trunc(ub = 1) ~ x, data = d, family = bernoulli())), "\n")

# quadrature = TRUE
set.seed(2)
dq <- data.frame(g = factor(rep(1:15, each = 8)), x = rnorm(120))
dq$y <- rpois(120, exp(0.3 + 0.4 * dq$x + rnorm(15)[dq$g]))
fq <- frm(y ~ x + (1 | g), family = poisson(), data = dq, quadrature = TRUE)
cat("quadrature = TRUE, residuals(type = 'osa'):",
    msg(residuals(fq, type = "osa")), "\n")
dqg <- dq
dqg$y <- 1 + 0.4 * dq$x + rnorm(15)[dq$g] + rnorm(120)
fqg <- tryCatch(frm(y ~ x + (1 | g), family = Gamma(link = "log"),
                    data = transform(dqg, y = exp(y / 3)), quadrature = TRUE),
                error = function(e) e)
if (!inherits(fqg, "error")) {
  cat("quadrature = TRUE Gamma, residuals(type = 'osa'):",
      msg(residuals(fqg, type = "osa")), "\n")
}

# glmmTMB's simulate() with weights on a poisson
suppressPackageStartupMessages(library(glmmTMB))
dw <- data.frame(x = rnorm(60))
dw$y <- rpois(60, exp(0.5 + 0.3 * dw$x))
dw$w <- rep(c(1, 3), 30)
gw <- glmmTMB(y ~ x, family = poisson, data = dw, weights = w)
ws <- character(0)
withCallingHandlers(invisible(simulate(gw, nsim = 1, seed = 1)),
                    warning = function(w) {
                      ws <<- c(ws, conditionMessage(w))
                      invokeRestart("muffleWarning")
                    })
cat("glmmTMB", as.character(packageVersion("glmmTMB")),
    "simulate() on a weighted poisson warns:",
    if (length(ws)) ws else "nothing", "\n")
cat("END review-osa-msgs\n")
