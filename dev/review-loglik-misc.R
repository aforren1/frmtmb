# Review of lane loglik: paths that re-resolve the prior after the fit
# (update(), a brms prior object), REML with a prior on sigma's
# coefficients checked on the tape, and logLik() conditions. Seed 601.
.libPaths(c("C:/Users/adf44/source/r/wt-loglik-lib",
            "C:/Users/adf44/source/r/rellib-r18",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("lib:", dirname(find.package("frmtmb")), format(packageVersion("frmtmb")),
    "\n")
tape_rel <- function(f0, fp) {
  stopifnot(identical(names(f0$opt$par), names(fp$opt$par)))
  v <- -f0$obj$fn(fp$opt$par)
  f0$obj$fn(f0$opt$par)
  abs(as.numeric(logLik(fp)) - v) / abs(v)
}
set.seed(601)
d <- data.frame(x = rnorm(100), z = rnorm(100),
                g = factor(rep(1:20, each = 5)))
d$y <- 1 + 0.5 * d$x + rnorm(20, 0, 0.7)[d$g] +
  rnorm(100, 0, exp(0.2 * d$z))

# REML, prior on sigma's coefficients: answered; is it right?
pr <- set_prior("normal(0, 0.1)", class = "b", dpar = "sigma") +
  set_prior("normal(0, 0.3)", class = "Intercept", dpar = "sigma")
fp <- frm(bf(y ~ x + (1 | g), sigma ~ z), data = d, REML = TRUE, prior = pr)
f0 <- frm(bf(y ~ x + (1 | g), sigma ~ z), data = d, REML = TRUE)
cat(sprintf("REML dpar sigma: logLik %.10f  rel to tape %.2e  comps read: %s\n",
            as.numeric(logLik(fp)), tape_rel(f0, fp),
            paste(unique(unlist(lapply(frmtmb:::resolve_prior_input(
              list(frame = fp$frame, spec = fp$spec), fp$prior)$entries,
              function(e) c(e$comp, e$offset$comp)))), collapse = ",")))

# a brms prior object
if (requireNamespace("brms", quietly = TRUE)) {
  bp <- brms::prior(normal(0, 0.2), class = b) +
    brms::prior(normal(0, 0.5), class = sd)
  fb <- frm(y ~ x + z + (1 | g), data = d, prior = bp)
  f0 <- frm(y ~ x + z + (1 | g), data = d)
  cat(sprintf("brmsprior: class(fit$prior) %s  logLik rel to tape %.2e\n",
              paste(class(fb$prior), collapse = "/"), tape_rel(f0, fb)))
}

# update() that drops a coefficient a prior names
fx <- frm(y ~ x + z + (1 | g), data = d,
          prior = set_prior("normal(0, 0.2)", class = "b", coef = "z") +
            set_prior("normal(0, 0.5)", class = "sd"))
msgs <- character(0)
fu <- withCallingHandlers(update(fx, . ~ . - z),
                          message = function(m) {
                            msgs <- c(msgs, conditionMessage(m))
                            invokeRestart("muffleMessage")
                          })
f0 <- frm(y ~ x + (1 | g), data = d)
cat(sprintf("update(-z): logLik rel to tape %.2e  log_prior %.6f\n",
            tape_rel(f0, fu), attr(logLik_penalized(fu), "log_prior")))
nc <- 0L
withCallingHandlers(for (i in 1:3) logLik(fu),
                    message = function(m) { nc <<- nc + 1L
                      invokeRestart("muffleMessage") },
                    warning = function(w) { nc <<- nc + 1L
                      invokeRestart("muffleWarning") })
cat("conditions from 3 logLik() calls after update():", nc, "\n")

# a fit saved and read back: the prior resolves from the stored frame
tf <- tempfile(fileext = ".rds")
saveRDS(fx, tf)
fr <- readRDS(tf)
cat(sprintf("readRDS: logLik equal %s\n",
            identical(as.numeric(logLik(fr)), as.numeric(logLik(fx)))))
unlink(tf)
cat("done\n")
