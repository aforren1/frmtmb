source("C:/Users/adf44/source/r/frmtmb-wt-predfix/dev/predfix-review2/r2-helpers.R")
# B1: an error raised inside the pre-fit's fit_assembled() must reach the
# caller unchanged. Injected after the frame is scaled, so it fires
# inside the muffling handler. Seed 506.
set.seed(506)
d <- data.frame(x = rnorm(100) * 1e-4); d$y <- 1 + 1e4 * d$x + rnorm(100)
suppressMessages(trace("build_objective", where = asNamespace("frmtmb"),
  tracer = quote(if (any(abs(frame$linpreds[[1]]$X[, 2]) > 0.5))
    stop("r2 injected pre-fit error")), print = FALSE))
e <- tryCatch(frm(y ~ x, data = d), error = function(e) e,
              warning = function(w) paste("WARNING instead:", conditionMessage(w)))
cat("class:", class(e), "\nmessage:", conditionMessage(e), "\n")
f <- frm(y ~ x, data = d, control = frmtmb_control(autoscale = FALSE))
cat("FALSE unaffected, logLik", as.numeric(logLik(f)), "\n")
