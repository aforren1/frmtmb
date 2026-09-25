source("C:/Users/adf44/source/r/frmtmb-wt-predfix/dev/predfix-prelude.R")
set.seed(1); n <- 200; x <- rnorm(n) * 1e-4 + 5e-4
d <- data.frame(x = x, y = 1 + 2000 * (x - 5e-4) + rnorm(n))
for (pr in list(set_prior("", class = "Intercept", ub = 0.5), set_prior("", class = "b", coef = "x", lb = 3e4))) {
  msgs <- character(0)
  f <- withCallingHandlers(suppressWarnings(frm(bf(y ~ x), family = gaussian(), data = d, prior = pr, control = frmtmb_control(verbose = TRUE))),
    message = function(m) { msgs <<- c(msgs, conditionMessage(m)); invokeRestart("muffleMessage") })
  cat(grep("pre-fit", msgs, value = TRUE), sep = "\n"); cat("engaged", !is.null(f$par_units), "\n")
}
