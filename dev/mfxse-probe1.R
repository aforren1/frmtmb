# Probe: which interop methods marginaleffects calls, with which dots, and the
# layout of a mixed fit's joint covariance.
source("C:/Users/adf44/source/r/frmtmb-wt-mfxse/dev/mfxse-lib.R")
set.seed(7411)
n <- 200
d <- data.frame(x = rnorm(n), g = factor(sample(1:8, n, TRUE)))
d$y <- rpois(n, exp(0.5 + 0.4 * d$x + rnorm(8, 0, 0.8)[d$g]))
fit <- quiet(frm(y ~ x + (1 + x | g), family = poisson(), data = d))
ns <- asNamespace("frmtmb")
cat("estimates:", paste(names(fit$estimates), lengths(fit$estimates)), "\n")
jc <- ns$get_joint_cov(fit)
print(table(jc$names))
cat("outer par:", names(fit$obj$par), "\n")
cat("random idx len:", length(fit$obj$env$random), "\n")
print(str(fit$frame$re_blocks[[1]][c("covstruct", "theta_idx", "b_idx")]))
cat("b head:", head(fit$estimates$b), "\n")
log <- new.env()
log$calls <- list()
wrap <- function(nm) {
  f0 <- get(paste0(nm, ".frmtmb_fit"), envir = ns)
  function(model, ...) {
    dn <- names(list(...))
    log$calls[[length(log$calls) + 1L]] <- paste(nm, paste(dn, collapse = ","))
    f0(model, ...)
  }
}
for (nm in c("get_coef", "get_vcov")) {
  registerS3method(nm, "frmtmb_fit", wrap(nm),
                   envir = asNamespace("marginaleffects"))
}
run <- function(lab, e) {
  log$calls <- list()
  r <- quiet(e)
  cat("\n==", lab, "\n")
  print(table(unlist(log$calls)))
  print(head(as.data.frame(r)[, intersect(c("term", "estimate", "std.error"),
                                          names(r))], 4))
}
run("predictions", predictions(fit, newdata = d[1:3, ]))
run("predictions re_formula NA", predictions(fit, newdata = d[1:3, ],
                                             re_formula = NA))
run("avg_slopes", avg_slopes(fit, variables = "x"))
run("hypotheses", hypotheses(fit))
run("hypotheses b2", hypotheses(fit, hypothesis = "b2 = 0"))
run("avg_predictions by", avg_predictions(fit, by = "g"))
