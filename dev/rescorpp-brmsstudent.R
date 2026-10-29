## Lane rescorpp: brms 2.23.0 on ONE response of a Student-t rescor fit.
##   Rscript dev/rescorpp-brmsstudent.R
.libPaths("C:/Users/adf44/AppData/Local/R/win-library/4.6")
suppressMessages(library(brms))
cat("brms", format(packageVersion("brms")), "\n")
fd <- "C:/Users/adf44/source/r/frmtmb-wt-rescorpp/dev/rescorpp-brmsfits"
for (w in c("g3", "st")) {
  b <- readRDS(file.path(fd, paste0(w, ".rds")))
  cat("\n########", w, "\n")
  set.seed(1)
  p1 <- posterior_predict(b, resp = "y1")
  cat("posterior_predict(resp = 'y1') dim:", dim(p1), " dimnames NULL:",
      is.null(dimnames(p1)), "\n")
  prep <- brms::prepare_predictions(b, resp = "y1")
  cat("prep class:", class(prep), " family:", prep$family$family, "\n")
  cat("names(prep$dpars):", names(prep$dpars), "\n")
  if (w == "st") {
    nu <- brms:::get_dpar(prep, "nu", i = 1)
    cat("str(nu):\n"); str(nu)
    dr <- as_draws_df(b)
    cat("nu draws head:", head(dr$nu), "\n")
    d <- b$data
    L1 <- log_lik(b, resp = "y1")
    m1 <- outer(dr$b_y1_Intercept, rep(1, nrow(d))) + outer(dr$b_y1_x, d$x)
    yy <- matrix(d$y1, nrow(m1), nrow(d), byrow = TRUE)
    z <- (yy - m1) / dr$sigma_y1
    for (cand in list(list("nu draws", dr$nu), list("nu = 1", 1),
                      list("gaussian", Inf))) {
      ref <- if (is.infinite(cand[[2]][1])) {
        dnorm(z, log = TRUE) - log(dr$sigma_y1)
      } else dt(z, cand[[2]], log = TRUE) - log(dr$sigma_y1)
      cat("log_lik(resp='y1') vs", cand[[1]], ": max abs diff",
          max(abs(L1 - ref)), "\n")
    }
    # the spread of posterior_predict(resp = 'y1') against the joint draw's
    # y1 margin, at the same draws: a quantile ratio
    pj <- posterior_predict(b)[, , 1]
    ep <- posterior_epred(b, resp = "y1")
    q <- function(e) quantile(abs(e), c(0.5, 0.9, 0.99))
    cat("|pp - epred| quantiles, resp = 'y1':", q(p1 - ep), "\n")
    cat("|pp - epred| quantiles, joint y1   :", q(pj - ep), "\n")
  }
}
