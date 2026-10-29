## Lane rescorpp: what brms 2.23.0 does on rescor fits, measured.
##   Rscript dev/rescorpp-brmsbehav.R
.libPaths("C:/Users/adf44/AppData/Local/R/win-library/4.6")
suppressMessages(library(brms))
cat("brms", format(packageVersion("brms")), "\n")
options(width = 120)
fd <- "C:/Users/adf44/source/r/frmtmb-wt-rescorpp/dev/rescorpp-brmsfits"
show <- function(tag, expr) {
  cat("\n==", tag, "==\n")
  w <- character(0)
  r <- withCallingHandlers(tryCatch(expr, error = function(e) e),
                           warning = function(x) {
                             w <<- c(w, conditionMessage(x))
                             invokeRestart("muffleWarning")
                           })
  if (inherits(r, "error")) cat("ERROR:", conditionMessage(r), "\n")
  else print(r)
  for (m in w) cat("WARNING:", substr(m, 1, 300), "\n")
  invisible(r)
}
shape <- function(a) list(dim = dim(a), dimnames3 = dimnames(a)[[3]])
# the mean over rows of the across-draw correlation of two responses'
# predictive residuals (prediction minus expectation)
rcor <- function(pp, ep, j, k) {
  e <- pp - ep
  mean(vapply(seq_len(dim(e)[2]), function(i) cor(e[, i, j], e[, i, k]), 1))
}
for (w in c("g3", "gs", "st", "g2f")) {
  cat("\n\n########", w, "########\n")
  b <- readRDS(file.path(fd, paste0(w, ".rds")))
  d <- b$data
  set.seed(1)
  pp <- show("posterior_predict() shape", shape(posterior_predict(b)))
  show("posterior_predict(resp = 'y1') shape",
       shape(posterior_predict(b, resp = "y1")))
  show("posterior_epred() shape", shape(posterior_epred(b)))
  show("posterior_linpred() shape", shape(posterior_linpred(b)))
  show("predictive_error() shape", shape(predictive_error(b)))
  show("fitted() shape", shape(fitted(b)))
  show("predict() shape", shape(predict(b)))
  show("residuals() shape", shape(residuals(b)))
  show("predictive_interval()", dim(predictive_interval(b)))
  show("log_lik() dim", dim(log_lik(b)))
  show("log_lik(combine = FALSE) dim", dim(log_lik(b, combine = FALSE)))
  show("log_lik(resp = 'y1') dim", dim(log_lik(b, resp = "y1")))
  show("bayes_R2()", bayes_R2(b))
  show("pp_check() without resp", class(pp_check(b, ndraws = 5)))
  show("pp_check(resp = 'y1')", class(pp_check(b, resp = "y1", ndraws = 5)))
  show("conditional_effects(method = 'posterior_predict') names",
       names(conditional_effects(b, method = "posterior_predict")))
  ## marginal or joint: log_lik(resp = 'y1') against the univariate density
  dr <- as_draws_df(b)
  L1 <- log_lik(b, resp = "y1")
  if (w %in% c("g3", "g2f")) {
    m1 <- outer(dr$b_y1_Intercept, rep(1, nrow(d))) + outer(dr$b_y1_x, d$x)
    ref <- dnorm(matrix(d$y1, nrow(m1), nrow(d), byrow = TRUE), m1,
                 dr$sigma_y1, log = TRUE)
    show("log_lik(resp = 'y1') vs univariate normal: max abs diff",
         max(abs(L1 - ref)))
  }
  if (w == "st") {
    m1 <- outer(dr$b_y1_Intercept, rep(1, nrow(d))) + outer(dr$b_y1_x, d$x)
    yy <- matrix(d$y1, nrow(m1), nrow(d), byrow = TRUE)
    ref <- dt((yy - m1) / dr$sigma_y1, dr$nu, log = TRUE) - log(dr$sigma_y1)
    show("student log_lik(resp = 'y1') vs univariate t(nu): max abs diff",
         max(abs(L1 - ref)))
  }
  if (w == "g3") {
    set.seed(2)
    p12 <- show("posterior_predict(resp = c('y1','y2')) shape",
                shape(posterior_predict(b, resp = c("y1", "y2"))))
    pp12 <- posterior_predict(b, resp = c("y1", "y2"))
    ep12 <- posterior_epred(b, resp = c("y1", "y2"))
    show("resp = c('y1','y2'): residual cor y1,y2", rcor(pp12, ep12, 1, 2))
    show("log_lik(resp = c('y1','y2')) dim",
         dim(log_lik(b, resp = c("y1", "y2"))))
  }
  set.seed(3)
  P <- posterior_predict(b)
  E <- posterior_epred(b)
  if (length(dim(P)) == 3L) {
    K <- dim(P)[3]
    for (j in 1:(K - 1)) for (k in (j + 1):K) {
      nm <- paste0("rescor__", dimnames(P)[[3]][j], "__",
                   dimnames(P)[[3]][k])
      rpost <- if (nm %in% names(dr)) mean(dr[[nm]]) else 0
      show(paste("residual cor of predictive draws", j, k,
                 "(mean over rows) | rescor posterior mean"),
           c(rcor(P, E, j, k), rpost))
    }
  }
  if (w == "gs") {
    # Sigma varies by row: the spread of y1's predictive residual per row
    # against that row's sigma posterior mean
    s1 <- posterior_epred(b, dpar = "sigma", resp = "y1")
    sdr <- apply(P[, , 1] - E[, , 1], 2, sd)
    show("gs: sd of y1 residual / mean sigma_y1, rows 1..5 and range",
         list(head(sdr / sqrt(colMeans(s1^2)), 5),
              range(sdr / sqrt(colMeans(s1^2)))))
  }
}
