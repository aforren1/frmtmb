## Review rescorpp, short check: resp on a univariate model, both arms.
##   Rscript dev/review-rescorpp-uniresp.R <lane|base>
source("C:/Users/adf44/source/r/frmtmb-wt-rescorpp/dev/review-rescorpp-common.R")
set.seed(5)
d1 <- data.frame(x = rnorm(80)); d1$y <- 1 + d1$x + rnorm(80)
ds1 <- quiet(frm_sample(frm(bf(y ~ x), data = d1), chains = 1, iter = 300,
                        refresh = 0, seed = 3))
tq <- function(e) tryCatch(paste(dim(suppressMessages(e)), collapse = "x"),
                           error = function(e) paste("ERROR:", substr(conditionMessage(e), 1, 90)))
for (rs in list("y", "foo", c("y", "y"), character(0))) {
  cat(sprintf("resp=%-5s pp %s | epred %s | log_lik %s | fitted %s\n", paste(rs, collapse = ","),
              tq(posterior_predict(ds1, resp = rs, ndraws = 5)),
              tq(posterior_epred(ds1, resp = rs, ndraws = 5)),
              tq(log_lik(ds1, resp = rs, ndraws = 5)),
              tq(fitted(ds1, resp = rs, ndraws = 5))))
}
