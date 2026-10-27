# Reviewer of lane loglik: what brms 2.23.0 does with addition terms
# beside a residual covariance term (stancode only, no compile), and
# whether its log_lik() reads weights there. Seed 92.
# Output: dev/review-loglik-log/brmsadd.txt
.libPaths("C:/Users/adf44/AppData/Local/R/win-library/4.6")
suppressPackageStartupMessages(library(brms))
set.seed(92)
d <- data.frame(g = rep(1:5, each = 4), t = rep(1:4, 5))
d$y <- rnorm(20)
d$w <- runif(20, 0.5, 2)
d$se <- runif(20, 0.1, 0.3)
d$cc <- 0L
d$lb <- -10
d$ym <- d$y
d$ym[3] <- NA
fs <- list(y | weights(w) ~ ar(t, g, cov = TRUE),
           y | se(se, sigma = TRUE) ~ ar(t, g, cov = TRUE),
           y | cens(cc) ~ ar(t, g, cov = TRUE),
           y | trunc(lb = lb) ~ ar(t, g, cov = TRUE),
           ym | mi() ~ ar(t, g, cov = TRUE),
           y | weights(w) ~ cosy(t, g),
           y | se(se, sigma = TRUE) ~ unstr(t, g),
           y | weights(w) ~ fcor(M))
M <- diag(20)
for (f in fs) {
  r <- tryCatch({
    sc <- make_stancode(f, data = d, data2 = list(M = M))
    lk <- grep("target \\+=|_lpdf\\(Y", strsplit(sc, "\n")[[1]],
               value = TRUE)
    paste("ACCEPTED:", paste(trimws(lk), collapse = " | "))
  }, error = function(e) paste("REFUSED:", conditionMessage(e)))
  cat(sprintf("%-46s %s\n", deparse1(f), substr(r, 1, 220)))
}
cat("DONE\n")
