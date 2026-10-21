# Fix round 2, N5: what brms 2.23.0 does with a binomial or
# beta_binomial row right-censored at (or above) its trials, an event
# of probability zero (P(Y > n) = 0), and with a count left-censored at
# 0 (P(Y <= 0) = P(Y = 0), possible). brm() with one short chain; the
# log density at a point through rstan::log_prob().
.libPaths(c("C:/Users/adf44/source/r/rellib-r10",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(brms))
set.seed(1)
n <- 40
d <- data.frame(y = rbinom(n, 5, 0.4), tr = 5, cc = 0)
try_brm <- function(lab, d, fam) {
  msgs <- character()
  r <- withCallingHandlers(
    tryCatch(brm(y | trials(tr) + cens(cc) ~ 1, data = d, family = fam,
                 chains = 1, iter = 40, warmup = 20, refresh = 0,
                 seed = 1, silent = 2),
             error = function(e) e),
    warning = function(w) {
      msgs <<- c(msgs, conditionMessage(w))
      invokeRestart("muffleWarning")
    },
    message = function(m) {
      msgs <<- c(msgs, conditionMessage(m))
      invokeRestart("muffleMessage")
    })
  if (inherits(r, "error")) {
    cat(sprintf("%-48s brm ERROR: %s\n", lab, substr(conditionMessage(r), 1, 200)))
  } else {
    up <- rstan::get_num_upars(r$fit)
    lp <- tryCatch(rstan::log_prob(r$fit, numeric(up), adjust_transform = FALSE),
                   error = function(e) paste("log_prob ERROR", conditionMessage(e)))
    cat(sprintf("%-48s brm ran; log_prob at 0: %s; draws %d\n", lab,
                format(lp), tryCatch(nrow(as.matrix(r)), error = function(e) 0L)))
  }
  if (length(msgs)) cat("   messages:", substr(paste(unique(msgs), collapse = " | "), 1, 300), "\n")
}
d1 <- d
d1$y[1] <- 5
d1$cc[1] <- 1
try_brm("binomial, right-censored at its trials", d1, binomial())
try_brm("beta_binomial, right-censored at its trials", d1, beta_binomial())
d2 <- d
d2$y[1] <- 0
d2$cc[1] <- -1
try_brm("binomial, left-censored at 0", d2, binomial())
cat("END\n")
