# Fix round 2: which of the shifted-mean and bimodal fits warns, and
# what each one's logLik is against (1 | g).
.libPaths(c("C:/Users/adf44/source/r/wt-spatial-lib",
            "C:/Users/adf44/source/r/rellib-r11",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
set.seed(1)
d <- data.frame(g = factor(rep(1:20, each = 5)), x = rnorm(100))
d$y <- 1 + d$x + rnorm(20)[d$g] + rnorm(100, 0, 0.5)
vc <- function(theta, blk) {
  matrix(exp(2 * theta[1]), 1, 1, dimnames = list(blk$cnms, blk$cnms))
}
for (m in c(600, 50, 5)) {
  nm <- paste0("shift", m)
  frmtmb_register_covstruct(nm, lpdf = local({mm <- m
    function(b, theta, blk) sum(RTMB::dnorm(b, mm, exp(theta[1]), log = TRUE))
  }), vcov = vc, npar = 1L, sd_idx = 1L)
  w <- character()
  f <- withCallingHandlers(frm(stats::as.formula(paste0("y ~ x + ", nm,
                                                        "(1 | g)")),
                               data = d),
                           warning = function(e) {
                             w <<- c(w, substr(conditionMessage(e), 1, 60))
                             invokeRestart("muffleWarning")
                           })
  cat(nm, format(as.numeric(logLik(f)), digits = 12), length(w), w, "\n")
}
cat("(1|g)", format(as.numeric(logLik(frm(y ~ x + (1 | g), data = d))),
                    digits = 12), "\n")
frmtmb_register_covstruct("twomodes", lpdf = function(b, theta, blk) {
  sd <- exp(theta[1])
  sum(log(0.5 * exp(RTMB::dnorm(b, -1, sd, log = TRUE)) +
            0.5 * exp(RTMB::dnorm(b, 1, sd, log = TRUE))))
}, vcov = function(theta, blk) {
  matrix(1 + exp(2 * theta[1]), 1, 1, dimnames = list(blk$cnms, blk$cnms))
}, npar = 1L, sd_idx = 1L)
w <- character()
f <- withCallingHandlers(frm(y ~ x + twomodes(1 | g), data = d),
                         warning = function(e) {
                           w <<- c(w, substr(conditionMessage(e), 1, 60))
                           invokeRestart("muffleWarning")
                         })
cat("twomodes", format(as.numeric(logLik(f)), digits = 12), length(w), w,
    "\n")
# data drawn from the bimodal structure itself
set.seed(3)
G <- 40
u <- sample(c(-1, 1), G, TRUE) + rnorm(G, 0, 0.3)
d2 <- data.frame(g = factor(rep(seq_len(G), each = 8)), x = rnorm(8 * G))
d2$y <- 1 + d2$x + u[d2$g] + rnorm(8 * G, 0, 0.3)
w <- character()
f <- withCallingHandlers(frm(y ~ x + twomodes(1 | g), data = d2),
                         warning = function(e) {
                           w <<- c(w, substr(conditionMessage(e), 1, 60))
                           invokeRestart("muffleWarning")
                         })
cat("twomodes on its own data", format(as.numeric(logLik(f)), digits = 12),
    "sd", exp(f$estimates$theta), length(w), w, "\n")
