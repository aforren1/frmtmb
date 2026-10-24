# Release review of 0.74.0, a cross-lane case no lane ran: hmm()
# (frmtmb.latent, which lane spatial changed to read core's field-term
# list) with a car() field (frmtmb.spatial) in a transition parameter,
# under init = "stationary", which needs constant transitions. Also
# mu ~ car(): must fit, and agree with 0.73.0.
#   Rscript dev/review-rel074-hmmcar.R old|new
a <- commandArgs(trailingOnly = TRUE)
user <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
if (a[1] == "old") {
  .libPaths(c("C:/Users/adf44/source/r/rellib-r11", user))
  suppressMessages({library(frmtmb); library(frmtmb.latent)})
} else {
  .libPaths(c("C:/Users/adf44/source/r/rellib-r12", user))
  suppressMessages({library(frmtmb); library(frmtmb.latent)
    library(frmtmb.spatial)})
}
cat("frmtmb", format(packageVersion("frmtmb")), "latent",
    format(packageVersion("frmtmb.latent")), "\n")
g <- expand.grid(r = 1:4, c = 1:4)
n <- nrow(g)
W <- matrix(0, n, n)
for (i in 1:n) for (j in 1:n) {
  if (abs(g$r[i] - g$r[j]) + abs(g$c[i] - g$c[j]) == 1) W[i, j] <- 1
}
dimnames(W) <- list(paste0("L", 1:n), paste0("L", 1:n))
set.seed(4006)
nt <- 30
dd <- data.frame(loc = factor(rep(paste0("L", 1:n), each = nt),
                              levels = paste0("L", 1:n)),
                 t = rep(1:nt, n))
s <- integer(nrow(dd))
for (i in seq_len(nrow(dd))) {
  s[i] <- if (dd$t[i] == 1) sample(1:2, 1) else
    if (runif(1) < 0.85) s[i - 1] else 3L - s[i - 1]
}
u <- rnorm(n, 0, 0.5)
dd$y <- c(0, 3)[s] + u[as.integer(dd$loc)] + rnorm(nrow(dd), 0, 0.6)
run <- function(label, expr) {
  r <- tryCatch({
    f <- withCallingHandlers(expr, warning = function(w) {
      cat("   warning:", substr(conditionMessage(w), 1, 120), "\n")
      invokeRestart("muffleWarning")})
    sprintf("FIT logLik %.10f df %d", as.numeric(logLik(f)),
            attr(logLik(f), "df"))
  }, error = function(e) paste("ERROR", substr(gsub("\\s+", " ",
                                                    conditionMessage(e)), 1, 220)))
  cat(sprintf("%-36s %s\n", label, r))
}
run("mu ~ car, init estimated",
    frm(bf(y ~ 1 + car(W, gr = loc)),
        family = hmm(K = 2, gaussian(), time = t, group = loc),
        data = dd, data2 = list(W = W)))
run("tr12 ~ car, init stationary",
    frm(bf(y ~ 1, tr12 ~ 1 + car(W, gr = loc)),
        family = hmm(K = 2, gaussian(), time = t, group = loc,
                     init = "stationary"),
        data = dd, data2 = list(W = W)))
run("tr12 ~ (1 | loc), init stationary",
    frm(bf(y ~ 1, tr12 ~ 1 + (1 | loc)),
        family = hmm(K = 2, gaussian(), time = t, group = loc,
                     init = "stationary"),
        data = dd))
