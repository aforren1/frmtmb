# Reviewer, item 2.3, claim 2: the 8.099-unit cold start, `diagnose()`
# printing "No convergence problems detected" on it, and `frm_allfit()`
# agreeing on the wrong answer across four optimizers.
#
# Built from the TEST FILE's own d4_data(), not from the lane's probe
# script, so that what is checked is the construction that ships in the
# suite. The two should be the same draw; whether they are is one of the
# things this script reports.
#
#   Rscript dev/rev-latent-repro81.R

source("C:/Users/adf44/source/r/frmtmb-wt-latent/dev/rev-latent-env.R")
suppressPackageStartupMessages({
  library(frmtmb)
  library(frmtmb.latent)
})
rev_env_report()

d4_data <- function() {
  set.seed(2026)
  K <- 2L; N <- 25L; Tg <- 30L
  G <- matrix(c(0.85, 0.15, 0.20, 0.80), 2, 2, byrow = TRUE)
  mu <- c(0, 3); sg <- c(0.6, 0.6); sd_b <- c(0.7, 0.5)
  stat <- local({
    A <- rbind(t(diag(K) - G), 1)
    drop(qr.solve(A, c(rep(0, K), 1)))
  })
  b <- cbind(stats::rnorm(N, 0, sd_b[1L]), stats::rnorm(N, 0, sd_b[2L]))
  d <- do.call(rbind, lapply(seq_len(N), function(g) {
    s <- integer(Tg)
    s[1L] <- sample.int(K, 1L, prob = stat)
    for (t in seq_len(Tg - 1L)) s[t + 1L] <- sample.int(K, 1L, prob = G[s[t], ])
    data.frame(ID = g, t = seq_len(Tg),
               y = stats::rnorm(Tg, mu[s] + b[g, s], sg[s]))
  }))
  d$gf <- factor(d$ID)
  d
}

dd <- d4_data()
cat("\nrows", nrow(dd), " sum(y)", format(sum(dd$y), digits = 12), "\n")

form <- bf(y ~ 1 + (1 | gf), mu2 ~ 1 + (1 | gf))
fam <- hmm(K = 2, gaussian(), time = t, group = ID, init = "stationary")

t0 <- Sys.time()
fit <- frm(form, family = fam, data = dd)
cat("cold fit seconds:", round(as.numeric(difftime(Sys.time(), t0,
                                                  units = "secs")), 2), "\n")
ll <- as.numeric(logLik(fit))
cat("cold logLik (16 sig)   :", format(ll, digits = 16), "\n")
cat("claim                  : -1096.09575602\n")
cat("relative gap to claim  :",
    format(abs(ll - (-1096.09575602)) / abs(ll), digits = 3), "\n")

cat("\n---- diagnose(), verbatim, NOT quiet ----\n")
dgo <- utils::capture.output(dg <- frmtmb::diagnose(fit))
cat(paste(dgo, collapse = "\n"), "\n")
cat("---- end diagnose() ----\n")
cat("prints 'No convergence problems detected':",
    any(grepl("No convergence problems detected", dgo)), "\n")
cat("convergence:", dg$convergence, " pdHess:", isTRUE(dg$pdHess),
    " nbadse:", length(dg$bad_se), " nflat:", length(dg$flat), "\n")
cat("max_grad (absolute):", format(dg$max_grad, digits = 6),
    "  claim 7e-4; diagnose()'s own clean test is max_grad < 1e-3\n")
cat("max_grad / |logLik| :", format(dg$max_grad / abs(ll), digits = 3),
    "  claim 6.4e-07\n")
cat("margin to diagnose()'s 1e-3 threshold, as a ratio:",
    format(dg$max_grad / 1e-3, digits = 4), "\n")

## the better optimum, restarted at probe D4's recorded point
tpl <- par_template(form, family = fam, data = dd)
st <- tpl
st$beta[["mu1_(Intercept)"]] <- -0.185185
st$beta[["mu2_(Intercept)"]] <- 3.121877
st$betad[["sigma1_(Intercept)"]] <- log(0.610455)
st$betad[["sigma2_(Intercept)"]] <- log(0.602054)
st$betad[["tr12_(Intercept)"]] <- log(0.173438 / 0.826562)
st$betad[["tr22_(Intercept)"]] <- log(0.828240 / 0.171760)
st$theta <- log(c(0.645297, 0.427590))
fit2 <- frm(form, family = fam, data = dd, start = st)
ll2 <- as.numeric(logLik(fit2))
cat("\nwarm logLik            :", format(ll2, digits = 16), "\n")
cat("claim                  : -1087.99646521\n")
cat("gap this run           :", format(ll2 - ll, digits = 10),
    "  claim 8.09929\n")
dg2 <- frmtmb::diagnose(fit2, quiet = TRUE)
cat("warm: convergence", dg2$convergence, " pdHess", isTRUE(dg2$pdHess),
    " max_grad", format(dg2$max_grad, digits = 3), "\n")

## ---- frm_allfit(): four optimizers, one start ----------------------
cat("\n---- frm_allfit() on the COLD fit ----\n")
t1 <- Sys.time()
af <- frmtmb::frm_allfit(fit)
cat("allfit seconds:", round(as.numeric(difftime(Sys.time(), t1,
                                                units = "secs")), 2), "\n")
print(af)
tb <- af$table %||% af[["table"]]
if (!is.null(tb)) {
  print(tb, row.names = FALSE, digits = 12)
  lls <- suppressWarnings(as.numeric(tb$logLik))
  ok <- is.finite(lls)
  cat("\noptimizers reporting a logLik:", sum(ok), "of", nrow(tb), "\n")
  cat("logLik spread over them      :",
      format(diff(range(lls[ok])), digits = 6), "  claim 7.84e-07\n")
  cat("worst |logLik - cold|        :",
      format(max(abs(lls[ok] - ll)), digits = 6), "\n")
  cat("any optimizer within 1 unit of the BETTER optimum:",
      any(abs(lls[ok] - ll2) < 1), "\n")
}
cat("\nnames on the allfit object:", paste(names(af), collapse = " "), "\n")
saveRDS(list(ll_cold = ll, ll_warm = ll2),
        "C:/Users/adf44/source/r/frmtmb-wt-latent/dev/rev-latent-repro81.rds")
