# Reviewer: equations between hmm() transition cells (frmtmb.latent from
# rellib-r3 over the lane's frmtmb). With tr12 = "tr22" both rows of the
# transition matrix are the same, so under init = "stationary" the chain
# is i.i.d. and the model is the two-component mixture.
.libPaths(c("C:/Users/adf44/source/r/wt-formula2-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(frmtmb); library(frmtmb.latent)})
cat("latent from", find.package("frmtmb.latent"), "frmtmb from",
    find.package("frmtmb"), "\n")
q <- function(expr) {
  tryCatch(suppressWarnings(suppressMessages(expr)),
           error = function(e) structure(conditionMessage(e), class = "ERR"))
}
set.seed(4001)
G <- 20; Tn <- 25
z <- unlist(lapply(seq_len(G), function(g) {
  s <- integer(Tn); s[1] <- sample(1:2, 1)
  for (t in 2:Tn) s[t] <- if (runif(1) < 0.8) s[t - 1] else 3 - s[t - 1]
  s
}))
dd <- data.frame(id = rep(seq_len(G), each = Tn), t = rep(seq_len(Tn), G))
dd$y <- rnorm(nrow(dd), c(0, 3)[z], 0.7)
h1 <- q(frm(bf(y ~ 1, tr12 = "tr22"),
            family = hmm(K = 2, gaussian(), time = t, group = id,
                         init = "stationary"), data = dd))
mx <- q(frm(bf(y ~ 1), family = mixture(gaussian(), gaussian()), data = dd))
if (inherits(h1, "ERR")) cat("hmm tr12 = tr22:", h1, "\n") else {
  cat("hmm tr12 = tr22 logLik", as.numeric(logLik(h1)), "df",
      attr(logLik(h1), "df"), "\n")
  cat("mixture          logLik", as.numeric(logLik(mx)), "df",
      attr(logLik(mx), "df"), "\n")
  cat("rel diff", (as.numeric(logLik(h1)) - as.numeric(logLik(mx))) /
        abs(as.numeric(logLik(mx))), "\n")
  cat("variables(h1):", format(q(variables(h1))), "
")
  cat("summary(h1):", class(q(summary(h1))), format(q(summary(h1)))[1], "
")
  cat("fixef(h1):", format(q(rownames(fixef(h1)))), "
")
  cat("hypothesis(h1):", class(q(hypothesis(h1, "tr12_Intercept = 0"))), "
")
}
h0 <- q(frm(bf(y ~ 1), family = hmm(K = 2, gaussian(), time = t,
                                    group = id), data = dd))
h2 <- q(frm(bf(y ~ 1, sigma1 = "sigma2"),
            family = hmm(K = 2, gaussian(), time = t, group = id),
            data = dd))
if (inherits(h2, "ERR")) cat("hmm sigma1 = sigma2:", h2, "\n") else {
  cat("hmm free df", attr(logLik(h0), "df"), "sigma-equated df",
      attr(logLik(h2), "df"), "logLik", as.numeric(logLik(h0)),
      as.numeric(logLik(h2)), "\n")
  cat("variables(h2):", format(q(variables(h2))), "
")
}
cat("DONE\n")
