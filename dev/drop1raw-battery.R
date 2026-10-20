# Does the data a fit holds reproduce the fit? For each model: fit on
# `dd` in the global environment, remove `dd`, then update(fit) with no
# arguments must give the same log-likelihood (bitwise, same process),
# and nobs, and drop1() must run where the model has a droppable term.
# Arg: base | lane. Log: dev/drop1raw-out/battery-<arm>.txt (seed 11).
arm <- commandArgs(TRUE)[1]
.libPaths(c(if (arm == "lane") "C:/Users/adf44/source/r/wt-drop1raw-lib",
            "C:/Users/adf44/source/r/rellib-r9",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({
  library(frmtmb)
  library(frmtmb.latent)
})
cat("frmtmb", format(packageVersion("frmtmb")), "from",
    dirname(find.package("frmtmb")), "\nfrmtmb.latent from",
    dirname(find.package("frmtmb.latent")), "\n")
q <- function(e) suppressWarnings(suppressMessages(e))
mk <- function(n = 160) {
  set.seed(11)
  d <- data.frame(x = rnorm(n), x2 = rexp(n) + 0.1,
                  g = factor(sample(letters[1:6], n, TRUE)),
                  g2 = factor(sample(letters[1:6], n, TRUE)),
                  w = runif(n, 0.5, 2),
                  sb = rep(c(TRUE, TRUE, FALSE), length.out = n),
                  nt = 10L, expo = runif(n, 1, 3),
                  o = factor(sample(c("lo", "mid", "hi"), n, TRUE),
                             levels = c("lo", "mid", "hi"), ordered = TRUE))
  d$t <- ave(seq_len(n), d$g, FUN = seq_along)
  d$z <- 1 + d$x + log(d$x2) + rnorm(6)[d$g] + rnorm(n)
  d$zp <- 2 * exp(0.3 * d$x) + rnorm(n, 0, 0.3)
  d$cc <- ifelse(seq_len(n) %% 5 == 0, "interval", "none")
  d$z2 <- d$z + 1
  d$cnt <- rbinom(n, 10, plogis(0.2 + 0.5 * d$x))
  d$pc <- rpois(n, d$expo * exp(0.2 + 0.4 * d$x))
  d$sdx <- 0.2
  d$xm <- d$x + rnorm(n, 0, 0.2)
  d$cat <- factor(sample(c("a", "b", "c"), n, TRUE))
  d$zmi <- d$z
  d$zmi[c(5, 9, 33)] <- NA
  d
}
cases <- list(
  log = list(bf(z ~ x + log(x2))),
  smooth = list(bf(z ~ s(x2, k = 5))),
  re_slope = list(bf(z ~ x + (1 + x | g))),
  mo = list(bf(z ~ mo(o) + log(x2))),
  me = list(bf(z ~ me(xm, sdx))),
  mi = list(bf(zmi | mi() ~ x + log(x2))),
  cens = list(bf(z | cens(cc, z2) ~ x + log(x2))),
  trunc = list(bf(z | trunc(lb = -20) ~ x + log(x2))),
  trials = list(bf(cnt | trials(nt) ~ x + log(x2)), binomial()),
  rate = list(bf(pc | rate(expo) ~ x + log(x2)), poisson()),
  weights = list(bf(z | weights(w) ~ x + log(x2))),
  subset = list(bf(z | subset(sb) ~ x + log(x2))),
  offset = list(bf(z ~ x + offset(log(x2)))),
  sigma = list(bf(z ~ x + log(x2), sigma ~ log(x2))),
  scale_resp = list(bf(scale(z) ~ x + poly(x2, 2))),
  nl = list(bf(zp ~ a * exp(b * x), a ~ 1, b ~ 1, nl = TRUE)),
  categorical = list(bf(cat ~ log(x2)), categorical()),
  cumulative = list(bf(o ~ log(x2) + x), cumulative()),
  ar = list(bf(z ~ x + log(x2) + ar(time = t, gr = g))),
  mm = list(bf(z ~ x + (1 | mm(g, g2)))),
  mv = list(mvbf(bf(z ~ x + log(x2)), bf(cnt ~ log(x2)) + poisson()))
)
dd <- mk()
cat(sprintf("%-12s %-10s %-8s %s\n", "case", "update", "same_ll", "drop1"))
for (nm in names(cases)) {
  cs <- cases[[nm]]
  fam <- if (length(cs) > 1L) cs[[2L]] else gaussian()
  dd <- mk()
  fit <- tryCatch(q(if (nm == "mv") frm(cs[[1L]], data = dd) else
                      frm(cs[[1L]], family = fam, data = dd)),
                  error = function(e) e)
  if (inherits(fit, "error")) {
    cat(sprintf("%-12s FIT ERROR: %s\n", nm, substr(conditionMessage(fit),
                                                    1, 100)))
    next
  }
  rm(dd)
  u <- tryCatch(q(update(fit)), error = function(e) e)
  same <- if (inherits(u, "error")) NA else
    identical(as.numeric(logLik(u)), as.numeric(logLik(fit))) &&
      identical(nobs(u), nobs(fit))
  d1 <- if (nm %in% c("nl", "mv", "smooth", "me", "offset")) "n/a" else
    tryCatch(paste(rownames(q(drop1(fit)))[-1L], collapse = ","),
             error = function(e) paste("ERROR:", substr(conditionMessage(e),
                                                        1, 70)))
  cat(sprintf("%-12s %-10s %-8s %s\n", nm,
              if (inherits(u, "error")) "ERROR" else "ok", same, d1))
  if (inherits(u, "error")) {
    cat("   ", substr(conditionMessage(u), 1, 110), "\n")
  } else if (!isTRUE(same)) {
    cat(sprintf("    logLik %.12f vs %.12f, nobs %d vs %d\n",
                as.numeric(logLik(u)), as.numeric(logLik(fit)),
                nobs(u), nobs(fit)))
  }
}

cat("\n== a variable found from the formula environment ==\n")
dd <- mk()
xe <- dd$x2 * 2
fit <- q(frm(bf(z ~ x + log(xe)), data = dd))
rm(dd, xe)
u <- tryCatch(q(update(fit)), error = function(e) e)
cat("update after rm(dd, xe):", if (inherits(u, "error"))
  conditionMessage(u) else identical(as.numeric(logLik(u)),
                                     as.numeric(logLik(fit))), "\n")

cat("\n== hmm_starts() after the data are removed ==\n")
set.seed(4)
G <- matrix(c(0.9, 0.1, 0.2, 0.8), 2, 2, byrow = TRUE)
hd <- do.call(rbind, lapply(1:12, function(id) {
  s <- integer(20); s[1] <- 1L
  for (tt in 2:20) s[tt] <- sample.int(2, 1, prob = G[s[tt - 1], ])
  data.frame(id = id, t = 1:20, x2 = rexp(20) + 0.1,
             y = rnorm(20, c(0, 3)[s], 0.6))
}))
hd$y <- hd$y + 0.5 * log(hd$x2)
hfit <- q(frm(bf(y ~ log(x2)),
              family = hmm(K = 2, gaussian(), time = t, group = id),
              data = hd))
r1 <- tryCatch(q(hmm_starts(hfit, n = 2, seed = 1)), error = function(e) e)
cat("data present:", if (inherits(r1, "error")) substr(conditionMessage(r1),
    1, 120) else paste(r1$table$status, collapse = ","), "\n")
rm(hd)
r2 <- tryCatch(q(hmm_starts(hfit, n = 2, seed = 1)), error = function(e) e)
cat("data removed:", if (inherits(r2, "error")) substr(conditionMessage(r2),
    1, 120) else paste(r2$table$status, collapse = ","), "\n")
if (!inherits(r1, "error") && !inherits(r2, "error")) {
  cat("same table logLik:", identical(r1$table$logLik, r2$table$logLik), "\n")
}
