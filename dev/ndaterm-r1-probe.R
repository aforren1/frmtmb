# Lane ndaterm, fix round 1: class-wise simulators (mixture, grouped
# mixture, hmm) with per-row trials(), and zero_inflated_binomial at the
# zi boundary (review item 8). The constructions are the reviewer's
# (dev/ndaterm-rev-probe2.R, probe3.R, probe6.R), seed 7.
# Usage: Rscript dev/ndaterm-r1-probe.R <arm>, arm mine or r9.
# Log: dev/ndaterm-out/r1-probe-<arm>.txt
args <- commandArgs(TRUE)
arm <- if (length(args)) args[1] else "mine"
base <- "C:/Users/adf44/source/r"
.libPaths(c(if (arm == "mine") file.path(base, "wt-ndaterm-lib"),
            file.path(base, "rellib-r9"),
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(frmtmb); library(frmtmb.latent)})
cat("arm", arm, find.package("frmtmb"), find.package("frmtmb.latent"),
    "\n")
show <- function(label, x) {
  cat(sprintf("%-50s %s\n", label,
              paste(format(x, digits = 5), collapse = "  ")))
}

cat("\n== mixture(binomial, binomial), N in {2, 50}\n")
set.seed(7)
n <- 200
d <- data.frame(x = rnorm(n), N = sample(c(2, 50), n, TRUE))
cls <- rbinom(n, 1, 0.5)
d$y <- rbinom(n, d$N, ifelse(cls == 1, plogis(-1 + 0.5 * d$x),
                             plogis(1 + 0.5 * d$x)))
fmx <- frm(bf(y | trials(N) ~ 1, mu1 ~ x, mu2 ~ x),
           family = mixture(binomial(), binomial()), data = d)
s <- as.matrix(simulate(fmx, nsim = 200, seed = 1))
show("simulate: share of draws > N", mean(s > d$N))
show("simulate: mean draw on N = 2 rows", mean(s[d$N == 2, ]))
show("fitted: mean on N = 2 rows", mean(fitted(fmx)[d$N == 2, 1]))
nd <- data.frame(x = c(0, 0), N = c(1000, 2))
set.seed(3)
p <- predict(fmx, newdata = nd, summary = FALSE, ndraws = 400)
show("predict(newdata) col means", colMeans(p))
show("predict(newdata) col max", apply(p, 2, max))
show("fitted(newdata)", fitted(fmx, newdata = nd)[, "Estimate"])

cat("\n== mixture(binomial, binomial, groups = ~id)\n")
set.seed(7)
ng <- 40; m <- 6
dg <- data.frame(id = rep(seq_len(ng), each = m), x = rnorm(ng * m),
                 N = sample(c(2, 50), ng * m, TRUE))
cl <- rep(rbinom(ng, 1, 0.5), each = m)
dg$y <- rbinom(nrow(dg), dg$N, ifelse(cl == 1, 0.2, 0.8))
fg <- frm(bf(y | trials(N) ~ 1),
          family = mixture(binomial(), binomial(), groups = ~id), data = dg)
s <- as.matrix(simulate(fg, nsim = 100, seed = 1))
show("simulate: share of draws > N", mean(s > dg$N))
show("simulate: mean draw on N = 2 rows", mean(s[dg$N == 2, ]))
show("fitted: mean on N = 2 rows", mean(fitted(fg)[dg$N == 2, 1]))

cat("\n== hmm(K = 2, binomial())\n")
set.seed(7)
d2 <- data.frame(id = rep(1:20, each = 30), t = rep(1:30, 20),
                 N = sample(c(2, 50), 600, TRUE))
st <- unlist(lapply(1:20, function(i) {
  s <- numeric(30); s[1] <- 1
  for (j in 2:30) s[j] <- if (runif(1) < 0.9) s[j - 1] else 3 - s[j - 1]
  s
}))
d2$y <- rbinom(600, d2$N, c(0.2, 0.8)[st])
fh <- frm(bf(y | trials(N) ~ 1),
          family = hmm(K = 2, binomial(), time = t, group = id), data = d2)
s <- as.matrix(simulate(fh, nsim = 50, seed = 1))
show("simulate: share of draws > N", mean(s > d2$N))
show("simulate: mean draw on N = 2 rows", mean(s[d2$N == 2, ]))
show("fitted: mean on N = 2 rows", mean(fitted(fh)[d2$N == 2, 1]))

cat("\n== zero_inflated_binomial, no inflation in the data\n")
set.seed(7)
dz <- data.frame(x = rnorm(n), N = sample(5:15, n, TRUE))
dz$y <- rbinom(n, dz$N, plogis(0.3 + 0.5 * dz$x))
fz <- frm(bf(y | trials(N) ~ x), family = zero_inflated_binomial(),
          data = dz)
show("zi estimate (link), its se",
     c(fz$opt$par[grep("zi", names(fz$opt$par))][1],
       sqrt(diag(vcov(fz, full = TRUE)))[grep("zi",
         rownames(vcov(fz, full = TRUE)))][1]))
print(fixef(fz))
set.seed(1)
p <- predict(fz, ndraws = 400)
show("predict() mean of Estimate / fitted() mean",
     c(mean(p[, "Estimate"]), mean(fitted(fz)[, "Estimate"])))
set.seed(1)
p0 <- predict(fz, ndraws = 400, propagate_error = FALSE)
show("predict(propagate_error = FALSE) mean", mean(p0[, "Estimate"]))
s <- as.matrix(simulate(fz, nsim = 400, seed = 1))
show("simulate() mean", mean(s))
