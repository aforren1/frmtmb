# Lane catmulti, fix round 2: ordinary fits keep 0.71.0's logLik() df,
# AIC and BIC exactly. The same fits on the lane build and on rellib-r9,
# with transformed terms; arg 1: "lane" or "base". The categorical
# reference spelling, refcat = "b", exists on the lane only and is
# compared with base's categorical() (the same likelihood).
# Output: dev/catmulti-out/df-sweep-<arm>.txt; compare with diff.
arm <- commandArgs(TRUE)[1]
.libPaths(c(if (arm == "lane") "C:/Users/adf44/source/r/wt-catmulti-lib",
            "C:/Users/adf44/source/r/rellib-r9",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
set.seed(31)
n <- 400
d <- data.frame(x = rnorm(n), x2 = exp(rnorm(n)),
                g = factor(sample(letters[1:10], n, TRUE)))
eta <- 0.3 + 0.5 * d$x - 0.4 * log(d$x2)
d$yg <- eta + rnorm(n)
d$yp <- rpois(n, exp(eta))
d$tr <- sample(5:15, n, TRUE)
d$yb <- rbinom(n, d$tr, plogis(eta))
d$yG <- rgamma(n, shape = 3, rate = 3 / exp(eta))
d$yo <- factor(cut(eta + rlogis(n), c(-Inf, -0.5, 0.5, Inf), labels = FALSE),
               ordered = TRUE)
E <- cbind(0, eta, -eta)
P <- exp(E) / rowSums(exp(E))
d$yc <- factor(apply(P, 1, function(p) sample(c("a", "b", "c"), 1, prob = p)))
Yc <- t(vapply(seq_len(n), function(i) as.numeric(rmultinom(1, 8, P[i, ])),
               numeric(3)))
d$Yc <- Yc
d$N <- 8L
G <- matrix(rgamma(3 * n, shape = 5 * P), n)
d$Yd <- G / rowSums(G)
rhs <- "poly(x, 2) + log(x2) + scale(x2) + I(x^2) + (1 | g)"
cases <- list(
  gaussian = list(paste("yg ~", rhs), gaussian()),
  poisson = list(paste("yp ~", rhs), poisson()),
  binomial = list(paste("yb | trials(tr) ~", rhs), binomial()),
  Gamma = list(paste("yG ~", rhs), Gamma(link = "log")),
  negbinomial = list(paste("yp ~", rhs), negbinomial()),
  cumulative = list(paste("yo ~", rhs), cumulative()),
  categorical = list(paste("yc ~", rhs), categorical()),
  multinomial = list(paste("Yc | trials(N) ~", rhs), multinomial(K = 3)),
  dirichlet = list(paste("Yd ~", "poly(x, 2) + log(x2)"), dirichlet()))
if (arm == "lane") {
  cases$categorical_refb <- list(paste("yc ~", rhs),
                                 categorical(refcat = "b"))
}
for (nm in names(cases)) {
  cs <- cases[[nm]]
  f <- suppressWarnings(suppressMessages(
    frm(bf(stats::as.formula(cs[[1]])), family = cs[[2]], data = d)))
  ll <- logLik(f)
  cat(sprintf("%-18s df %d  AIC %.10f  BIC %.10f\n",
              if (nm == "categorical_refb") "categorical" else nm,
              attr(ll, "df"), AIC(f), BIC(f)))
}
