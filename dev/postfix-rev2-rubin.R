# Reviewer of lane postfix, round 2: fitted() of frm_multiple() on the
# link scale, against Rubin's rules by hand over each imputation's
# fitted(scale = "linear"), and its bounds inside the support.
.libPaths(c("C:/Users/adf44/source/r/wt-postfix-lib",
            "C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(frmtmb); library(mice)})
cat("lib:", find.package("frmtmb"), "\n")
rel <- function(a, b) max(abs(a - b) / pmax(abs(b), 1e-300))
hand <- function(fm, newdata, linkinv) {
  m <- fm$m
  per <- lapply(fm$fits, fitted, newdata = newdata, scale = "linear")
  Q <- sapply(per, function(p) p[, "Estimate"])
  U <- sapply(per, function(p) p[, "Est.Error"]^2)
  B <- apply(Q, 1, var); Tt <- rowMeans(U) + (1 + 1/m) * B
  lam <- (1 + 1/m) * B / Tt
  dfcom <- df.residual(fm$fits[[1]])
  dfo <- (m - 1) / lam^2
  dfb <- (dfcom + 1) / (dfcom + 3) * dfcom * (1 - lam)
  df <- dfo * dfb / (dfo + dfb)
  q <- rowMeans(Q)
  list(lin = cbind(q, sqrt(Tt), q + qt(0.025, df) * sqrt(Tt),
                   q + qt(0.975, df) * sqrt(Tt)),
       resp = cbind(linkinv(q), NA, linkinv(q + qt(0.025, df) * sqrt(Tt)),
                    linkinv(q + qt(0.975, df) * sqrt(Tt))))
}
## poisson (1 | g), as in dev/postfix-rev-rubin.R section 1
set.seed(41)
n <- 200; G <- 15
g <- factor(sample(seq_len(G), n, TRUE))
u <- rnorm(G, 0, 0.5)
x <- rnorm(n); z <- rnorm(n)
y <- rpois(n, exp(0.3 + 0.4 * x - 0.3 * z + u[g]))
x[sample(n, 50)] <- NA
imp <- mice(data.frame(y, x, z, g), m = 5, seed = 42, print = FALSE,
            method = c("", "pmm", "", ""))
fm <- frm_multiple(y ~ x + z + (1 | g), family = poisson(), data = imp)
for (nd in list(NULL, data.frame(x = c(-3, 0, 3), z = c(3, 0, -3),
                                 g = factor(c(1, 2, 3))))) {
  lbl <- if (is.null(nd)) "in sample" else "newdata"
  h <- hand(fm, nd, exp)
  fr <- fitted(fm, newdata = nd)
  fl <- fitted(fm, newdata = nd, scale = "linear")
  cat("poisson", lbl, ": linear vs hand", rel(fl[, c(1, 2, 3, 4)], h$lin),
      " response Estimate", rel(fr[, 1], h$resp[, 1]), " Q2.5",
      rel(fr[, 3], h$resp[, 3]), " Q97.5", rel(fr[, 4], h$resp[, 4]), "\n")
  cat("   min Q2.5", min(fr[, 3]), " Est.Error / (mu * se_eta) range",
      range(fr[, 2] / (exp(h$lin[, 1]) * h$lin[, 2])), "\n")
}
## bernoulli, as in dev/postfix-rev-rubin.R section 2
set.seed(43)
n <- 40
x <- rnorm(n); y <- rbinom(n, 1, plogis(2 + 2 * x))
x[sample(n, 10)] <- NA
imp2 <- mice(data.frame(y, x), m = 5, seed = 44, print = FALSE)
fb <- frm_multiple(y ~ x, family = bernoulli(), data = imp2)
nd <- data.frame(x = c(-4, -2, 0, 1.5, 2.5, 5))
h <- hand(fb, nd, plogis)
fr <- fitted(fb, newdata = nd)
print(fr)
cat("bernoulli newdata: Estimate", rel(fr[, 1], h$resp[, 1]), " Q2.5",
    rel(fr[, 3], h$resp[, 3]), " Q97.5", rel(fr[, 4], h$resp[, 4]), "\n")
fi <- fitted(fb)
cat("bernoulli in sample bounds range", range(fi[, 3:4]), "\n")
ce <- conditional_effects(fb, "x")$x
fc <- fitted(fb, newdata = data.frame(x = ce$x))
cat("ce vs fitted on the ce grid: est", rel(fc[, 1], ce$estimate__),
    " lower", rel(fc[, 3], ce$lower__), " upper", rel(fc[, 4], ce$upper__), "\n")
## cumulative: category probabilities stay in [0, 1]
set.seed(1)
dd <- data.frame(x = rnorm(60))
dd$y <- 1 + dd$x + rnorm(60)
imo <- lapply(1:3, function(i) {
  d <- dd; d$o <- cut(d$y + rnorm(60, 0, .3), c(-Inf, 0, 1, 2, Inf),
                      ordered_result = TRUE); d })
fmo <- frm_multiple(bf(o ~ x), family = cumulative(), data = imo)
fo <- fitted(fmo, newdata = data.frame(x = c(-3, 0, 3)))
cat("cumulative bounds range", range(fo[, 3:4, ]), " probs sum",
    range(apply(fo[, 1, ], 1, sum)), "\n")
## gaussian identity link: response = linear
gi <- lapply(1:3, function(i) { d <- dd; d$y <- d$y + rnorm(60, 0, .1); d })
fg <- frm_multiple(bf(y ~ x), family = gaussian(), data = gi)
cat("gaussian response vs linear", rel(fitted(fg)[, 1:4],
                                       fitted(fg, scale = "linear")[, 1:4]), "\n")
