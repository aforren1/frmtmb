# Reviewer of lane ordmix: post-fit outputs of an ordinal mixture,
# checked by hand. Data seed 20261099, n = 400.
.libPaths(c("C:/Users/adf44/source/r/wt-ordmix-lib",
            "C:/Users/adf44/source/r/rellib-r5",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({
  library(frmtmb)
})
cat("frmtmb", find.package("frmtmb"), "\n")
set.seed(20261099)
n <- 400
d <- data.frame(x = rnorm(n), z = rnorm(n),
                g = factor(sample(c("a", "b"), n, TRUE)))
cls <- rbinom(n, 1, plogis(-0.3 + 0.8 * d$z))
lat <- ifelse(cls == 1, 1.5 * d$x + 1, -0.8 * d$x - 1) + rlogis(n)
d$y <- 1L + (lat > -1.5) + (lat > 0) + (lat > 1.5)

## 1. fitted() against the theta-weighted sum computed by hand
f1 <- frm(bf(y ~ x, theta1 ~ z), family = mixture(cumulative("probit"),
                                                   sratio()), data = d)
fx <- fixef(f1)[, "Estimate"]
print(round(fx, 4))
t1 <- fx[paste0("mu1_Intercept[", 1:3, "]")]
t2 <- fx[paste0("mu2_Intercept[", 1:3, "]")]
e1 <- fx[["mu1_x"]] * d$x
e2 <- fx[["mu2_x"]] * d$x
th <- plogis(fx[["theta1_Intercept"]] + fx[["theta1_z"]] * d$z)
p1 <- t(sapply(seq_len(n), function(i) diff(c(0, pnorm(t1 - e1[i]), 1))))
h2 <- function(i) plogis(t2 - e2[i])
p2 <- t(sapply(seq_len(n), function(i) {
  h <- h2(i)
  c(h[1], (1 - h[1]) * h[2], (1 - h[1]) * (1 - h[2]) * h[3],
    prod(1 - h))
}))
Ph <- th * p1 + (1 - th) * p2
Pf <- fitted(f1)[, "Estimate", ]
cat(sprintf("1 fitted vs hand: max|diff| %.3g; max|rowSums(fitted) - 1| %.3g; max|rowSums(hand) - 1| %.3g\n",
            max(abs(unname(Pf) - Ph)), max(abs(rowSums(Pf) - 1)),
            max(abs(rowSums(Ph) - 1))))
# the density before the normalization ord_mix_probs() applies
fam <- f1$spec$responses[[1]]$family
cat("  dimnames(fitted):", paste(dimnames(fitted(f1))[[3]], collapse = " "),
    "\n")

## 2. simulate(): category counts against the fitted probabilities
S <- 500
sm <- simulate(f1, nsim = S, seed = 1)
Y <- as.matrix(sm)
cat("  simulate() dims", dim(Y), "\n")
obs <- sapply(1:4, function(k) sum(Y == k))
ex <- colSums(Ph) * S
sdv <- sqrt(colSums(Ph * (1 - Ph)) * S)
cat(sprintf("2 simulate nsim=%d: category z = %s\n", S,
            paste(sprintf("%.2f", (obs - ex) / sdv), collapse = " ")))
# per decile of P(Y = 1): each decile's share of category 1
dec <- cut(Ph[, 1], quantile(Ph[, 1], 0:10 / 10), include.lowest = TRUE)
zd <- sapply(levels(dec), function(l) {
  r <- dec == l
  (sum(Y[r, ] == 1) - S * sum(Ph[r, 1])) /
    sqrt(S * sum(Ph[r, 1] * (1 - Ph[r, 1])))
})
cat(sprintf("  per-decile z for category 1: max |z| %.2f over 10\n",
            max(abs(zd))))

## 3. predict()
pr <- predict(f1)
cat("3 predict() class", class(pr), "dims", paste(dim(pr), collapse = "x"),
    "\n")
print(head(pr, 3))

## 4. newdata with an unseen threshold group
f2 <- frm(bf(y | thres(gr = g) ~ x), family = mixture(cumulative(),
                                                       sratio()), data = d)
nd <- data.frame(x = c(0, 1), g = factor(c("a", "c")))
r4 <- tryCatch(fitted(f2, newdata = nd)[, "Estimate", ],
               error = function(e) paste("ERROR:", conditionMessage(e)))
cat("4 mixture, unseen group 'c':\n")
print(r4)
fp <- frm(bf(y | thres(gr = g) ~ x), family = cumulative(), data = d)
r4p <- tryCatch(fitted(fp, newdata = nd)[, "Estimate", ],
                error = function(e) paste("ERROR:", conditionMessage(e)))
cat("  plain cumulative, unseen group 'c':\n")
print(r4p)
r4s <- tryCatch(simulate(f2, newdata = nd, nsim = 2, seed = 1),
                error = function(e) paste("ERROR:", conditionMessage(e)))
cat("  mixture simulate(newdata) unseen group:\n")
print(r4s)

## 5. conditional_effects(): default and categorical = TRUE
ce0 <- withCallingHandlers(conditional_effects(f1, effects = "x"),
                           warning = function(w) {
                             cat("  CE default WARNING:",
                                 conditionMessage(w), "\n")
                             invokeRestart("muffleWarning")
                           })
c0 <- ce0[[1]]
cat("5 CE default columns:", paste(names(c0), collapse = " "), "\n")
cat("  CE default has cats__:", "cats__" %in% names(c0), " rows", nrow(c0),
    "\n")
ce1 <- conditional_effects(f1, effects = "x", categorical = TRUE)[[1]]
cat("  CE categorical rows", nrow(ce1), "cats", levels(ce1$cats__), "\n")
# brms's categorical display at a grid point: posterior_epred at the
# grid with z at its mean. Here: fitted() at the same newdata
gx <- unique(ce1$x)
ndc <- data.frame(x = gx, z = mean(d$z))
Pg <- fitted(f1, newdata = ndc)[, "Estimate", ]
mx <- max(sapply(1:4, function(k) {
  max(abs(ce1$estimate__[ce1$cats__ == as.character(k)] - Pg[, k]))
}))
cat(sprintf("  CE categorical vs fitted(newdata) at the grid: max|diff| %.3g\n",
            mx))
# plain cumulative's default CE, for the convention
cp <- withCallingHandlers(conditional_effects(fp, effects = "x")[[1]],
                          warning = function(w) {
                            cat("  plain CE default WARNING:",
                                conditionMessage(w), "\n")
                            invokeRestart("muffleWarning")
                          })
cat("  plain cumulative CE default has cats__:", "cats__" %in% names(cp),
    "\n")
if (!"cats__" %in% names(c0)) {
  # brms's default: the expected category sum_k k p_k, with a warning
  Pd <- fitted(f1, newdata = data.frame(x = unique(c0$x),
                                        z = mean(d$z)))[, "Estimate", ]
  cat(sprintf("  CE default vs sum_k k p_k: max|diff| %.3g\n",
              max(abs(c0$estimate__ - as.vector(Pd %*% (1:4))))))
}

## 6. residuals and OSA refused by name
for (ty in c("response", "pearson", "osa")) {
  r <- tryCatch({
    residuals(f1, type = ty)
    "returned"
  }, error = function(e) paste("ERROR:", substr(conditionMessage(e), 1, 150)))
  cat("6 residuals type =", ty, ":", r, "\n")
}
