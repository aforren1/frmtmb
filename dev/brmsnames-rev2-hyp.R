## Reviewer recheck, BLOCKER 1 on frm_multiple() and on spellings the pin
## does not reach: dots and digits in names, functions, powers, decimals,
## named hypothesis vectors, an sd_ of an interaction group, an r_ name
## whose level holds brms's renamed characters. Hand values come from
## fixef(), varcorr_matrices() and the draws matrix, never from a parser.
##   Rscript dev/brmsnames-rev2-hyp.R
## Data seed 52 (dev/brmsnames-rev2-data.R), second imputation seed 91.
.libPaths(c("C:/Users/adf44/source/r/brmsnames-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
q <- function(e) suppressWarnings(suppressMessages(e))
try1 <- function(e) {
  tryCatch(e, error = function(err) paste("ERROR:", conditionMessage(err)))
}
q(library(frmtmb)); q(library(frmtmb.sample))
source("dev/brmsnames-rev2-data.R")
d <- rev2_data()
d$x.1 <- d$x + rnorm(nrow(d), 0, 0.5)
d2 <- d
set.seed(91)
d2$y <- d2$y + rnorm(nrow(d2), 0, 0.3)
f <- bf(y ~ x.1 + x * f + (1 | g:h))
fm <- q(frm_multiple(f + gaussian(), data = list(d, d2)))
fit1 <- q(frm(f, family = gaussian(), data = d))
fit2 <- q(frm(f, family = gaussian(), data = d2))
fe <- function(fit) {
  v <- unlist(fixef(fit))
  names(v) <- sub("^mu[.]", "", names(v))
  v
}
pool <- (fe(fit1) + fe(fit2)) / 2
cat("variables(fm):", try1(variables(fm)), "\n")
hs <- c(
  a = "x:fe:f > 0",
  b = "x.1 - 0.5 = 0",
  c = "exp(x:fcMd) > 1",
  d = "fe:f^2 = 0",
  e = "2 * x + .5 * x:fcMd = 0"
)
hand <- c(pool[["x:fe:f"]], pool[["x.1"]] - 0.5, exp(pool[["x:fc-d"]]) - 1,
          pool[["fe:f"]]^2, 2 * pool[["x"]] + 0.5 * pool[["x:fc-d"]])
r <- try1(q(hypothesis(fm, hs))$hypothesis)
if (is.data.frame(r)) {
  cat("frm_multiple, pooled point estimate against the mean of the two",
      "fits:\n")
  print(data.frame(h = hs, label = r$Hypothesis, est = r$Estimate,
                   hand = hand, diff = r$Estimate - hand))
} else cat("frm_multiple:", r, "\n")

cat("\n== the fit, sd_ of an interaction group ==\n")
sdg <- sqrt(varcorr_matrices(fit1)[[1]][1, 1])
r <- try1(q(hypothesis(fit1, "sd_g:h__Intercept = 0", class = NULL))$
            hypothesis$Estimate)
cat("hand", sdg, " hypothesis", r, "\n")
r <- try1(q(hypothesis(fit1, "Intercept = 0", class = "sd", group = "g:h"))$
            hypothesis$Estimate)
cat("class = sd, group = g:h:", r, "\n")

cat("\n== draws, r_ names with renamed levels ==\n")
set.seed(3)
ds <- q(frm_sample(fit1, chains = 1, iter = 100, refresh = 0, seed = 3))
cn <- grep("^r_", colnames(ds$draws), value = TRUE)
cat("first r_ names:", head(cn, 3), "\n")
h <- paste(cn[1], "-", cn[2], "= 0")
r <- try1(q(hypothesis(ds, h, class = NULL))$hypothesis$Estimate)
cat(h, ": hand", mean(ds$draws[, cn[1]] - ds$draws[, cn[2]]),
    " hypothesis", r, "\n")
r <- try1(q(hypothesis(ds, "Intercept > 0", scope = "ranef",
                       group = "g:h"))$hypothesis)
if (is.data.frame(r)) {
  cat("scope ranef group g:h: rows", nrow(r), " first group",
      as.character(r$Group[1]), " est", r$Estimate[1], " hand",
      mean(ds$draws[, cn[1]]), "\n")
} else cat("scope ranef:", r, "\n")
