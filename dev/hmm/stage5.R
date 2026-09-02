# Stage 5: post-processing - forward-backward, Viterbi, fitted,
# residuals, simulate, and the refusals around them.
suppressMessages(pkgload::load_all("C:/Users/adf44/source/r/frmtmb-wt-hmm",
                                   quiet = TRUE))
source("C:/Users/adf44/source/r/frmtmb-wt-hmm/dev/hmm/hmm-common.R")
D <- "C:/Users/adf44/source/r/frmtmb-wt-hmm/dev/hmm/"
b1 <- readRDS(paste0(D, "probeB1.rds"))
dd <- b1$dat

fb <- frm(bf(y ~ 1),
          family = hmm(K = 2, gaussian(), time = t, group = g,
                       init = "estimated", trans = ~x), data = dd)

P <- hmm_probs(fb)
cat("dim:", dim(P), " rowSums range:", range(rowSums(P)), "\n")
cat("local decoding accuracy :", mean(max.col(P) == dd$state), "\n")
v <- hmm_viterbi(fb)
cat("global decoding accuracy:", mean(v == dd$state), "\n")
cat("probe C reference: 98.7% local\n")

## --- forward-backward against a brute-force per-group reference ------
e <- unlist(fixef(fb))
mu <- c(e[["mu1.(Intercept)"]], e[["mu2.(Intercept)"]])
sg <- exp(c(e[["sigma1.(Intercept)"]], e[["sigma2.(Intercept)"]]))
B <- matrix(0, 2, 2); Bx <- matrix(0, 2, 2)
B[1, 2] <- e[["tr12.(Intercept)"]]; Bx[1, 2] <- e[["tr12.x"]]
B[2, 2] <- e[["tr22.(Intercept)"]]; Bx[2, 2] <- e[["tr22.x"]]
dl <- softmax0(fb$estimates[["hmm_ldel"]])
lpm <- lpmat_gauss(dd$y, mu, sg)
rows <- hmm_seq_index(dd$g, dd$t)
Pref <- matrix(NA_real_, nrow(dd), 2)
Vref <- integer(nrow(dd))
for (r in rows) {
  Tl <- length(r)
  A <- matrix(0, Tl, 2); Bm <- matrix(0, Tl, 2)
  A[1, ] <- log(dl) + lpm[r[1], ]
  for (s in 2:Tl) {
    Gp <- log(tpm_tv_num(B, Bx, dd$x[r[s - 1L]], 2L))
    for (j in 1:2) {
      z <- Gp[, j] + A[s - 1L, ]
      A[s, j] <- max(z) + log(sum(exp(z - max(z))))
    }
    A[s, ] <- A[s, ] + lpm[r[s], ]
  }
  Bm[Tl, ] <- 0
  for (s in (Tl - 1):1) {
    Gp <- log(tpm_tv_num(B, Bx, dd$x[r[s]], 2L))
    for (i in 1:2) {
      z <- Gp[i, ] + lpm[r[s + 1L], ] + Bm[s + 1L, ]
      Bm[s, i] <- max(z) + log(sum(exp(z - max(z))))
    }
  }
  L <- A + Bm
  Pref[r, ] <- exp(L - apply(L, 1, max)) / rowSums(exp(L - apply(L, 1, max)))
  # Viterbi
  Dl <- matrix(0, Tl, 2); Ptr <- matrix(0L, Tl, 2)
  Dl[1, ] <- log(dl) + lpm[r[1], ]
  for (s in 2:Tl) {
    Gp <- log(tpm_tv_num(B, Bx, dd$x[r[s - 1L]], 2L))
    for (j in 1:2) {
      z <- Gp[, j] + Dl[s - 1L, ]
      Ptr[s, j] <- which.max(z); Dl[s, j] <- max(z)
    }
    Dl[s, ] <- Dl[s, ] + lpm[r[s], ]
  }
  path <- integer(Tl); path[Tl] <- which.max(Dl[Tl, ])
  for (s in (Tl - 1):1) path[s] <- Ptr[s + 1L, path[s + 1L]]
  Vref[r] <- path
}
cat("\nhmm_probs vs brute force  : max abs diff",
    max(abs(P - Pref)), "\n")
cat("hmm_viterbi vs brute force: disagreements", sum(v != Vref), "\n")

## --- fitted / residuals / predict -------------------------------------
fv <- fitted(fb)
emean <- as.numeric(Pref %*% mu)
cat("\nfitted vs occupancy-weighted mean: max abs diff",
    max(abs(fv - emean)), "\n")
cat("cor(fitted, y):", round(cor(fv, dd$y), 4),
    " (probe C rung-1 reported a CONSTANT fitted, cor undefined)\n")
cat("range(fitted):", round(range(fv), 4),
    " range(y):", round(range(dd$y), 4), "\n")
cat("predict(type='response') identical to fitted():",
    isTRUE(all.equal(as.numeric(predict(fb, type = "response")),
                     as.numeric(fv))), "\n")
r1 <- residuals(fb)
cat("residuals response == y - fitted:",
    isTRUE(all.equal(as.numeric(r1), dd$y - as.numeric(fv))), "\n")
r2 <- residuals(fb, type = "pearson")
cat("pearson sd:", round(sd(r2), 4), "\n")

## --- simulate ----------------------------------------------------------
set.seed(7)
s1 <- simulate(fb, nsim = 2)
cat("\nsimulate dims:", dim(s1), " mean/sd:", round(mean(s1[[1]]), 3),
    round(sd(s1[[1]]), 3), " data:", round(mean(dd$y), 3),
    round(sd(dd$y), 3), "\n")

## --- refusals -----------------------------------------------------------
for (ex in list(
  quote(predict(fb, type = "response", se.fit = TRUE)),
  quote(predict(fb, type = "response", newdata = dd)),
  quote(residuals(fb, type = "osa")),
  quote(residuals(fb, type = "deviance")),
  quote(conditional_effects(fb)))) {
  r <- try(eval(ex), silent = TRUE)
  cat("\n", deparse(ex), "->\n  ",
      substr(conditionMessage(attr(r, "condition")), 1, 120), "...\n")
}
cat("\npredict(dpar='mu2') head:",
    round(head(as.numeric(predict(fb, dpar = "mu2")), 3), 5), "\n")
cat("predict(dpar='tr12') head:",
    round(head(as.numeric(predict(fb, dpar = "tr12")), 3), 5), "\n")
