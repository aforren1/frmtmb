# Part C experiments: the 41-draw ridge tally with a candidate change
# assigned into the installed namespace (no install).
#   CODE8_ARM=before|after Rscript dev/code8-ridge-exp.R <variant> <label>
# variants: none, retake (fd_hessian keeps the retaken columns whenever a
# column's rise is above the band)
a <- commandArgs(TRUE)
variant <- a[1]; lab <- a[2]
source("dev/code8-common.R")
options(width = 200)
ns <- asNamespace("frmtmb")
if (variant == "retake") {
  fdh <- function(par, gr, tau_band = 1e-2, tau_step = 1e-4, eps = 0.001) {
    n <- length(par)
    eps <- rep_len(eps, n)
    dpar <- par
    H <- matrix(0, n, n)
    for (i in seq_len(n)) {
      dpar[i] <- dpar[i] + eps[i]
      df1 <- as.numeric(gr(dpar))
      dpar[i] <- dpar[i] - 2 * eps[i]
      df2 <- as.numeric(gr(dpar))
      H[, i] <- (1 * (df1 - df2)) / (2 * eps[i] * 1 * 1)
      dpar[i] <- dpar[i] + eps[i]
    }
    d <- diag(H)
    redo <- which(is.finite(d) & d > 0 & 0.5 * d * eps^2 > tau_band)
    for (i in redo) {
      h <- sqrt(2 * tau_step / d[i])
      dpar[i] <- par[i] + h
      df1 <- as.numeric(gr(dpar))
      dpar[i] <- par[i] - h
      df2 <- as.numeric(gr(dpar))
      H[, i] <- (df1 - df2) / (2 * h)
      dpar[i] <- par[i]
    }
    E <- abs(H - t(H)) / 2
    H <- 0.5 * (H + t(H))
    dimnames(H) <- list(names(par), names(par))
    list(H = H, E = E)
  }
  environment(fdh) <- ns
  assignInNamespace("fd_hessian", fdh, "frmtmb")
}
set.seed(77)
dn <- data.frame(g = factor(rep(seq_len(8), each = 10)))
dn$x <- rnorm(nrow(dn))
u <- rnorm(8, 0, 0.6)[dn$g]
dn$yn <- 1 + 0.3 * dn$x + exp(0.3 + u)^0.8 + rnorm(nrow(dn), 0, 0.3)
y0 <- dn$yn
form <- bf(yn ~ c0 + exp(a)^k, c0 ~ 1 + x, a ~ 1 + (1 | g), k ~ 1, nl = TRUE)
rows <- list()
for (j in -20:20) {
  dn$yn <- y0 * (1 + 2e-15 * j)
  r <- rfit(form, data = dn)
  f <- r$f
  nm <- fns$outer_par_names(f)
  sdr <- quiet(fns$sdr_of(f))
  se <- sqrt(diag(as.matrix(sdr$cov.fixed)))
  lost <- sdr$se_lost
  rows[[length(rows) + 1L]] <- data.frame(
    j = j, obj = f$opt$objective, a = f$opt$par[3], k = f$opt$par[4],
    se_a = se[3], se_c0 = se[1], se_sigma = se[5],
    lost = paste(sprintf("%s:%s", names(lost), lost), collapse = " "),
    conds = substr(gsub("Standard errors are not available for", "SE-NA",
                        cond_txt(r)), 1, 90))
}
X <- do.call(rbind, rows)
write.table(X, sprintf("dev/code8-out/ridge-exp-%s.tsv", lab), sep = "\t",
            quote = FALSE, row.names = FALSE)
X$obj <- sprintf("%.8f", X$obj)
print(X[, c("j", "obj", "a", "k", "se_a", "se_c0", "se_sigma", "lost")],
      right = FALSE, digits = 4)
cat("LOST SETS\n")
tb <- table(X$lost)
for (i in seq_along(tb)) cat(tb[[i]], ":", names(tb)[i], "\n")
cat("CONDS\n")
tb <- table(X$conds)
for (i in seq_along(tb)) cat(tb[[i]], ":", names(tb)[i], "\n")
cat("DONE\n")
