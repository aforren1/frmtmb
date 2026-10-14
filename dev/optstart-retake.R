# Lane optstart: what the retake of fd_hessian() does on a matrix that
# is singular by construction (test-se-tier.R's curved ridge) and on
# R1's oversized steps (brms_missings fit_imp2 with 0.70.0's units).
# Per case: eigenvalues of the first pass and of the retake, the noise
# (spectral norm of E, the asymmetry), and the columns retaken.
#   Rscript dev/optstart-retake.R <lib>
.libPaths(c(commandArgs(TRUE)[1], "C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
cat("lib", dirname(find.package("frmtmb")), "\n")
show <- function(label, p, gr) {
  n <- length(p)
  H1 <- matrix(0, n, n)
  for (i in seq_len(n)) {
    a <- p; b <- p
    a[i] <- a[i] + 1e-3; b[i] <- b[i] - 1e-3
    H1[, i] <- (as.numeric(gr(a)) - as.numeric(gr(b))) / 2e-3
  }
  tau <- 0.5 * diag(H1) * 1e-6
  fh <- ns$fd_hessian(p, gr)
  e1 <- eigen((H1 + t(H1)) / 2, symmetric = TRUE, only.values = TRUE)$values
  e2 <- eigen(fh$H, symmetric = TRUE, only.values = TRUE)$values
  cat(sprintf("%s\n  tau above 1e-2: %s\n  first pass eigen %s | noise %.3g\n  returned eigen %s | noise %.3g\n",
              label, paste(which(tau > 1e-2), collapse = ","),
              paste(signif(range(e1), 4), collapse = " to "),
              norm(abs(H1 - t(H1)) / 2, "2"),
              paste(signif(range(e2), 4), collapse = " to "),
              norm(fh$E, "2")))
}
set.seed(77)
dn <- data.frame(g = factor(rep(seq_len(8), each = 10)))
dn$x <- rnorm(nrow(dn))
u <- rnorm(8, 0, 0.6)[dn$g]
dn$yn <- 1 + 0.3 * dn$x + exp(0.3 + u)^0.8 + rnorm(nrow(dn), 0, 0.3)
f <- suppressMessages(suppressWarnings(
  frm(bf(yn ~ c0 + exp(a)^k, c0 ~ 1 + x, a ~ 1 + (1 | g), k ~ 1,
         nl = TRUE), data = dn)))
uu <- f$par_units %||% rep(1, length(f$opt$par))
show(paste("ridge, units", paste(signif(uu, 3), collapse = " ")),
     f$opt$par / uu, function(q) as.numeric(f$obj$gr(q * uu)) * uu)
data("nhanes", package = "mice")
f2 <- suppressWarnings(suppressMessages(frm(bf(bmi | mi() ~ age * mi(chl)) +
  bf(chl | mi() ~ age) + set_rescor(FALSE), data = nhanes)))
u8 <- c(8, 8, 8, 8, 128, 128, 1, 1)
show("fit_imp2, 0.70.0's units", f2$opt$par / u8,
     function(q) as.numeric(f2$obj$gr(q * u8)) * u8)
