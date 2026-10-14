# Lane ci701: what the standard-error check says at the boundary-
# correlation fit of test-review-v29.R (set.seed(61)) when it runs.
# On Windows the fit ends with a convergence warning, so the check never
# runs; on the macOS runner it ran. Here se_fit_analysis(force = TRUE) is
# asked at each fit's point, and the outer Hessian's eigen-structure is
# printed, for data moved by a few ulps (seed 1000 + k, as
# dev/ci701-bcor-probe.R).
#
#   Rscript dev/ci701-bcor-force.R <lib|base> <k0> <k1>
args <- commandArgs(TRUE)
lib <- args[1]; k0 <- as.integer(args[2]); k1 <- as.integer(args[3])
base <- "C:/Users/adf44/source/r/rellib-r8"
ul <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
.libPaths(if (lib == "base") c(base, ul) else c(lib, base, ul))
suppressMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), as.character(packageVersion("frmtmb")),
    "\n")

set.seed(61)
ng <- 12
u <- stats::rnorm(ng, 0, 0.8)
per <- 5
dd <- data.frame(id = factor(rep(paste0("a", seq_len(ng)), each = per)),
                 x = stats::rnorm(ng * per))
dd$y1 <- 1 + 0.5 * dd$x + u[as.integer(dd$id)] +
  stats::rnorm(nrow(dd), 0, 0.4)
dd$y2 <- -1 + 0.2 * dd$x + 0.6 * u[as.integer(dd$id)] +
  stats::rnorm(nrow(dd), 0, 0.3)
d0 <- dd
mdata <- function(k) {
  dd <- d0
  if (k > 0) {
    set.seed(1000 + k)
    for (v in c("x", "y1", "y2")) {
      dd[[v]] <- dd[[v]] * (1 + sample(0:4, nrow(dd), TRUE) * 2^-52)
    }
  }
  dd
}
f <- mvbf(bf(y1 ~ x + (1 | q | id)) + gaussian(),
          bf(y2 ~ x + (1 | q | id)) + gaussian())
for (k in k0:k1) {
  fit <- suppressWarnings(frm(f, data = mdata(k)))
  g <- drop(fit$obj$gr(fit$opt$par))
  an <- frmtmb:::se_fit_analysis(fit, force = TRUE)
  H <- fit$cache$hessian_fixed$H
  nm <- names(fit$opt$par)
  th <- which(nm == "theta")
  cat(sprintf("RES k=%d code=%d gmax=%.3g th=%s tier=%s lost=%s\n", k,
              fit$opt$convergence, max(abs(g)),
              paste(sprintf("%.5g", fit$estimates$theta), collapse = ","),
              an$tier,
              if (length(an$lost)) {
                paste(names(an$lost), an$lost, sep = ":", collapse = ",")
              } else "none"))
  if (!is.null(H)) {
    D <- 1 / sqrt(abs(diag(H)))
    Hu <- H * outer(D, D)
    e <- eigen(Hu, symmetric = TRUE)
    cat("   diagH(theta):", sprintf("%.4g", diag(H)[th]), "\n")
    cat("   eig(unit):", sprintf("%.4g", e$values), "\n")
    v <- e$vectors[, length(e$values)]
    cat("   smallest vec on theta:", sprintf("%.4f", v[th]), "\n")
    sh <- an$shown
    if (!is.null(sh)) {
      cat("   shown se(theta):", sprintf("%.4g", sqrt(diag(sh))[th]), "\n")
    }
  }
}
