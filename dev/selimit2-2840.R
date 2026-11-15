# Lane selimit2, backlog items 2840 and 2888 (rerun on 0.82.0): the
# boundary correlation of test-review-v29.R's `(1 | q | id)` fit on
# ulp-moved data (dev/ci701-rev-tstart.R's mdata(), seeds 8000 + k), from
# the default start. Per fit: code, convergence warnings, the correlation
# parameter t, the lost set and the two log-sd standard errors, which
# the rank-one model at the boundary gives as 0.217708 and 0.222476
# (dev/ci701-findings.md).
#   Rscript dev/selimit2-2840.R <lib or base> <k0> <k1>
a <- commandArgs(trailingOnly = TRUE)
.libPaths(c(if (a[1] != "base") a[1], "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
cat("lib", dirname(find.package("frmtmb")), "blas", La_library(), "\n")
k0 <- as.integer(a[2]); k1 <- as.integer(a[3])
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
    set.seed(8000 + k)
    for (v in c("x", "y1", "y2")) {
      dd[[v]] <- dd[[v]] * (1 + sample(0:4, nrow(dd), TRUE) * 2^-52)
    }
  }
  dd
}
f <- mvbf(bf(y1 ~ x + (1 | q | id)) + gaussian(),
          bf(y2 ~ x + (1 | q | id)) + gaussian())
for (k in k0:k1) {
  w <- character(0)
  fit <- withCallingHandlers(frm(f, data = mdata(k)),
    warning = function(x) {
      w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning")
    }, message = function(m) invokeRestart("muffleMessage"))
  th <- fit$estimates$theta
  se <- suppressWarnings(sqrt(diag(vcov(fit, full = TRUE))))
  sth <- se[names(se) == "theta" | grepl("^sd_|^theta", names(se))]
  lost <- ns$sdr_of(fit)$se_lost
  conv <- any(grepl("Optimizer did not report|Large maximum", w))
  if (k == k0) cat(paste0("  W: ", substr(w, 1, 300), "\n"), sep = "")
  cat(sprintf("k %2d code %d conv %d t %10.4g lost %-22s se_th %s ll %.10f\n",
              k, fit$opt$convergence, conv, th[length(th)],
              paste(names(lost), lost, sep = ":", collapse = ","),
              paste(sprintf("%.5f", sth), collapse = ","),
              as.numeric(logLik(fit))))
}
