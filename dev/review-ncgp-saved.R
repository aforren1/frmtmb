# Item 3: a saved fit read in a fresh process, tape pointers dead. Args: save|read
source("C:/Users/adf44/source/r/frmtmb-wt-ncgp/dev/review-ncgp-lib.R")
rev_libs("lane")
mode <- commandArgs(TRUE)[1]
pth <- file.path(REV_OUT, "saved-fit.rds"); pth2 <- file.path(REV_OUT, "saved-vals.rds")
rd <- function(fit, nd) list(fit = suppressWarnings(fitted(fit)), new = suppressWarnings(fitted(fit, newdata = nd)), vc = as.matrix(vcov(fit, full = TRUE)),
                             ci = suppressWarnings(confint(fit)), ran = unlist(ranef(fit)))
nd <- data.frame(x = c(1.11, 2.22, 7.77))
if (mode == "save") {
  set.seed(11); n <- 200; x <- sort(runif(n, 0, 8)); d <- data.frame(x = x, y = sin(x) + rnorm(n, sd = 0.3))
  fit <- suppressWarnings(frm(y ~ gp(x), d))
  saveRDS(fit, pth); saveRDS(rd(fit, nd), pth2)
  cat("file MB", file.size(pth) / 1e6, "\n")
} else {
  fit <- readRDS(pth); v0 <- readRDS(pth2)
  t0 <- proc.time()[[3]]
  v1 <- rd(fit, nd)
  for (k in names(v0)) cat(sprintf("RSA %-4s max rel diff %.2e\n", k, max(abs(as.numeric(v1[[k]]) - as.numeric(v0[[k]])) / pmax(abs(as.numeric(v0[[k]])), 1e-8), na.rm = TRUE)))
  cat(sprintf("RSA readers %.1f s, compress %d, matmul %s\n", proc.time()[[3]] - t0, TMB::config(DLL = "RTMB")$tmbad.sparse_hessian_compress, RTMB::TapeConfig()$matmul))
  fit$obj$retape()
  e <- fit$obj$env
  cat(sprintf("RSA retape: fn %.8f vs objective %.8f; gr max %.2e\n", fit$obj$fn(fit$opt$par), fit$opt$objective, max(abs(fit$obj$gr(fit$opt$par)))))
  fresh <- tryCatch(suppressWarnings(update(fit)), error = function(e) e)
  cat("RSA update:", if (inherits(fresh, "error")) conditionMessage(fresh) else format(logLik(fresh)), "\n")
  cat(sprintf("RSA peak working set %s MB\n", system2("powershell", c("-NoProfile", "-Command", shQuote(sprintf("(Get-Process -Id %d).PeakWorkingSet64/1MB", Sys.getpid()))), stdout = TRUE)[1]))
}
cat("DONE\n")
