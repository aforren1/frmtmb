# Short check of lane ci082 after fix round 2.
#   Rscript dev/reviews/ci082-short-r2.R
# Part A: for each structure, the rows varcorr_fixed_rows() treats as log
#   sds (is_cor FALSE) against the "sd" rows varcorr_trans_rows() reports
#   for that block, and the fixed flags.
# Part B: can an estimated log-sd row be called fixed? A registered
#   structure whose sd is sqrt(c + exp(2 theta)), at theta = -30, where
#   exp(2 theta) is absorbed: c = 1 (log sd exactly 0) and c = 4 (log sd
#   log 2), with the default start (0) and with a start at -30.
.libPaths(c("C:/Users/adf44/source/r/wt-ci082-lib",
            "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", dirname(system.file(package = "frmtmb")), "\n")
ns <- asNamespace("frmtmb")
cat("probe has is_cor:",
    any(grepl("is_cor", deparse(ns$varcorr_fixed_rows))), "\n")

q <- function(expr) suppressWarnings(suppressMessages(expr))
set.seed(2)
G <- 40L
d <- expand.grid(t = factor(c("a", "b", "c")), rep = 1:3,
                 g = factor(sprintf("g%02d", seq_len(G))))
u <- matrix(rnorm(3 * G), G) %*% chol(0.5^abs(outer(1:3, 1:3, "-")))
d$y <- 1 + u[cbind(as.integer(d$g), as.integer(d$t))] +
  rnorm(nrow(d), 0, 0.4)
d$tn <- num_factor(rep(c(0, 1, 3), length.out = nrow(d)))
d$x <- rnorm(nrow(d))

# registered: sd-only (diag-like, dim 3, not in uncorrelated_structures)
# and one fixed correlation (xcs at 0.5)
frmtmb_register_covstruct(
  "xdg",
  lpdf = function(b, theta, blk) {
    dim(b) <- c(3L, length(b) %/% 3L)
    sum(RTMB::dnorm(b, 0, exp(theta), log = TRUE))
  },
  vcov = function(theta, blk) {
    V <- diag(exp(2 * theta), 3L)
    dimnames(V) <- list(blk[["cnms"]], blk[["cnms"]])
    V
  },
  npar = 3L, sd_idx = 1:3)
xc <- matrix(0.5, 3, 3) + diag(0.5, 3)
frmtmb_register_covstruct(
  "xcs",
  lpdf = function(b, theta, blk) {
    dim(b) <- c(3L, length(b) %/% 3L)
    sum(RTMB::dmvnorm(t(b), 0, exp(2 * theta[1]) * xc, log = TRUE))
  },
  vcov = function(theta, blk) {
    V <- exp(2 * theta[1]) * xc
    dimnames(V) <- list(blk[["cnms"]], blk[["cnms"]])
    V
  },
  npar = 1L, sd_idx = 1L)

forms <- list(
  us = y ~ 1 + (0 + t | g), cs = y ~ 1 + cs(0 + t | g),
  homcs = y ~ 1 + homcs(0 + t | g), ar1 = y ~ 1 + ar1(0 + t | g),
  toep = y ~ 1 + toep(0 + t | g), diag = y ~ 1 + diag(0 + t | g),
  ou = y ~ 1 + ou(tn + 0 | g), us1 = y ~ 1 + (1 | g),
  usx = y ~ x + (1 + x | g), xdg = y ~ 1 + xdg(0 + t | g),
  xcs = y ~ 1 + xcs(0 + t | g))
for (nm in names(forms)) {
  fit <- tryCatch(q(frm(forms[[nm]], family = gaussian(), data = d)),
                  error = function(e) {
                    cat("FIT ERROR", nm, conditionMessage(e), "\n")
                    NULL
                  })
  if (is.null(fit)) next
  tr <- ns$varcorr_trans_rows(fit)
  bk <- fit$frame[["re_blocks"]][[1]]
  t0 <- fit$estimates[["theta"]][bk[["theta_idx"]]]
  e <- ns$cs_entry(bk[["covstruct"]])
  n_sd_probe <- NROW(e$vcov(t0, bk))
  fx <- ns$varcorr_fixed_rows(bk, t0)
  n_sd_tr <- sum(tr$type == "sd")
  cat(sprintf(paste("RESULT %-5s covstruct %-5s dim %d has_cor %s rows %d",
                    "sd(probe) %d sd(trans_rows) %d match %s fixed %s\n"),
              nm, bk[["covstruct"]], bk[["dim"]], ns$block_has_cor(bk),
              nrow(tr), n_sd_probe, n_sd_tr,
              n_sd_probe == n_sd_tr && length(fx) == nrow(tr),
              paste(as.integer(fx), collapse = "")))
}

for (cc in c(1, 4)) for (st in c(NA, -30)) {
  nm <- sprintf("xab%d%s", cc, if (is.na(st)) "d" else "s")
  args <- list(
    nm,
    lpdf = function(b, theta, blk) 0,
    vcov = local({
      cc <- cc
      function(theta, blk) {
        V <- diag(cc + exp(2 * theta[1]), 2L)
        dimnames(V) <- list(blk[["cnms"]], blk[["cnms"]])
        V
      }
    }),
    npar = 1L, sd_idx = 1L)
  if (!is.na(st)) args$start <- function(blk) -30
  do.call(frmtmb_register_covstruct, args)
  bk <- list(covstruct = nm, dim = 2L, cnms = c("a", "b"), theta_idx = 1L)
  fx <- ns$varcorr_fixed_rows(bk, -30)
  cat(sprintf("RESULT absorbed c %g start %s log sd %.17g fixed %s\n", cc,
              if (is.na(st)) "default" else "-30",
              log(sqrt(cc + exp(-60))), paste(as.integer(fx), collapse = "")))
}
