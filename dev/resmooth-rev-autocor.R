# Reviewer, claim 1, the autocor(cov = FALSE) branch. The reference is
# Monte Carlo over the same joint covariance, with the differenced set
# taken from the frame rather than from smooth_b_idx().
ARM <- if (identical(Sys.getenv("REVLIB"), "base")) "base" else "lane"
LIB <- if (ARM == "base") "C:/Users/adf44/source/r/rellib-r3" else
  "C:/Users/adf44/source/r/wt-resmooth-lib"
.libPaths(c(LIB, "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true")
suppressMessages(library(frmtmb))
cat("ARM:", ARM, "| frmtmb from:", find.package("frmtmb"), "\n")
ns <- asNamespace("frmtmb")
gt <- function(nm) get(nm, envir = ns)
NDRAW <- 6000L

set.seed(77)
nt <- 8; nsub <- 30
d <- expand.grid(t = seq_len(nt), subj = factor(seq_len(nsub)))
d$x <- stats::runif(nrow(d), -2, 2)
e <- as.vector(apply(matrix(stats::rnorm(nt * nsub), nt, nsub), 2,
                     function(z) as.vector(stats::filter(z, 0.6,
                                                         "recursive"))))
d$y <- sin(1.5 * d$x) + stats::rnorm(nsub, 0, 0.6)[d$subj] + 0.5 * e
fit <- suppressWarnings(frm(bf(y ~ s(x, k = 8) + (1 | subj) +
                                 ar(t, subj)), family = gaussian(),
                           data = d))
cat("blocks:", paste(vapply(fit$frame[["re_blocks"]],
                            function(b) b[["covstruct"]], ""),
                     collapse = ","), "| n_b:",
    length(fit$estimates[["b"]]), "\n")
cat("autocor is cond:",
    gt("autocor_is_cond")(fit$frame[["autocor"]][["y"]]), "\n")

ref_b <- function(fit, keep_re) {
  bl <- fit$frame[["re_blocks"]]
  is_sm <- vapply(bl, function(b) {
    b[["covstruct"]] %in% c("smooth", "gp", "hsgp")
  }, NA)
  take <- if (keep_re) rep(TRUE, length(bl)) else is_sm
  sort(unique(unlist(lapply(bl[take], `[[`, "b_idx"))))
}
mc_sd <- function(fit, f, idx, ndraw = NDRAW, seed = 13) {
  map <- gt("fit_draw_space")(fit)$map
  v0 <- gt("fit_outer_vector")(fit, map)
  p <- length(v0); b0 <- fit$estimates[["b"]]
  V <- gt("fd_joint_cov")(fit, map, idx)
  V <- (V + t(V)) / 2
  L <- t(chol(V + diag(1e-10, nrow(V))))
  set.seed(seed)
  m0 <- as.vector(f(fit))
  out <- matrix(NA_real_, ndraw, length(m0))
  for (i in seq_len(ndraw)) {
    z <- as.vector(L %*% stats::rnorm(nrow(V)))
    g <- gt("fit_set_outer")(fit, v0 + z[seq_len(p)], map)
    g$estimates[["b"]][idx] <- b0[idx] + z[p + seq_along(idx)]
    g$cache <- new.env(parent = emptyenv())
    out[i, ] <- as.vector(f(g))
  }
  apply(out, 2, stats::sd)
}

rows <- c(1, 4, 9, 30, 100)
for (rf in list(NA, NULL)) {
  lab <- if (is.null(rf)) "NULL" else "NA"
  got <- frm_linpred(fit, re_formula = rf, se.fit = TRUE)
  se <- got$se.fit[rows]
  f <- function(x) {
    frm_linpred(x, re_formula = rf, se.fit = FALSE)
  }
  idx <- ref_b(fit, is.null(rf))
  mc <- mc_sd(fit, f, idx)[rows]
  cat(sprintf("re_formula = %-4s b in reference: %3d\n", lab, length(idx)))
  cat("  shipped se.fit:", sprintf("%.5f", se), "\n")
  cat("  Monte Carlo sd:", sprintf("%.5f", mc), "\n")
  cat(sprintf("  ratio: min %.4f med %.4f max %.4f\n", min(se / mc),
              stats::median(se / mc), max(se / mc)))
}
saveRDS(lapply(list(NA, NULL), function(rf) {
  frm_linpred(fit, re_formula = rf, se.fit = TRUE)$se.fit
}), file.path("dev", paste0("resmooth-rev-autocor-", ARM, ".rds")))
cat("DONE\n")
