# Reviewer re-check, item (b), the 0.64.0 arm. 0.64.0 has no
# smooth_b_idx(), so its fitted_point_se() configuration is replicated
# here as that release spells it: nothing at re_formula = NA, and
# re_governed_b() intersected with re_used_b() when re_formula keeps the
# group effects.
.libPaths(c("C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true")
suppressMessages(library(frmtmb))
cat("ARM: base | frmtmb:", find.package("frmtmb"), "\n")
ns <- asNamespace("frmtmb")
gt <- function(nm) get(nm, envir = ns)
cat("smooth_b_idx present:",
    exists("smooth_b_idx", envir = ns, inherits = FALSE),
    "| re_batch_try present:",
    exists("re_batch_try", envir = ns, inherits = FALSE), "\n\n")

base_cfg <- function(fit, nd, rf) {
  b_idx <- if (gt("re_form_keeps")(rf)) {
    gov <- gt("re_governed_b")(fit)
    used <- gt("re_used_b")(fit, nd, NULL, FALSE)
    if (is.null(used)) gov else intersect(gov, used)
  }
  bt <- if (length(b_idx)) {
    gt("re_b_batches")(fit, nd, NULL, FALSE, b_idx)
  }
  list(b_idx = b_idx, bt = bt)
}
counted <- function(fit, f, b_idx, bt) {
  n <- 0L
  g <- function(x) { n <<- n + 1L; f(x) }
  list(se = as.vector(gt("fit_fd_se")(fit, g, b_idx = b_idx, b_batch = bt)),
       evals = n)
}
report <- function(lab, fit, nd, rf) {
  f <- function(x) gt("fitted_point")(x, nd, rf, "response", NULL, NULL,
                                      FALSE)
  cfg <- base_cfg(fit, nd, rf)
  sh <- counted(fit, f, cfg$b_idx, cfg$bt)
  ft <- as.vector(fitted(fit, newdata = nd, re_formula = rf)[, "Est.Error", ])
  cat(sprintf(paste0("%-34s %-9s evals %-5d b %-4d batches %-5s | equals",
                     " fitted() %s\n"),
              lab, if (is.null(rf)) "NULL" else "NA", sh$evals,
              length(cfg$b_idx %||% integer(0)),
              if (is.null(cfg$bt)) "NULL" else as.character(length(cfg$bt)),
              identical(sh$se, ft)))
}
`%||%` <- function(a, b) if (is.null(a)) b else a
ord <- function(v) {
  factor(cut(v, c(-Inf, -0.8, 0.8, Inf), labels = FALSE), ordered = TRUE)
}

set.seed(41)
nD <- 160
dG <- data.frame(x = stats::runif(nD, -2, 2))
dG$y <- ord(1.1 * sin(2 * dG$x) + stats::rlogis(nD))
fG <- suppressWarnings(frm(bf(y ~ gp(x)), family = cumulative(), data = dG))
cat("gp(x) n_b =", length(fG$estimates[["b"]]), "\n")
report("A gp(x), newdata 3 NEW positions", fG,
       data.frame(x = c(-1.234, 0.777, 1.611)), NA)
report("B gp(x), in sample", fG, NULL, NA)
set.seed(31)
ngD <- 40; perD <- 20
dD <- data.frame(g = factor(rep(seq_len(ngD), each = perD)),
                 x = stats::runif(ngD * perD, -2, 2))
dD$y <- ord(sin(1.5 * dD$x) + stats::rnorm(ngD, 0, 0.8)[dD$g] +
              stats::rlogis(nrow(dD)))
fD <- suppressWarnings(frm(bf(y ~ s(x, k = 8) + (1 | g)),
                          family = cumulative(), data = dD))
report("D s(x,k=8)+(1|g) 40 lv, in sample", fD, NULL, NULL)
report("E s(x,k=8)+(1|g) 40 lv, in sample", fD, NULL, NA)
cat("DONE\n")
