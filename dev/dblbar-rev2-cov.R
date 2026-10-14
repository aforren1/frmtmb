# Reviewer round 2, lane dblbar: the gr(cov =)/gr(prec =) structure
# check now in the frame (within_struct). Multi-column cases must be
# refused, one-column cases must fit, and the fits are compared with the
# base where the base was right. Data: dev/dblbar-rev-data.R (seed
# 20261013) plus tim. Usage: Rscript ... base|lane
# Output: dev/dblbar-rev2-log/cov-<arm>.txt and .rds
arm <- commandArgs(TRUE)[1]
lib <- if (identical(arm, "base")) character(0) else
  "C:/Users/adf44/source/r/wt-dblbar-lib"
.libPaths(c(lib, "C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
wt <- "C:/Users/adf44/source/r/frmtmb-wt-dblbar"
source(file.path(wt, "dev/dblbar-rev-data.R"))
suppressPackageStartupMessages(library(frmtmb))
dir.create(file.path(wt, "dev/dblbar-rev2-log"), showWarnings = FALSE)
sink(file.path(wt, sprintf("dev/dblbar-rev2-log/cov-%s.txt", arm)),
     split = TRUE)
cat("lib:", find.package("frmtmb"), "\n")
dd <- dblbar_data()
set.seed(99)
dd$tim <- factor(sample(sprintf("t%d", 1:4), nrow(dd), TRUE))
A <- diag(24) * 0.5 + 0.5
dimnames(A) <- list(levels(dd$g), levels(dd$g))
Q <- solve(A)
cases <- c(
  # several columns: every one must be refused on the lane
  m_dbar_f = "y ~ f + (1 + f || gr(g, cov = A))",
  m_corF_f = "y ~ f + (1 + f | gr(g, cov = A, cor = FALSE))",
  m_diag_x = "y ~ x + diag(1 + x | gr(g, cov = A))",
  m_diag_f = "y ~ f + diag(0 + f | gr(g, cov = A))",
  m_hd_f = "y ~ f + homdiag(0 + f | gr(g, cov = A))",
  m_hd_poly = "y ~ poly(x, 2) + homdiag(0 + poly(x, 2) | gr(g, cov = A))",
  m_dbar_poly = "y ~ poly(x, 2) + (1 + poly(x, 2) || gr(g, cov = A))",
  m_dbar_fh = "y ~ f:h + (0 + f:h || gr(g, cov = A))",
  m_diag_fh = "y ~ f:h + diag(0 + f:h | gr(g, cov = A))",
  m_hd_fh = "y ~ f:h + homdiag(0 + f:h | gr(g, cov = A))",
  m_hd_xz = "y ~ x * z + homdiag(0 + x * z | gr(g, cov = A))",
  m_ar1 = "y ~ 1 + ar1(0 + tim | gr(g, cov = A))",
  m_cs = "y ~ 1 + cs(0 + tim | gr(g, cov = A))",
  m_toep = "y ~ 1 + toep(0 + tim | gr(g, cov = A))",
  m_rr = "y ~ 1 + rr(0 + tim | gr(g, cov = A), d = 1)",
  m_prec_f = "y ~ f + (1 + f || gr(g, prec = Q))",
  m_prec_hd = "y ~ f + homdiag(0 + f | gr(g, prec = Q))",
  # the correlated known-covariance term: must fit, as on the base
  u_x = "y ~ x + (1 + x | gr(g, cov = A))",
  u_f = "y ~ f + (1 + f | gr(g, cov = A))",
  u_prec_f = "y ~ f + (1 + f | gr(g, prec = Q))",
  # one column: must fit, and equal the plain one-column gr(cov =) term
  o_us_x = "y ~ x + (0 + x | gr(g, cov = A))",
  o_hd_x = "y ~ x + homdiag(0 + x | gr(g, cov = A))",
  o_diag_x = "y ~ x + diag(0 + x | gr(g, cov = A))",
  o_ar1_x = "y ~ x + ar1(0 + x | gr(g, cov = A))",
  o_us_1 = "y ~ x + (1 | gr(g, cov = A))",
  o_hd_1 = "y ~ x + homdiag(1 | gr(g, cov = A))",
  o_diag_1 = "y ~ x + diag(1 | gr(g, cov = A))",
  o_prec_1 = "y ~ x + (1 | gr(g, prec = Q))",
  o_prec_hd1 = "y ~ x + homdiag(1 | gr(g, prec = Q))",
  o_dbar_x = "y ~ x + (1 + x || gr(g, cov = A))",
  o_split = "y ~ x + (1 | gr(g, cov = A)) + (0 + x | gr(g, cov = A))",
  o_dbar_prec = "y ~ x + (1 + x || gr(g, prec = Q))",
  o_dbar_mo = "y ~ mo(xo) + (1 + mo(xo) || gr(g, cov = A))"
)
res <- list()
for (nm in names(cases)) {
  fit <- tryCatch(suppressWarnings(suppressMessages(
    frm(bf(stats::as.formula(cases[[nm]])), family = gaussian(), data = dd,
        data2 = list(A = A, Q = Q)))), error = function(e) e)
  if (inherits(fit, "error")) {
    res[[nm]] <- list(error = conditionMessage(fit))
    cat(sprintf("%-11s %-64s ERROR %s\n", nm, cases[[nm]],
                substr(gsub("\n", " ", conditionMessage(fit)), 1, 100)))
    next
  }
  bl <- vapply(fit$frame$re_blocks, function(b) {
    paste0(b$covstruct, "/", b$dim)
  }, "")
  res[[nm]] <- list(ll = as.numeric(logLik(fit)), blocks = bl,
                    est = unlist(fit$estimates),
                    vars = grep("^(sd|cor)_", variables(fit), value = TRUE))
  cat(sprintf("%-11s %-64s ok %s ll %.10f %s\n", nm, cases[[nm]],
              paste(bl, collapse = ","), res[[nm]]$ll,
              paste(res[[nm]]$vars, collapse = " ")))
}
saveRDS(res, file.path(wt, sprintf("dev/dblbar-rev2-log/cov-%s.rds", arm)))
sink()
