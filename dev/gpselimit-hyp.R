# Lane gpselimit: the white-limit message's hypothesis() strings, evaluated.
# Fits: sweep row 3 (1-D), and two-dimensional fields (isotropic and
# anisotropic, exp_quad and matern32) on a rough surface with a small k,
# which end at the white limit. For each: the conditions, the lost reasons,
# and hypothesis() of every string the message gives.
# Usage: Rscript gpselimit-hyp.R <lane|base>
source("C:/Users/adf44/source/r/frmtmb-wt-gpselimit/dev/gpselimit-lib.R")
a <- commandArgs(TRUE)
gs_libs(a[1])
one <- function(lab, ...) {
  q <- gs_fit(...)
  fit <- q$fit
  if (inherits(fit, "error")) {
    cat("HYP", lab, "ERROR", conditionMessage(fit), "\n"); return()
  }
  ss <- gs_ses(fit)
  cat(sprintf("HYP %s | lost [%s] | warn [%s] | msg [%s]\n", lab,
              paste(paste0(names(ss$lost), "=", ss$lost), collapse = ","),
              paste(unique(q$wclass), collapse = ";"),
              paste(unique(q$mclass), collapse = ";")))
  bk <- Filter(function(b) b$covstruct == "hsgp", fit$frame$re_blocks)
  for (b in bk) {
    th <- fit$estimates$theta[b$theta_idx]
    cat("   block", b$term_label, "sd", format(exp(th[1]), digits = 3),
        "kind", paste(fit$cache$se_boundary_kind, collapse = ","),
        "ls/spacing",
        format(exp(th[-1]) / b$gp_spacing, digits = 3), "\n")
  }
  m <- q$msg[q$mclass == "frmtmb_gp_basis"]
  if (!length(m)) return(invisible())
  hs <- regmatches(m, gregexpr("\"[^\"]*\"", m))[[1]]
  hs <- gsub("\"", "", hs)
  if (!length(hs)) cat("   no hypothesis string in the message\n")
  for (h in hs) {
    r <- tryCatch(hypothesis(fit, h, class = "")$hypothesis,
                  error = function(e) conditionMessage(e))
    if (is.character(r)) {
      cat("   hypothesis", h, "ERROR", r, "\n")
    } else {
      cat("   hypothesis", h, "est", format(r$Estimate, digits = 5),
          "se", format(r$Est.Error, digits = 5), "\n")
    }
  }
}
s <- gs_sweep_row(3)
one("row3", bf(y ~ gp(x, k = s$k)), data = s$d, family = gaussian())
# two-dimensional fields drawn from the white limit itself: the basis of a
# first fit, coefficients iid N(0, 1), noise 0.1
ns <- asNamespace("frmtmb")
set.seed(5150)
n <- 150
d2 <- data.frame(x1 = runif(n), x2 = runif(n), y = rnorm(n))
for (cv in c("exp_quad", "matern32")) for (iso in c(TRUE, FALSE)) {
  form <- bf(y ~ gp(x1, x2, k = 4, iso = iso, cov = cv))
  f0 <- suppressWarnings(suppressMessages(frm(form, data = d2,
                                              family = gaussian())))
  bk <- Filter(function(b) b$covstruct == "hsgp", f0$frame$re_blocks)[[1]]
  lp <- f0$frame$linpreds[[ns$gp_block_lp(bk, f0$frame)]]
  Z <- as.matrix(ns$lp_Z(f0, lp))[, bk$c_idx, drop = FALSE]
  for (s in 1:3) {
    set.seed(5160 + s)
    dd <- d2
    dd$y <- as.vector(Z %*% rnorm(ncol(Z))) + 0.1 * rnorm(n)
    one(sprintf("2d %s iso %s seed %d", cv, iso, s), form, data = dd,
        family = gaussian())
  }
}
