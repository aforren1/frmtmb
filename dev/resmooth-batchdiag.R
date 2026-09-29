# Lane wt-resmooth, nits round. WHY batching is abandoned on the three
# fits the reviewer timed, counted rather than reasoned about: per block,
# whether it is positionwise, whether a row loads more than one of its
# columns at a position, and how many coefficients each route ends up
# differencing.
#   Rscript dev/resmooth-batchdiag.R > dev/resmooth-batchdiag.txt
.libPaths(c("C:/Users/adf44/source/r/wt-resmooth-lib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true")
suppressMessages(library(frmtmb))
cat("frmtmb from:", find.package("frmtmb"), "\n")

diag_one <- function(lab, fit, newdata, re_formula) {
  cat("==", lab, "| re_formula =",
      if (is.null(re_formula)) "NULL" else "NA", "\n")
  b_idx <- {
    want <- frmtmb:::smooth_b_idx(fit)
    if (frmtmb:::re_form_keeps(re_formula)) {
      want <- sort(unique(c(want, frmtmb:::re_governed_b(fit))))
    }
    if (!length(want)) NULL else {
      used <- frmtmb:::re_used_b(fit, newdata, NULL, FALSE)
      if (is.null(used)) want else intersect(want, used)
    }
  }
  cat("   b_idx length:", length(b_idx), "of",
      length(fit$estimates[["b"]]), "\n")
  S <- frmtmb:::re_row_support(fit, newdata, NULL, FALSE)
  cat("   re_row_support:", if (is.null(S)) "NULL" else
    paste(dim(S), collapse = "x"), "\n")
  for (i in seq_along(fit$frame[["re_blocks"]])) {
    bk <- fit$frame[["re_blocks"]][[i]]
    ci <- bk[["c_idx"]]
    bi <- bk[["b_idx"]]
    if (!any(bi %in% b_idx)) {
      cat(sprintf("   block %d %-7s not in b_idx\n", i, bk[["covstruct"]]))
      next
    }
    pw <- frmtmb:::block_b_positionwise(bk)
    D <- max(1L, bk[["dim"]])
    mx <- NA
    if (!is.null(S) && pw && length(ci) == length(bi)) {
      for (k in seq_len(D)) {
        at <- seq.int(k, length(ci), by = D)
        keep <- at[bi[at] %in% b_idx]
        if (!length(keep)) next
        mx <- max(mx, max(Matrix::rowSums(S[, ci[keep], drop = FALSE])),
                  na.rm = TRUE)
      }
    }
    cat(sprintf(paste0("   block %d %-7s dim %d, %d coefs, positionwise ",
                       "%-5s, max row load %s\n"),
                i, bk[["covstruct"]], D, length(ci), pw,
                if (is.na(mx)) "n/a" else format(mx)))
  }
  bt <- frmtmb:::re_b_batches(fit, newdata, NULL, FALSE, b_idx)
  cat("   re_b_batches:",
      if (is.null(bt)) "NULL (batching abandoned for the whole fit)" else
        paste(length(bt), "batches covering",
              length(unique(unlist(lapply(bt, `[[`, "idx")))),
              "coefficients"), "\n")
}

set.seed(41)
nD <- 160
dD <- data.frame(x = stats::runif(nD, -2, 2))
latD <- 1.1 * sin(2 * dD$x) + stats::rlogis(nD)
dD$y <- factor(cut(latD, c(-Inf, -0.8, 0.8, Inf), labels = FALSE),
               ordered = TRUE)
fD <- suppressWarnings(frm(bf(y ~ gp(x)), family = cumulative(), data = dD))
diag_one("cumulative gp(x), n_b = 160", fD,
         data.frame(x = c(-1.5, 0, 1.5)), NA)

set.seed(31)
ng <- 40
per <- 20
dB <- data.frame(g = factor(rep(seq_len(ng), each = per)),
                 x = stats::runif(ng * per, -2, 2))
latB <- sin(1.5 * dB$x) + stats::rnorm(ng, 0, 0.8)[dB$g] +
  stats::rlogis(nrow(dB))
dB$y <- factor(cut(latB, c(-Inf, -0.8, 0.8, Inf), labels = FALSE),
               ordered = TRUE)
fB <- suppressWarnings(frm(bf(y ~ s(x, k = 8) + (1 | g)),
                           family = cumulative(), data = dB))
diag_one("cumulative s(x,k=8)+(1|g), 40 levels", fB, NULL, NULL)
diag_one("cumulative s(x,k=8)+(1|g), 40 levels", fB, NULL, NA)

set.seed(59)
dH <- data.frame(g = factor(rep(seq_len(ng), each = per)),
                 x = stats::runif(ng * per, -2, 2))
latH <- 0.8 * dH$x + stats::rnorm(ng, 0, 0.7)[dH$g] +
  stats::rlogis(nrow(dH))
dH$y <- factor(cut(latH, c(-Inf, -0.8, 0.8, Inf), labels = FALSE),
               ordered = TRUE)
fH <- suppressWarnings(frm(bf(y ~ x + (1 | g)), family = cumulative(),
                           data = dH))
diag_one("CONTROL cumulative x+(1|g), 40 levels", fH, NULL, NULL)
