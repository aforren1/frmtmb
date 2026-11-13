# Sweep of lane ppdraw: for every core family with a zero-inflation,
# hurdle or mixture parameter, predict()'s mean against fitted()'s at an
# interior fit and at a boundary fit. z is the difference in units of the
# Monte Carlo standard error of the draws' mean (each replicate is one
# parameter draw, so the replicates' row means are independent); for a
# category-valued response it is the largest |z| over the category shares.
# Usage: Rscript dev/ppdraw-sweep.R <lib or "base"> <out.tsv>
args <- commandArgs(TRUE)
lib <- if (args[1] != "base") args[1] else NULL
.libPaths(c(lib, "C:/Users/adf44/source/r/rellib-r20",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), "\n")
out <- args[2]
B <- 1000L
n <- 300L

source("dev/ppdraw-sweep-gen.R")
rows <- list()
for (fam in fams) for (case in c("interior", "boundary")) {
  d <- gen(fam, case)
  s <- spec(fam)
  conds <- character(0)
  keep <- function(c) conds <<- c(conds, substr(gsub("\n", " ",
                                    conditionMessage(c)), 1, 70))
  fit <- tryCatch(withCallingHandlers(
    frm(s[[1]], family = s[[2]], data = d),
    warning = function(w) { keep(w); invokeRestart("muffleWarning") },
    message = function(m) { keep(m); invokeRestart("muffleMessage") }),
    error = function(e) { keep(e); NULL })
  r <- list(family = fam, case = case, fit_said = paste(conds, collapse = " / "),
            held = "", z = NA_real_, pred = NA_real_, fitted = NA_real_,
            na_share = NA_real_, z_pe0 = NA_real_, pred_said = "")
  if (!is.null(fit)) {
    psaid <- character(0)
    held <- character(0)
    set.seed(5)
    p <- tryCatch(withCallingHandlers(
      predict(fit, ndraws = B, summary = FALSE),
      frmtmb_draw_held = function(m) {
        held <<- c(held, regmatches(conditionMessage(m),
                                    gregexpr("`[^`]+`",
                                             conditionMessage(m)))[[1]])
        invokeRestart("muffleMessage")
      },
      warning = function(w) {
        psaid <<- c(psaid, substr(conditionMessage(w), 1, 70))
        invokeRestart("muffleWarning")
      }), error = function(e) {
        psaid <<- c(psaid, substr(conditionMessage(e), 1, 70))
        NULL
      })
    set.seed(5)
    p0 <- tryCatch(suppressWarnings(predict(fit, ndraws = B, summary = FALSE,
                                            propagate_error = FALSE)),
                   error = function(e) NULL)
    fe <- suppressWarnings(fitted(fit))
    r$held <- paste(unique(held), collapse = " ")
    r$pred_said <- paste(unique(psaid), collapse = " / ")
    if (!is.null(p)) {
      r$na_share <- mean(is.na(p))
      if (length(dim(fe)) == 3L) {
        K <- dim(fe)[3]
        code0 <- if (fam == "hurdle_cumulative") 0L else 1L
        zf <- function(pp) {
          sh <- t(apply(pp, 2L, function(col) {
            col <- col[!is.na(col)]
            tabulate(as.integer(col) + 1L - code0, K) / length(col)
          }))
          P <- fe[, "Estimate", ]
          z <- (sh - P) / sqrt(pmax(P * (1 - P), 1e-12) / B)
          max(abs(z))
        }
        r$z <- zf(p)
        if (!is.null(p0)) r$z_pe0 <- zf(p0)
      } else {
        zf <- function(pp) {
          m <- rowMeans(pp, na.rm = TRUE)
          (mean(m) - mean(fe[, "Estimate"])) / (stats::sd(m) / sqrt(B))
        }
        r$pred <- mean(rowMeans(p, na.rm = TRUE))
        r$fitted <- mean(fe[, "Estimate"])
        r$z <- zf(p)
        if (!is.null(p0)) r$z_pe0 <- zf(p0)
      }
    }
  }
  rows[[length(rows) + 1L]] <- as.data.frame(r)
  cat(sprintf("%-28s %-8s z=%8.2f z_pe0=%7.2f na=%.4f held=%s\n", fam, case,
              r$z, r$z_pe0, r$na_share, r$held))
}
tab <- do.call(rbind, rows)
utils::write.table(tab, out, sep = "\t", row.names = FALSE, quote = FALSE)
cat("done\n")
