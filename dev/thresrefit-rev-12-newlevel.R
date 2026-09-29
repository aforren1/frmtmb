## REVIEW: influence(data = ) with a thres(gr = ) level the fit never saw,
## on BOTH builds, so the change in FAILURE MODE is measured and not
## inferred. Same for a response that reaches a higher category.
lib <- Sys.getenv("FRMTMB_LIB")
lib <- if (identical(lib, "base")) "C:/Users/adf44/source/r/rellib-r3" else
  "C:/Users/adf44/source/r/wt-thresrefit-lib"
.libPaths(c(lib, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("## lib =", lib, "\n\n")

say <- function(tag, expr) {
  cat("==", tag, "\n")
  out <- tryCatch(expr, error = function(e) e)
  if (inherits(out, "condition")) {
    cat("    ERROR (", class(out)[1L], "):",
        gsub("\n", " ", conditionMessage(out)), "\n\n")
  } else {
    cat("    no error;"); print(out); cat("\n")
  }
}

set.seed(1103)
k <- 96
g <- factor(rep(c("a", "b", "c"), length.out = k))
xg <- stats::rnorm(k)
tau <- list(a = c(-0.6, 0.6, 2.4), b = c(-0.5, 1.2), c = c(-0.3, 0.9))
yg <- vapply(seq_len(k), function(i) {
  1L + sum(stats::runif(1) >
             stats::plogis(tau[[as.character(g[i])]] - 0.5 * xg[i]))
}, 1L)
dg <- data.frame(x = xg, g = g, y = yg)
fg <- suppressWarnings(frm(bf(y | thres(gr = g) ~ x), family = cumulative(),
                          data = dg))
cat("fitted nthres =",
    paste(fg$spec$responses$y$family[["thres"]][["nthres"]], collapse = "/"),
    " nraw =", length(fg$estimates[["tau_raw"]]), "\n")
dnew <- dg
dnew$g <- as.character(dnew$g)
dnew$g[1:8] <- "d"
dnew$g <- factor(dnew$g)
dnew$y[1:8] <- pmin(dnew$y[1:8], 3L)
say("influence(data = ) with a NEW thres(gr = ) level 'd'",
    {
      inf <- suppressWarnings(influence(fg, data = dnew, force = TRUE))
      c(nrow = nrow(inf$fixed), ncol = ncol(inf$fixed),
        na = sum(is.na(inf$fixed)),
        warns = 0)
    })
## count the warnings separately: a silent all-NA table is the concern
w <- withCallingHandlers({
  n <- 0L
  inf <- tryCatch(influence(fg, data = dnew, force = TRUE),
                  error = function(e) NULL)
  n
}, warning = function(x) {
  cat("    WARNING seen:", conditionMessage(x), "\n")
  invokeRestart("muffleWarning")
})
cat("warnings printed above (none means the all-NA table is silent)\n")
cat("DONE rev-12\n")
