# Shared header for lane reslot's probes. The first command-line argument
# picks the library stack: "before" reads rellib-r14 (frmtmb 0.76.0),
# "after" puts the lane library in front of it.
reslot_arm <- commandArgs(trailingOnly = TRUE)[1L]
if (is.na(reslot_arm)) reslot_arm <- "after"
reslot_lib <- "C:/Users/adf44/source/r/wt-reslot-lib"
reslot_base <- c("C:/Users/adf44/source/r/rellib-r14",
                 "C:/Users/adf44/AppData/Local/R/win-library/4.6")
.libPaths(if (identical(reslot_arm, "before")) reslot_base else
  c(reslot_lib, reslot_base))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win",
           FRMTMB_STAN_CACHE =
             "C:/Users/adf44/source/r/frmtmb-wt-reslot/dev/stan-cache")
suppressPackageStartupMessages(library(frmtmb))
cat("arm:", reslot_arm, " lib:", dirname(find.package("frmtmb")),
    " frmtmb", format(packageVersion("frmtmb")), "\n")

# the resolved entries of a prior, one line each: component, positions,
# density kind, scale
show_entries <- function(ri) {
  for (e in ri$entries) {
    cat(sprintf("  entry %-8s %-10s %-8s %s\n", e$comp,
                paste(e$idx, collapse = ","), e$dist$kind, e$scale))
  }
  if (length(ri$lower)) cat("  lower:", paste(names(ri$lower), ri$lower,
                                              sep = "=", collapse = " "),
                            "\n")
  if (length(ri$upper)) cat("  upper:", paste(names(ri$upper), ri$upper,
                                              sep = "=", collapse = " "),
                            "\n")
}

try_show <- function(label, expr) {
  cat("--", label, "\n")
  r <- tryCatch(expr, error = function(e) e)
  if (inherits(r, "error")) {
    cat("  ERROR:", conditionMessage(r), "\n")
  } else if (is.list(r) && !is.null(r$entries)) {
    show_entries(r)
  } else if (!is.null(r)) {
    print(r)
  }
  invisible(r)
}

re_data <- function(seed = 1L, G = 8L, n = 10L) {
  set.seed(seed)
  g <- factor(rep(seq_len(G), each = n))
  x <- rnorm(G * n)
  z <- rnorm(G * n)
  b0 <- rnorm(G, 0, 0.7)[g]
  b1 <- rnorm(G, 0, 0.4)[g]
  y <- 1 + 0.5 * x - 0.3 * z + b0 + b1 * x + rnorm(G * n, 0, 0.6)
  y2 <- 0.5 * y + rnorm(G * n)
  t <- rep(seq_len(n), G)
  data.frame(y, y2, x, z, g, t)
}
