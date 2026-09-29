# Two gaps of the gradient trip-wire, constructed rather than assumed.
#
#   V1  a stall the ABSOLUTE gradient cannot see, because a badly scaled
#       column makes the gradient tiny while the log likelihood left on
#       the table is large. autoscale_small_sd()'s comment records this
#       from the other side.
#   V2  a flat ridge where the headroom is above grad_tol while the
#       gradient is below it.
#
# Both are gaps the check had BEFORE this lane as well, so they bound
# what the new criterion does and does not buy.
#
#   Rscript dev/gradcheck-07-gaps.R base|lane

args <- commandArgs(trailingOnly = TRUE)
which_lib <- if (length(args)) args[[1L]] else "base"
source("dev/gradcheck-helpers.R")
gc_libpaths(which_lib)
library(frmtmb)
cat("library:", which_lib, " frmtmb", format(packageVersion("frmtmb")),
    "\n\n")

rows <- list()
add <- function(r) {
  rows[[length(rows) + 1L]] <<- r
  gc_print(r)
}
hr <- function(s) cat("\n== ", s, "\n", sep = "")

hr("V1. a column spread far BELOW one, autoscale off")
set.seed(501)
n <- 2000
for (sp in c(1e-3, 1e-5, 1e-7)) {
  dd <- data.frame(x = rnorm(n))
  dd$xs <- dd$x * sp
  dd$y <- rnorm(n, 1 + 2 * dd$x, 1)
  good <- frm(bf(y ~ x), family = gaussian(), data = dd)
  cp <- gc_catch(frm(bf(y ~ xs), family = gaussian(), data = dd,
                     control = frmtmb_control(autoscale = FALSE)))
  f <- cp$value
  add(gc_row(paste0("spread ", sp, ", autoscale off"), 501, f, cp$warnings))
  cat("  logLik shortfall against the well-scaled fit:",
      format(as.numeric(logLik(good)) - as.numeric(logLik(f)),
             digits = 6), "\n")
  cat("  warnings:", if (!length(cp$warnings)) "none" else "", "\n")
  for (w in cp$warnings) cat("   -", substr(w, 1, 100), "\n")
}

hr("V2. collinearity deep enough to hide the headroom")
for (eps in c(1e-4, 1e-5, 1e-6, 1e-7, 1e-8)) {
  set.seed(502)
  n <- 1000
  x1 <- rnorm(n)
  d4 <- data.frame(x1 = x1, x2 = x1 + rnorm(n, 0, eps))
  d4$y <- rnorm(n, 1 + 2 * d4$x1, 1)
  cp <- gc_catch(frm(bf(y ~ x1 + x2), family = gaussian(), data = d4,
                     control = frmtmb_control(restarts = 0)))
  f <- cp$value
  if (inherits(f, "gc_error")) {
    cat("  ERROR at eps =", eps, ":", substr(as.character(f), 1, 90), "\n")
    next
  }
  add(gc_row(paste0("collinear eps = ", eps), 502, f, cp$warnings))
}

saveRDS(rows, file.path("dev", paste0("gradcheck-07-", which_lib, ".rds")))
cat("\nrows:", length(rows), "  warned:",
    sum(vapply(rows, function(r) isTRUE(r$warned), TRUE)), "\n")
