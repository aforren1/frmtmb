# Item 4836: skew_normal with (1 | g) under quadrature, started at the
# stall, ends 5.5e-2 short (dev/quadseam-seed6.R, ported). Prints the
# stop, the escape and settle records, and the Laplace arm.
#   Rscript dev/quadseam2-seed6.R <libpath|base> [seeds]
a <- commandArgs(TRUE)
.libPaths(c(if (a[1] != "base") a[1], "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("lib:", dirname(find.package("frmtmb")), "\n")
seeds <- if (length(a) > 1) eval(parse(text = a[2])) else 6
for (seed in seeds) {
  set.seed(seed); n <- 200
  xs <- -abs(rnorm(n)) * 3
  y <- xs + (abs(rnorm(n)) - sqrt(2 / pi)) * 1.5 + rnorm(n, 0, 0.3)
  dd <- data.frame(y = y, xs = xs, g = factor(rep(1:20, length.out = n)))
  m3 <- mean((y - mean(y))^3) / sd(y)^3
  bad <- list(betad = c(log(sd(y)), 2 * sign(m3) + 0.5 * m3))
  fo <- bf(y ~ xs + (1 | g), sigma ~ 1, alpha ~ 1)
  for (q in c(FALSE, TRUE)) for (st in list(NULL, bad)) {
    w <- character()
    t0 <- proc.time()[["elapsed"]]
    f <- withCallingHandlers(suppressMessages(frm(fo, family = skew_normal(),
           data = dd, quadrature = q, start = st)),
           warning = function(x) {
             w <<- c(w, class(x)[1]); invokeRestart("muffleWarning")
           })
    t1 <- proc.time()[["elapsed"]] - t0
    cat(sprintf(paste0("SEED %d quad %-5s start %-7s ll %.9f alpha %.4f ",
                       "conv %d secs %.1f esc %s settle %s polish %s ",
                       "warn %s\n"), seed, q,
                if (is.null(st)) "default" else "stall",
                as.numeric(logLik(f)), f$estimates$betad[2],
                f$opt$convergence, t1,
                paste(format(f$opt$stationary_escape, digits = 4),
                      collapse = "/"),
                paste(format(f$opt$quad_settle, digits = 3), collapse = "/"),
                paste(format(f$opt$polish, digits = 3), collapse = "/"),
                if (length(w)) paste(w, collapse = ",") else "-"))
  }
}
cat("END\n")
