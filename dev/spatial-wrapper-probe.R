# Probe: what an unknown covariance wrapper around a bar term does,
# foo(1 | g), on the library in SPATIAL_LIBS.
libs <- strsplit(Sys.getenv("SPATIAL_LIBS"), ";", fixed = TRUE)[[1]]
.libPaths(c(libs, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
set.seed(1)
d <- data.frame(g = gl(10, 5), x = rnorm(50))
d$y <- rnorm(50) + rnorm(10)[d$g]
for (fo in c("y ~ x + foo(1 | g)", "y ~ x + xcs(0 + x | g)", "y ~ x + us_t(1 | g)", "y ~ x + diag(1 | g)",
             "y ~ x + (1 | g)")) {
  r <- tryCatch({
    f <- frm(bf(stats::as.formula(fo)), family = gaussian(), data = d)
    paste("FIT covstruct =", f$frame$re_blocks[[1]]$covstruct,
          "logLik =", format(as.numeric(logLik(f)), digits = 10))
  }, error = function(e) paste("ERR", conditionMessage(e)))
  cat(fo, "->", r, "\n")
}
