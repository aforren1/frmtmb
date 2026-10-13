# Lane optrap, punch round 1 (m6): huber() seed 19 of
# dev/optrap-rev-scale.R at scale 1. Step 1 (lib = rellib-r7) saves the
# point where r7 stopped at false convergence; step 2 (lib = the lane)
# fits at scale 1000, then walks the objective on the segment from the lane's optimum
# to r7's point, to see whether the higher log likelihood is an optimum
# or the edge of a jump in the Laplace objective.
#   Rscript dev/optrap-p1-huber19.R <lib> <r7|lane>
a <- commandArgs(TRUE)
.libPaths(c(a[1], "C:/Users/adf44/source/r/rellib-r7",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
set.seed(19)
d <- data.frame(g = factor(rep(1:20, each = 5)), x = rnorm(100))
b <- rnorm(20)[d$g]
d$y <- 1 + 0.5 * d$x + 0.5 * b + rt(100, 3)
sc <- if (a[2] == "r7") 1 else 1000
d$y <- d$y * sc
f <- suppressWarnings(frm(bf(y ~ x + (1 | g)), family = huber(), data = d))
cat("lib", dirname(find.package("frmtmb")), "code", f$opt$convergence,
    "logLik", format(as.numeric(logLik(f)), digits = 10), "\n")
pp <- f$opt$par
print(pp, digits = 8)
rds <- "dev/optrap-log/p1-huber19-r7.rds"
if (a[2] == "r7") {
  saveRDS(pp, rds)
} else {
  # r7's point at scale 1, moved to scale 1000: beta times sc, the log
  # sigma and the log sd plus log(sc)
  p7 <- readRDS(rds)
  p7[1:2] <- p7[1:2] * sc
  p7[3:4] <- p7[3:4] + log(sc)
  fn <- f$obj$fn
  gr <- f$obj$gr
  cat("\nfn at the lane optimum and r7's stop:",
      format(-fn(pp), digits = 10), format(-fn(p7), digits = 10), "\n")
  cat("max |gr| at each:", max(abs(gr(pp))), max(abs(gr(p7))), "\n")
  t <- seq(0, 1, length.out = 201)
  v <- vapply(t, function(s) -fn(pp + s * (p7 - pp)), 0)
  dv <- diff(v)
  cat("\nsegment: 201 points; steps of log likelihood\n")
  cat("median |step|", signif(median(abs(dv)), 3), " largest 5 |step| at t:\n")
  o <- order(-abs(dv))[1:5]
  print(data.frame(t = t[o], step = dv[o]), digits = 4, row.names = FALSE)
  # A jump is a step many times the step its neighbors take.
  # A local optimum is a sign change in the slope.
  cat("max along segment", format(max(v), digits = 10), "at t", t[which.max(v)],
      "\n")
}
invisible()
