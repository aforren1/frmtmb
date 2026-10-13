# Reviewer of lane optrap, consolidation check: are the Newton steps on
# fits with random effects deterministic across core types? Degree 6
# and 7 raw polynomials with (1 | g) (y x 1) and the y ~ big + (1 | g)
# design, seeds 1 to 10; objective in hex and the steps. Run pinned to P
# cores and to E cores and compare the two logs.
.libPaths(c("C:/Users/adf44/source/r/wt-optrap-lib",
            "C:/Users/adf44/source/r/rellib-r7",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
for (deg in 6:7) for (s in 1:10) {
  set.seed(1000 * deg + s)
  d <- data.frame(x = runif(200, 1, 2), g = factor(rep(1:20, 10)))
  d$y <- sin(3 * d$x) + rnorm(20, 0, 0.3)[d$g] + rnorm(200, 0, 0.2)
  for (k in 1:deg) d[[paste0("x", k)]] <- d$x^k
  fo <- as.formula(paste("y ~", paste0("x", 1:deg, collapse = " + "),
                         "+ (1 | g)"))
  f <- suppressWarnings(suppressMessages(frm(fo, data = d)))
  cat(sprintf("deg%d %d obj %a steps %s\n", deg, s, f$opt$objective,
              f$opt$polish[["steps"]] %||% 0))
}
for (s in 1:10) {
  set.seed(s)
  d <- data.frame(g = factor(rep(1:40, each = 10)), big = 1e5 + runif(400))
  d$y2 <- 0.5 * (d$big - 1e5) + rnorm(40, 0, 0.5)[d$g] + rnorm(400)
  f <- suppressWarnings(suppressMessages(frm(y2 ~ big + (1 | g), data = d)))
  cat(sprintf("big %d obj %a steps %s\n", s, f$opt$objective,
              f$opt$polish[["steps"]] %||% 0))
}
