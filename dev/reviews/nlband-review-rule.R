# Review of lane nlband, check 1 (rule level, no fits): how often
# nl_pole_holes() flags a hole on designs the lane did not simulate.
# Walls are those of a / (x + b) (w = -x) unless the design says 1/x
# (a / (1 + b * x), walls b = -1 / x). Seed 20261117, 200 draws a design.
# Rscript dev/reviews/nlband-review-rule.R
source("dev/nlband-common.R")
holes <- ns$nl_pole_holes
set.seed(20261117)
R <- 200
walls_of <- function(x, kind = "x") {
  w <- if (kind == "x") -x else -1 / x[x != 0]
  sort(unique(w))
}
rate <- function(tag, gen, kind = "x") {
  nh <- vapply(seq_len(R), function(i) NROW(holes(walls_of(gen(), kind))),
               1)
  cat(sprintf("RULE %-34s | any %.3f mean %.2f max %d\n", tag,
              mean(nh > 0), mean(nh), max(nh)))
}
# designs with no hole a pole could need
rate("dose5 x8 reps", function() rep(c(0, 1, 3, 10, 30), each = 8))
rate("dose12 x30 reps", function() rep(1:12, each = 30))
rate("dose2fold 24 lev x4", function() rep(0.001 * 2^(0:23), each = 4))
rate("dose2fold 30 lev x4", function() rep(0.001 * 2^(0:29), each = 4))
rate("dose3fold 25 lev x4", function() rep(0.001 * 3^(0:24), each = 4))
rate("dose3fold 30 lev x4", function() rep(0.001 * 3^(0:29), each = 4))
rate("grid 1:60 x3", function() rep(1:60, each = 3))
rate("grid 1:60 jitter 0.01", function() rep(1:60, each = 3) +
       runif(180, -0.01, 0.01))
rate("logunif 6 dec n400", function() 10^runif(400, -3, 3))
rate("logunif 6 dec n2000", function() 10^runif(2000, -3, 3))
rate("lnorm sdlog 2 n400", function() rlnorm(400, 0, 2))
rate("U(0,4) n400 walls -1/x", function() runif(400, 0, 4), "inv")
rate("U(0.1,5) n400 walls -1/x", function() runif(400, 0.1, 5), "inv")
rate("U(-2,2) n400 walls -1/x", function() runif(400, -2, 2), "inv")
rate("visits 10 x40 jit 0.2", function() rep(1:10, each = 40) +
       runif(400, -0.2, 0.2))
rate("clust3 n400", function() c(runif(133, 0, 1), runif(133, 2, 3),
                                runif(134, 4, 5)))
rate("rounded 0.1 U(0,4) n400", function() round(runif(400, 0, 4), 1))
rate("rounded 0.05 U(0,4) n400 hole", function() {
  x <- round(runif(400, 0, 4), 2); x[abs(x - 2) > 0.3]
})
# detection of the band hole (|x - 2| < 0.3 removed, n 400) with k rows
# put back inside it
for (k in 0:5) {
  rate(sprintf("band hole, %d rows inside", k), function() {
    x <- runif(400, 0, 4)
    c(x[abs(x - 2) > 0.3], runif(k, 1.7, 2.3))
  })
}
# a hole near the edge: k rows on one side (pole 0.25, hole 0.1 to 0.4)
for (k in c(3, 6, 9, 10, 12, 20)) {
  rate(sprintf("edge hole, %d rows below", k), function() {
    c(runif(k, 0, 0.1), runif(400 - k, 0.4, 4))
  })
}
# the hole of half width 0.3 at n 100 and 50 (fewer rows)
for (n in c(50, 100, 200)) {
  rate(sprintf("band hole n %d", n), function() {
    x <- runif(n, 0, 4); x[abs(x - 2) > 0.3]
  })
}
# visits: which gap does the rule take? truth pole between visits 5 and 6
pick <- vapply(seq_len(R), function(i) {
  x <- rep(1:10, each = 40) + runif(400, -0.2, 0.2)
  h <- holes(walls_of(x))
  if (is.null(h)) return(NA_real_)
  any(h[, "lo"] < -5.5 & h[, "hi"] > -5.5)
}, NA)
cat(sprintf("VISITS hole 5|6 among the starts on %.3f of draws (NA %d)\n",
            mean(pick, na.rm = TRUE), sum(is.na(pick))))
cat("DONE\n")
