# Lane optima, item 4b: would a one-evaluation screen find F11 seed 8's
# better optimum? The default fit's objective against the objective
# with one smoothing SD's log moved to -10 (everything else held, the
# inner problem re-solved), per SD; and the same screen on seeds 1:40
# of the same design, to count how often it would fire and gain.
#   Rscript dev/optima-f11-screen.R
.libPaths(c("C:/Users/adf44/source/r/wt-optima-lib",
            "C:/Users/adf44/source/r/rellib-r6",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
src <- readLines("dev/optima-f11.R")
eval(parse(text = src[grep("^gs <- ", src):(grep("^d <- gs", src) - 1L)]))
fo <- bf(y ~ s(x1), sigma ~ s(x2))
rows <- list()
for (seed in 1:40) {
  d <- gs(6, seed)
  f <- suppressWarnings(frm(fo, data = d))
  p <- f$opt$par
  it <- which(names(p) == "theta")
  sc <- vapply(it, function(i) {
    q <- p
    q[i] <- -10
    f$obj$fn(q)
  }, 0)
  fire <- which(sc < f$opt$objective)
  gain <- 0
  if (length(fire)) {
    st <- f$estimates
    th <- st$theta
    th[fire] <- -10
    g <- suppressWarnings(frm(fo, data = d, start = list(theta = th)))
    gain <- as.numeric(logLik(g)) - as.numeric(logLik(f))
  }
  rows[[seed]] <- data.frame(seed = seed, theta = paste(round(p[it], 2),
                                                        collapse = " "),
                             fire = length(fire), gain = gain)
}
X <- do.call(rbind, rows)
print(X)
cat("screen fired on", sum(X$fire > 0), "of", nrow(X), "seeds; gains:",
    paste(signif(X$gain[X$fire > 0], 3), collapse = " "), "\n")
