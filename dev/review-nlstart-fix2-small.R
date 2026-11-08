# Review of lane nlstart, second short check: the two-stage start sizes.
# (a) `b1 * exp(b2 * xs)` and `b1 * (1 - exp(-b2 * xs))` with xs = x * 1e-6;
# (b) a coefficient inside the curve and multiplying it: `b1 * exp(-b1 * xp)`
#     and `b1 * (1 - exp(-b1 * b2 * x))`.
# Default fit against the best of 16 random starts and the truth, se = FALSE.
# Usage: REV_ARM=after|before Rscript dev/review-nlstart-fix2-small.R [ref]
do_ref <- isTRUE(commandArgs(trailingOnly = TRUE)[1] == "ref")
source("dev/review-nlstart-common.R")
F <- list(gaussian = gaussian(), poisson = poisson(), bernoulli = bernoulli(),
          cumulative = cumulative())
B <- list(
  exp6 = list(f = "b1 * exp(b2 * xs)", p = "b1 + b2", tr = c(1.5, 6e5),
              sc = c(2, 6e5), lat = function(d) 1.5 * exp(0.6 * d$x)),
  sat6 = list(f = "b1 * (1 - exp(-b2 * xs))", p = "b1 + b2", tr = c(-1.5, -6e5),
              sc = c(2, 6e5), lat = function(d) 1.5 * exp(0.6 * d$x)),
  self = list(f = "b1 * exp(-b1 * xp)", p = "b1", tr = 2, sc = 3,
              lat = function(d) 2 * exp(-2 * d$xp)),
  selfprod = list(f = "b1 * (1 - exp(-b1 * b2 * x))", p = "b1 + b2",
                  tr = c(-1.5, 0.4), sc = c(2, 1),
                  lat = function(d) 1.5 * exp(0.6 * d$x)))
cat("ROWHEAD body fam seed ll conv nwarn multi ref est\n")
for (bn in names(B)) for (fam in names(F)) for (seed in 1:5) {
  set.seed(seed * 37 + 11)
  n <- 300
  d <- data.frame(x = round(rnorm(n), 3), xp = round(runif(n, 0, 2), 3))
  d$xs <- d$x * 1e-6
  lat <- B[[bn]]$lat(d)
  # a body on a link scale: centre a positive mean for the counts
  d$y <- switch(fam, gaussian = lat + rnorm(n, 0, 0.3),
                poisson = rpois(n, pmax(lat, 0.05)),
                bernoulli = rbinom(n, 1, plogis(lat - 1)),
                cumulative = ordered(cut(lat + rlogis(n),
                                         c(-Inf, 0.3, 1, 2, Inf), labels = 1:4)))
  fm <- bf(as.formula(paste("y ~", B[[bn]]$f)),
           as.formula(paste(B[[bn]]$p, "~ 1")), nl = TRUE)
  r <- tt(frm(fm, family = F[[fam]], data = d))
  best <- NA
  if (do_ref) {
    nm <- paste0(strsplit(B[[bn]]$p, " \\+ ")[[1]], "_(Intercept)")
    set.seed(seed)
    sts <- c(lapply(1:16, function(i) rnorm(length(nm)) * B[[bn]]$sc),
             list(B[[bn]]$tr))
    best <- max(vapply(sts, function(st) llv(tt(frm(fm, family = F[[fam]],
      data = d, se = FALSE, start = list(beta = setNames(st, nm))))), 0),
      na.rm = TRUE)
  }
  cat(sprintf("ROW %s %s %d %.6f %s %d %d %.6f %s\n", bn, fam, seed, llv(r),
              if (iserr(r)) NA else r$opt$convergence, length(attr(r, "w")),
              as.integer(any(attr(r, "mc") == "frmtmb_nl_multistart")), best,
              if (iserr(r)) gsub("\\s+", "_", substr(r$msg, 1, 50)) else
                paste(format(r$estimates$beta, digits = 4), collapse = ",")))
  flush.console()
}
cat("DONE\n")
