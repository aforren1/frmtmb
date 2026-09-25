# Punch round 1, item 4: the newdata design is now built once per call
# instead of once per replicate. The draws must not change. Saves
# simulate(newdata = ) on a grid of designs and re_formula values, from
# the build installed in the lane library, and compares two saves.
#   Rscript dev/simnewdata-cache-bitwise.R before
#   Rscript dev/simnewdata-cache-bitwise.R after
#   Rscript dev/simnewdata-cache-bitwise.R compare
arm <- commandArgs(trailingOnly = TRUE)[1]
OUT <- "dev/simnewdata-log"
if (identical(arm, "compare")) {
  a <- readRDS(file.path(OUT, "cache-before.rds"))
  b <- readRDS(file.path(OUT, "cache-after.rds"))
  stopifnot(identical(names(a), names(b)))
  same <- vapply(names(a), function(k) identical(a[[k]], b[[k]]), NA)
  for (k in names(a)) cat(sprintf("%-40s %s\n", k,
                                  if (same[[k]]) "identical" else "DIFFERS"))
  cat("identical", sum(same), "of", length(same), "\n")
  quit(save = "no")
}
source("dev/simnewdata-prelude.R")
suppressMessages(library(frmtmb))
set.seed(101)
n <- 240
d <- data.frame(x = stats::rnorm(n), z = stats::runif(n),
                g = factor(rep(seq_len(12), each = 20)),
                h = factor(rep(seq_len(8), 30)))
u <- stats::rnorm(12, 0, 0.5)[d$g]
eta <- 0.3 + 0.5 * d$x + u
d$yg <- stats::rnorm(n, eta, 0.7)
d$yp <- stats::rpois(n, exp(eta))
d$nt <- sample(5:15, n, replace = TRUE)
d$yb <- stats::rbinom(n, d$nt, stats::plogis(eta))
d$yo <- factor(cut(eta + stats::rlogis(n), c(-Inf, -0.5, 0.5, 1.5, Inf),
                   labels = FALSE), ordered = TRUE)
d$yk <- factor(sample(c("a", "b", "c"), n, replace = TRUE))
d$ys <- sin(3 * d$z) + stats::rnorm(n, 0, 0.3)
d$ym <- ifelse(stats::runif(n) < 0.4, stats::rnorm(n, -2), stats::rnorm(n, 2))
d$ystu <- eta + stats::rt(n, 4)
d$time <- rep(1:20, length.out = n)
designs <- list(
  gaussian = list(bf(yg ~ x + (1 + x | g) + (1 | h)), gaussian()),
  poisson = list(bf(yp ~ x + offset(z) + (1 | g)), poisson()),
  binomial = list(bf(yb | trials(nt) ~ x + (1 | g)), binomial()),
  cumulative_cs = list(bf(yo ~ cs(x)), sratio()),
  categorical = list(bf(yk ~ x), categorical()),
  smooth = list(bf(ys ~ s(z) + (1 | g)), gaussian()),
  distributional = list(bf(yg ~ x + (1 | g), sigma ~ x), gaussian()),
  mixture = list(bf(ym ~ 1), mixture(gaussian(), gaussian())),
  student_ar = list(bf(ystu ~ x + ar(time, gr = g, cov = TRUE)), student())
)
nd <- d[c(1:30, 200:215), ]
nd_new <- nd
nd_new$g <- factor(rep(c("1", "new"), length.out = nrow(nd)))
res <- list()
for (nm in names(designs)) {
  fit <- suppressWarnings(frm(designs[[nm]][[1]], family = designs[[nm]][[2]],
                              data = d))
  go <- function(...) {
    tryCatch(simulate(fit, nsim = 5, seed = 3, ...),
             error = function(e) paste("ERROR:", conditionMessage(e)))
  }
  res[[paste(nm, "NULL")]] <- go(newdata = nd)
  res[[paste(nm, "NA")]] <- go(newdata = nd, re_formula = NA)
  res[[paste(nm, "allow")]] <- go(newdata = nd_new, allow_new_levels = TRUE)
  res[[paste(nm, "NA new")]] <- go(newdata = nd_new, re_formula = NA)
  if (nm == "gaussian") {
    res[[paste(nm, "partial")]] <- go(newdata = nd, re_formula = ~ (1 | g))
  }
}
saveRDS(res, file.path(OUT, paste0("cache-", arm, ".rds")))
cat("saved", length(res), "\n")
