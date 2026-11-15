# Review of lane quadseam2, attack 2 and 3 cost: fit and vcov seconds,
# median of 3, on quadrature fits the lane's new paths may touch. On the
# lane, the seconds inside quad_sd_up_escape() and se_edge_limit() are
# measured by wrapping them.
#   Rscript quadseam2-review-time.R <base|lib>
a <- commandArgs(TRUE)
if (a[1] == "lib") a[1] <- "C:/Users/adf44/source/r/wt-quadseam2-lib"
.libPaths(c(if (a[1] != "base") a[1], "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("lib:", dirname(find.package("frmtmb")), "\n")
ns <- asNamespace("frmtmb")
acc <- new.env()
acc$t <- c(esc = 0, lim = 0)
wrap <- function(nm, key) {
  if (!exists(nm, ns)) return()
  f0 <- get(nm, ns)
  unlockBinding(nm, ns)
  assign(nm, function(...) {
    t0 <- proc.time()[["elapsed"]]
    on.exit(acc$t[[key]] <- acc$t[[key]] + proc.time()[["elapsed"]] - t0)
    f0(...)
  }, envir = ns)
}
wrap("quad_sd_up_escape", "esc")
wrap("se_edge_limit", "lim")
sn_draw <- function(n, alpha) {
  dl <- alpha / sqrt(1 + alpha^2)
  dl * abs(rnorm(n)) + sqrt(1 - dl^2) * rnorm(n)
}
mk <- list(
  bernsep3 = function() {
    set.seed(3)
    d <- data.frame(g = factor(rep(1:20, each = 9)),
                    f = factor(rep(c("a", "b", "c"), 60)), x = rnorm(180))
    d$y <- rbinom(180, 1, plogis(0.2 + 0.5 * d$x + rnorm(20, 0, 0.7)[d$g]))
    d$y[d$f == "c"] <- 1L
    list(fo = y ~ f + x + (1 | g), fam = bernoulli(), d = d)
  },
  poisrun = function() {
    set.seed(5)
    d <- data.frame(g = factor(rep(1:20, each = 9)),
                    f = factor(rep(c("a", "b", "c"), 60)), x = rnorm(180))
    d$y <- rpois(180, exp(0.5 + 0.3 * d$x + rnorm(20, 0, 0.5)[d$g]))
    d$y[d$f == "c"] <- 0L
    list(fo = y ~ f + x + (1 | g), fam = poisson(), d = d)
  },
  poisok = function() {
    set.seed(5)
    d <- data.frame(g = factor(rep(1:20, each = 9)),
                    f = factor(rep(c("a", "b", "c"), 60)), x = rnorm(180))
    d$y <- rpois(180, exp(0.5 + 0.3 * d$x + rnorm(20, 0, 0.5)[d$g]))
    list(fo = y ~ f + x + (1 | g), fam = poisson(), d = d)
  },
  bernsd0 = function() {
    set.seed(1)
    d <- data.frame(x = rnorm(200), g = factor(rep(1:20, length.out = 200)))
    d$y <- rbinom(200, 1, plogis(0.2 + 0.6 * d$x))
    list(fo = y ~ x + (1 | g), fam = bernoulli(), d = d)
  },
  snquad = function() {
    set.seed(1)
    d <- data.frame(x = rnorm(200), g = factor(rep(1:20, length.out = 200)))
    d$y <- 1 + 0.5 * d$x + rnorm(20, 0, 0.6)[d$g] + sn_draw(200, 4)
    list(fo = y ~ x + (1 | g), fam = skew_normal(), d = d)
  }
)
for (nm in names(mk)) {
  cs <- mk[[nm]]()
  tf <- tv <- te <- tl <- numeric(3)
  for (r in 1:3) {
    acc$t[] <- 0
    t0 <- proc.time()[["elapsed"]]
    f <- suppressWarnings(suppressMessages(
      frm(cs$fo, family = cs$fam, data = cs$d, quadrature = TRUE)))
    t1 <- proc.time()[["elapsed"]]
    v <- suppressWarnings(suppressMessages(vcov(f)))
    t2 <- proc.time()[["elapsed"]]
    tf[r] <- t1 - t0; tv[r] <- t2 - t1
    te[r] <- acc$t[["esc"]]; tl[r] <- acc$t[["lim"]]
  }
  cat(sprintf("TIME %-9s fit %.2f vcov %.2f | in sd escape %.2f in limit rule %.2f | sd %.3g ll %.9f\n",
              nm, median(tf), median(tv), median(te), median(tl),
              exp(f$estimates$theta[1]), as.numeric(logLik(f))))
}
cat("END\n")
