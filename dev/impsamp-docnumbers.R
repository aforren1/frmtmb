# The numbers the docs quote, on the fixed build. Log:
# dev/impsamp-out/docnumbers.txt
.libPaths(c("C:/Users/adf44/source/r/wt-impsamp-lib",
            "C:/Users/adf44/source/r/rellib-r16",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
sq <- function(expr) suppressWarnings(suppressMessages(expr))
# grid: the rise and the rounds, generated
d <- do.call(rbind, lapply(Sys.glob("dev/impsamp-out/grid-v4-*.tsv"),
                           utils::read.delim))
tol <- pmax(3 * d$mcse, 0.1)
cat("grid v4 fits", nrow(d), "ok", sum(d$status == "ok"), "\n")
cat("rise: min", format(min(d$rise, na.rm = TRUE), digits = 3), "median",
    format(median(d$rise, na.rm = TRUE), digits = 3), "max",
    format(max(d$rise, na.rm = TRUE), digits = 3), "\n")
cat("max rise / tolerance:", format(max(d$rise / tol, na.rm = TRUE), digits = 3),
    "\n")
print(table(rounds = d$rounds, N = d$N))
print(quantile(d$ess_min[d$status == "ok"], c(0, .01, .05, .25, .5)))
# the package's designs
set.seed(11)
ng <- 60L
per <- 8L
g <- rep(seq_len(ng), each = per)
x <- rnorm(ng * per)
b <- MASS::mvrnorm(ng, c(0, 0), matrix(c(1, 0.35, 0.35, 0.49), 2))
dp <- data.frame(y = rbinom(ng * per, 1, plogis(-0.5 + 0.8 * x + b[g, 1] +
                                                  b[g, 2] * x)),
                 x = x, g = factor(g))
fp <- sq(frm(bf(y ~ x + (x | g)) + bernoulli(), data = dp, importance = 1000L))
cat("probe 1000 draws: ess_min", fp$importance$ess_min, "logLik",
    as.numeric(logLik(fp)), "mcse", fp$importance$mcse, "rounds",
    fp$importance$rounds, "\n")
set.seed(601)
ng <- 100L
per <- 4L
g <- factor(rep(seq_len(ng), each = per))
x <- rnorm(ng * per)
u <- rnorm(ng, 0, 1.2)
ds <- data.frame(y = rbinom(ng * per, 1, plogis(-0.5 + 0.7 * x + u[g])),
                 x = x, g = g)
fs <- sq(frm(bf(y ~ x + (1 | g)) + bernoulli(), data = ds, importance = 1000L))
cat("scalar 1000 draws: ess_min", fs$importance$ess_min, "\n")
set.seed(7)
ng <- 40
per <- 5
gg <- rep(seq_len(ng), each = per)
xx <- rnorm(ng * per)
bb <- cbind(rnorm(ng, 0, 1.5), rnorm(ng, 0, 1))
bd <- data.frame(y = rbinom(ng * per, 1, plogis(-0.5 + 0.8 * xx + bb[gg, 1] +
                                                  bb[gg, 2] * xx)),
                 x = xx, g = factor(gg))
fv <- sq(frm(bf(y ~ x + (x | g)) + bernoulli(), data = bd, importance = 500))
cat("vignette design 500: ess_min", fv$importance$ess_min, "mcse",
    fv$importance$mcse, "\n")
for (s in c(1L, 2L, 3L, 7L, 11L, 101L, 601L)) {
  set.seed(s)
  ng <- 40L
  g <- factor(rep(seq_len(ng), each = 3L))
  x <- rnorm(ng * 3L)
  u <- rnorm(ng, 0, 2)
  dh <- data.frame(y = rbinom(ng * 3L, 1, plogis(-0.5 + 0.5 * x + u[g])),
                   x = x, g = g)
  fh <- sq(frm(bf(y ~ x + (1 | g)) + bernoulli(), data = dh, importance = 1000L))
  cat("forty of three, seed", s, "ess_min", fh$importance$ess_min, "\n")
}
s <- readLines("tests/testthat/test-importance.R")
a <- grep("^imp_dpar_data <- function", s)
b <- a + which(s[(a + 1L):length(s)] == "}")[1L]
eval(parse(text = s[a:b]))
dd <- imp_dpar_data()
for (N in c(500L, 2000L)) {
  fd <- sq(frm(bf(y ~ x + (1 + x | g), sigma ~ 1 + (1 | g)), family = gaussian(),
               data = dd, importance = N))
  cat("mu slope + sigma block,", N, "draws: ess_min", fd$importance$ess_min,
      "\n")
}
# cost per draw on the 480-row probe
fr <- fp$frame
lap_fit <- sq(frm(bf(y ~ x + (x | g)) + bernoulli(), data = dp))
tpl <- lap_fit$obj$env$parList(lap_fit$opt$par)
lap <- RTMB::MakeADFun(ns$build_objective(fr), tpl, random = "b",
                       map = fr[["map"]], silent = TRUE)
lay <- ns$imp_layout(fr[["re_blocks"]])
gm <- ns$imp_group_map(fr, lay)
at <- lap_fit$opt$par
clock <- function(f) {
  best <- Inf
  for (r in 1:3) {
    k <- 0L
    t0 <- proc.time()[["elapsed"]]
    repeat {
      f(k)
      k <- k + 1L
      el <- proc.time()[["elapsed"]] - t0
      if (el > 1.2) break
    }
    best <- min(best, el / k)
  }
  best
}
for (N in c(200L, 1000L)) {
  plan <- ns$imp_plan(lap, fr, lay, at, N, 1L)
  io <- ns$build_importance_objective(fr, lay, gm, plan)
  invisible(gc(reset = TRUE))
  m0 <- sum(gc()[, 2L])
  obj <- RTMB::MakeADFun(io$fn, ns$imp_template(tpl, "b"), silent = TRUE)
  tg <- clock(function(k) obj$gr(at + (k %% 2L) * 1e-7))
  cat(sprintf("probe N=%d: gr %.4f s, %.3f ms per draw\n", N, tg, 1e3 * tg / N))
}
tg <- clock(function(k) lap$gr(at + (k %% 2L) * 1e-7))
cat(sprintf("Laplace gr %.4f s\n", tg))
