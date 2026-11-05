# What the rewritten tests can assert, on the fixed build (scratch log
# dev/impsamp-out/tune.txt): the ESS along a far profile bound; a real
# fit that warns on ESS; a real fit capped at one round; the numDeriv
# gap on the sigma-block design.
.libPaths(c("C:/Users/adf44/source/r/wt-impsamp-lib",
            "C:/Users/adf44/source/r/rellib-r16",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
source("dev/impsamp-agh.R")
# profile design
set.seed(11)
ng <- 40L
per <- 4L
g <- factor(rep(seq_len(ng), each = per))
x <- rnorm(ng * per)
u <- rnorm(ng, 0, 1.1)
dd <- data.frame(y = rbinom(ng * per, 1, plogis(-0.4 + 0.7 * x + u[g])),
                 x = x, g = g)
fi <- frm(bf(y ~ x + (1 | g)) + bernoulli(), data = dd, importance = 500L)
prop <- ns$imp_frozen_proposal(fi)
cat("ess at fit", fi$importance$ess_min, "from prop",
    min(ns$imp_ess_at(prop, fi$opt$par)$ess), "\n")
i <- which(names(fi$opt$par) == "theta")
for (d in c(-3, -2, -1, -0.5, 0.5, 1, 2, 3)) {
  at <- fi$opt$par
  at[i] <- at[i] + d
  cat(sprintf("theta_1 %+.1f: min ESS %.3f\n", d,
              min(ns$imp_ess_at(prop, at)$ess)))
}
ci <- withCallingHandlers(confint(fi, parm = "theta_1", method = "profile"),
                          warning = function(w) {
                            cat("WARN:", conditionMessage(w), "\n")
                            invokeRestart("muffleWarning")
                          })
print(ci)
cat("x bound:\n")
print(confint(fi, parm = "x", method = "profile"))
qd <- frm(bf(y ~ x + (1 | g)) + bernoulli(), data = dd, quadrature = TRUE)
print(confint(qd, parm = "theta_1", method = "profile"))
# a real ESS warning
set.seed(5)
ng <- 30L
g <- factor(rep(seq_len(ng), each = 2L))
x <- rnorm(2L * ng)
du <- data.frame(y = rbinom(2L * ng, 1, 0.4), x = x, g = g)
r <- withCallingHandlers(frm(bf(y ~ x + (x | g)) + bernoulli(), data = du,
                             importance = 500L),
                         warning = function(w) {
                           cat("WARN:", substr(conditionMessage(w), 1, 200), "\n")
                           invokeRestart("muffleWarning")
                         })
cat("ess_min", r$importance$ess_min, "\n")
print(sort(r$importance$ess)[1:3])
# capped at one round
sc <- function() {
  set.seed(601)
  ng <- 100L
  per <- 4L
  g <- factor(rep(seq_len(ng), each = per))
  x <- rnorm(ng * per)
  u <- rnorm(ng, 0, 1.2)
  data.frame(y = rbinom(ng * per, 1, plogis(-0.5 + 0.7 * x + u[g])), x = x,
             g = g)
}
for (rr in 1:3) {
  w <- character(0)
  f <- withCallingHandlers(frm(bf(y ~ x + (1 | g)) + bernoulli(), data = sc(),
                               importance = 200L,
                               control = frmtmb_control(importance_rounds = rr)),
                           warning = function(e) {
                             w <<- c(w, conditionMessage(e))
                             invokeRestart("muffleWarning")
                           })
  cat("rounds cap", rr, "capped", f$importance$capped, "moves",
      format(f$importance$moves, digits = 3), "grads",
      format(f$importance$grads, digits = 3), "|", substr(w, 1, 90), "\n")
}
# numDeriv on the sigma block design
set.seed(21)
ng <- 12L
per <- 10L
g <- factor(rep(seq_len(ng), each = per))
x <- rnorm(ng * per)
d <- data.frame(y = rnorm(ng * per, 1 + 0.5 * x + rnorm(ng)[g], 0.7), x = x, g = g)
s <- readLines("tests/testthat/test-importance.R")
a <- grep("^imp_dpar_data <- function", s)
b <- a + which(s[(a + 1L):length(s)] == "}")[1L]
eval(parse(text = s[a:b]))
a <- grep("^imp_parts <- function", s)
b <- a + which(s[(a + 1L):length(s)] == "}")[1L]
eval(parse(text = s[a:b]))
dd <- imp_dpar_data()
p <- imp_parts(bf(y ~ x + (1 + x | g), sigma ~ 1 + (1 | g)) + gaussian(), dd)
plan <- ns$imp_plan(p$lap, p$frame, p$lay, p$opt$par, 500L, 1L)
io <- ns$build_importance_objective(p$frame, p$lay, p$gmap, plan)
otpl <- ns$imp_template(p$tpl, "b")
obj <- RTMB::MakeADFun(io$fn, otpl, silent = TRUE)
for (shift in c(0, 0.3)) {
  at <- p$opt$par + shift
  nd1 <- numDeriv::grad(obj$fn, at)
  nd2 <- numDeriv::grad(obj$fn, at, method.args = list(eps = 1e-5, d = 1e-3))
  cd <- sapply(seq_along(at), function(j) {
    e <- replace(numeric(length(at)), j, 1e-5)
    (obj$fn(at + e) - obj$fn(at - e)) / 2e-5
  })
  ad <- as.numeric(obj$gr(at))
  cat("shift", shift, "\n ad   ", format(ad, digits = 7), "\n nd   ",
      format(nd1, digits = 7), "\n nd2  ", format(nd2, digits = 7),
      "\n cd   ", format(cd, digits = 7), "\n max rel ad-nd",
      format(max(abs(ad - nd1) / pmax(abs(ad), 1)), digits = 3),
      " ad-cd", format(max(abs(ad - cd) / pmax(abs(ad), 1)), digits = 3),
      " nd-nd2", format(max(abs(nd1 - nd2) / pmax(abs(ad), 1)), digits = 3), "\n")
}
