# Reviewer: cost of quadrature = TRUE by nesting depth and size. One
# configuration per process; run under dev/fuzzfix-rev-capped.ps1 for the
# peak working set and a memory cap.
#   Rscript dev/fuzzfix-rev-quad-cost.R <lib|base> <single|nested|nested3>
#           <ng> <m> [bypass]
# ng outer levels, 2 inner levels each (2 x 2 for nested3), m rows per
# innermost cell. "bypass" switches off quad_check_depth() in this
# process only, to measure what the refusal prevents.
a <- commandArgs(TRUE)
base <- "C:/Users/adf44/source/r/rellib-r10"
ul <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
.libPaths(if (a[1] == "base") c(base, ul) else c(a[1], base, ul))
suppressMessages(library(frmtmb))
cat("lib:", format(packageVersion("frmtmb")), find.package("frmtmb"), "\n")
st <- a[2]; ng <- as.integer(a[3]); m <- as.integer(a[4])
if (length(a) > 4 && a[5] == "bypass") {
  utils::assignInNamespace("quad_check_depth", function(frame) invisible(),
                           "frmtmb")
}
set.seed(7)
d <- expand.grid(r = seq_len(m), gc = c("p", "q"), gb = c("i", "ii"),
                 ga = factor(seq_len(ng)))
if (st != "nested3") d <- d[d$gc == "p", ]
d$x <- rnorm(nrow(d))
ab <- interaction(d$ga, d$gb, drop = TRUE)
d$y <- rbinom(nrow(d), 1, plogis(0.3 + 0.5 * d$x +
                                   rnorm(ng, 0, 0.8)[d$ga] +
                                   rnorm(nlevels(ab), 0, 0.6)[ab]))
re <- c(single = "(1 | ga)", nested = "(1 | ga/gb)",
        nested3 = "(1 | ga/gb/gc)")[[st]]
fo <- bf(stats::as.formula(paste("y ~ x +", re)))
tm <- function(expr) {
  t0 <- proc.time()[[3]]
  v <- tryCatch(suppressWarnings(suppressMessages(expr)),
                error = function(e) e)
  list(v = v, t = proc.time()[[3]] - t0)
}
lap <- tm(frm(fo, family = bernoulli(), data = d))
qd <- tm(frm(fo, family = bernoulli(), data = d, quadrature = TRUE))
if (inherits(qd$v, "error")) {
  cat(sprintf("COST %s ng %d m %d rows %d laplace %.2fs quad ERROR %s\n", st,
              ng, m, nrow(d), lap$t, conditionMessage(qd$v)))
  quit(status = 0)
}
sm <- tm(summary(qd$v))
cat(sprintf(paste0("COST %s ng %d m %d rows %d nb %d laplace %.2fs ",
                   "quad %.2fs summary %.2fs left_random %d ",
                   "logLik quad %.8f laplace %.8f\n"),
            st, ng, m, nrow(d), length(qd$v$estimates$b), lap$t, qd$t, sm$t,
            length(qd$v$obj$env$random), as.numeric(logLik(qd$v)),
            as.numeric(logLik(lap$v))))
