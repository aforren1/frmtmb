# Reviewer: cost of nested quadrature for a censored Beta response on the
# trial merge (cdfs adds the Beta CDF). One size per process; run under
# dev/fuzzfix-rev-capped.ps1.
#   Rscript dev/fuzzfix-rev-betacens.R <lib> <ng> <nested|single>
a <- commandArgs(TRUE)
.libPaths(c(a[1], "C:/Users/adf44/source/r/rellib-r10",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("lib:", format(packageVersion("frmtmb")), find.package("frmtmb"), "\n")
ng <- as.integer(a[2]); m <- 3
set.seed(5)
d <- expand.grid(r = seq_len(m), gb = factor(c("i", "ii")),
                 ga = factor(seq_len(ng)))
d$x <- rnorm(nrow(d))
ab <- interaction(d$ga, d$gb, drop = TRUE)
mu <- plogis(0.2 + 0.4 * d$x + rnorm(ng, 0, 0.5)[d$ga] +
               rnorm(nlevels(ab), 0, 0.4)[ab])
y <- rbeta(nrow(d), mu * 8, (1 - mu) * 8)
d$cc <- ifelse(y > 0.75, "right", "none")
d$y <- pmin(y, 0.75)
re <- if (a[3] == "nested") "(1 | ga/gb)" else "(1 | ga)"
fo <- bf(stats::as.formula(paste("y | cens(cc) ~ x +", re)))
t0 <- proc.time()[[3]]
f <- tryCatch(suppressWarnings(suppressMessages(
  frm(fo, family = Beta(), data = d, quadrature = TRUE))),
  error = function(e) e)
cat(sprintf("BETACENS %s ng %d rows %d censored %d time %.1fs %s\n", a[3], ng,
            nrow(d), sum(d$cc == "right"), proc.time()[[3]] - t0,
            if (inherits(f, "error")) paste("ERROR", conditionMessage(f))
            else sprintf("logLik %.6f left_random %d", as.numeric(logLik(f)),
                         length(f$obj$env$random))))
