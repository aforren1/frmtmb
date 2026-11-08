# Review of lane latperf, attacks 2 and 3: time, peak memory and saved size of
# the readers of design "two" (poisson, ar(cov = TRUE) on one series, plus
# (1 | id) with 50 levels; the generator of dev/latperf-prof.R, seed 1), and
# what the saved fit holds, part by part.
# Usage: Rscript dev/review-latperf-scale.R <n> [nosel]
#   nosel  TMB's solveSubset() mocked away: every fitted row is solved
#          directly (the fallback the lane names)
args <- commandArgs(TRUE)
n <- as.integer(args[1])
nosel <- "nosel" %in% args
.libPaths(c("C:/Users/adf44/source/r/wt-latperf-lib",
            "C:/Users/adf44/source/r/rellib-r18",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
cat(sprintf("lib: %s frmtmb %s\n", dirname(find.package("frmtmb")),
            format(packageVersion("frmtmb"))))
cat(sprintf("n %d nosel %s\n", n, nosel))
if (nosel) {
  unlockBinding("jcs_selinv", ns)
  assign("jcs_selinv", function(S) NULL, envir = ns)
}
peak_ws <- function() {
  out <- system2("powershell", c("-NoProfile", "-Command",
    sprintf("(Get-Process -Id %d).PeakWorkingSet64", Sys.getpid())),
    stdout = TRUE)
  as.numeric(out[length(out)]) / 2^20
}
heap_mb <- function() sum(gc()[, 6])
mb <- function(x) length(serialize(x, NULL)) / 2^20
set.seed(1)
d <- data.frame(t = seq_len(n), g = factor(1), x = rnorm(n),
                id = factor(sample.int(50, n, TRUE)))
e <- as.numeric(arima.sim(list(ar = .5), n, sd = .5))
u <- rnorm(50, 0, 0.3)
d$y <- rpois(n, exp(.5 + .3 * d$x + e + u[d$id]))
invisible(gc(reset = TRUE))
tf <- system.time(fit <- suppressWarnings(
  frm(y ~ x + (1 | id) + ar(t, g, cov = TRUE), family = poisson(),
      data = d)))[["elapsed"]]
cat(sprintf("RESULT fit time %.2f heap %.0f ws %.0f saved %.2f\n", tf,
            heap_mb(), peak_ws(), mb(fit)))
ws0 <- peak_ws()
invisible(gc(reset = TRUE))
tt <- system.time(v <- suppressWarnings(fitted(fit)))[["elapsed"]]
cat(sprintf(paste("RESULT fitted time %.2f heap %.0f ws %.0f (ws before %.0f)",
                  "se_na %d se_first %.10g\n"), tt, heap_mb(), peak_ws(), ws0,
            sum(is.na(v[, "Est.Error"])), v[1, "Est.Error"]))
invisible(gc(reset = TRUE))
tl <- system.time(l <- frm_linpred(fit, se.fit = TRUE))[["elapsed"]]
cat(sprintf("RESULT linpred time %.2f heap %.0f ws %.0f\n", tl, heap_mb(),
            peak_ws()))
S <- fit$cache$Vjoint$V
cat(sprintf(paste("RESULT saved total %.2f | Vjoint %.2f: L %.2f J %.2f",
                  "Vf %.3f keys %.2f | Qjoint %.2f | sdr %.2f\n"),
            mb(fit), mb(fit$cache$Vjoint), mb(S$L), mb(S$J), mb(S$Vf),
            if (is.list(S$memo$keys)) mb(S$memo$keys) else 0,
            mb(fit$cache$Qjoint), mb(frmtmb:::sdr_of(fit))))
cat("cache:", paste(ls(fit$cache), collapse = ","), "\n")
# the lookup is a memo: dropped and rebuilt, what does it cost
if (!nosel) {
  rm("keys", envir = S$memo)
  cat(sprintf("RESULT saved without keys %.2f\n", mb(fit)))
  tk <- system.time(frmtmb:::jcs_selinv_keys(S))[["elapsed"]]
  cat(sprintf("RESULT rebuild keys time %.2f\n", tk))
  # a reloaded fit, its cache as saved
  f2 <- unserialize(serialize(fit, NULL))
  t2 <- system.time(v2 <- suppressWarnings(fitted(f2)))[["elapsed"]]
  cat(sprintf("RESULT reloaded fitted time %.2f identical %s\n", t2,
              identical(v2, v)))
}
cat(sprintf("RESULT end ws %.0f\n", peak_ws()))
