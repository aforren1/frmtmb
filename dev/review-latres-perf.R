# Attack 6: the performance claims, every post-fit step clocked, peak memory of
# the process. One shape per process. Run alone, three repeats of the fit.
# Usage: Rscript review-latres-perf.R <G> <k> <term> <outfile>
source("C:/Users/adf44/source/r/frmtmb-wt-latres/dev/review-latres-common.R")
args <- commandArgs(TRUE)
G <- as.integer(args[[1]])
k <- as.integer(args[[2]])
term <- args[[3]]
outf <- args[[4]]
set.seed(20261102)
n <- G * k
d <- data.frame(t = rep(seq_len(k), G), g = factor(rep(seq_len(G), each = k)))
d$x <- rnorm(n)
d$ft <- factor(d$t)
e <- numeric(n)
for (i in seq_len(G)) {
  e[(i - 1) * k + seq_len(k)] <- as.numeric(arima.sim(list(ar = .5), k, sd = .5))
}
d$y <- rpois(n, exp(.5 + .3 * d$x + e))
peak <- function() {
  x <- system(sprintf('powershell -NoProfile -Command "(Get-Process -Id %d).PeakWorkingSet64"', Sys.getpid()), intern = TRUE)
  as.numeric(x[1]) / 1e9
}
form <- as.formula(paste("y ~ x +", term))
clock <- function(lab, expr, rep = 1) {
  ts <- numeric(rep)
  val <- NULL
  for (i in seq_len(rep)) {
    t0 <- proc.time()[[3]]
    val <- tryCatch(suppressWarnings(eval.parent(substitute(expr))), error = function(e) structure(conditionMessage(e), class = "err"))
    ts[i] <- proc.time()[[3]] - t0
  }
  r <- data.frame(G = G, k = k, n = n, term = term, step = lab, secs = median(ts), min = min(ts), max = max(ts), peakGB = peak(),
                  note = if (inherits(val, "err")) substr(val, 1, 80) else "")
  write.table(r, outf, sep = "\t", row.names = FALSE, col.names = !file.exists(outf), append = file.exists(outf), quote = FALSE)
  cat(sprintf("%-24s %8.2f s  [%.2f, %.2f]  peak %.2f GB %s\n", lab, r$secs, r$min, r$max, r$peakGB, r$note))
  invisible(val)
}
f <- clock("frm", frm(form, family = poisson(), data = d), rep = 3)
if (inherits(f, "err")) quit("no")
env <- f$obj$env
H <- tryCatch(env$spHess(f$opt$par, random = TRUE), error = function(e) NULL)
if (is.null(H)) H <- tryCatch(env$spHess(env$last.par.best, random = TRUE), error = function(e) NULL)
nnzH <- if (!is.null(H)) length(H@x) else NA
cat("inner dim", length(env$random), " stored Hessian entries", nnzH, " logLik", as.numeric(logLik(f)), " code", f$opt$convergence, "\n")
write.table(data.frame(G = G, k = k, n = n, term = term, step = "INFO", secs = nnzH, min = length(env$random), max = NA, peakGB = NA, note = paste0("ll=", round(as.numeric(logLik(f)), 3), " code=", f$opt$convergence)), outf, sep = "\t", row.names = FALSE, col.names = FALSE, append = TRUE, quote = FALSE)
clock("fn", f$obj$fn(f$opt$par), rep = 5)
clock("gr", f$obj$gr(f$opt$par), rep = 5)
clock("vcov full", vcov(f, full = TRUE))
clock("summary", summary(f))
clock("confint_varcorr", confint_varcorr(f))
if (n <= 2000) {
  clock("fitted", fitted(f))
  clock("frm_linpred se", frm_linpred(f, se.fit = TRUE))
  clock("residuals", residuals(f))
}
clock("frm_linpred", frm_linpred(f))
clock("frm_linpred response", frm_linpred(f, type = "response"))
clock("simulate 1", simulate(f, nsim = 1))
clock("predict 20", predict(f, ndraws = 20))
nd <- d[seq_len(min(n, max(k, 2000))), ]
if (n <= 10000 && nrow(nd) <= 2000) clock("fitted newdata", fitted(f, newdata = nd))
clock("predict newdata 20", predict(f, newdata = nd, ndraws = 20))
clock("simulate newdata 1", simulate(f, nsim = 1, newdata = nd))
clock("predict newdata all 5", predict(f, newdata = d, ndraws = 5))
clock("saveRDS", {tmp <- tempfile(); saveRDS(f, tmp); file.size(tmp)})
cat("rds MB", file.size(tmp) / 1e6, "\n")
clock("update", update(f, . ~ . + I(x^2)))
cat("object.size frame MB", as.numeric(object.size(f$frame)) / 1e6, " aux MB", as.numeric(object.size(f$frame$re_blocks[[1]]$aux)) / 1e6, "\n")
