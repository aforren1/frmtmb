# Lane ordstart, item 4: how many starts. The rule of mixture_starts()
# (a candidate with a step-function component is disqualified; a proper
# default is replaced only by a gain above 1e-3; a step default by the
# best step-free candidate) applied to the default plus the first n
# jittered runs of dev/ordstart-starts.R's 40-start output, for several
# n. Per n: the seeds whose choice has the logLik of the n = 40 choice
# (within 1e-3), the seeds whose answer moved from the default, and the
# seeds left with a step-function answer.
#   Rscript dev/ordstart-count.R file.tsv [file.tsv ...]
args <- commandArgs(trailingOnly = TRUE)
X <- do.call(rbind, lapply(args, utils::read.delim))
X$step <- pmax(X$sharpen1, X$sharpen2, na.rm = TRUE) > -0.1
choose <- function(S) {
  S <- S[!is.na(S$ll), ]
  d0 <- S[S$run == 0, ]
  ok <- S[!S$step, ]
  if (!nrow(ok)) return(d0)
  j <- ok[which.max(ok$ll), ]
  if (d0$step || j$ll - d0$ll > 1e-3) j else d0
}
keys <- unique(paste(X$family, X$seed))
ref <- vapply(keys, function(k) {
  choose(X[paste(X$family, X$seed) == k, ])$ll
}, 0)
for (n in c(0, 1, 2, 3, 5, 8, 10, 15, 20, 30, 40)) {
  ch <- lapply(keys, function(k) {
    S <- X[paste(X$family, X$seed) == k & X$run <= n, ]
    choose(S)
  })
  ll <- vapply(ch, `[[`, 0, "ll")
  ll0 <- vapply(keys, function(k) {
    X$ll[paste(X$family, X$seed) == k & X$run == 0]
  }, 0)
  cat(sprintf(paste("n %2d: at the n = 40 choice %2d of %d; moved from the",
                    "default %2d; step answer %2d; mean shortfall %.4f\n"),
              n, sum(abs(ll - ref) <= 1e-3), length(keys),
              sum(abs(ll - ll0) > 1e-3),
              sum(vapply(ch, `[[`, NA, "step")), mean(ref - ll)))
}
