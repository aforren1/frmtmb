# Prototype of a simplex plateau escape, measured against the exact
# profile maximum of dev/nanse-mo-profile.R.
#
# At a saturated mo() simplex (a weight below `sat`) the softmax's
# gradient vanishes, so the optimizer stops wherever the weight got
# small, whether or not the likelihood is maximal there. The escape
# moves the simplex a step `eps` toward each vertex in turn, in weight
# space, where the move is visible to the likelihood, and reoptimizes
# from the best such point when it improves the objective; repeated
# until no vertex direction improves.
#   Rscript dev/nanse-mo-escape.R [lib] [seeds] [profile.tsv]
args <- commandArgs(trailingOnly = TRUE)
lib <- if (length(args)) args[1] else "C:/Users/adf44/source/r/rellib-r5"
seeds <- if (length(args) > 1) eval(parse(text = args[2])) else 1:200
ptsv <- if (length(args) > 2) args[3] else "dev/nanse-log/mo-profile-base.tsv"
.libPaths(unique(c(lib, "C:/Users/adf44/source/r/rellib-r5",
                   "C:/Users/adf44/AppData/Local/R/win-library/4.6")))
suppressMessages(library(frmtmb))
P <- utils::read.delim(ptsv)
mk <- function(s) {
  set.seed(s)
  lev <- c("below_20", "20_to_40", "40_to_100", "greater_100")
  income <- factor(sample(lev, 100, TRUE), levels = lev, ordered = TRUE)
  ls <- c(30, 60, 70, 75)[income] + rnorm(100, sd = 7)
  d <- data.frame(income, ls)
  d$age <- rnorm(100, mean = 40, sd = 10)
  d
}
escape <- function(obj, par, sat = 1e-6, eps = 0.01, rounds = 10) {
  nm <- names(par)
  blocks <- split(seq_along(par), nm)[grepl("^zeta", unique(nm))]
  blocks <- split(seq_along(par), nm)
  blocks <- blocks[grepl("^zeta", names(blocks))]
  best <- par
  fbest <- obj$fn(par)
  n_re <- 0L
  for (r in seq_len(rounds)) {
    cand <- NULL
    fc <- fbest
    for (b in blocks) {
      w <- exp(c(0, best[b]))
      w <- w / sum(w)
      if (min(w) > sat) next
      for (j in seq_along(w)) {
        e <- numeric(length(w))
        e[j] <- 1
        w2 <- (1 - eps) * w + eps * e
        p <- best
        p[b] <- log(w2[-1] / w2[1])
        f2 <- obj$fn(p)
        if (is.finite(f2) && f2 < fc - 1e-8) {
          fc <- f2
          cand <- p
        }
      }
    }
    if (is.null(cand)) break
    o <- stats::nlminb(cand, obj$fn, obj$gr,
                       control = list(eval.max = 1000, iter.max = 1000))
    n_re <- n_re + 1L
    if (o$objective < fbest - 1e-10) {
      best <- stats::setNames(o$par, nm)
      fbest <- o$objective
    } else break
  }
  list(par = best, objective = fbest, rounds = n_re)
}
rows <- list()
for (s in seeds) {
  f <- suppressWarnings(frm(ls ~ mo(income) * age, data = mk(s)))
  t0 <- proc.time()[[3]]
  e <- escape(f$obj, f$opt$par)
  dt <- proc.time()[[3]] - t0
  ll_prof <- P$ll_profile[P$seed == s]
  rows[[length(rows) + 1L]] <- data.frame(
    seed = s, ll_fit = -f$opt$objective, ll_escape = -e$objective,
    ll_profile = ll_prof, gap_fit = ll_prof + f$opt$objective,
    gap_escape = ll_prof + e$objective, rounds = e$rounds, sec = dt)
}
X <- do.call(rbind, rows)
utils::write.table(X, "dev/nanse-log/mo-escape.tsv", sep = "\t",
                   quote = FALSE, row.names = FALSE)
cat("seeds", nrow(X), "\n")
for (v in c("gap_fit", "gap_escape")) {
  cat(v, ": > 1e-3", sum(X[[v]] > 1e-3), " > 1e-2", sum(X[[v]] > 1e-2),
      " > 0.1", sum(X[[v]] > 0.1), " max", format(max(X[[v]]), digits = 3),
      " min", format(min(X[[v]]), digits = 3), "\n")
}
cat("escape rounds used:\n")
print(table(X$rounds))
cat("seconds per escape: median", median(X$sec), "max", max(X$sec), "\n")
print(X[X$gap_escape > 1e-3, ])
