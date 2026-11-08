# Fix round, B1: the cost of item 1 as shipped (the body start and the second
# descent from 0.79.1's start) against 0.79.1's path in the same build
# (gp_body_start() patched to return NULL, so neither runs), interleaved in
# one process, order alternating by row, on the review's two 96-row sets.
# The rows are dev/fitedge-r2-gpfresh.R's, read from that file. Args: arm.
arm <- commandArgs(TRUE)[1]
source("C:/Users/adf44/source/r/frmtmb-wt-fitedge/dev/fitedge-lib.R")
fe_libs(arm)
ns <- asNamespace("frmtmb")
orig <- get("gp_body_start", envir = ns)
use <- function(on) {
  unlockBinding("gp_body_start", ns)
  assign("gp_body_start", if (on) orig else function(...) NULL, envir = ns)
  lockBinding("gp_body_start", ns)
}
one <- function(on, sp, d) {
  use(on)
  t0 <- proc.time()[[3]]
  q <- fe_fit(sp$nl, data = d, family = sp$fam)
  list(t = proc.time()[[3]] - t0,
       ll = if (inherits(q$fit, "error")) NA else as.numeric(logLik(q$fit)))
}
src <- readLines(file.path(FE_WT, "dev/fitedge-r2-gpfresh.R"))
i <- 0L
for (set in c("fresh", "lane96")) {
  # the set's rows, data and spellings, as the sweep builds them
  cut <- grep("^# theta starts in each design", src)
  env <- new.env()
  env$set <- set
  env$FE_WT <- FE_WT
  env$arm <- "x"
  body_src <- src[grep("^if \\(set == \"fresh\"\\) \\{", src):(cut - 1L)]
  eval(parse(text = body_src), env)
  on_t <- off_t <- 0
  better <- worse <- 0L
  for (k in seq_len(nrow(env$rows))) {
    r <- env$rows[k, ]
    d <- env$getd(r)
    sp <- env$getsp(r)
    i <- i + 1L
    if (i %% 2L) {
      a <- one(TRUE, sp, d); b <- one(FALSE, sp, d)
    } else {
      b <- one(FALSE, sp, d); a <- one(TRUE, sp, d)
    }
    on_t <- on_t + a$t
    off_t <- off_t + b$t
    better <- better + isTRUE(a$ll > b$ll + 1e-3)
    worse <- worse + isTRUE(a$ll < b$ll - 1e-3)
  }
  cat(sprintf("COST2 %s rows %d: shipped %.1f s, 0.79.1 path %.1f s, ratio %.2f; shipped above %d, below %d\n",
              set, nrow(env$rows), on_t, off_t, on_t / off_t, better, worse))
}
cat("DONE\n")
