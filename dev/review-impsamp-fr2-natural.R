# Fix round 2, item 3. Usage: REV_ARM=... Rscript ... <mode>
# same: the lane's 7 designs, data seeds 1:12, draws 2 to 64 (420 fits).
# sep: the slope near separation alone, for the arm comparison.
# search: seeds 13:36, draws 2 to 16, and rounds capped at 1 and 2, plus harder designs.
source("C:/Users/adf44/source/r/frmtmb-wt-impsamp/dev/review-impsamp-common.R")
arm <- Sys.getenv("REV_ARM", "after")
set_libs(arm)
suppressPackageStartupMessages(library(frmtmb))
cat("lib", find.package("frmtmb"), "\n")
mode <- commandArgs(TRUE)[[1]]
src <- readLines("C:/Users/adf44/source/r/frmtmb-wt-impsamp/dev/impsamp-natural.R")
i <- grep("^mk <- list\\(", src)
j <- grep("^n <- c\\(fires", src)
mk <- eval(parse(text = sub("^mk <- ", "", paste(src[i:(j - 1L)], collapse = "\n"))))
mk$sep_tight <- function(s) {
  set.seed(s); g <- factor(rep(1:12, each = 3L)); x <- rnorm(36)
  list(d = data.frame(y = as.integer(2 * x + rnorm(12, 0, 3)[g] > 0), x = x, g = g),
       f = bf(y ~ x + (x | g)), fam = bernoulli())
}
mk$cauchy_two <- function(s) {
  set.seed(s); g <- factor(rep(1:15, each = 2L)); x <- rnorm(30)
  list(d = data.frame(y = 0.5 * x + rnorm(15, 0, 2)[g] + rt(30, 1), x = x, g = g),
       f = bf(y ~ x + (1 | g), nu = 2), fam = student())
}
mk$pois_big <- function(s) {
  set.seed(s); g <- factor(rep(1:10, each = 2L)); x <- rnorm(20)
  list(d = data.frame(y = rpois(20, exp(2 + rnorm(10, 0, 3)[g])), x = x, g = g),
       f = bf(y ~ x + (x | g)), fam = poisson())
}
grid <- switch(mode,
  same = expand.grid(nm = names(mk)[1:7], s = 1:12, N = c(2L, 4L, 8L, 16L, 64L), R = 5L, stringsAsFactors = FALSE),
  sep = expand.grid(nm = "bern_slope_sep", s = 1:12, N = c(2L, 4L, 8L, 16L, 64L), R = 5L, stringsAsFactors = FALSE),
  search = rbind(
    expand.grid(nm = names(mk), s = 13:30, N = c(2L, 8L, 32L), R = 5L, stringsAsFactors = FALSE),
    expand.grid(nm = names(mk), s = 1:12, N = c(2L, 8L, 32L), R = c(1L, 2L), stringsAsFactors = FALSE)))
rows <- vector("list", nrow(grid))
for (k in seq_len(nrow(grid))) {
  g <- grid[k, ]
  m <- mk[[g$nm]](g$s)
  w <- character(0)
  r <- withCallingHandlers(
    tryCatch(frm(m$f, family = m$fam, data = m$d, importance = g$N,
                 control = frmtmb_control(importance_rounds = g$R)), error = function(e) e),
    warning = function(x) {
      w <<- c(w, substr(conditionMessage(x), 1, 90))
      invokeRestart("muffleWarning")
    })
  e <- list(design = g$nm, s = g$s, N = g$N, R = g$R, warn = w)
  if (inherits(r, "error")) {
    e$msg <- conditionMessage(r)
  } else {
    e[c("rounds", "moves", "capped", "ess", "mcse", "ll", "conv")] <-
      list(r$importance$rounds, r$importance$moves, r$importance$capped, r$importance$ess_min,
           r$importance$mcse, as.numeric(logLik(r)), r$opt$convergence)
  }
  rows[[k]] <- e
}
fire <- vapply(rows, function(e) !is.null(e$msg) && grepl("did not settle", e$msg, fixed = TRUE), TRUE)
oth <- vapply(rows, function(e) !is.null(e$msg) && !grepl("did not settle", e$msg, fixed = TRUE), TRUE)
cap <- vapply(rows, function(e) isTRUE(e$capped), TRUE)
cat(sprintf("%s %s: fits %d, fires %d, other stops %d, capped %d\n", arm, mode, length(rows), sum(fire), sum(oth), sum(cap)))
tb <- table(vapply(rows, function(e) e$design, ""), ifelse(fire, "fire", ifelse(oth, "other", ifelse(cap, "capped", "ok"))))
print(tb)
for (e in rows[fire]) cat(sprintf("FIRE %s s=%d N=%d R=%d: %s\n", e$design, e$s, e$N, e$R, substr(e$msg, 1, 520)))
for (e in rows[oth][seq_len(min(5, sum(oth)))]) cat(sprintf("OTHER %s s=%d N=%d: %s\n", e$design, e$s, e$N, substr(e$msg, 1, 160)))
if (any(cap)) {
  sp <- vapply(rows[cap], function(e) (max(e$moves) - min(e$moves)) / mean(e$moves), 0)
  cat("capped: spread of moves", format(range(sp), digits = 3), " without a rounds warning",
      sum(vapply(rows[cap], function(e) !any(grepl("rounds", e$warn)), TRUE)), "\n")
}
