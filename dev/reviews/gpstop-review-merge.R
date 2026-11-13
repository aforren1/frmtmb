# Review of lane gpstop: does lane nlstop's nl_limit_message() fire on a
# fit that the relinearization improved? nlstop's R/fit-restart.R is
# sourced into an environment under frmtmb's namespace (lane library) and
# called on each fit of both arms. Rows: E wig 120 2 and L wig 120 2
# (dev/review-fitedge-designs.R), the sweep's pow_a gaussian n40 s2 and the
# gate's gp(by =) row.
source("C:/Users/adf44/source/r/frmtmb-wt-gpstop/dev/reviews/gpstop-review-lib.R")
rv_libs("lane")
rv_instrument()
source(file.path(RV_WT, "dev/review-fitedge-designs.R"))
nlx <- new.env(parent = asNamespace("frmtmb"))
sys.source("C:/Users/adf44/source/r/frmtmb-wt-nlstop/R/fit-restart.R",
           envir = nlx)
mk_pow <- function() {
  bodies <- c("exp_a", "exp_b", "logis_a", "mm_V", "mm_K", "pow_a", "pow_b")
  n <- 40L; seed <- 2L
  set.seed(81100 + 1000 * 1 + 100 * seed + n + 10 * match("pow_a", bodies))
  x <- sort(round(runif(n, 0, 10), 3))
  z <- round(rnorm(n), 3)
  w <- round(runif(n, 0.2, 5), 3)
  f <- 2 * sin(3 * x) + x
  fs <- f / max(abs(f))
  y <- (2 + fs) * w^0.5 + 0.1 * rnorm(n)
  data.frame(x, z, w, y)
}
dE <- rv_mk("E", 2, 120L, 2)
dby <- dE
dby$y <- ifelse(dby$fby == "p", dby$y, 0.6 * dby$y)
rows <- list(
  E120w2 = list(rv_spell$E$nl, rv_mk("E", 2, 120L, 2)),
  L120w2 = list(rv_spell$L$nl, rv_mk("L", 2, 120L, 2)),
  pow_a_n40_s2 = list(bf(y ~ a * w^b, a ~ gp(x), b ~ 1, nl = TRUE), mk_pow()),
  by = list(bf(y ~ a * exp(b * z), a ~ gp(x, by = fby), b ~ 1, nl = TRUE),
            dby))
for (nm in names(rows)) for (arm in c("lane", "norelin")) {
  q <- rv_run(nm, arm, rows[[nm]][[1]], data = rows[[nm]][[2]])
  if (inherits(q$fit, "error")) next
  msg <- tryCatch(nlx$nl_limit_message(q$fit),
                  error = function(e) paste("ERROR", conditionMessage(e)))
  cat(sprintf("LIMIT %s %s %s\n", nm, arm,
              if (is.null(msg)) "NULL" else substr(msg, 1, 160)))
  flush.console()
}
cat("DONE\n")
