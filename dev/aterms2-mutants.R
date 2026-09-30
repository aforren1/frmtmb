# Each new code path against a mutant that removes it, to show the test
# that covers it fails without it. A mutant replaces one internal
# function in the installed namespace for the session only.
# Output: dev/aterms2-log-mutants.txt
.libPaths(c("C:/Users/adf44/source/r/wt-aterms2-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
mutate <- function(name, f) {
  old <- get(name, envir = ns)
  environment(f) <- ns
  unlockBinding(name, ns)
  assign(name, f, envir = ns)
  lockBinding(name, ns)
  invisible(old)
}
restore <- function(name, old) {
  unlockBinding(name, ns)
  assign(name, old, envir = ns)
  lockBinding(name, ns)
}
report <- function(label, expr) {
  r <- tryCatch(expr, error = function(e) paste("ERROR:", conditionMessage(e)))
  cat(sprintf("%-58s %s\n", label, paste(format(r), collapse = " ")))
}

# --- rate() in predict(newdata =): has_rate() routes the exposure
set.seed(24)
n <- 200
d <- data.frame(x = rnorm(n), time = runif(n, 0.5, 4))
d$y <- rpois(n, exp(0.2 + 0.5 * d$x) * d$time)
f <- frm(y | rate(time) ~ x, data = d, family = poisson())
nd <- d[1:2, ]
nd$time <- c(1, 1000)
ratio <- function() {
  set.seed(1)
  pr <- predict(f, newdata = nd, ndraws = 400)
  mu <- fitted(f, newdata = nd, dpar = "mu")[, "Estimate"]
  pr[2, "Estimate"] / mu[2]
}
report("rate predict ratio, lane (test wants > 900)", ratio())
old <- mutate("has_rate", function(rspec) FALSE)
report("rate predict ratio, has_rate() -> FALSE", ratio())
report("rate fitted/mu at time 1000, has_rate() -> FALSE",
       fitted(f, newdata = nd)[2, "Estimate"] /
         fitted(f, newdata = nd, dpar = "mu")[2, "Estimate"])
restore("has_rate", old)

# --- subset(): the NA rule
set.seed(25)
n <- 120
ds <- data.frame(x = rnorm(n), z = rnorm(n),
                 s1 = rep(c(TRUE, FALSE), n / 2),
                 s2 = c(rep(TRUE, 90), rep(FALSE, 30)))
ds$y1 <- 1 + ds$x + rnorm(n)
ds$y2 <- rpois(n, exp(0.5 - 0.3 * ds$z))
ds$y2[!ds$s2] <- NA
ds$z[!ds$s2] <- NA
nobs_of <- function() {
  suppressMessages(frm(bf(y1 | subset(s1) ~ x) + gaussian() +
                         bf(y2 | subset(s2) ~ z) + poisson(), data = ds,
                       dry_run = "frame"))$n_obs
}
report("subset n_obs, lane (test wants 120)", nobs_of())
old <- mutate("subset_na_rows", function(spec, mf, resp_cols, exempt) {
  Reduce(`|`, lapply(setdiff(names(mf), exempt), function(cn) {
    is.na(mf[[cn]])
  }), rep(FALSE, nrow(mf)))
})
report("subset n_obs, NA rule -> every NA drops its row", nobs_of())
restore("subset_na_rows", old)

# --- subset(): the per-response rows
old <- mutate("frame_rows", function(mf, rows) mf)
report("subset logLik with frame_rows() a no-op",
       as.numeric(logLik(frm(bf(y1 | subset(s1) ~ x) + gaussian() +
                               bf(y2 | subset(s2) ~ z) + poisson(),
                             data = ds))))
restore("frame_rows", old)

# --- mi(x, idx = ): the row map
set.seed(26)
n <- 120
dm <- data.frame(g1 = sample(seq(1, n - 1, 2), n, TRUE), g2 = seq_len(n),
                 s = rep(c(TRUE, FALSE), n / 2), w = rnorm(n))
dm$x <- rnorm(n)
dm$y <- 1 + 0.5 * dm$x[match(dm$g1, dm$g2)] + 0.3 * dm$w + rnorm(n, sd = 0.5)
bm <- bf(y ~ mi(x, idx = g1) + w) + bf(x | mi() + index(g2) + subset(s) ~ 1)
dm$xm <- dm$x[match(dm$g1, dm$g2)]
ref <- as.numeric(logLik(frm(y ~ xm + w, data = dm, family = gaussian()))) +
  as.numeric(logLik(frm(x ~ 1, data = dm[dm$s, ], family = gaussian())))
report("mi idx logLik minus separate fits, lane",
       as.numeric(logLik(frm(bm, data = dm, family = gaussian()))) - ref)
old <- mutate("mi_idx_rows", function(ent, vn, resp, tgt, mf, index_vals,
                                      sub_rows) NULL)
report("mi idx logLik minus separate fits, idxl -> NULL",
       suppressWarnings(as.numeric(logLik(frm(bm, data = dm,
                                               family = gaussian()))) - ref))
restore("mi_idx_rows", old)
