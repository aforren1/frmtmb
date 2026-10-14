# Reviewer of lane illse, attack 7: the mo() split direction. Seeds 15
# and 19 of `ls ~ mo(inc)` (the 0.70.0 consolidation's generator, an
# unobserved interior category), and a guard-absent design where every
# category is observed. Per fit: code, warnings, zeta SEs, the split
# directions se_mo_null() adds (lane), and the objective along them.
# Saves fixef/vcov/fitted for the bitwise comparison of the guard-absent
# fits.
#   REV_ARM=lane|base Rscript dev/illse-rev-mo.R <out.rds>
source("C:/Users/adf44/source/r/frmtmb-wt-illse/dev/illse-rev-mp.R")
out <- commandArgs(TRUE)[1]
gen_gap <- function(seed) {
  set.seed(seed)
  inc <- sample(c(0L, 1L, 3L), 100, TRUE)
  ls <- c(30, 60, 70, 75)[inc + 1L] + stats::rnorm(100, sd = 7)
  data.frame(inc, ls, age = stats::rnorm(100, mean = 40, sd = 10))
}
gen_full <- function(seed) {
  set.seed(seed)
  inc <- sample(0:3, 100, TRUE)
  ls <- c(30, 60, 70, 75)[inc + 1L] + stats::rnorm(100, sd = 7)
  data.frame(inc, ls, age = stats::rnorm(100, mean = 40, sd = 10))
}
res <- list()
run <- function(tag, fo, d) {
  r <- cap(frm(fo, data = d))
  f <- r$v
  sdr <- suppressWarnings(ns$sdr_of(f))
  ss <- summary(sdr, "fixed")
  z <- ss[grepl("^zeta", rownames(ss)), 2]
  gaps <- if (exists("se_mo_gaps", ns)) length(ns$se_mo_gaps(f)) else NA
  flat <- ""
  if (exists("se_mo_null", ns) && isTRUE(gaps > 0)) {
    N <- ns$se_mo_null(f)
    p <- f$opt$par
    f0 <- as.numeric(f$obj$fn(p))
    for (j in seq_len(ncol(N))) {
      dj <- N[, j] / sqrt(sum(N[, j]^2))
      dv <- vapply(c(1e-3, 1e-2, 1e-1), function(t)
        as.numeric(f$obj$fn(p + t * dj)) - f0, 0)
      flat <- paste0(flat, sprintf(" dir%d: f(p+t d)-f(p) at t=1e-3,1e-2,1e-1: %s", j,
                                   paste(sprintf("%.2e", dv), collapse = " ")))
    }
    f$obj$fn(p)
  }
  cat(sprintf("%-22s code %d | gaps %s | lost %s | zeta SE %s | warn %s\n%s\n",
              tag, f$opt$convergence, gaps,
              paste(names(sdr$se_lost), collapse = ","),
              paste(sprintf("%.4g", z), collapse = " "),
              if (length(r$w)) substr(paste(r$w, collapse = " / "), 1, 160) else "-",
              flat))
  res[[tag]] <<- list(code = f$opt$convergence, par = f$opt$par,
                      vcov = suppressWarnings(vcov(f, full = TRUE)),
                      fitted = suppressWarnings(fitted(f)[, "Est.Error"]),
                      warn = r$w)
}
for (s in c(15, 19)) run(sprintf("gap mo(inc) s%d", s), ls ~ mo(inc), gen_gap(s))
run("gap mo(inc)+age s19", ls ~ mo(inc) + age, gen_gap(19))
run("gap mo(inc)+age s26", ls ~ mo(inc) + age, gen_gap(26))
for (s in 1:8) run(sprintf("full mo(inc) s%d", s), ls ~ mo(inc) + age, gen_full(s))
saveRDS(res, out)
