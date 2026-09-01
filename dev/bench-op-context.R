## Context for the optimParallel question: where does InstEval time go,
## and what do the existing levers (profile, sparse_x, optimizer choice)
## do? Plus lme4::lmer on the same models.
LIB <- "C:/Users/adf44/AppData/Local/Temp/1/claude/c--Users-adf44-source-r-frmtmb/529b6e73-d28f-46aa-a279-7dbeeb58fd4f/scratchpad/lib"
.libPaths(c(LIB, .libPaths()))
suppressPackageStartupMessages({library(frmtmb); library(lme4)})
d <- lme4::InstEval
REPS <- as.integer(Sys.getenv("BENCH_REPS", "3"))
# NB: substitute/eval, not force() - a promise evaluates once and every
# later repetition would time nothing.
med_time <- function(expr) {
  e <- substitute(expr); pf <- parent.frame()
  ts <- numeric(REPS)
  for (i in seq_len(REPS)) {
    gc(); t <- proc.time()[["elapsed"]]
    v <- suppressWarnings(eval(e, pf))
    ts[i] <- proc.time()[["elapsed"]] - t
  }
  list(t = median(ts), all = ts, v = v)
}
forms <- list(service = y ~ service + (1 | s) + (1 | d),
              plain   = y ~ 1 + (1 | s) + (1 | d))

cat("### frm defaults vs options vs lmer (median of", REPS, ") ###\n")
for (fnm in names(forms)) {
  f <- forms[[fnm]]
  cat("\n-- ", fnm, " --\n", sep = "")
  a <- med_time(frm(f, data = d, family = gaussian(), se = FALSE))
  cat(sprintf("frm default (nlminb, restarts=1, se=FALSE) : %6.2f s  obj %.4f\n",
              a$t, a$v$opt$objective))
  b <- med_time(frm(f, data = d, family = gaussian(), se = TRUE))
  cat(sprintf("frm default + se = TRUE                    : %6.2f s\n", b$t))
  cc <- med_time(frm(f, data = d, family = gaussian(), se = FALSE,
                     control = frmtmb_control(optimizer = "optim")))
  cat(sprintf("frm optimizer = 'optim' (L-BFGS-B)         : %6.2f s  obj %.4f\n",
              cc$t, cc$v$opt$objective))
  e <- med_time(frm(f, data = d, family = gaussian(), se = FALSE,
                    control = frmtmb_control(profile = TRUE)))
  cat(sprintf("frm profile = TRUE                         : %6.2f s  obj %.4f\n",
              e$t, e$v$opt$objective))
  g <- med_time(frm(f, data = d, family = gaussian(), se = FALSE,
                    control = frmtmb_control(sparse_x = TRUE)))
  cat(sprintf("frm sparse_x = TRUE                        : %6.2f s  obj %.4f\n",
              g$t, g$v$opt$objective))
  h <- med_time(frm(f, data = d, family = gaussian(), REML = TRUE,
                    se = FALSE))
  cat(sprintf("frm REML = TRUE                            : %6.2f s\n", h$t))
  l1 <- med_time(lmer(f, data = d, REML = FALSE))
  cat(sprintf("lme4::lmer(REML = FALSE)                   : %6.2f s  dev %.4f\n",
              l1$t, deviance(l1$v) / 2))
  l2 <- med_time(lmer(f, data = d, REML = TRUE))
  cat(sprintf("lme4::lmer(REML = TRUE)                    : %6.2f s\n", l2$t))
  cat("frm ML objective vs lmer -logLik(ML):",
      format(a$v$opt$objective, digits = 10), "vs",
      format(-as.numeric(logLik(l1$v)), digits = 10), "\n")
}

## --- where does an obj$fn call go? ----------------------------------
cat("\n### Rprof attribution of one full fit (service model) ###\n")
pf <- tempfile()
Rprof(pf, interval = 0.005, line.profiling = FALSE)
fit <- suppressWarnings(frm(forms$service, data = d, family = gaussian(),
                            se = TRUE))
Rprof(NULL)
sp <- summaryRprof(pf)
cat("\nby.self (top 15):\n")
print(head(sp$by.self, 15))
cat("\nsampling time:", sp$sampling.time, "s\n")

cat("\n### fn/gr are pure C++ after taping: R-level profile of 20 fn calls ###\n")
obj <- fit$obj; p <- fit$opt$par
pf2 <- tempfile()
Rprof(pf2, interval = 0.005)
for (i in 1:20) obj$fn(p + rnorm(length(p), 0, 0.02))
Rprof(NULL)
print(head(summaryRprof(pf2)$by.self, 10))

cat("\n### model size ###\n")
cat("outer parameters:", length(obj$par), "\n")
cat("inner (random) parameters:", length(obj$env$random), "\n")
h <- obj$env$spHess(random = TRUE)
cat("inner Hessian:", nrow(h), "x", ncol(h), ", nnz =", length(h@x), "\n")
