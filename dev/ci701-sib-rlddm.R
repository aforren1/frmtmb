# Lane ci701: sibling probe. frmtmb.learn test-stan-identity.R "rlddm
# reproduces Stan with a correlated block on all four" fits a 4 x 4
# correlated block whose correlation parameter ran to 2413 on this
# machine (dev/ci701-run1.R trace), past the clamp, where the outer
# gradient is noise. The test allows only "false convergence (8)" and
# the large-gradient warning. Here the fit is repeated with rt moved by
# a few ulps (seed 5000 + k) and every warning and message counted.
#
#   Rscript dev/ci701-sib-rlddm.R <lib|base> <k0> <k1>
args <- commandArgs(TRUE)
lib <- args[1]; k0 <- as.integer(args[2]); k1 <- as.integer(args[3])
base <- "C:/Users/adf44/source/r/rellib-r8"
ul <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
.libPaths(if (lib == "base") c(base, ul) else c(lib, base, ul))
suppressMessages(library(frmtmb))
suppressMessages(library(frmtmb.learn))
cat("lib:", find.package("frmtmb"), as.character(packageVersion("frmtmb")),
    "\n")
set.seed(204)
ns <- 20L
fam <- rlddm(subject = id, trial = trial)
d <- frm_task_design("bandit2arm", n_subject = ns, n_trial = 120,
                     seed = 204)
R <- matrix(c(1.0,  0.4,  0.0,  0.0,
              0.4,  1.0, -0.3,  0.0,
              0.0, -0.3,  1.0,  0.3,
              0.0,  0.0,  0.3,  1.0), 4L, 4L, byrow = TRUE)
U <- matrix(stats::rnorm(4L * ns), ns, 4L) %*% chol(R)
U <- sweep(U, 2L, c(0.5, 1.0, 0.2, 0.15), "*")
i <- as.integer(d$id)
d0 <- frm_task_simulate(
  fam, d,
  pars = list(alpha = stats::plogis(stats::qlogis(0.35) + U[i, 1L]),
              drift = 2.5 + U[i, 2L], bs = 1.5 * exp(U[i, 3L]),
              ndt = 0.25 * exp(U[i, 4L]), bias = 0.5),
  seed = 204)[[1L]]
allowed <- c("Optimizer did not report convergence: false convergence (8)",
             "Large maximum absolute gradient at the optimum")
for (k in k0:k1) {
  d <- d0
  if (k > 0) {
    set.seed(5000 + k)
    d$rt <- d$rt * (1 + sample(0:4, nrow(d), TRUE) * 2^-52)
  }
  w <- character(0); bm <- 0L
  fit <- withCallingHandlers(
    frmtmb::frm(
      frmtmb::bf(rt | dec(choice) + reward(pay1, pay2) +
                   ndt_group(id) ~ 1 + (1 | p | id),
                 drift ~ 1 + (1 | p | id), bs ~ 1 + (1 | p | id),
                 ndt ~ 1 + (1 | p | id), bias = 0.5),
      family = fam, data = d),
    warning = function(cnd) {
      w <<- c(w, conditionMessage(cnd)); invokeRestart("muffleWarning")
    },
    frmtmb_boundary_fit = function(m) {
      bm <<- bm + 1L; invokeRestart("muffleMessage")
    })
  esc <- w[!vapply(w, function(m) {
    any(vapply(allowed, function(a) grepl(a, m, fixed = TRUE), NA))
  }, NA)]
  th <- fit$estimates$theta
  cat(sprintf("RES k=%d code=%d msg=%s max|t|=%.1f ll=%.6f nwarn=%d bmsg=%d escaped=%d\n",
              k, fit$opt$convergence, gsub(" ", "_", fit$opt$message),
              max(abs(th[-(1:4)])), as.numeric(logLik(fit)), length(w), bm,
              length(esc)))
  for (m in unique(esc)) cat("   ESC:", substr(m, 1, 160), "\n")
}
