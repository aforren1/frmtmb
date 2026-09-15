## The ?cross_wishart example, timed in its own process against the
## library named on the command line, with a CONTROL beside it.
##
## The lane's R CMD check reported the examples NOTE at 4.69 + 0.59 s
## of CPU against a 5 s threshold and blamed machine load. CPU time is
## supposed to be load independent, so the claim needs the control:
## `calib` is a fixed arithmetic block that does the same work in every
## process, so its CPU time measures how fast this machine is RIGHT NOW
## and any inflation in the example has to show up there too.
args <- commandArgs(trailingOnly = TRUE)
lib <- args[1L]
tag <- args[2L]
rounds <- as.integer(args[3L])
.libPaths(c(lib, "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
suppressMessages(library(frmtmb.coupling))
cat(tag, "using", getNamespaceInfo("frmtmb.coupling", "path"), "\n")

calib <- function() {
  x <- matrix(seq_len(400L * 400L) / 1e5, 400L, 400L)
  s <- 0
  for (i in 1:6) s <- s + sum(x %*% x)
  s
}
ex <- function() {
  set.seed(4)
  src <- rnorm(4096)
  xs <- frm_cross_spectrum(src + rnorm(4096), 0.9 * src + rnorm(4096),
                           segments = 16)
  xs <- xs[xs$freq < 0.1, ]
  fit <- frmtmb::frm(
    frmtmb::bf(w11 | vreal(w22, w12r, w12i) + vint(n) ~ 1,
               pow2 ~ 1, coh ~ 1, phase ~ 1),
    family = cross_wishart(), data = xs)
  frm_coherence(fit)[1, ]
  fit
}
fit <- NULL
for (i in seq_len(rounds)) {
  tc <- system.time(calib())
  t <- system.time(fit <- ex())
  cat(sprintf(paste0("RESULT %s round %d ex_cpu %.2f ex_user %.2f ",
                     "ex_sys %.2f ex_elapsed %.2f calib_cpu %.3f ",
                     "calib_elapsed %.3f\n"),
              tag, i, t[["user.self"]] + t[["sys.self"]],
              t[["user.self"]], t[["sys.self"]], t[["elapsed"]],
              tc[["user.self"]] + tc[["sys.self"]], tc[["elapsed"]]))
}
cat(sprintf("WORK %s objective %.10f iters %s evals %s npar %d\n", tag,
            as.numeric(fit$opt$objective),
            paste(fit$opt$iterations, collapse = ","),
            paste(fit$opt$evaluations, collapse = ","),
            length(fit$opt$par)))
