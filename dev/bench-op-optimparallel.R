## optimParallel vs the built-in optimizers on InstEval.
##
## Modes:
##   nlminb   package default, sequential, exact AD gradient
##   optim    package "optim" (L-BFGS-B), sequential, exact AD gradient
##   opP-cold optimParallel + exact gradient, 2 workers, cluster and
##            worker tapes built inside the fit (honest end to end)
##   opP-warm same, but cluster + tapes prepared once beforehand
##   opP-fd   optimParallel native mode: gr = NULL, central differences,
##            2p+1 tasks on 2p+1 workers
##
## Every mode is instrumented so "time inside fn/gr" is separable from
## optimizer and communication overhead. Run this script alone.
LIB <- "C:/Users/adf44/AppData/Local/Temp/1/claude/c--Users-adf44-source-r-frmtmb/529b6e73-d28f-46aa-a279-7dbeeb58fd4f/scratchpad/lib"
.libPaths(c(LIB, .libPaths()))
suppressPackageStartupMessages({
  library(frmtmb); library(parallel); library(optimParallel)
})
d <- lme4::InstEval
form <- y ~ service + (1 | s) + (1 | d)
REPS <- as.integer(Sys.getenv("BENCH_REPS", "3"))

tm <- function(expr) {
  t <- proc.time()[["elapsed"]]; v <- force(expr)
  list(t = proc.time()[["elapsed"]] - t, v = v)
}

## --- worker side ----------------------------------------------------
W_build_tape <- function(frame) {
  nll <- frmtmb:::build_objective(frame)
  tpl <- frmtmb:::make_start(frame, NULL)
  RTMB::MakeADFun(nll, tpl, random = "b", map = frame$map, silent = TRUE)
}
worker_init <- bquote({
  .libPaths(c(.(LIB), .libPaths()))
  suppressPackageStartupMessages({library(RTMB); library(frmtmb)})
  invisible(TRUE)
})
prep_cluster <- function(ncl, frame) {
  cl <- makePSOCKcluster(ncl)
  invisible(clusterCall(cl, eval, worker_init, .GlobalEnv))
  clusterExport(cl, c("W_frame", "W_build_tape"), envir = .GlobalEnv)
  invisible(clusterEvalQ(cl, {
    wobj <<- W_build_tape(W_frame); invisible(TRUE)
  }))
  cl
}
# these closures are serialized with .GlobalEnv by reference, so on the
# worker `wobj` resolves to the worker's own rebuilt tape
W_fn <- function(par) as.numeric(wobj$fn(par))
W_gr <- function(par) as.numeric(wobj$gr(par))
environment(W_fn) <- .GlobalEnv
environment(W_gr) <- .GlobalEnv

## --- instrumented sequential optimizers -----------------------------
E <- new.env()
reset <- function() {
  E$n_fn <- 0L; E$n_gr <- 0L; E$t_fn <- 0; E$t_gr <- 0
  E$rounds <- 0L; E$t_round <- 0; E$setup <- 0
}
wrap <- function(f, which) function(p) {
  t <- proc.time()[["elapsed"]]; v <- f(p)
  dt <- proc.time()[["elapsed"]] - t
  if (which == "fn") { E$n_fn <- E$n_fn + 1L; E$t_fn <- E$t_fn + dt }
  else { E$n_gr <- E$n_gr + 1L; E$t_gr <- E$t_gr + dt }
  v
}
opt_nlminb <- function(par, fn, gr, lower, upper, control) {
  reset()
  r <- stats::nlminb(par, wrap(fn, "fn"), wrap(gr, "gr"),
                     control = control, lower = lower, upper = upper)
  E$iters <- r$iterations
  list(par = r$par, objective = r$objective,
       convergence = r$convergence, message = r$message)
}
opt_lbfgsb <- function(par, fn, gr, lower, upper, control) {
  reset()
  ctl <- control[names(control) %in% c("maxit", "factr", "pgtol")]
  if (is.null(ctl$maxit)) ctl$maxit <- 1000
  r <- stats::optim(par, wrap(fn, "fn"), wrap(gr, "gr"),
                    method = "L-BFGS-B", lower = lower, upper = upper,
                    control = ctl)
  E$iters <- r$counts[["function"]]
  list(par = r$par, objective = r$value, convergence = r$convergence,
       message = r$message %||% "")
}
`%||%` <- function(a, b) if (is.null(a)) b else a

## --- optimParallel optimizers ---------------------------------------
make_opP <- function(ncl, use_gr = TRUE, forward = FALSE, cl = NULL) {
  function(par, fn, gr, lower, upper, control) {
    reset()
    own <- is.null(cl)
    if (own) {
      s <- tm(prep_cluster(ncl, W_frame)); cl_use <- s$v; E$setup <- s$t
    } else cl_use <- cl
    ctl <- control[names(control) %in% c("maxit", "factr", "pgtol")]
    if (is.null(ctl$maxit)) ctl$maxit <- 1000
    t0 <- proc.time()[["elapsed"]]
    r <- optimParallel::optimParallel(
      par, W_fn, if (use_gr) W_gr else NULL,
      lower = lower, upper = upper, control = ctl,
      parallel = list(cl = cl_use, forward = forward, loginfo = TRUE))
    E$t_round <- proc.time()[["elapsed"]] - t0
    E$rounds <- nrow(r$loginfo)
    E$iters <- r$counts[["function"]]
    if (own) stopCluster(cl_use)
    list(par = r$par, objective = r$value,
         convergence = r$convergence, message = r$message %||% "")
  }
}

## --- driver ---------------------------------------------------------
W_frame <<- frm(form, data = d, family = gaussian(), dry_run = "frame")
n_outer <- 5L

run_mode <- function(label, optimizer) {
  out <- vector("list", REPS)
  for (i in seq_len(REPS)) {
    gc()
    r <- tm(suppressWarnings(frm(
      form, data = d, family = gaussian(), se = FALSE,
      control = frmtmb_control(optimizer = optimizer, restarts = 0))))
    mg <- max(abs(r$v$obj$gr(r$v$opt$par)))
    out[[i]] <- c(wall = r$t, obj = r$v$opt$objective,
                  maxgrad = mg, conv = r$v$opt$convergence,
                  n_fn = E$n_fn, n_gr = E$n_gr, t_fn = E$t_fn,
                  t_gr = E$t_gr, rounds = E$rounds,
                  t_opt = E$t_round, setup = E$setup, iters = E$iters)
  }
  m <- do.call(rbind, out)
  cat("\n=== ", label, " ===\n", sep = "")
  print(round(m, 3))
  med <- apply(m, 2, median)
  cat("median: ")
  cat(paste(sprintf("%s=%.3f", names(med), med), collapse = "  "), "\n")
  invisible(med)
}

res <- list()
res$nlminb <- run_mode("nlminb (sequential, exact gr)", opt_nlminb)
res$optim  <- run_mode("optim L-BFGS-B (sequential, exact gr)", opt_lbfgsb)
res$opP_cold <- run_mode("optimParallel exact gr, 2 workers, cold cluster",
                         make_opP(2L, use_gr = TRUE))

warm_cl <- prep_cluster(2L, W_frame)
res$opP_warm <- run_mode("optimParallel exact gr, 2 workers, warm cluster",
                         make_opP(2L, use_gr = TRUE, cl = warm_cl))
stopCluster(warm_cl)

nfd <- 2L * n_outer + 1L
cat("\nfinite-difference mode uses", nfd, "workers (2p+1, p =", n_outer, ")\n")
res$opP_fd <- run_mode(
  sprintf("optimParallel finite-difference (central), %d workers", nfd),
  make_opP(nfd, use_gr = FALSE, forward = FALSE))

nfw <- n_outer + 1L
res$opP_fwd <- run_mode(
  sprintf("optimParallel finite-difference (forward), %d workers", nfw),
  make_opP(nfw, use_gr = FALSE, forward = TRUE))

saveRDS(res, file.path(dirname(LIB), "optimparallel.rds"))
cat("\n--- summary (median seconds) ---\n")
print(round(do.call(rbind, res), 3))
