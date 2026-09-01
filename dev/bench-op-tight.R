## Follow-up to bench-op-optimparallel.R.
##
## L-BFGS-B stopped at max|grad| 0.165 where nlminb reached 0.002, so the
## L-BFGS-B timings above are partly bought with looser convergence.
## Re-time sequential and parallel L-BFGS-B at a tight factr, and repeat
## the head-to-head on the plain (no service) model.
LIB <- "C:/Users/adf44/AppData/Local/Temp/1/claude/c--Users-adf44-source-r-frmtmb/529b6e73-d28f-46aa-a279-7dbeeb58fd4f/scratchpad/lib"
.libPaths(c(LIB, .libPaths()))
suppressPackageStartupMessages({
  library(frmtmb); library(parallel); library(optimParallel)
})
d <- lme4::InstEval
REPS <- as.integer(Sys.getenv("BENCH_REPS", "3"))
`%||%` <- function(a, b) if (is.null(a)) b else a
tm <- function(expr) {
  t <- proc.time()[["elapsed"]]; v <- force(expr)
  list(t = proc.time()[["elapsed"]] - t, v = v)
}

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
prep_cluster <- function(ncl) {
  cl <- makePSOCKcluster(ncl)
  invisible(clusterCall(cl, eval, worker_init, .GlobalEnv))
  clusterExport(cl, c("W_frame", "W_build_tape"), envir = .GlobalEnv)
  invisible(clusterEvalQ(cl, {wobj <<- W_build_tape(W_frame); invisible(TRUE)}))
  cl
}
W_fn <- function(par) as.numeric(wobj$fn(par))
W_gr <- function(par) as.numeric(wobj$gr(par))
environment(W_fn) <- .GlobalEnv
environment(W_gr) <- .GlobalEnv

E <- new.env()
reset <- function() {
  E$n_fn <- 0L; E$n_gr <- 0L; E$t_fn <- 0; E$t_gr <- 0
  E$rounds <- 0L; E$t_opt <- 0
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
make_lbfgsb <- function(factr) function(par, fn, gr, lower, upper, control) {
  reset()
  r <- stats::optim(par, wrap(fn, "fn"), wrap(gr, "gr"),
                    method = "L-BFGS-B", lower = lower, upper = upper,
                    control = list(maxit = 1000, factr = factr))
  E$iters <- r$counts[["function"]]
  list(par = r$par, objective = r$value, convergence = r$convergence,
       message = r$message %||% "")
}
make_opP <- function(cl, factr) function(par, fn, gr, lower, upper, control) {
  reset()
  t0 <- proc.time()[["elapsed"]]
  r <- optimParallel::optimParallel(
    par, W_fn, W_gr, lower = lower, upper = upper,
    control = list(maxit = 1000, factr = factr),
    parallel = list(cl = cl, forward = FALSE, loginfo = TRUE))
  E$t_opt <- proc.time()[["elapsed"]] - t0
  E$rounds <- nrow(r$loginfo); E$iters <- r$counts[["function"]]
  list(par = r$par, objective = r$value, convergence = r$convergence,
       message = r$message %||% "")
}

run_mode <- function(label, optimizer, form) {
  out <- vector("list", REPS)
  for (i in seq_len(REPS)) {
    gc()
    r <- tm(suppressWarnings(frm(
      form, data = d, family = gaussian(), se = FALSE,
      control = frmtmb_control(optimizer = optimizer, restarts = 0))))
    out[[i]] <- c(wall = r$t, obj = r$v$opt$objective,
                  maxgrad = max(abs(r$v$obj$gr(r$v$opt$par))),
                  n_fn = E$n_fn, n_gr = E$n_gr, t_fn = E$t_fn,
                  t_gr = E$t_gr, rounds = E$rounds, t_opt = E$t_opt,
                  iters = E$iters)
  }
  m <- apply(do.call(rbind, out), 2, median)
  cat(sprintf("%-42s wall %6.2f  obj %.4f  max|g| %8.4f  fn %3.0f gr %3.0f  t_fn %5.2f t_gr %5.2f  rounds %3.0f t_opt %6.2f  iters %3.0f\n",
              label, m[["wall"]], m[["obj"]], m[["maxgrad"]],
              m[["n_fn"]], m[["n_gr"]], m[["t_fn"]], m[["t_gr"]],
              m[["rounds"]], m[["t_opt"]], m[["iters"]]))
  invisible(m)
}

forms <- list(service = y ~ service + (1 | s) + (1 | d),
              plain   = y ~ 1 + (1 | s) + (1 | d))
res <- list()
for (fnm in names(forms)) {
  form <- forms[[fnm]]
  W_frame <<- frm(form, data = d, family = gaussian(), dry_run = "frame")
  cat("\n########## ", fnm, " ##########\n", sep = "")
  cl <- prep_cluster(2L)
  res[[paste0(fnm, ".nlminb")]] <-
    run_mode("nlminb (default optCtrl)", opt_nlminb, form)
  for (fc in c(1e7, 1e2)) {
    res[[paste0(fnm, ".seq.", fc)]] <-
      run_mode(sprintf("L-BFGS-B seq, factr=%g", fc), make_lbfgsb(fc), form)
    res[[paste0(fnm, ".opP.", fc)]] <-
      run_mode(sprintf("optimParallel(2w) warm, factr=%g", fc),
               make_opP(cl, fc), form)
  }
  stopCluster(cl)
}
saveRDS(res, file.path(dirname(LIB), "tight.rds"))
