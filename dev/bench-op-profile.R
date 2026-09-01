## Where does "optimize" go? A custom optimizer that wraps fn/gr in
## timing counters and then calls nlminb exactly as the package does.
## Answers: evaluation time vs optimizer overhead, and the true nlminb
## iteration/evaluation counts (restarts = 0 so nothing is masked).
.libPaths(c("C:/Users/adf44/AppData/Local/Temp/1/claude/c--Users-adf44-source-r-frmtmb/529b6e73-d28f-46aa-a279-7dbeeb58fd4f/scratchpad/lib",
            .libPaths()))
suppressPackageStartupMessages(library(frmtmb))
d <- lme4::InstEval

counter_env <- new.env()

timing_nlminb <- function(par, fn, gr, lower, upper, control) {
  e <- counter_env
  e$n_fn <- 0L; e$n_gr <- 0L; e$t_fn <- 0; e$t_gr <- 0
  e$fn_times <- numeric(0); e$gr_times <- numeric(0)
  wfn <- function(p) {
    t <- proc.time()[["elapsed"]]
    v <- fn(p)
    dt <- proc.time()[["elapsed"]] - t
    e$n_fn <- e$n_fn + 1L; e$t_fn <- e$t_fn + dt
    e$fn_times <- c(e$fn_times, dt)
    v
  }
  wgr <- function(p) {
    t <- proc.time()[["elapsed"]]
    v <- gr(p)
    dt <- proc.time()[["elapsed"]] - t
    e$n_gr <- e$n_gr + 1L; e$t_gr <- e$t_gr + dt
    e$gr_times <- c(e$gr_times, dt)
    v
  }
  t0 <- proc.time()[["elapsed"]]
  r <- stats::nlminb(par, wfn, wgr, control = control,
                     lower = lower, upper = upper)
  e$wall <- proc.time()[["elapsed"]] - t0
  e$iterations <- r$iterations
  e$evaluations <- r$evaluations
  list(par = r$par, objective = r$objective,
       convergence = r$convergence, message = r$message)
}

forms <- list(
  service = y ~ service + (1 | s) + (1 | d),
  plain   = y ~ 1 + (1 | s) + (1 | d))

for (nm in names(forms)) {
  cat("\n=== ", nm, " ===\n", sep = "")
  gc()
  fit <- suppressWarnings(frm(
    forms[[nm]], data = d, family = gaussian(), se = FALSE,
    control = frmtmb_control(optimizer = timing_nlminb, restarts = 0)))
  e <- counter_env
  cat(sprintf("nlminb iterations: %d   evaluations: fn %d, gr %d\n",
              e$iterations, e$evaluations[[1]], e$evaluations[[2]]))
  cat(sprintf("wrapper counts:    fn %d, gr %d\n", e$n_fn, e$n_gr))
  cat(sprintf("optimizer wall:    %.2f s\n", e$wall))
  cat(sprintf("  time in fn:      %.2f s (%.1f%%), median call %.4f s\n",
              e$t_fn, 100 * e$t_fn / e$wall, median(e$fn_times)))
  cat(sprintf("  time in gr:      %.2f s (%.1f%%), median call %.4f s\n",
              e$t_gr, 100 * e$t_gr / e$wall, median(e$gr_times)))
  cat(sprintf("  nlminb overhead: %.2f s (%.1f%%)\n",
              e$wall - e$t_fn - e$t_gr,
              100 * (e$wall - e$t_fn - e$t_gr) / e$wall))
  cat("first 6 fn call times:", paste(round(head(e$fn_times, 6), 3),
                                      collapse = " "), "\n")
  cat("last 6 fn call times: ", paste(round(tail(e$fn_times, 6), 3),
                                      collapse = " "), "\n")
  cat("first 6 gr call times:", paste(round(head(e$gr_times, 6), 3),
                                      collapse = " "), "\n")
  cat("last 6 gr call times: ", paste(round(tail(e$gr_times, 6), 3),
                                      collapse = " "), "\n")
  cat("objective:", format(fit$opt$objective, digits = 10),
      " convergence:", fit$opt$convergence, "\n")
  cat("outer par at optimum:\n"); print(round(fit$opt$par, 5))
  saveRDS(as.list(e), sprintf(
    "C:/Users/adf44/AppData/Local/Temp/1/claude/c--Users-adf44-source-r-frmtmb/529b6e73-d28f-46aa-a279-7dbeeb58fd4f/scratchpad/prof-%s.rds",
    nm))
}
