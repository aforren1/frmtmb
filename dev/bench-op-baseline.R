## Baseline: InstEval fits with frm(), stage timings from verbose = 1,
## nlminb iteration/evaluation counts, and per-call cost of obj$fn and
## obj$gr at the optimum.
##
## Run alone (timings on this machine are noisy with other R processes).
.libPaths(c("C:/Users/adf44/AppData/Local/Temp/1/claude/c--Users-adf44-source-r-frmtmb/529b6e73-d28f-46aa-a279-7dbeeb58fd4f/scratchpad/lib",
            .libPaths()))
suppressPackageStartupMessages(library(frmtmb))
data(InstEval, package = "lme4")
d <- lme4::InstEval

# stage lines go through message(); capture and parse them
capture_stages <- function(expr) {
  lines <- character(0)
  val <- withCallingHandlers(
    force(expr),
    message = function(m) {
      lines <<- c(lines, conditionMessage(m))
      invokeRestart("muffleMessage")
    })
  st <- regmatches(lines, regexec(
    "^frmtmb: ([a-z ]+[a-z0-9]*) \\[([0-9.]+)s\\]", lines))
  st <- Filter(function(x) length(x) == 3L, st)
  nm <- vapply(st, `[`, "", 2L)
  tm <- as.numeric(vapply(st, `[`, "", 3L))
  list(value = val, stages = stats::setNames(tm, nm), lines = lines)
}

run_one <- function(form, se = TRUE) {
  gc()
  t0 <- proc.time()[["elapsed"]]
  r <- capture_stages(
    frm(form, data = d, family = gaussian(), se = se,
        control = frmtmb_control(verbose = 1)))
  wall <- proc.time()[["elapsed"]] - t0
  c(r$stages, total_wall = wall)
}

forms <- list(
  service = y ~ service + (1 | s) + (1 | d),
  plain   = y ~ 1 + (1 | s) + (1 | d))

res <- list()
for (nm in names(forms)) {
  reps <- replicate(3, run_one(forms[[nm]]), simplify = FALSE)
  keys <- unique(unlist(lapply(reps, names)))
  m <- vapply(keys, function(k)
    stats::median(vapply(reps, function(r) unname(r[k]), 0)), 0)
  res[[nm]] <- m
  cat("\n== ", nm, " (median of 3) ==\n", sep = "")
  print(round(m, 3))
  cat("all reps:\n")
  print(round(do.call(rbind, lapply(reps, function(r) r[keys])), 3))
}

## --- optimizer counts + per-call fn/gr cost -------------------------
for (nm in names(forms)) {
  cat("\n### ", nm, " optimizer detail ###\n", sep = "")
  fit <- frm(forms[[nm]], data = d, family = gaussian(), se = FALSE)
  cat("iterations:", fit$opt$iterations, "\n")
  cat("evaluations (function, gradient):",
      paste(fit$opt$evaluations, collapse = ", "), "\n")
  cat("objective:", format(fit$opt$objective, digits = 10), "\n")
  cat("n outer:", length(fit$obj$par), " n inner:",
      length(fit$obj$env$random), "\n")
  obj <- fit$obj
  p <- fit$opt$par

  # (a) repeated calls at exactly the optimum: the inner Newton solve
  #     warm-starts at its own solution, so this is a lower bound
  tfn <- replicate(20, system.time(obj$fn(p))[["elapsed"]])
  tgr <- replicate(20, system.time(obj$gr(p))[["elapsed"]])
  # (b) jittered parameters: each call does a real inner solve, which is
  #     what an optimizer step actually pays
  set.seed(1)
  jf <- replicate(20, {
    q <- p + stats::rnorm(length(p), 0, 0.02)
    system.time(obj$fn(q))[["elapsed"]]
  })
  set.seed(2)
  jg <- replicate(20, {
    q <- p + stats::rnorm(length(p), 0, 0.02)
    system.time(obj$gr(q))[["elapsed"]]
  })
  cat(sprintf("fn at optimum : median %.4f s (min %.4f, max %.4f)\n",
              median(tfn), min(tfn), max(tfn)))
  cat(sprintf("gr at optimum : median %.4f s (min %.4f, max %.4f)\n",
              median(tgr), min(tgr), max(tgr)))
  cat(sprintf("fn jittered   : median %.4f s\n", median(jf)))
  cat(sprintf("gr jittered   : median %.4f s\n", median(jg)))
  ev <- fit$opt$evaluations
  cat(sprintf("implied eval time (jittered): %.2f s  (%d fn x %.4f + %d gr x %.4f)\n",
              ev[[1]] * median(jf) + ev[[2]] * median(jg),
              ev[[1]], median(jf), ev[[2]], median(jg)))
  cat(sprintf("implied eval time (at optimum): %.2f s\n",
              ev[[1]] * median(tfn) + ev[[2]] * median(tgr)))
}

saveRDS(res, "C:/Users/adf44/AppData/Local/Temp/1/claude/c--Users-adf44-source-r-frmtmb/529b6e73-d28f-46aa-a279-7dbeeb58fd4f/scratchpad/baseline.rds")
