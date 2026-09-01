## Does an RTMB tape survive the trip to a PSOCK worker?
## (1) locate the external pointers inside obj, (2) clusterExport the obj
## and show what the pointers look like on the worker, (3) call obj$fn
## there and record the failure, (4) time the only viable alternative:
## ship the frame and rebuild the tape in the worker.
LIB <- "C:/Users/adf44/AppData/Local/Temp/1/claude/c--Users-adf44-source-r-frmtmb/529b6e73-d28f-46aa-a279-7dbeeb58fd4f/scratchpad/lib"
.libPaths(c(LIB, .libPaths()))
suppressPackageStartupMessages({library(frmtmb); library(parallel)})
d <- lme4::InstEval
form <- y ~ service + (1 | s) + (1 | d)

worker_init <- bquote({
  .libPaths(c(.(LIB), .libPaths()))
  suppressPackageStartupMessages({library(RTMB); library(frmtmb)})
  invisible(TRUE)
})

## --- build a frame and a tape locally -------------------------------
t0 <- proc.time()[["elapsed"]]
frame <- frm(form, data = d, family = gaussian(), dry_run = "frame")
t_frame <- proc.time()[["elapsed"]] - t0
cat(sprintf("frame build: %.2f s, serialized size %.1f MB\n",
            t_frame, length(serialize(frame, NULL)) / 2^20))

build_tape <- function(frame) {
  nll <- frmtmb:::build_objective(frame)
  tpl <- frmtmb:::make_start(frame, NULL)
  RTMB::MakeADFun(nll, tpl, random = "b", map = frame$map, silent = TRUE)
}
t0 <- proc.time()[["elapsed"]]
obj <- build_tape(frame)
t_tape <- proc.time()[["elapsed"]] - t0
cat(sprintf("tape build : %.2f s\n", t_tape))
cat(sprintf("obj serialized size: %.1f MB\n",
            length(serialize(obj, NULL)) / 2^20))

## --- (1) locate every external pointer ------------------------------
find_ep <- function(x, path = "obj", depth = 0L, out = list()) {
  if (depth > 4L) return(out)
  if (typeof(x) == "externalptr") {
    out[[path]] <- utils::capture.output(print(x))[1]
    return(out)
  }
  if (is.environment(x)) {
    for (nm in ls(x, all.names = TRUE)) {
      out <- find_ep(tryCatch(get(nm, x), error = function(e) NULL),
                     paste0(path, "$", nm), depth + 1L, out)
    }
    return(out)
  }
  if (is.list(x) && length(x)) {
    nms <- names(x) %||% rep("", length(x))
    for (i in seq_along(x)) {
      p <- if (nzchar(nms[i])) paste0(path, "$", nms[i]) else
        paste0(path, "[[", i, "]]")
      out <- find_ep(x[[i]], p, depth + 1L, out)
    }
  }
  out
}
`%||%` <- function(a, b) if (is.null(a)) b else a
eps <- find_ep(obj)
cat("\nexternal pointers reachable from obj (depth <= 4):\n")
for (nm in names(eps)) cat("  ", nm, " -> ", eps[[nm]], "\n", sep = "")

## --- (2)/(3) export the tape and use it on a worker ------------------
cl <- makePSOCKcluster(1)
invisible(clusterCall(cl, eval, worker_init, .GlobalEnv))
clusterExport(cl, c("obj", "find_ep"), envir = environment())
cat("\n--- pointers as seen ON THE WORKER after clusterExport ---\n")
rp <- clusterEvalQ(cl, {
  `%||%` <- function(a, b) if (is.null(a)) b else a
  find_ep(obj)
})[[1]]
for (nm in names(rp)) cat("  ", nm, " -> ", rp[[nm]], "\n", sep = "")

cat("\n--- obj$fn(obj$par) ON THE WORKER ---\n")
res <- tryCatch(clusterEvalQ(cl, try(obj$fn(obj$par), silent = TRUE))[[1]],
                error = function(e) structure(conditionMessage(e),
                                              class = "worker_died"))
if (inherits(res, "worker_died")) {
  cat("worker process DIED. master-side error:\n  ", res, "\n", sep = "")
} else {
  cat("class:", paste(class(res), collapse = "/"), "\n")
  print(res)
}
try(stopCluster(cl), silent = TRUE)

cat("\n--- same call locally, for contrast ---\n")
print(obj$fn(obj$par))

## --- (4) ship the frame, rebuild the tape in the worker --------------
cat("\n--- ship frame, rebuild tape in worker ---\n")
t0 <- proc.time()[["elapsed"]]
cl <- makePSOCKcluster(2)
t_start <- proc.time()[["elapsed"]] - t0
t0 <- proc.time()[["elapsed"]]
invisible(clusterCall(cl, eval, worker_init, .GlobalEnv))
t_lib <- proc.time()[["elapsed"]] - t0
t0 <- proc.time()[["elapsed"]]
clusterExport(cl, c("frame", "build_tape"), envir = environment())
t_ship <- proc.time()[["elapsed"]] - t0
t0 <- proc.time()[["elapsed"]]
invisible(clusterEvalQ(cl, {wobj <<- build_tape(frame); TRUE}))
t_rebuild <- proc.time()[["elapsed"]] - t0
remote2 <- clusterEvalQ(cl, wobj$fn(wobj$par))[[1]]
cat(sprintf("PSOCK 2-worker startup   : %.2f s\n", t_start))
cat(sprintf("library(frmtmb) on both  : %.2f s\n", t_lib))
cat(sprintf("ship frame to 2 workers  : %.2f s\n", t_ship))
cat(sprintf("rebuild tape on both     : %.2f s\n", t_rebuild))
cat(sprintf("total per-fit setup      : %.2f s\n",
            t_start + t_lib + t_ship + t_rebuild))
cat("remote fn after rebuild:", format(remote2, digits = 10), "\n")
cat("local  fn             :", format(obj$fn(obj$par), digits = 10), "\n")
stopCluster(cl)
