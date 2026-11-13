# Review of lane gpstop: shared helpers. Arms in ONE process on the lane
# library: "lane" (as installed) and "norelin" (gp_body_relin_applies()
# returns FALSE, which is 0.81.0's route). Counters wrap gp_body_relin()
# and gp_second_descent().
RV_WT <- "C:/Users/adf44/source/r/frmtmb-wt-gpstop"
RV_LIB <- "C:/Users/adf44/source/r/wt-gpstop-lib"
RV_BASE <- "C:/Users/adf44/source/r/rellib-r20"
RV_USER <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
RV_OUT <- file.path(RV_WT, "dev/reviews/gpstop-review-out")
dir.create(RV_OUT, showWarnings = FALSE)
rv_libs <- function(arm = "lane") {
  .libPaths(c(if (arm == "lane") RV_LIB, RV_BASE, RV_USER))
  Sys.setenv(NOT_CRAN = "true")
  suppressPackageStartupMessages(library(frmtmb))
  cat("lib:", dirname(system.file(package = "frmtmb")), "frmtmb",
      as.character(packageVersion("frmtmb")), "\n")
}
rv_put <- function(name, f) {
  ns <- asNamespace("frmtmb")
  environment(f) <- ns
  unlockBinding(name, ns); assign(name, f, envir = ns); lockBinding(name, ns)
}
CNT <- new.env()
CNT$relin <- 0L
CNT$relin_hit <- 0L
rv_instrument <- function() {
  ns <- asNamespace("frmtmb")
  if (!exists("gp_body_relin", ns)) return(invisible())
  CNT$real_relin <- get("gp_body_relin", ns)
  CNT$real_applies <- get("gp_body_relin_applies", ns)
  wrap <- function(obj, frame, par) {
    .GlobalEnv$CNT$relin <- .GlobalEnv$CNT$relin + 1L
    p <- .GlobalEnv$CNT$real_relin(obj, frame, par)
    if (!is.null(p)) .GlobalEnv$CNT$relin_hit <- .GlobalEnv$CNT$relin_hit + 1L
    p
  }
  rv_put("gp_body_relin", wrap)
}
rv_arm <- function(arm) {
  if (arm == "norelin") {
    rv_put("gp_body_relin_applies", function(...) FALSE)
  } else if (!is.null(CNT$real_applies)) {
    rv_put("gp_body_relin_applies", CNT$real_applies)
  }
}
rv_fit <- function(...) {
  w <- character(0)
  f <- tryCatch(withCallingHandlers(frm(...), warning = function(e) {
    w <<- c(w, paste0(class(e)[1], ": ", substr(conditionMessage(e), 1, 90)))
    invokeRestart("muffleWarning")
  }, message = function(m) invokeRestart("muffleMessage")),
  error = function(e) e)
  list(fit = f, warn = w)
}
# max |gradient| and the smallest eigenvalue of the outer Hessian (finite
# differences of the tape's gradient) at the fit's optimum
rv_kkt <- function(fit) {
  p <- fit$opt$par
  g <- tryCatch(drop(fit$obj$gr(p)), error = function(e) NA)
  H <- tryCatch(stats::optimHess(p, fit$obj$fn, fit$obj$gr),
                error = function(e) NULL)
  ev <- if (is.null(H) || !all(is.finite(H))) NA else
    min(eigen((H + t(H)) / 2, symmetric = TRUE, only.values = TRUE)$values)
  tryCatch(fit$obj$fn(p), error = function(e) NULL)
  c(gmax = max(abs(g)), evmin = ev)
}
rv_row <- function(tag, arm, q, wall, n_relin, n_hit) {
  f <- q$fit
  if (inherits(f, "error")) {
    cat(sprintf("ROW %s %s ERROR %s\n", tag, arm, conditionMessage(f)))
    return(invisible())
  }
  th <- f$opt$par[names(f$opt$par) == "theta"]
  k <- rv_kkt(f)
  left <- any(c(".frmtmb_gp_relin", ".frmtmb_gp_alt") %in%
                ls(f$obj$env, all.names = TRUE))
  gs <- f$opt$gp_second
  cat(sprintf(paste0("ROW %s %s ll %.6f conv %d nwarn %d wall %.2f relin %d ",
                     "hit %d runs %s gmax %.2e evmin %.3e left %s theta %s",
                     " | %s\n"),
              tag, arm, as.numeric(logLik(f)), f$opt$convergence,
              length(q$warn), wall, n_relin, n_hit,
              if (is.null(gs)) "NA" else format(gs[["runs"]]),
              k[["gmax"]], k[["evmin"]], left,
              paste(sprintf("%.3f", th), collapse = ","),
              paste(unique(sub(":.*", "", q$warn)), collapse = ";")))
  flush.console()
}
rv_run <- function(tag, arm, ...) {
  rv_arm(arm)
  CNT$relin <- 0L
  CNT$relin_hit <- 0L
  t0 <- proc.time()[[3]]
  q <- rv_fit(...)
  wall <- proc.time()[[3]] - t0
  rv_row(tag, arm, q, wall, CNT$relin, CNT$relin_hit)
  invisible(q)
}
