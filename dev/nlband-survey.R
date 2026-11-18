# Lane nlband: run one of lane nlstop2's survey scripts (fresh rows, scan
# rows, limit battery) on this lane's arms, and count what the pole
# restart costs. The script's own setup line is replaced by
# dev/nlband-common.R. Prints, besides the script's ROW lines, one ESC line
# per call of nl_pole_escape() (after arm: the hole and past starts it
# ran) and a COST line at the end: the objective and gradient evaluations
# of every optimize_obj() run, in all and inside nl_pole_escape().
# NLBAND_ARM=before|after Rscript dev/nlband-survey.R <script> [args]
argv <- commandArgs(trailingOnly = TRUE)
script <- argv[1]
rest <- argv[-1]
source("dev/nlband-common.R")
llp <- function(r) {
  if (iserr(r) || is.null(r$opt)) return(NA_real_)
  v <- suppressWarnings(as.numeric(logLik_penalized(r)))
  if (is.finite(v)) v else NA_real_
}
fresh_defs <- function(env = parent.frame()) {
  src <- parse("dev/nlstop-fresh.R")
  keep <- vapply(src, function(e) {
    is.call(e) && identical(e[[1L]], as.name("<-")) &&
      as.character(e[[2L]])[1] %in% c("F", "pos", "B", "gen", "form")
  }, NA)
  for (e in src[keep]) eval(e, env)
  invisible(NULL)
}
cost <- new.env()
cost$ev <- 0
cost$ev_esc <- 0
cost$in_esc <- FALSE
cost$n_esc <- 0L
cost$holes <- 0L
cost$past <- 0L
suppressMessages(trace("optimize_obj", exit = quote({
  .nlb <- get("cost", envir = globalenv())
  .nlv <- returnValue()$evals
  if (is.numeric(.nlv)) {
    .nlb$ev <- .nlb$ev + sum(.nlv)
    if (.nlb$in_esc) .nlb$ev_esc <- .nlb$ev_esc + sum(.nlv)
  }
}), print = FALSE, where = ns))
suppressMessages(trace("nl_pole_escape", tracer = quote({
  .nlb <- get("cost", envir = globalenv()); .nlb$in_esc <- TRUE
}), exit = quote({
  .nlb <- get("cost", envir = globalenv())
  .nlb$in_esc <- FALSE
  .nlb$n_esc <- .nlb$n_esc + 1L
  if (exists("checks", inherits = FALSE)) {
    .nlk <- vapply(checks, function(ck) ck$kind %||% "past", "")
    .nlb$holes <- .nlb$holes + sum(.nlk == "hole")
    .nlb$past <- .nlb$past + sum(.nlk == "past")
    cat("ESC holes", sum(.nlk == "hole"), "past", sum(.nlk == "past"), "\n")
  }
}), print = FALSE, where = ns))
`%||%` <- function(a, b) if (is.null(a)) b else a
lines <- readLines(script)
lines <- sub('^source[(]"dev/nlstop2-common.R"[)]', "invisible(NULL)", lines)
commandArgs <- function(trailingOnly = FALSE) rest
.nlb_t0 <- proc.time()[["elapsed"]]
eval(parse(text = lines), envir = globalenv())
cat(sprintf(paste0("COST script %s args %s arm %s | escapes %d holes %d ",
                   "past %d | evals %d in escape %d | secs %.1f\n"),
            script, paste(rest, collapse = " "), arm, cost$n_esc,
            cost$holes, cost$past, as.integer(cost$ev),
            as.integer(cost$ev_esc), proc.time()[["elapsed"]] - .nlb_t0))
