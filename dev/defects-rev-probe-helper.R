# Reviewer of lane defects: probe helper, loaded AFTER helper-brms-suite.R
# in a scratch copy of the test directory (named helper-zz-probe.R there).
# For every assertion that HELD, it re-evaluates each all()/any()/anyNA()
# operand inside it and logs the operand's length and the hollow flags,
# so that vacuous forms the harness does not recognize are listed.
.defrev_log <- Sys.getenv("DEFREV_PROBE_LOG")
.defrev_subs <- function(e, acc = list()) {
  if (!is.call(e)) return(acc)
  h <- e[[1]]
  if (is.name(h) && as.character(h) %in% c("all", "any", "anyNA")) {
    acc <- c(acc, list(e))
  }
  for (a in as.list(e)[-1]) {
    if (is.name(a) && !nzchar(as.character(a))) next
    acc <- .defrev_subs(a, acc)
  }
  acc
}
brms_port_orig <- brms_port
brms_port <- function(id, verdict, reason, code) {
  call <- substitute(code)
  env <- parent.frame()
  res <- eval(as.call(list(quote(brms_port_orig), id, verdict, reason,
                           call)), env)
  subs <- .defrev_subs(call)
  for (s in subs) {
    x <- tryCatch(suppressWarnings(suppressMessages(eval(s[[2]], env))),
                  error = function(e) structure(NA, err = TRUE))
    cat(paste(c(id, verdict, isTRUE(res$held), isTRUE(res$vacuous),
                as.character(s[[1]]), length(x),
                gsub("[\t\n]", " ", substr(deparse1(call), 1, 150))),
              collapse = "\t"), "\n", sep = "", file = .defrev_log,
        append = TRUE)
  }
  invisible(res)
}
