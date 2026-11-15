# Item 3164: the positive control of run-session.R's extension-state
# snapshot. Takes ext_state() and ext_state_leaks() out of the runner's
# source, loads frmtmb.ode and frmtmb.eam, and changes their session
# state the way a fit would; the change must be printed, and a second
# snapshot with nothing changed must print nothing.
#   Rscript dev/hygiene-extstate.R
.libPaths(c("C:/Users/adf44/source/r/wt-hygiene-lib",
            "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
ex <- parse("dev/release/run-session.R")[[1]][[2]]
env <- new.env()
for (e in as.list(ex)[-1]) {
  if (is.call(e) && identical(e[[1]], as.name("<-")) &&
      as.character(e[[2]]) %in% c("ext_state", "ext_state_leaks")) {
    eval(e, env)
  }
}
stopifnot(is.function(env$ext_state), is.function(env$ext_state_leaks))
for (q in c("frmtmb.ode", "frmtmb.eam", "frmtmb.sample")) loadNamespace(q)
s0 <- env$ext_state()
cat("environments read:", paste(names(s0), collapse = " "), "\n")
cat("-- unchanged:\n")
env$ext_state_leaks("none.R", s0, env$ext_state())
assign("probe-key", TRUE, envir = asNamespace("frmtmb.ode")$ode_warned)
cat("-- after setting a key of ode_warned:\n")
env$ext_state_leaks("probe.R", s0, env$ext_state())
rm("probe-key", envir = asNamespace("frmtmb.ode")$ode_warned)
