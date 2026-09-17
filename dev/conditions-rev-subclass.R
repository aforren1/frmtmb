# Reviewer, lane wt-conditions: where the run-time subclass rule gives a
# wrong or missing class, what frm_stop() does with a condition object,
# and what a raise costs against stop().
#   Rscript dev/conditions-rev-subclass.R
.libPaths(c("C:/Users/adf44/source/r/conditions-lib",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
rscript <- "C:/Program Files/R/R-4.6.1/bin/Rscript.exe"
cls <- function(expr) {
  e <- tryCatch(expr, condition = identity)
  paste(class(e), collapse = "/")
}
show <- function(label, x) cat(sprintf("%-46s %s\n", label, x))

cat("== 1. extension loaded, frmtmb NOT attached (fresh process)\n")
code <- paste0(
  ".libPaths(c('C:/Users/adf44/source/r/conditions-lib', ",
  "'C:/Users/adf44/source/r/pinlib', ",
  "'C:/Users/adf44/AppData/Local/R/win-library/4.6'));",
  "e <- tryCatch(frmtmb.eam::wiener()$links[['ndt']]$linkinv(0), ",
  "error = identity);",
  "cat(class(e), '|', conditionMessage(e), '|', ",
  "'frmtmb' %in% .packages(), '\\n')")
cat(system2(rscript, c("-e", shQuote(code)), stdout = TRUE,
            stderr = TRUE), sep = "\n")

suppressPackageStartupMessages({library(frmtmb); library(frmtmb.eam)})
eam <- asNamespace("frmtmb.eam")
in_eam <- function(f) { environment(f) <- eam; f }

cat("\n== 2. helper reached from inside an extension function\n")
show("direct", cls(in_eam(function() frm_stop("x"))()))
show("anonymous closure in lapply",
     cls(in_eam(function() lapply(1, function(i) frm_stop("x")))()))
show("tryCatch handler",
     cls(in_eam(function() tryCatch(stop("a"),
                                    error = function(e) frm_stop("x")))()))
show("do.call(frm_stop, ...)",
     cls(in_eam(function() do.call(frm_stop, list("x")))()))
show("do.call(\"frm_stop\", ...)",
     cls(in_eam(function() do.call("frm_stop", list("x")))()))
show("lapply(msgs, frm_stop)  [helper as a value]",
     cls(in_eam(function() lapply("x", frm_stop))()))
show("Map(frm_warning, msgs)  [helper as a value]",
     cls(in_eam(function() Map(frm_warning, "x"))()))
show("vapply(msgs, frm_message, NULL)",
     cls(in_eam(function() vapply("x", frm_message, NULL))()))
show("Reduce(function(a, b) frm_stop(b), ...)",
     cls(in_eam(function() Reduce(function(a, b) frm_stop(b), 1:2))()))
show("local({ frm_stop() })",
     cls(in_eam(function() local(frm_stop("x")))()))
show("eval(quote(frm_stop()), new.env())",
     cls(in_eam(function() eval(quote(frm_stop("x")), new.env()))()))
show("eval(quote(frm_stop()), globalenv())",
     cls(in_eam(function() eval(quote(frm_stop("x")), globalenv()))()))
pe <- new.env(); assign(".packageName", "frmtmb.eam", envir = pe)
f <- function() frm_stop("x"); environment(f) <- pe
show("env carrying .packageName (sys.source style)", cls(f()))
# a core helper called by an extension raises the plain class, by design
show("core validator called from extension",
     cls(in_eam(function() frmtmb:::frm_condition_class("error",
                                                        environment()))()))

cat("\n== 3. frm_stop() given a condition object, as stop(e) accepts\n")
e0 <- simpleError("boom", call = quote(f(1)))
a <- tryCatch(stop(e0), error = identity)
b <- tryCatch(frm_stop(e0), error = identity)
show("stop(e) message", conditionMessage(a))
show("frm_stop(e) message", conditionMessage(b))
w0 <- simpleWarning("careful")
show("warning(w) message",
     conditionMessage(tryCatch(warning(w0), warning = identity)))
show("frm_warning(w) message",
     conditionMessage(tryCatch(frm_warning(w0), warning = identity)))

cat("\n== 4. cost of one raise and catch, stop() against frm_stop()\n")
# interleaved arms in one process, blocks grown past 1.2 s, minimum of
# 7 rounds, with a control arm made from the same code that must be 1.0
fs <- function() stop("bad x: ", 3L)
ff <- function() frm_stop("bad x: ", 3L)
fe <- in_eam(function() frm_stop("bad x: ", 3L))
fw <- function() warning("w")
fwf <- function() frm_warning("w")
arm <- function(g, n) {
  t0 <- proc.time()[[3L]]
  for (i in seq_len(n)) tryCatch(g(), condition = function(e) NULL)
  proc.time()[[3L]] - t0
}
n <- 2000L
while (arm(fs, n) < 1.2) n <- n * 2L
cat("iterations per block", n, "\n")
tab <- replicate(7L, c(stop = arm(fs, n), stop_ctl = arm(fs, n),
                       frm_stop = arm(ff, n), frm_stop_ext = arm(fe, n),
                       warning = arm(fw, n), frm_warning = arm(fwf, n)))
mins <- apply(tab, 1L, min)
print(round(mins / n * 1e6, 1))   # microseconds per raise
cat(sprintf("ratio control %.3f  frm_stop/stop %.3f  ext/stop %.3f  frm_warning/warning %.3f\n",
            mins[["stop_ctl"]] / mins[["stop"]],
            mins[["frm_stop"]] / mins[["stop"]],
            mins[["frm_stop_ext"]] / mins[["stop"]],
            mins[["frm_warning"]] / mins[["warning"]]))
