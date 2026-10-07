# Reviewer of the 0.69.0 consolidation: where does fuzz seed 20379118
# (Gamma, REML, ar1, mo(xo) * z, shape ~ 1 + z, weights) die on the
# merged build under OpenBLAS? Fits the spec with variants of the
# control and traces the optima stages.
#
#   Rscript dev/relrev069-fuzzdiag.R <lib or "base"> [base lib]
a <- commandArgs(trailingOnly = TRUE)
base <- if (length(a) >= 2) a[2] else "C:/Users/adf44/source/r/rellib-r7"
.libPaths(c(if (!identical(a[1], "base")) a[1], base,
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(testthat); library(frmtmb)})
root <- "C:/Users/adf44/source/r/frmtmb-wt-release"
env <- new.env(parent = asNamespace("frmtmb"))
sys.source(file.path(root, "tests/testthat/helper-fuzz.R"), envir = env)
cat("lib:", as.character(packageVersion("frmtmb")),
    dirname(find.package("frmtmb")), "| R:", R.home(), "| threads:",
    Sys.getenv("OPENBLAS_NUM_THREADS"), "\n")
sp <- list(family = "Gamma", aterm = "weights", re = "ar1",
           special = "mo_int", dpar = "dpar_x", mode = "reml",
           op = "confint", seed = 20379118L)
d <- env$fuzz_data(sp)
bform <- eval(parse(text = env$fuzz_bf_text(sp)))
ns <- asNamespace("frmtmb")
has <- function(f) exists(f, envir = ns, inherits = FALSE)
for (f in c("mo_search", "escape_stationary", "nlminb_best_par",
            "optimize_obj")) {
  if (has(f)) {
    trace(f, where = ns, print = FALSE,
          tracer = bquote(cat("  ENTER", .(f), "\n")),
          exit = bquote(cat("  EXIT", .(f), "\n")))
  }
}
if (has("nlminb_best_par")) {
  trace("nlminb_best_par", where = ns, print = FALSE,
        exit = quote(cat("    best_par: rejected_last =",
                         isTRUE(returnValue()$rejected_last),
                         "conv", returnValue()$convergence,
                         "obj", format(returnValue()$objective, digits = 12),
                         "msg", returnValue()$message, "\n")))
}
one <- function(lab, ctl) {
  cat("==", lab, "\n")
  r <- tryCatch(withCallingHandlers(
    frm(bform, data = d, REML = TRUE, control = ctl),
    warning = function(w) {
      cat("  WARN:", substr(conditionMessage(w), 1, 160), "\n")
      invokeRestart("muffleWarning")
    }, message = function(m) {
      cat("  MSG:", substr(conditionMessage(m), 1, 160), "\n")
      invokeRestart("muffleMessage")
    }), error = function(e) e)
  if (inherits(r, "error")) {
    cat("  ERROR:", substr(conditionMessage(r), 1, 200), "\n")
  } else {
    cat("  OK code", r$opt$convergence, "objective",
        sprintf("%.10f", r$opt$objective), "\n")
    print(round(r$opt$par, 4))
    if (!is.null(r$opt$mo_search)) str(r$opt$mo_search)
  }
  invisible(r)
}
one("default", frmtmb_control())
if ("mo_search" %in% names(formals(frmtmb_control))) {
  one("mo_search = FALSE", frmtmb_control(mo_search = FALSE))
}
one("optimizer optim", frmtmb_control(optimizer = "optim"))

# What the restart starts from: the objective and gradient at the start
# point of every optimizer run (fn and gr are optimize_obj's counters).
trace("run_optimizer", where = ns, print = FALSE,
      tracer = quote({
        v <- tryCatch(fn(par), error = function(e) NA)
        g <- tryCatch(gr(par), error = function(e) NA)
        cat("  run_optimizer start: fn", format(v, digits = 12),
            "| non-finite gr", sum(!is.finite(g)), "of", length(g), "\n")
      }))
one("default, start of each run traced", frmtmb_control())
untrace("run_optimizer", where = ns)
# Arm: the 0.68.1 softmax chart for mo(), everything else 0.69.0
if (has("mo_simplex")) {
  ms <- get("mo_simplex", ns)
  ms2 <- ms
  body(ms2) <- bquote({
    if (is.null(chart)) chart <- "softmax"
    .(body(ms))
  })
  unlockBinding("mo_simplex", ns)
  assign("mo_simplex", ms2, envir = ns)
  one("softmax chart (0.68.1's), mo_search = FALSE",
      frmtmb_control(mo_search = FALSE))
  assign("mo_simplex", ms, envir = ns)
}
# Arm: nlminb_best_par as identity (0.68.1's par), stereographic chart
if (has("nlminb_best_par")) {
  untrace("nlminb_best_par", where = ns)
  nb <- get("nlminb_best_par", ns)
  unlockBinding("nlminb_best_par", ns)
  assign("nlminb_best_par", function(res, fnw) res, envir = ns)
  one("nlminb_best_par identity", frmtmb_control())
  assign("nlminb_best_par", nb, envir = ns)
}
