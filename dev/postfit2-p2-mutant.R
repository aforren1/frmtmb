# Punch round 2: each fix is pinned by a test. A mutant undoes one fix
# in the process (no install); the named test file must then fail, and
# it passes without the mutant (dev/postfit2-log/p2-*-lane.txt).
#   Rscript dev/postfit2-p2-mutant.R <mutant> <pkg> <file>
# Mutants:
#   byall    a gr(g, by = f) block is read by every row (round 1)
#   mmskip   an mm() block gets no plan: not held, no new level drawn
#   mmgv     ce_group_vars() without the mm() member columns (round 1)
#   rekey    re_formula left out of the boot reuse key (round 1)
#   nosmooth the smooth-level refusal removed (round 1)
args <- commandArgs(TRUE)
mutant <- args[1]
pkg <- args[2]
file <- args[3]
.libPaths(c("C:/Users/adf44/source/r/wt-postfit2-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true")
suppressMessages({
  library(testthat)
  library(frmtmb)
  if (pkg != "frmtmb") library(pkg, character.only = TRUE)
})
swap <- function(e, from, to, count) {
  if (identical(e, from)) {
    count$n <- count$n + 1L
    return(to)
  }
  if (is.call(e)) {
    for (i in seq_along(e)) {
      if (!is.null(e[[i]]) && !identical(e[[i]], quote(expr = ))) {
        e[[i]] <- swap(e[[i]], from, to, count)
      }
    }
  }
  e
}
put <- function(name, fn) {
  ns <- asNamespace("frmtmb")
  environment(fn) <- ns
  utils::assignInNamespace(name, fn, ns = "frmtmb")
  # an importing extension holds its own copy of an exported binding
  if (isNamespaceLoaded("frmtmb.sample")) {
    imp <- parent.env(asNamespace("frmtmb.sample"))
    if (exists(name, envir = imp, inherits = FALSE)) {
      unlockBinding(name, imp)
      assign(name, fn, envir = imp)
    }
  }
}
mutate <- function(name, from, to) {
  fn <- get(name, envir = asNamespace("frmtmb"))
  count <- new.env()
  count$n <- 0L
  body(fn) <- swap(body(fn), from, to, count)
  # a mutant that changed nothing would pass for the wrong reason
  if (count$n < 1L) stop("mutant ", mutant, ": no replacement in ", name)
  cat("mutant", mutant, ":", count$n, "replacement(s) in", name, "\n")
  put(name, fn)
}
switch(mutant,
  byall = put("ce_by_reads", function(bk, nd) rep(TRUE, nrow(nd))),
  mmskip = mutate("ce_level_plan", quote(length(mmv) > 0L), quote(FALSE)),
  mmgv = mutate("ce_group_vars", quote(unique(c(g, mm))), quote(g)),
  rekey = mutate("ce_boot_key", quote(ce_re_key(re_form)), quote("NA")),
  nosmooth = put("ce_check_smooth_levels", function(...) invisible(NULL)),
  none = NULL,
  stop("unknown mutant ", mutant))
dir <- if (pkg == "frmtmb") {
  "C:/Users/adf44/source/r/frmtmb-wt-postfit2/tests/testthat"
} else {
  file.path("C:/Users/adf44/source/r/frmtmb-wt-postfit2/extensions", pkg,
            "tests/testthat")
}
r <- test_file(file.path(dir, file), package = pkg, env = test_env(pkg),
               reporter = ListReporter$new(), stop_on_failure = FALSE)
df <- as.data.frame(r)
cat(sprintf("RESULT %s %s [%s]: failed=%d error=%d passed=%d\n", pkg, file,
            mutant, sum(df$failed), sum(df$error), sum(df$passed)))
print(df[df$failed > 0 | df$error, c("test", "failed", "error")])
