# Review of lane fitedge: one test file under one in-process mutant of the
# lane's build. Args: mutant, test file name. The mutant is applied to the
# namespace after library(frmtmb) and before test_file().
local({
.libPaths(c("C:/Users/adf44/source/r/wt-fitedge-lib",
            "C:/Users/adf44/source/r/rellib-r18",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win",
           NOT_CRAN = "true", FRMTMB_BRMS_FIT_TESTS = "true",
           FRMTMB_STAN_CACHE =
             "C:/Users/adf44/source/r/frmtmb-wt-fitedge/dev/stan-cache")
suppressMessages(library(testthat))
suppressMessages(library(frmtmb))
a <- commandArgs(trailingOnly = TRUE)
mut <- a[1]
f <- file.path("C:/Users/adf44/source/r/frmtmb-wt-fitedge/tests/testthat", a[2])
cat("lib: frmtmb", as.character(utils::packageVersion("frmtmb")),
    dirname(find.package("frmtmb")), "mutant", mut, "\n")
ns <- asNamespace("frmtmb")
put <- function(nm, val) {
  unlockBinding(nm, ns)
  assign(nm, val, envir = ns)
  lockBinding(nm, ns)
}
# text substitution in one function; stops unless it matched exactly once
sub1 <- function(nm, from, to) {
  fn <- get(nm, envir = ns)
  txt <- paste(deparse(fn, width.cutoff = 500L), collapse = "\n")
  k <- lengths(regmatches(txt, gregexpr(from, txt, fixed = TRUE)))
  if (k != 1L) stop("mutant ", mut, ": '", from, "' matched ", k, " times")
  g <- eval(parse(text = sub(from, to, txt, fixed = TRUE)))
  environment(g) <- environment(fn)
  put(nm, g)
}
switch(mut,
  control = NULL,
  # item 1: no body start at all (0.79.1's path)
  nogn = put("gp_body_start", function(...) NULL),
  # item 1: one Gauss-Newton step, not two
  gn1 = sub1("gp_body_start", "for (step in 1:2)", "for (step in 1:1)"),
  # item 1: the body's other coefficients left out of the linearization
  noothers = sub1("gp_body_start", "for (k in seq_len(NROW(others)))",
                  "for (k in integer(0))"),
  # item 3: the condition gate at 1e12, and lifted
  kappa12 = put("gp_marginal_kappa_max", 1e12),
  nokappa = put("gp_marginal_kappa_max", Inf),
  # item 3: the value check lifted
  novalue = put("gp_marginal_value_tol", Inf),
  # item 4: the refusal back on every mi() response
  miback = sub1("latres_refuse_laplace",
                "length(mi_ml_frame(frame)[[\"mi_map\"]][[rn]][[\"rows\"]])",
                "TRUE"),
  # item 4: the refusal never fires
  minever = sub1("latres_refuse_laplace",
                 "length(mi_ml_frame(frame)[[\"mi_map\"]][[rn]][[\"rows\"]])",
                 "FALSE"),
  # item 2: no ridge pass
  noridge = put("nl_all_bodies", function(frame) list()),
  # item 2: the ridge pass from the kept values
  ridgekept = sub1("nl_fit_off_saddle", "ridge <- length(fresh) > 0L",
                   "ridge <- FALSE"),
  stop("unknown mutant ", mut))
res <- test_file(f, package = "frmtmb", env = testthat::test_env("frmtmb"),
                 reporter = "silent")
for (t in res) {
  for (x in t$results) {
    k <- if (inherits(x, "expectation_failure")) "FAIL" else
      if (inherits(x, "expectation_error")) "ERROR" else
        if (inherits(x, "expectation_warning")) "WARN" else NA
    if (!is.na(k)) {
      cat("DETAIL ", k, " [", t$test, "] ",
          substr(gsub("\n", " | ", conditionMessage(x)), 1, 300), "\n",
          sep = "")
    }
  }
}
r <- as.data.frame(res)
cat("RESULT ", mut, " ", basename(f), " pass=", sum(r$passed), " fail=",
    sum(r$failed), " err=", sum(r$error), " skip=", sum(r$skipped),
    " warn=", sum(r$warning), "\n", sep = "")
})
