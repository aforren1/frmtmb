# Why does the mocked METHOD run? The test file's comment says it is
# the caller-environment leg of UseMethod(). That explanation predicts
# the mock is invisible to a caller outside the namespace, and the
# first run showed it is NOT. Settle which mechanism it is, because a
# wrong reason in a comment is what the next test author will rely on.
.libPaths(c("C:/Users/adf44/source/r/sgrev-lib",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
q <- function(e) suppressWarnings(suppressMessages(e))
q(library(frmtmb.sample)); q(library(testthat))
arg <- commandArgs(trailingOnly = TRUE)
if (length(arg) && identical(arg[[1]], "owner")) q(library(rstantools))

ns <- asNamespace("frmtmb.sample")
tbl <- get(".__S3MethodsTable__.", envir = ns, inherits = FALSE)
rt <- if (isNamespaceLoaded("rstantools"))
  get(".__S3MethodsTable__.", envir = asNamespace("rstantools"),
      inherits = FALSE) else NULL
fd <- structure(list(), class = "frmtmb_draws")
gen <- get("log_lik", envir = globalenv())
cat("rstantools loaded: ", isNamespaceLoaded("rstantools"),
    "   generic reached from globalenv is defined in: ",
    environmentName(environment(gen)), "\n", sep = "")

before_tbl <- get("log_lik.frmtmb_draws", envir = tbl, inherits = FALSE)
before_ns <- get("log_lik.frmtmb_draws", envir = ns, inherits = FALSE)
local_mocked_bindings(
  log_lik.frmtmb_draws = function(object, ...) "MOCK RAN",
  .package = "frmtmb.sample")
after_ns <- get("log_lik.frmtmb_draws", envir = ns, inherits = FALSE)
after_tbl <- get("log_lik.frmtmb_draws", envir = tbl, inherits = FALSE)
cat("namespace binding changed by the mock: ",
    !identical(before_ns, after_ns), "\n")
cat("frmtmb.sample's method TABLE entry changed too: ",
    !identical(before_tbl, after_tbl), "\n")
if (!is.null(rt)) {
  rt_e <- tryCatch(get("log_lik.frmtmb_draws", envir = rt,
                       inherits = FALSE), error = function(e) NULL)
  cat("rstantools' table entry is the mock: ",
      identical(rt_e, after_ns), "\n")
}
g <- function(x) log_lik(x)   # caller is globalenv
cat("called from globalenv: ",
    tryCatch(g(fd), error = function(e)
      substr(conditionMessage(e), 1, 70)), "\n")
h <- function(x) log_lik(x); environment(h) <- ns
cat("called from inside the namespace: ",
    tryCatch(h(fd), error = function(e)
      substr(conditionMessage(e), 1, 70)), "\n")
cat("DONE\n")
