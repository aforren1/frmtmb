# Item 5. The active-binding mocking trap: the mechanism, the old
# spelling seen failing, the new spelling seen doing the right thing
# for the right reason, and the size of the exposed surface.
.libPaths(c("C:/Users/adf44/source/r/sgrev-lib",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
q <- function(e) suppressWarnings(suppressMessages(e))

cat("== 1. the mechanism, on a throwaway environment ==\n")
e <- new.env()
calls <- 0
makeActiveBinding("z", function() { calls <<- calls + 1; "OWNER" }, e)
cat("get:      ", get("z", envir = e), "  binding calls so far ",
    calls, "\n", sep = "")
r <- tryCatch({ assign("z", function(a) a, envir = e); "assigned" },
              error = function(err) conditionMessage(err))
cat("assign:   ", r, "  binding calls now ", calls, "\n", sep = "")
# the same with a zero-argument binding function, which is what
# frm_bind_generic() installs
e2 <- new.env()
makeActiveBinding("z2", function() "OWNER", e2)
r2 <- tryCatch({ assign("z2", 1, envir = e2); "assigned" },
               error = function(err) conditionMessage(err))
cat("assign to a NO-ARG active binding: ", r2, "\n", sep = "")
cat("  (an assignment CALLS the binding function with the value; a",
    " function of no arguments therefore errors)\n\n")

cat("== 2. the size of the surface ==\n")
q(library(frmtmb)); q(library(frmtmb.sample))
for (p in c("frmtmb", "frmtmb.sample")) {
  ns <- asNamespace(p)
  nmz <- ls(ns, all.names = TRUE)
  act <- nmz[vapply(nmz, function(n)
    tryCatch(bindingIsActive(n, ns), error = function(e) FALSE), NA)]
  cat(p, ": ", length(act), " active bindings\n", sep = "")
  cat("   ", paste(sort(act), collapse = " "), "\n")
}

cat("\n== 3. the old mock spelling, seen failing ==\n")
q(library(testthat))
old <- tryCatch({
  testthat::local_mocked_bindings(
    log_lik = function(x, ...) 1, .package = "frmtmb.sample")
  "mock installed"
}, error = function(err) conditionMessage(err))
cat("local_mocked_bindings(log_lik = ) on the FIX build: ", old, "\n")

cat("\n== 4. the new spelling: is the mocked METHOD what runs? ==\n")
# The test's claim is that UseMethod() finds a method in the CALLING
# environment, which is frmtmb.sample's namespace, before the method
# table. Assert the claim directly rather than trusting a green test.
ns <- asNamespace("frmtmb.sample")
tbl <- get(".__S3MethodsTable__.", envir = ns, inherits = FALSE)
cat("log_lik.frmtmb_draws exists in the namespace: ",
    exists("log_lik.frmtmb_draws", envir = ns, inherits = FALSE), "\n")
cat("log_lik.frmtmb_draws exists in this table:    ",
    exists("log_lik.frmtmb_draws", envir = tbl, inherits = FALSE), "\n")
cat("the namespace copy and the table copy are the same object: ",
    identical(get("log_lik.frmtmb_draws", envir = ns, inherits = FALSE),
              get("log_lik.frmtmb_draws", envir = tbl,
                  inherits = FALSE)), "\n")
fd <- structure(list(), class = "frmtmb_draws")
probe <- function() {
  # a caller INSIDE frmtmb.sample's namespace, the way loo_matrix() is
  f <- function(x) log_lik(x)
  environment(f) <- ns
  f(fd)
}
testthat::local_mocked_bindings(
  log_lik.frmtmb_draws = function(object, ...) "MOCK RAN",
  .package = "frmtmb.sample")
cat("mocked method, called from inside the namespace: ",
    tryCatch(probe(), error = function(e) conditionMessage(e)), "\n")
g <- function(x) log_lik(x)   # caller is globalenv
cat("mocked method, called from the GLOBAL environment: ",
    tryCatch(g(fd), error = function(e) substr(conditionMessage(e), 1, 60)),
    "\n")
cat("  (the caller-environment leg is what makes the mock work, so a",
    " mock of this shape only covers callers in the same namespace)\n")
cat("DONE\n")
