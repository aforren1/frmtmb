# Does UseMethod's "look in the calling environment" leg reach the
# SEARCH PATH, or only the caller's own frame? The answer decides
# whether an owner that EXPORTS a `.default` can capture dispatch from
# a rival generic. Measured, not read.
.libPaths(c("C:/Users/adf44/source/r/sgrev-lib",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(loo))
cat("loo exports loo_moment_match.default: ",
    "loo_moment_match.default" %in%
      parseNamespaceFile("loo", dirname(system.file(package = "loo")))$exports,
    "\n")
cat("visible from globalenv: ",
    exists("loo_moment_match.default", envir = globalenv()), "\n")

# a rival generic in a namespace of its own, with NO method table
e <- new.env(parent = globalenv())
assign(".__S3MethodsTable__.", new.env(), envir = e)
rival <- function(x, ...) UseMethod("loo_moment_match")
environment(rival) <- e
obj <- structure(list(), class = "zzz")
cat("rival generic, class zzz, caller = globalenv: ")
print(tryCatch(rival(obj), error = function(err) conditionMessage(err)))

# the same call from inside a function whose environment is globalenv
f <- function() rival(obj)
cat("rival generic, caller = a function frame: ")
print(tryCatch(f(), error = function(err) conditionMessage(err)))

# and with a .default defined in the CALLER's own frame
g <- function() {
  loo_moment_match.default <- function(x, ...) "LOCAL DEFAULT RAN"
  rival(obj)
}
cat("with a .default in the caller's own frame: ")
print(tryCatch(g(), error = function(err) conditionMessage(err)))
