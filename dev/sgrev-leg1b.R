# Second half of the caller-env question, on a plain base generic so
# nothing about the rival construction can be blamed.
.libPaths(c("C:/Users/adf44/source/r/sgrev-lib",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
obj <- structure(1, class = "sgrevcls")

# 1. method at top level, generic called from top level
format.sgrevcls <- function(x, ...) "TOPLEVEL"
cat("1 top-level method, top-level call: ", format(obj), "\n")

# 2. method at top level, generic called from inside a function whose
#    frame is NOT globalenv
f <- function(z) format(z)
cat("2 top-level method, called inside a function: ", f(obj), "\n")

# 3. method attached on the search path (a package-like environment
#    spliced in), generic called from top level
env <- new.env()
assign("format.sgrevcls2", function(x, ...) "SEARCHPATH", envir = env)
attach(env, name = "sgrevtest", warn.conflicts = FALSE)
obj2 <- structure(1, class = "sgrevcls2")
cat("3 method on the search path, top-level call: ",
    tryCatch(format(obj2), error = function(e) conditionMessage(e)), "\n")
cat("3 method on the search path, called in a function: ",
    tryCatch(f(obj2), error = function(e) conditionMessage(e)), "\n")
detach("sgrevtest")
