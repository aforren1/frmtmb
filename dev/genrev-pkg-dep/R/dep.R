# The owner of the generic as seen from INSIDE this namespace, which is
# what a third package's own code reaches.
genrev_owner <- function(g) {
  f <- tryCatch(get(g), error = function(e) NULL)
  if (!is.function(f)) return("<absent>")
  e <- environment(f)
  if (is.null(e)) return("base")
  n <- environmentName(topenv(e))
  if (nzchar(n)) n else "<anon>"
}

# Whether a method for `cls` is reachable through the generic this
# namespace sees.
genrev_lookup <- function(g, cls) {
  f <- tryCatch(get(g), error = function(e) NULL)
  if (!is.function(f)) return("<absent>")
  genv <- environment(f)
  tbl <- tryCatch(get(".__S3MethodsTable__.", envir = genv),
                  error = function(e) NULL)
  nm <- paste0(g, ".", cls)
  if (!is.null(tbl) && exists(nm, envir = tbl, inherits = FALSE)) {
    m <- get(nm, envir = tbl)
    return(environmentName(topenv(environment(m))))
  }
  "<none>"
}
