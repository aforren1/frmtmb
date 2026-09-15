zap <- function(x, ...) UseMethod("zap")
zap.mine <- function(x, ...) "ADOPTER-METHOD"

.onLoad <- function(libname, pkgname) {
  ns <- asNamespace(pkgname)
  fallback <- get("zap", envir = ns, inherits = FALSE)
  rm(list = "zap", envir = ns)
  makeActiveBinding("zap", function() {
    if (isNamespaceLoaded("genrevown")) {
      g <- tryCatch(getExportedValue("genrevown", "zap"),
                    error = function(e) NULL)
      if (is.function(g)) return(g)
    }
    fallback
  }, ns)
  invisible()
}
