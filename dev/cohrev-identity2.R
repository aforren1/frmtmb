## Same as cohrev-identity.R, but deparsing WITHOUT srcrefs. The two
## libraries were installed with different keep.source settings, so the
## first pass compared comments and whitespace rather than code.
args <- commandArgs(trailingOnly = TRUE)
.libPaths(c(args[1L], "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb.coupling))
ns <- asNamespace("frmtmb.coupling")
nms <- sort(ls(ns, all.names = TRUE))
strip <- function(o) {
  if (is.function(o)) {
    attr(o, "srcref") <- NULL
    environment(o) <- globalenv()
    body(o) <- removeSource(body(o))
  }
  o
}
h <- vapply(nms, function(n) {
  paste(deparse(strip(get(n, envir = ns)),
                control = c("keepInteger", "keepNA", "niceNames")),
        collapse = "\n")
}, character(1))
bc <- vapply(nms, function(n) {
  o <- get(n, envir = ns)
  if (is.function(o)) {
    class(.Internal(bodyCode(o)))[1L]
  } else "not a function"
}, character(1))
saveRDS(list(nms = nms, h = h, bc = bc), args[2L])
cat("objects:", length(nms), " byte-compiled:",
    sum(bc == "bytecode"), "\n")
