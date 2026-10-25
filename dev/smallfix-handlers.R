# Every error handler in frmtmb's R/ that can swallow the refusal for a
# part of the model this session cannot use (class
# frmtmb_extension_missing). Review finding 2 of
# dev/reviews/2026-10-24-smallfix.md: fixef() wrapped vcov() in
# tryCatch(error = function(e) NULL) and answered NA.
#
# The call graph of the installed namespace (codetools::findGlobals) is
# walked backwards from the functions that raise the refusal, to find
# every function that can reach one. A tryCatch() or try() in R/ is
# listed when its expression calls such a function and its error
# handler does not re-signal the refusal (keep_missing(), or
# a handler that names the class).
#   Rscript dev/smallfix-handlers.R <lib;lib>
# Log: dev/smallfix-out/handlers.txt
a <- commandArgs(trailingOnly = TRUE)
.libPaths(c(strsplit(a[1], ";", fixed = TRUE)[[1]],
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
ns <- asNamespace("frmtmb")
fns <- Filter(function(n) is.function(get(n, ns)), ls(ns, all.names = TRUE))
calls <- lapply(fns, function(n) {
  intersect(codetools::findGlobals(get(n, ns), merge = FALSE)$functions, fns)
})
names(calls) <- fns
roots <- c("require_extension", "require_model_packages",
           "require_structure_packages", "stop_extension_parts", "cs_entry")
reach <- roots
repeat {
  more <- names(Filter(function(cl) any(cl %in% reach), calls))
  more <- setdiff(more, reach)
  if (!length(more)) break
  reach <- c(reach, more)
}
# S3 generics dispatch to methods the graph cannot see by name: a call
# to vcov() reaches vcov.frmtmb_fit()
gen_of <- unique(sub("[.]frmtmb_fit$", "", grep("[.]frmtmb_fit$", reach,
                                                value = TRUE)))
reach_names <- unique(c(reach, gen_of))
cat("functions that can reach the refusal:", length(reach), "\n")

found <- list()
for (f in list.files("R", "[.]R$", full.names = TRUE)) {
  ex <- parse(f, keep.source = TRUE)
  sd <- getParseData(ex)
  walk <- function(e, owner) {
    if (!is.call(e)) return(invisible())
    h <- e[[1L]]
    hn <- if (is.name(h)) as.character(h) else ""
    if (hn %in% c("tryCatch", "try")) {
      nms <- names(e) %||% rep("", length(e))
      body <- e[[2L]]
      handler <- if (hn == "try") "try()" else
        paste(deparse(e[[which(nms == "error")[1L]]]), collapse = " ")
      if (hn == "try" || any(nms == "error")) {
        used <- tryCatch(all.names(body), error = function(err) character())
        hit <- intersect(used, reach_names)
        keeps <- grepl("keep_missing|frmtmb_extension_missing",
                       handler)
        if (length(hit) && !keeps) {
          found[[length(found) + 1L]] <<- data.frame(
            file = f, fun = owner, calls = paste(head(hit, 4), collapse = ","),
            handler = substr(handler, 1, 60))
        }
      }
    }
    for (i in seq_along(e)) {
      ai <- tryCatch(if (is.call(e[[i]])) e[[i]], error = function(err) NULL)
      if (!is.null(ai)) walk(ai, owner)
    }
  }
  for (k in seq_along(ex)) {
    e <- ex[[k]]
    owner <- if (is.call(e) && identical(e[[1L]], as.name("<-"))) {
      deparse(e[[2L]])
    } else "?"
    walk(e, owner)
  }
}
`%||%` <- function(a, b) if (is.null(a)) b else a
out <- do.call(rbind, found)
cat("handlers that can swallow it:", nrow(out %||% data.frame()), "\n")
if (!is.null(out)) print(out, right = FALSE, row.names = FALSE)
