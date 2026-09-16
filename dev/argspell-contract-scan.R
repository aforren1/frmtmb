## The scan itself, shared by dev/argspell-contract.R (which prints the
## table) and tests/testthat/test-arg-refusal.R (which asserts it is
## empty). One implementation, so the test cannot drift from the table.
##
## Sourced, not a package function: it has to run against the INSTALLED
## package in the test and against a working tree in the script.

# Walks a call tree without ever BINDING an element. `x[i, ]` holds the
# empty symbol, and assigning that to a variable raises "argument is
# missing, with no default"; comparing the one-element sublist does not
# force it. This walker is the third attempt at the problem and the
# first that survives base R.
frm_walk_calls <- function(e, fn) {
  if (!is.call(e)) return(invisible(NULL))
  fn(e)
  l <- as.list(e)
  for (k in seq_along(l)) {
    if (identical(l[k], list(quote(expr = )))) next
    frm_walk_calls(l[[k]], fn)
  }
  invisible(NULL)
}

# The head of a call, as a plain name: `f(...)`, `pkg::f(...)` and
# `pkg:::f(...)` all report `f`.
frm_call_head <- function(e) {
  h <- e[[1L]]
  if (is.name(h)) return(as.character(h))
  if (is.call(h) && is.name(h[[1L]]) &&
        as.character(h[[1L]]) %in% c("::", ":::") && length(h) == 3L) {
    return(as.character(h[[3L]]))
  }
  ""
}

# Every S3 generic the packages register a dots-refusing method on, with
# the union of the names those methods can already take.
frm_dots_generics <- function(pkgs) {
  out <- list()
  for (pkg in pkgs) {
    reg <- parseNamespaceFile(
      pkg, dirname(system.file(package = pkg)))$S3methods
    for (i in seq_len(nrow(reg))) {
      gen <- reg[i, 1]
      fun <- if (is.na(reg[i, 3])) paste(reg[i, 1], reg[i, 2], sep = ".")
      else reg[i, 3]
      o <- tryCatch(get(fun, envir = asNamespace(pkg)),
                    error = function(e) NULL)
      if (!is.function(o) || !"..." %in% names(formals(o))) next
      # a method that forwards its dots cannot be broken by an extra
      # name reaching it, so only the REFUSING ones matter here
      if (!any(grepl("frm_check_dots", deparse(body(o)), fixed = TRUE))) {
        next
      }
      out[[gen]] <- unique(c(out[[gen]],
                             setdiff(names(formals(o)), "...")))
    }
  }
  out
}

# `contract` is the generic-to-names table the package ships. A hit is a
# named argument that is neither a formal of our method nor in the
# table, passed at a real call site to that generic.
frm_contract_scan <- function(pkgs, scan_pkgs, contract = NULL) {
  gens <- frm_dots_generics(pkgs)
  if (is.null(contract)) contract <- list()
  rows <- list()
  for (pk in scan_pkgs) {
    ns <- asNamespace(pk)
    for (nm in ls(ns, all.names = TRUE)) {
      o <- tryCatch(get(nm, envir = ns), error = function(e) NULL)
      if (!is.function(o)) next
      b <- tryCatch(body(o), error = function(e) NULL)
      if (is.null(b) || !is.call(b)) next
      tryCatch(frm_walk_calls(b, function(e) {
        g <- frm_call_head(e)
        if (!nzchar(g) || is.null(gens[[g]])) return(invisible(NULL))
        # A call to `plot()` inside `plot.dendrogram()` dispatches on a
        # dendrogram's own parts, not on whatever the caller was handed,
        # so it says nothing about the contract our method must honor.
        if (startsWith(nm, paste0(g, "."))) return(invisible(NULL))
        nn <- names(e)
        if (is.null(nn)) return(invisible(NULL))
        nn <- nn[nzchar(nn)]
        bad <- setdiff(nn, c(gens[[g]], contract[[g]], "..."))
        for (a in bad) {
          rows[[length(rows) + 1L]] <<- data.frame(
            generic = g, arg = a, where = paste0(pk, "::", nm),
            stringsAsFactors = FALSE)
        }
        invisible(NULL)
      }), error = function(e) NULL)
    }
  }
  hits <- if (length(rows)) do.call(rbind, rows) else
    data.frame(generic = character(), arg = character(),
               where = character(), stringsAsFactors = FALSE)
  hits <- hits[!duplicated(hits[, c("generic", "arg")]), , drop = FALSE]
  list(n_generics = length(gens), generics = names(gens), hits = hits)
}
