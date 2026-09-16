## The exemption list, measured the way it has to be.
##
## The first measurement (dev/argspell-exempt-evidence.R) counted call
## sites that forward dots into a generic ANYWHERE, including sites that
## can never reach a frmtmb_fit: `emmeans::emm_basis.gamm` calling
## `emm_basis(object$gam, ...)` dispatches on a gam. That is weaker than
## the sentence around it claimed.
##
## What a guard would actually refuse is a call site passing a NAMED
## argument our method has no formal for. That is what this counts.
##
## Run: Rscript dev/argspell-exempt2.R
LIB <- "C:/Users/adf44/source/r/argspell-lib"
.libPaths(c(LIB,
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))

methods <- c(
  emm_basis = "emm_basis.frmtmb_fit",
  recover_data = "recover_data.frmtmb_fit",
  get_coef = "get_coef.frmtmb_fit",
  get_predict = "get_predict.frmtmb_fit",
  get_vcov = "get_vcov.frmtmb_fit",
  set_coef = "set_coef.frmtmb_fit",
  get_varcov = "get_varcov.frmtmb_fit",
  get_parameters = "get_parameters.frmtmb_fit",
  find_formula = "find_formula.frmtmb_fit",
  find_random = "find_random.frmtmb_fit",
  find_statistic = "find_statistic.frmtmb_fit",
  link_function = "link_function.frmtmb_fit",
  link_inverse = "link_inverse.frmtmb_fit")

pkgs <- c("emmeans", "insight", "marginaleffects")

# Every call in `pkg` to `gen`, with the names it passes explicitly.
call_names <- function(pkg, gen) {
  if (!requireNamespace(pkg, quietly = TRUE)) return(list())
  ns <- asNamespace(pkg)
  out <- list()
  for (nm in ls(ns, all.names = TRUE)) {
    o <- tryCatch(get(nm, envir = ns), error = function(e) NULL)
    if (!is.function(o)) next
    b <- tryCatch(body(o), error = function(e) NULL)
    if (is.null(b)) next
    walk <- function(e) {
      if (!is.call(e)) return(invisible(NULL))
      h <- e[[1L]]
      hn <- if (is.name(h)) as.character(h) else
        if (is.call(h) && identical(as.character(h[[1L]]), "::")) {
          as.character(h[[3L]])
        } else ""
      if (identical(hn, gen)) {
        nn <- names(e)
        nn <- nn[nzchar(nn)]
        out[[length(out) + 1L]] <<- list(fun = paste0(pkg, "::", nm),
                                         args = nn)
      }
      parts <- as.list(e)
      for (i in seq_along(parts)) {
        el <- tryCatch(parts[[i]], error = function(err) NULL)
        if (is.null(el)) next
        if (is.name(el) && !nzchar(as.character(el))) next
        walk(el)
      }
      invisible(NULL)
    }
    tryCatch(walk(b), error = function(e) NULL)
  }
  out
}

cat("---- GENERATED: exemption evidence, re-measured ----\n")
cat("A call site COUNTS when it passes a name the method has no formal",
    "for.\n\n")
rows <- list()
for (i in seq_along(methods)) {
  gen <- names(methods)[i]
  meth <- methods[[i]]
  o <- tryCatch(get(meth, envir = asNamespace("frmtmb")),
                error = function(e) NULL)
  fo <- if (is.null(o)) character() else setdiff(names(formals(o)), "...")
  hits <- 0L
  example <- ""
  for (p in pkgs) {
    for (cl in call_names(p, gen)) {
      bad <- setdiff(cl$args, fo)
      if (length(bad)) {
        hits <- hits + 1L
        if (!nzchar(example)) {
          example <- paste0(cl$fun, " passes ",
                            paste(bad, collapse = ", "))
        }
      }
    }
  }
  reads <- !is.null(o) && "..." %in% names(formals(o)) &&
    any(all.names(body(o)) %in% c("...", "..1", "...length", "...names"))
  rows[[length(rows) + 1L]] <- data.frame(
    method = meth, formals = paste(fo, collapse = ", "),
    sites = hits, reads_dots = reads, example = example,
    stringsAsFactors = FALSE)
}
tab <- do.call(rbind, rows)
tab <- tab[order(-tab$sites, tab$method), ]
for (i in seq_len(nrow(tab))) {
  cat(sprintf("%-28s sites %3d  reads_dots %-5s\n",
              tab$method[i], tab$sites[i], tab$reads_dots[i]))
  if (nzchar(tab$example[i])) cat("    e.g. ", tab$example[i], "\n")
}
cat(sprintf("\nexemptions with measured need : %d\n", sum(tab$sites > 0)))
cat(sprintf("exemptions inert (reads dots) : %d\n", sum(tab$reads_dots)))
cat(sprintf("exemptions with NO evidence   : %d\n",
            sum(tab$sites == 0 & !tab$reads_dots)))
cat("---- END GENERATED: exemption evidence ----\n")
