# Fix round 2, item 9: installed code against the tree, all nine packages; dates.
.libPaths(c("C:/Users/adf44/source/r/wt-nlfix-lib2", "C:/Users/adf44/source/r/rellib-r16", "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
L <- "C:/Users/adf44/source/r/wt-nlfix-lib2"
for (p in c("frmtmb", paste0("frmtmb.", c("coupling", "eam", "latent", "learn", "ode", "sample", "spatial", "survival")))) {
  src <- if (p == "frmtmb") "R" else file.path("extensions", p, "R")
  suppressPackageStartupMessages(library(p, character.only = TRUE))
  ns <- asNamespace(p); e <- new.env(parent = ns)
  fs <- list.files(src, "[.]R$", full.names = TRUE)
  res <- tryCatch({ for (f in fs) sys.source(f, envir = e, keep.source = FALSE); "ok" }, error = function(er) conditionMessage(er))
  nm <- ls(e, all.names = TRUE); isf <- vapply(nm, function(x) is.function(e[[x]]), TRUE)
  dif <- character(0); miss <- character(0)
  for (x in nm[isf]) { g <- get0(x, envir = ns, inherits = FALSE); if (is.null(g)) miss <- c(miss, x) else if (!identical(paste(deparse(body(e[[x]])), collapse = ""), paste(deparse(body(g)), collapse = ""))) dif <- c(dif, x) }
  built <- utils::packageDescription(p, lib.loc = L)$Built
  srct <- max(file.info(fs)$mtime); inst <- file.info(file.path(L, p, "R", paste0(p, ".rdb")))$mtime
  cat(sprintf("%-17s from %s | source: %s, functions %d, differ %d, missing %d | newest source %s, installed %s, source older: %s\n", p, dirname(find.package(p)), res, sum(isf), length(dif), length(miss), format(srct, "%m-%d %H:%M:%S"), format(inst, "%m-%d %H:%M:%S"), srct < inst))
  if (length(dif)) print(head(dif))
}
