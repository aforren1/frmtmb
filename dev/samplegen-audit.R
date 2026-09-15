# Ownership audit for the 28 generics frmtmb.sample DEFINES.
#
# Who owns a name is decided in two steps, and grep is used for
# neither. parseNamespaceFile() over every installed package gives the
# candidates that EXPORT the name (nlme's multi-name export block is
# why this is not a grep). A package that exports a name need not
# DEFINE it, so each candidate is then loaded and the environment of
# its exported function says who defined it. That is the question that
# decides the target, because an S3 method registers into the table of
# the namespace where the GENERIC was defined.
#
#   Rscript dev/samplegen-audit.R
RR <- "C:/Users/adf44/source/r/rellib-r3"
UL <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
SYS <- file.path(R.home(), "library")
.libPaths(c(RR, "C:/Users/adf44/source/r/pinlib", UL, SYS))

sample_ns <- parseNamespaceFile("frmtmb.sample", RR)
suppressMessages(loadNamespace("frmtmb.sample"))
sns <- asNamespace("frmtmb.sample")
is_gen <- function(f) {
  is.function(f) && !is.null(body(f)) &&
    any(grepl("UseMethod", deparse(body(f)), fixed = TRUE))
}
home_of <- function(f) environmentName(topenv(environment(f)))
exps <- sample_ns$exports
gens <- Filter(function(n) is_gen(get0(n, envir = sns)), exps)
own <- gens[vapply(gens, function(n)
  home_of(get(n, envir = sns)) == "frmtmb.sample", NA)]

libs <- c(SYS, UL, "C:/Users/adf44/source/r/pinlib")
inst <- do.call(rbind, lapply(libs, function(l) {
  d <- list.dirs(l, recursive = FALSE, full.names = FALSE)
  d <- d[file.exists(file.path(l, d, "NAMESPACE"))]
  if (!length(d)) return(NULL)
  data.frame(pkg = d, lib = l, stringsAsFactors = FALSE)
}))
inst <- inst[!duplicated(inst$pkg), ]
inst <- inst[!grepl("^frmtmb", inst$pkg), ]
cat(sprintf("installed packages scanned        %d\n", nrow(inst)))

exporters <- list()
for (i in seq_len(nrow(inst))) {
  ex <- tryCatch(parseNamespaceFile(inst$pkg[i], inst$lib[i])$exports,
                 error = function(e) character())
  for (h in intersect(own, ex)) {
    exporters[[h]] <- c(exporters[[h]], inst$pkg[i])
  }
}
cands <- sort(unique(unlist(exporters)))
cat("candidate exporters               ",
    paste(cands, collapse = ", "), "\n")
for (p in cands) {
  ok <- tryCatch({ suppressMessages(loadNamespace(p)); TRUE },
                 error = function(e) FALSE)
  if (!ok) cat("  COULD NOT LOAD", p, "\n")
}

sig <- function(f) paste(names(formals(f)), collapse = ", ")
tab_has <- function(p, m) {
  ns <- asNamespace(p)
  tb <- get0(".__S3MethodsTable__.", envir = ns, inherits = FALSE)
  !is.null(tb) && exists(m, envir = tb, inherits = FALSE)
}

cat(sprintf("\nfrmtmb.sample exports %d S3 generics, defines %d itself\n",
            length(gens), length(own)))
cat("\n### per name: definers (generic formals), importers, brmsfit method\n")
rows <- list()
for (g in sort(own)) {
  defs <- character(); imps <- character(); nongen <- character()
  for (p in exporters[[g]]) {
    f <- tryCatch(getExportedValue(p, g), error = function(e) NULL)
    if (is.null(f)) next
    h <- home_of(f)
    if (!is_gen(f)) {
      nongen <- c(nongen, sprintf("%s(not a generic)", p))
    } else if (identical(h, p)) {
      defs <- c(defs, sprintf("%s(%s)", p, sig(f)))
    } else {
      imps <- c(imps, sprintf("%s<-%s", p, h))
    }
  }
  m <- paste0(g, ".brmsfit")
  where <- cands[vapply(cands, tab_has, NA, m)]
  cat(sprintf("\n%s   frmtmb.sample(%s)\n", g, sig(get(g, envir = sns))))
  cat("  defines:  ", paste(defs, collapse = "  "), "\n")
  if (length(imps)) cat("  imports:  ", paste(imps, collapse = "  "), "\n")
  if (length(nongen)) cat("  other:    ", paste(nongen, collapse = "  "), "\n")
  cat("  brmsfit method in table of:",
      if (length(where)) paste(where, collapse = ", ") else "(none)", "\n")
  rows[[g]] <- data.frame(gen = g,
    definers = paste(sub("[(].*", "", defs), collapse = " "),
    brmsfit = paste(where, collapse = " "), stringsAsFactors = FALSE)
}

cat("\n### frmtmb.sample's own S3method directives on these names\n")
m <- sample_ns$S3methods
for (g in sort(own)) {
  hit <- m[m[, 1] == g, , drop = FALSE]
  tg <- ifelse(is.na(hit[, 4]), "(own)", hit[, 4])
  cat(sprintf("  %-20s %s\n", g, paste(tg, collapse = " ")))
}
cat("\nDONE\n")
