# How much of the same defect is left in frmtmb.sample, and how much of
# it did core's fix reach for free?
#
# frmtmb.sample re-exports 21 generics from core (those follow core's
# binding, because its importFrom resolves at load time and frmtmb
# always loads first) and DEFINES a further set of its own, which core
# cannot reach.
LIB <- "C:/Users/adf44/source/r/generics-lib"
.libPaths(c(LIB, "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6",
            file.path(R.home(), "library")))
suppressMessages(library(frmtmb.sample))

ns <- asNamespace("frmtmb.sample")
exps <- parseNamespaceFile("frmtmb.sample",
                           dirname(find.package("frmtmb.sample")))$exports

is_gen <- function(nm, env) {
  f <- tryCatch(get(nm, envir = env), error = function(e) NULL)
  if (!is.function(f)) return(FALSE)
  any(grepl("UseMethod", deparse(body(f)), fixed = TRUE))
}
gens <- Filter(function(n) is_gen(n, ns), exps)

home <- vapply(gens, function(n) {
  f <- tryCatch(get(n, envir = ns), error = function(e) NULL)
  e <- environmentName(topenv(environment(f)))
  if (nzchar(e)) e else "?"
}, "")

libs <- c(file.path(R.home(), "library"),
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
inst <- unique(unlist(lapply(libs, function(l)
  list.dirs(l, recursive = FALSE, full.names = FALSE))))
inst <- setdiff(inst, c("frmtmb", grep("^frmtmb[.]", inst, value = TRUE)))
owners <- list()
for (p in inst) {
  ex <- tryCatch(parseNamespaceFile(p, dirname(find.package(p)))$exports,
                 error = function(e) character())
  for (h in intersect(gens, ex)) owners[[h]] <- c(owners[[h]], p)
}

own_gen <- names(home)[home == "frmtmb.sample"]
rival <- own_gen[vapply(own_gen, function(g) !is.null(owners[[g]]), NA)]

cat("```\n")
cat("== frmtmb.sample, dev/generics-audit4.R ==\n")
cat(sprintf("exported S3 generics                       %d\n", length(gens)))
cat(sprintf("  re-exported from core, adopted with core %d\n",
            sum(home != "frmtmb.sample")))
cat(sprintf("  defined by frmtmb.sample itself          %d\n",
            length(own_gen)))
cat(sprintf("  of those, a name another package owns    %d\n",
            length(rival)))
cat("\nthe names frmtmb.sample still owns that someone else owns:\n")
for (g in sort(rival)) {
  cat(sprintf("  %-22s %s\n", g, paste(owners[[g]], collapse = ", ")))
}
cat("\nhow those re-exports resolved in THIS process, where no optional owner was attached:\n")
tb <- table(home[home != "frmtmb.sample"])
cat("  ", paste(sprintf("%s=%d", names(tb), as.integer(tb)),
                collapse = " "), "\n")
cat("```\n")
