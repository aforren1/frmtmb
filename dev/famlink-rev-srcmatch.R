## Does the lane library hold the worktree's R/ sources? Every function
## defined in R/*.R is compared, by deparsed body and formals, with the
## installed namespace. Run: Rscript dev/famlink-rev-srcmatch.R
ARM <- "lane"
source("dev/famlink-rev-common.R")
ns <- asNamespace("frmtmb")
env <- new.env()
n_fun <- 0L; bad <- character()
for (f in list.files("R", pattern = "[.]R$", full.names = TRUE)) {
  ex <- parse(f, keep.source = FALSE)
  for (e in ex) {
    if (is.call(e) && (identical(e[[1]], as.name("<-")) ||
                       identical(e[[1]], as.name("="))) &&
        is.name(e[[2]]) && is.call(e[[3]]) &&
        identical(e[[3]][[1]], as.name("function"))) {
      nm <- as.character(e[[2]])
      if (!exists(nm, envir = ns, inherits = FALSE)) next
      got <- get(nm, envir = ns)
      if (!is.function(got)) next
      want <- eval(e[[3]], env)
      n_fun <- n_fun + 1L
      if (!identical(deparse(body(want)), deparse(body(got))) ||
          !identical(deparse(formals(want)), deparse(formals(got)))) {
        bad <- c(bad, paste0(basename(f), ":", nm))
      }
    }
  }
}
cat("functions compared:", n_fun, " mismatched:", length(bad), "\n")
print(head(bad, 30))
