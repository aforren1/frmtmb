## REVIEW: the lane measured 5 of the 13 exempt generics. Measure all
## 13: for each, find call sites in the owning packages that pass a
## NAMED argument our method has no formal for. That is the argument a
## guard would refuse, and it is what makes an exemption necessary.
.libPaths(c("C:/Users/adf44/source/r/asrev-lib",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
exempt <- c("emm_basis", "recover_data", "get_coef", "get_predict",
            "get_vcov", "set_coef", "get_varcov", "get_parameters",
            "find_formula", "find_random", "find_statistic",
            "link_function", "link_inverse")
owners <- c("emmeans", "insight", "marginaleffects", "datawizard",
            "parameters", "performance", "modelbased", "ggeffects")
owners <- Filter(function(p) requireNamespace(p, quietly = TRUE), owners)
cat("owners loadable:", paste(owners, collapse = ", "), "\n\n")
for (g in exempt) {
  ours <- tryCatch(get(paste0(g, ".frmtmb_fit"), envir = ns),
                   error = function(e) NULL)
  if (is.null(ours)) { cat(g, ": NO METHOD\n"); next }
  fo <- setdiff(names(formals(ours)), "...")
  hits <- character()
  for (pk in owners) {
    nsp <- asNamespace(pk)
    for (nm in ls(nsp, all.names = TRUE)) {
      o <- tryCatch(get(nm, envir = nsp), error = function(e) NULL)
      if (!is.function(o)) next
      src <- tryCatch(deparse(body(o)), error = function(e) character())
      ln <- grep(paste0("(^|[^._[:alnum:]])", g, "\\("), src, value = TRUE)
      for (l in ln) {
        # named arguments in the call text
        args <- regmatches(l, gregexpr("[A-Za-z._][A-Za-z._0-9]* *=", l))[[1]]
        args <- trimws(sub("=$", "", trimws(args)))
        extra <- setdiff(args, c(fo, "..."))
        # only names that are plausibly arguments of THIS call
        extra <- intersect(extra, c("component", "effects", "vcov", "type",
                                    "newdata", "verbose", "flatten",
                                    "split_nested", "summary", "data",
                                    "trms", "xlev", "grid", "coefs",
                                    "model", "parameters", "dpar"))
        if (length(extra)) {
          hits <- c(hits, sprintf("%s::%s | %s | would-refuse: %s",
                                  pk, nm, trimws(l),
                                  paste(extra, collapse = ", ")))
        }
      }
    }
  }
  cat(sprintf("%-16s ours(%s)  sites passing a name we lack: %d\n", g,
              paste(fo, collapse = ","), length(hits)))
  if (length(hits)) cat(paste0("    ", utils::head(hits, 3)), sep = "\n")
}
cat("\nDONE\n")
