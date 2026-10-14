# Release reviewer of 0.71.0: the exports frmtmb 0.71.0 adds over 0.70.0
# (and over 0.68.1), and which extension's R code names any of them as
# a word, which a floor below 0.71.0 would not provide.
exports_of <- function(lib) {
  e <- parseNamespaceFile("frmtmb", lib)
  unique(c(e$exports, names(e$S3methods)))
}
L <- "C:/Users/adf44/source/r/"
x9 <- exports_of(paste0(L, "rellib-r9"))
x8 <- exports_of(paste0(L, "rellib-r8"))
x6 <- exports_of(paste0(L, "rellib-r6"))
new8 <- setdiff(x9, x8)
new6 <- setdiff(x9, x6)
cat("exports new over 0.70.0:", length(new8), "\n ", paste(sort(new8),
                                                           collapse = ", "), "\n")
cat("removed since 0.70.0:", paste(setdiff(x8, x9), collapse = ", "), "\n")
root <- paste0(L, "frmtmb-wt-release/extensions/")
for (p in c("frmtmb.sample", "frmtmb.eam", "frmtmb.latent", "frmtmb.spline",
            "frmtmb.coupling", "frmtmb.learn", "frmtmb.ode")) {
  files <- list.files(file.path(root, p, "R"), "\\.R$", full.names = TRUE)
  src <- unlist(lapply(files, readLines, warn = FALSE))
  src <- sub("#.*$", "", src)
  defs <- sub("^([A-Za-z_.][A-Za-z0-9_.]*) *<- *function.*$", "\\1",
              grep("^[A-Za-z_.][A-Za-z0-9_.]* *<- *function", src,
                   value = TRUE))
  hit <- function(set) {
    set <- setdiff(set, defs)
    set[vapply(set, function(s) {
      any(grepl(paste0("(^|[^A-Za-z0-9_.])", gsub(".", "\\.", s, fixed = TRUE),
                       " *\\("), src))
    }, NA)]
  }
  cat(sprintf("%-16s calls exports new over 0.70.0: %s | over 0.68.1: %s\n", p,
              paste(hit(new8), collapse = ", "),
              paste(hit(new6), collapse = ", ")))
}
