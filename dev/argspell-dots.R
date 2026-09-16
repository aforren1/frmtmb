# Which registered S3 methods take `...` and never read it? Those are the
# methods that swallow a misspelled or unsupported argument in silence.
LIB <- "C:/Users/adf44/source/r/argspell-lib"
.libPaths(c(LIB,
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
pkgs <- c("frmtmb", "frmtmb.sample")
for (p in pkgs) suppressMessages(library(p, character.only = TRUE))

# Deparse rather than walk the AST: a call with an empty argument, such as
# `x[i, ]`, holds the empty symbol, and touching it raises "argument is
# missing" inside a recursive walker.
uses_dots <- function(b) {
  any(grepl("...", deparse(b), fixed = TRUE))
}

rows <- list()
for (p in pkgs) {
  nsf <- parseNamespaceFile(p, dirname(system.file(package = p)))
  s3 <- nsf$S3methods
  seen <- character()
  for (i in seq_len(nrow(s3))) {
    gen <- s3[i, 1]; cls <- s3[i, 2]
    fun <- if (!is.na(s3[i, 3])) s3[i, 3] else paste0(gen, ".", cls)
    key <- paste(p, fun)
    if (key %in% seen) next
    seen <- c(seen, key)
    obj <- tryCatch(get(fun, envir = asNamespace(p)), error = function(e) NULL)
    if (is.null(obj)) next
    has <- "..." %in% names(formals(obj))
    rows[[length(rows) + 1L]] <- data.frame(
      pkg = p, generic = gen, class = cls, fun = fun,
      has_dots = has,
      reads_dots = if (has) uses_dots(body(obj)) else NA,
      stringsAsFactors = FALSE)
  }
}
tab <- do.call(rbind, rows)
write.csv(tab, file.path("dev", "argspell-dots.csv"), row.names = FALSE)
cat(sprintf("registered methods (deduped): %d\n", nrow(tab)))
cat(sprintf("  with `...`            : %d\n", sum(tab$has_dots)))
cat(sprintf("  with `...`, never read: %d\n",
            sum(tab$has_dots & !tab$reads_dots)))
cat("\n== SWALLOWERS: take `...`, body never mentions it ==\n")
sw <- tab[tab$has_dots & !tab$reads_dots, ]
for (i in seq_len(nrow(sw))) cat(sprintf("%-14s %s\n", sw$pkg[i], sw$fun[i]))
cat("\n== pass-through: take `...` and use it ==\n")
pt <- tab[tab$has_dots & tab$reads_dots, ]
for (i in seq_len(nrow(pt))) cat(sprintf("%-14s %s\n", pt$pkg[i], pt$fun[i]))
