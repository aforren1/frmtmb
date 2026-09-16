# Formals audit: every registered S3 method in frmtmb and frmtmb.sample
# against the brms method for the same generic. Run against the BASE
# install in the shared reference library, which is what the recorded
# "before" column in dev/argspell-findings.md describes.
LIB <- "C:/Users/adf44/source/r/argspell-lib"
.libPaths(c(LIB,
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))

pkgs <- c("frmtmb", "frmtmb.sample")
for (p in pkgs) suppressMessages(library(p, character.only = TRUE))
suppressMessages(library(brms))

fmt_formals <- function(f) {
  if (is.null(f)) return(NA_character_)
  paste(names(formals(f)), collapse = ", ")
}

rows <- list()
for (p in pkgs) {
  nsf <- parseNamespaceFile(p, dirname(system.file(package = p)))
  s3 <- nsf$S3methods
  for (i in seq_len(nrow(s3))) {
    gen <- s3[i, 1]
    cls <- s3[i, 2]
    fun <- if (!is.na(s3[i, 3])) s3[i, 3] else paste0(gen, ".", cls)
    obj <- tryCatch(get(fun, envir = asNamespace(p)), error = function(e) NULL)
    brms_cls <- c("brmsfit", "brmsformula", "brmsterms", "brmsprior",
                  "mvbrmsformula", "brms_environment")
    bm <- NULL
    bname <- NA_character_
    for (bc in brms_cls) {
      cand <- paste0(gen, ".", bc)
      o <- tryCatch(get(cand, envir = asNamespace("brms")),
                    error = function(e) NULL)
      if (!is.null(o)) { bm <- o; bname <- cand; break }
    }
    if (is.null(bm)) {
      o <- tryCatch(get(gen, envir = asNamespace("brms")),
                    error = function(e) NULL)
      if (!is.null(o) && is.function(o)) { bm <- o; bname <- gen }
    }
    rows[[length(rows) + 1L]] <- data.frame(
      pkg = p, generic = gen, class = cls, fun = fun,
      ours = fmt_formals(obj),
      brms_fun = bname, theirs = fmt_formals(bm),
      stringsAsFactors = FALSE)
  }
}
tab <- do.call(rbind, rows)
tab <- tab[order(tab$pkg, tab$generic, tab$class), ]
write.csv(tab, file.path("dev", "argspell-formals.csv"), row.names = FALSE)

cat("== methods with a brms counterpart ==\n")
have <- tab[!is.na(tab$brms_fun), ]
for (i in seq_len(nrow(have))) {
  a <- strsplit(have$ours[i], ", ")[[1]]
  b <- strsplit(have$theirs[i], ", ")[[1]]
  only_b <- setdiff(b, c(a, "object", "x", "formula", "..."))
  only_a <- setdiff(a, c(b, "object", "x", "formula", "..."))
  cat(sprintf("%-14s %-28s vs %-24s\n  ours  : %s\n  brms  : %s\n",
              have$pkg[i], have$fun[i], have$brms_fun[i],
              have$ours[i], have$theirs[i]))
  cat(sprintf("  brms-only: %s\n  ours-only: %s\n\n",
              paste(only_b, collapse = ", "),
              paste(only_a, collapse = ", ")))
}
cat("== methods with NO brms counterpart ==\n")
non <- tab[is.na(tab$brms_fun), ]
for (i in seq_len(nrow(non))) {
  cat(sprintf("%-14s %-34s ours: %s\n", non$pkg[i], non$fun[i], non$ours[i]))
}

cat("\n== our methods whose formals are exactly (object/x, ...) ==\n")
bare <- tab[tab$ours %in% c("object, ...", "x, ...", "..."), ]
print(bare[, c("pkg", "fun", "ours", "brms_fun")], row.names = FALSE)
