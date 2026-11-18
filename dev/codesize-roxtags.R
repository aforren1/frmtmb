# Part of Q7 of the codesize lane: what the roxygen blocks that make no Rd
# page (blocks with @noRd) hold, by tag. Lines of a tag include its
# continuation lines up to the next tag.
#   Rscript dev/codesize-roxtags.R
# Reads dev/codesize-out/parsed.rds.
out <- "C:/Users/adf44/source/r/frmtmb-wt-codesize/dev/codesize-out"
P <- readRDS(file.path(out, "parsed.rds"))
tot <- c(x = 0)[0]
nblk <- 0
for (f in names(P$core)) {
  x <- P$core[[f]]$lines
  rox <- grepl("^\\s*#'", x)
  r <- rle(rox); e <- cumsum(r$lengths); s <- e - r$lengths + 1
  for (i in which(r$values)) {
    b <- x[s[i]:e[i]]
    if (!any(grepl("@noRd", b))) next
    nblk <- nblk + 1
    tag <- character(length(b)); cur <- "(description)"
    for (j in seq_along(b)) {
      m <- regmatches(b[j], regexpr("^\\s*#'\\s*@[A-Za-z]+", b[j]))
      if (length(m)) cur <- sub("^.*@", "@", m)
      if (grepl("^\\s*#'\\s*$", b[j])) tag[j] <- "(blank #')" else
        tag[j] <- cur
    }
    t <- table(tag)
    for (k in names(t)) tot[k] <- (if (k %in% names(tot)) tot[[k]] else 0) + t[[k]]
  }
}
tot <- sort(tot, decreasing = TRUE)
cat("BEGIN codesize-roxtags\n")
cat(sprintf("@noRd blocks: %d, lines %d\n", nblk, sum(tot)))
cat(sprintf("  %s: %d\n", names(tot), tot), sep = "")
cat("END codesize-roxtags\n")
