# Reviewer, item 8: every lane NEWS bullet in the release NEWS once.
# A bullet is matched by its first 60 characters after "* " (bold
# markers dropped), then by its longest shared 8-word run, so an edited
# bullet still finds its source.
#   Rscript dev/relrev-news.R > dev/relrev-log/news.txt
root <- "C:/Users/adf44/source/r"
bullets <- function(path, head_rx) {
  l <- readLines(path, warn = FALSE)
  s <- grep(head_rx, l)[1]
  e <- grep("^# ", l); e <- e[e > s][1] - 1L
  if (is.na(e)) e <- length(l)
  l <- l[s:e]
  b <- character(); cur <- NULL; sec <- ""; secs <- character()
  for (x in l) {
    if (grepl("^## ", x)) { if (!is.null(cur)) { b <- c(b, cur); secs <- c(secs, sec) }; cur <- NULL; sec <- sub("^## ", "", x); next }
    if (grepl("^\\* ", x)) { if (!is.null(cur)) { b <- c(b, cur); secs <- c(secs, sec) }; cur <- sub("^\\* ", "", x) }
    else if (!is.null(cur) && nzchar(trimws(x))) cur <- paste(cur, trimws(x))
    else if (!is.null(cur) && !nzchar(trimws(x))) { b <- c(b, cur); secs <- c(secs, sec); cur <- NULL }
  }
  if (!is.null(cur)) { b <- c(b, cur); secs <- c(secs, sec) }
  data.frame(text = gsub("\\*\\*", "", b), sec = secs, stringsAsFactors = FALSE)
}
words <- function(x) strsplit(gsub("[^A-Za-z0-9_()=.\"' -]", " ", x), " +")[[1]]
run8 <- function(a, b) {
  wa <- words(a); wb <- paste(words(b), collapse = " ")
  if (length(wa) < 8) return(grepl(paste(wa, collapse = " "), wb, fixed = TRUE))
  any(vapply(seq_len(length(wa) - 7), function(i)
    grepl(paste(wa[i:(i + 7)], collapse = " "), wb, fixed = TRUE), NA))
}
for (pkg in c("core", "sample")) {
  sub_ <- if (pkg == "core") "NEWS.md" else "extensions/frmtmb.sample/NEWS.md"
  hx <- if (pkg == "core") "^# frmtmb 0.68.0" else "^# frmtmb.sample 0.16.0"
  rel <- bullets(file.path(root, "frmtmb-wt-release", sub_), hx)
  cat("\n####", pkg, ": release bullets", nrow(rel), "; by section:\n")
  print(table(factor(rel$sec, levels = unique(rel$sec))))
  hit <- integer(nrow(rel))
  for (lane in c("fixes", "gpby", "ordmix", "vigport")) {
    p <- file.path(root, paste0("frmtmb-wt-", lane), sub_)
    if (!file.exists(p)) next
    lb <- tryCatch(bullets(p, "development version"), error = function(e) NULL)
    if (is.null(lb) || !nrow(lb)) { cat(lane, ": 0 bullets\n"); next }
    found <- vapply(lb$text, function(t) {
      w <- which(vapply(rel$text, function(r) run8(t, r), NA))
      if (length(w)) paste(w, collapse = ",") else "MISSING"
    }, "")
    cat(lane, ":", nrow(lb), "bullets; missing", sum(found == "MISSING"),
        "; matched to release bullets:", paste(found, collapse = " "), "\n")
    for (i in which(found == "MISSING")) cat("   MISSING:", substr(lb$text[i], 1, 110), "\n")
    for (w in unlist(strsplit(found[found != "MISSING"], ","))) hit[as.integer(w)] <- hit[as.integer(w)] + 1L
  }
  cat("release bullets matched by no lane bullet (consolidation's):", sum(hit == 0), "\n")
  for (i in which(hit == 0)) cat("   ", substr(rel$text[i], 1, 100), "\n")
  cat("release bullets matched by more than one lane bullet:", sum(hit > 1), "\n")
  for (i in which(hit > 1)) cat("   [", hit[i], "]", substr(rel$text[i], 1, 100), "\n")
}
