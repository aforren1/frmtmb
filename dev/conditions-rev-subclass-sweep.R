# Reviewer, lane wt-conditions (round 2): in each extension's own suite,
# how many caught frmtmb_error conditions lack that extension's subclass,
# and which of those name an extension function in their text.
#   Rscript dev/conditions-rev-subclass-sweep.R <lane sweep dir>
fs <- list.files(commandArgs(TRUE)[1L], "rds$", full.names = TRUE)
x <- do.call(rbind, lapply(fs, readRDS))
x <- x[grepl("^expect_error", x$key) & grepl("frmtmb_error", x$class), ]
x$pkg <- regmatches(x$file, regexpr("frmtmb[.][a-z]+|^tests", x$file))
for (p in setdiff(unique(x$pkg), "tests")) {
  own <- paste0(gsub(".", "_", p, fixed = TRUE), "_error")
  y <- x[x$pkg == p, ]
  has <- grepl(own, y$class)
  cat(sprintf("\n%-16s caught %3d  with %s %3d  without %3d  distinct without %d\n",
              p, nrow(y), own, sum(has), sum(!has),
              length(unique(y$message[!has]))))
  u <- unique(y$message[!has])
  cat(paste0("   - ", substr(gsub("\\s+", " ", u), 1, 100)), sep = "\n")
}

# which of those carry text written in the extension's own R sources:
# the first 30 characters of the message, searched as a literal in the
# concatenated sources (a lower bound: messages pasted from pieces are
# missed)
cat("\n== messages without the subclass whose opening text is in the extension's sources\n")
root <- "C:/Users/adf44/source/r/frmtmb-wt-conditions"
for (p in setdiff(unique(x$pkg), "tests")) {
  own <- paste0(gsub(".", "_", p, fixed = TRUE), "_error")
  y <- x[x$pkg == p & !grepl(own, x$class), ]
  u <- unique(y$message)
  src <- paste(unlist(lapply(list.files(file.path(root, "extensions", p, "R"),
                                        full.names = TRUE), readLines,
                             warn = FALSE)), collapse = " ")
  src <- gsub("[[:space:]]+", " ", src)
  core <- paste(unlist(lapply(list.files(file.path(root, "R"),
                                         full.names = TRUE), readLines,
                              warn = FALSE)), collapse = " ")
  core <- gsub("[[:space:]]+", " ", core)
  hit_ext <- vapply(u, function(m) grepl(substr(gsub("[[:space:]]+", " ", m), 1, 30),
                                         src, fixed = TRUE), NA)
  hit_core <- vapply(u, function(m) grepl(substr(gsub("[[:space:]]+", " ", m), 1, 30),
                                          core, fixed = TRUE), NA)
  cat(sprintf("%-16s distinct %2d  opening text in extension sources %2d  in core sources %2d  neither %2d\n",
              p, length(u), sum(hit_ext), sum(hit_core & !hit_ext),
              sum(!hit_ext & !hit_core)))
}
