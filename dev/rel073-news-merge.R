# Build the 0.73.0 NEWS sections from the lanes' development-version
# sections (read from each lane's worktree, which is not edited) and the
# consolidation's leads (dev/rel073-news-*-lead.md), regrouped:
# breaking changes, new features, bug fixes, then any other section.
# Within a section the bullets keep the merge order of the lanes
# (loodraws, cdfs, nb2dm, fuzzfix). The released sections below come
# from the HEAD blob (4854e425, 0.72.0), so the merge's conflict markers
# at the top of the file never reach the result, and the script is
# rerunnable. Adapted from dev/rel072-news-merge.R.
#
# A bullet marked "BREAKING:" moves to "Breaking changes" and loses the
# prefix. A lane's own lead paragraph is replaced by the consolidation's
# lead, which names the version. The sentences over 25 words are split
# (dev/rel073-news-split.R), every bullet paragraph of the new sections
# is rewrapped to 72 columns, and the new sections name no lane and no
# `dev/` path (the merge stops if one is left): the trail is
# dev/rel073-news-map.md, which dev/round-20261021.md pastes.
#
#   Rscript dev/rel073-news-merge.R
root <- "C:/Users/adf44/source/r"
rel <- Sys.getenv("REL073_ROOT", file.path(root, "frmtmb-wt-release"))
git_dir <- Sys.getenv("REL073_GIT", rel)
setwd(rel)
lanes <- strsplit(Sys.getenv("REL073_LANES", "loodraws cdfs nb2dm fuzzfix"),
                  " ")[[1]]

rd <- function(f) sub("\r$", "", readLines(f, warn = FALSE))
wr <- function(x, f) {
  con <- file(f, "wb")
  writeLines(x, con, sep = "\n")
  close(con)
}
dev_section <- function(f) {
  if (!file.exists(f)) return(NULL)
  x <- rd(f)
  h <- grep("^# ", x)
  if (!grepl("development version", x[h[1]])) return(NULL)
  x[(h[1] + 1):(h[2] - 1)]
}
split_secs <- function(v, loose = "New features") {
  h <- grep("^## ", v)
  lead <- if (length(h)) v[seq_len(h[1] - 1)] else v
  secs <- list()
  s0 <- grep("^\\* ", lead)
  if (length(s0)) {
    secs[[loose]] <- lead[s0[1]:length(lead)]
    lead <- lead[seq_len(s0[1] - 1)]
  }
  for (i in seq_along(h)) {
    end <- if (i < length(h)) h[i + 1] - 1 else length(v)
    body <- if (end > h[i]) v[(h[i] + 1):end] else character(0)
    k <- sub("^## ", "", v[h[i]])
    secs[[k]] <- c(secs[[k]], body)
  }
  list(lead = lead, secs = secs)
}
bullets <- function(body) {
  s <- grep("^\\* ", body)
  if (!length(s)) return(list())
  lapply(seq_along(s), function(i) {
    end <- if (i < length(s)) s[i + 1] - 1 else length(body)
    y <- body[s[i]:end]
    while (length(y) && y[length(y)] == "") y <- y[-length(y)]
    y
  })
}
unbreak <- function(b) {
  for (p in c("* **BREAKING: ", "* BREAKING: ")) {
    if (startsWith(b[1], p)) {
      r <- substring(b[1], nchar(p) + 1L)
      b[1] <- paste0(sub("BREAKING: $", "", p),
                     toupper(substr(r, 1, 1)), substr(r, 2, nchar(r)))
      return(b)
    }
  }
  NULL
}
newsmap <- list()
build <- function(pkg_dir, pkg, version, prev_version, lead,
                  order_first) {
  f <- file.path(pkg_dir, "NEWS.md")
  rel_f <- sub(paste0(rel, "/?"), "", f)
  sec <- list()
  count <- integer(0)
  for (l in lanes) {
    lf <- file.path(root, paste0("frmtmb-wt-", l), rel_f)
    d <- dev_section(lf)
    if (is.null(d)) next
    s <- split_secs(d)
    n <- 0L
    for (k in names(s$secs)) {
      for (b in bullets(s$secs[[k]])) {
        kk <- k
        u <- unbreak(b)
        if (!is.null(u)) {
          b <- u
          kk <- "Breaking changes"
          cat("  moved to Breaking changes (", l, "): ", b[1], "\n", sep = "")
        }
        sec[[kk]] <- c(sec[[kk]], list(b))
        n <- n + 1L
        hd <- gsub("[*]", "", sub("^[*] ", "", b[1]))
        newsmap[[length(newsmap) + 1L]] <<- c(pkg, kk, l,
          paste0("`dev/", l, "-findings.md`"), substr(hd, 1, 60))
      }
    }
    count[l] <- n
  }
  keys <- c(intersect(order_first, names(sec)),
            setdiff(names(sec), order_first))
  out <- c(paste0("# ", pkg, " ", version), "", lead, "")
  for (k in keys) {
    out <- c(out, paste("##", k), "")
    for (b in sec[[k]]) out <- c(out, b, "")
  }
  x <- sub("\r$", "", system2("git", c("-C", git_dir, "show",
                                       paste0("HEAD:", rel_f)),
                              stdout = TRUE))
  prev <- which(x == paste0("# ", pkg, " ", prev_version))
  stopifnot(length(prev) == 1L)
  res <- c(out, x[prev:length(x)])
  res <- res[!(c(FALSE, res[-1] == "" & res[-length(res)] == ""))]
  while (res[length(res)] == "") res <- res[-length(res)]
  stopifnot(!any(grepl("^(<<<<<<<|=======|>>>>>>>)", res)),
            !any(grepl("^# .*development version", res)),
            !any(grepl("frmtmb release|next frmtmb|next release",
                       res[seq_along(out)])))
  wr(res, f)
  cat(pkg, version, ": sections", paste(keys, collapse = " | "), "\n")
  cat("  bullets per lane:", paste(names(count), count, collapse = ", "),
      "| total:", sum(lengths(sec)), "\n")
  for (k in keys) cat("  ", k, ": ", length(sec[[k]]), "\n", sep = "")
}
order_first <- c("Breaking changes", "New features", "Bug fixes",
                 "Performance", "Documentation")
src <- function(nm) rd(file.path(rel, "dev", paste0("rel073-news-", nm, ".md")))
news_files <- c("NEWS.md", "extensions/frmtmb.sample/NEWS.md")
build(rel, "frmtmb", "0.73.0", "0.72.0", src("core-lead"), order_first)
build(file.path(rel, "extensions/frmtmb.sample"), "frmtmb.sample",
      "0.21.0", "0.20.0", src("sample-lead"), order_first)

top_of <- function(x) {
  h <- grep("^# ", x)
  list(a = h[1] + 1L, b = h[2] - 1L)
}
# the sentences over 25 words, split (dev/rel073-news-split.R); each
# edit is a fixed text matched once in the file's top section with any
# run of whitespace standing for a space, so it does not depend on
# where the lane wrapped its lines
ws_regex <- function(s) {
  parts <- strsplit(s, "\\s+")[[1]]
  esc <- gsub("([][{}()+*^$|\\\\?.])", "\\\\\\1", parts)
  paste(esc, collapse = "\\s+")
}
sp <- file.path(rel, "dev/rel073-news-split.R")
splits <- if (file.exists(sp)) source(sp)$value else list()
# then the release review's fixes (dev/rel073-news-review.R), on the
# split text
rv <- source(file.path(rel, "dev/rel073-news-review.R"))$value
for (nf in names(rv)) splits[[nf]] <- c(splits[[nf]], rv[[nf]])
for (nf in names(splits)) {
  fp <- file.path(rel, nf)
  x <- rd(fp)
  t <- top_of(x)
  top <- paste(x[t$a:t$b], collapse = "\n")
  for (e in splits[[nf]]) {
    re <- ws_regex(e[1])
    k <- sum(gregexpr(re, top, perl = TRUE)[[1]] > 0)
    if (k != 1L) stop(nf, ": split, expected 1 match, found ", k, ": ",
                      substr(e[1], 1, 60))
    top <- sub(re, gsub("\\", "\\\\", gsub("\\s+", " ", e[2]), fixed = TRUE),
               top, perl = TRUE)
  }
  orig <- x[t$a:t$b]
  trail <- length(orig) - max(which(nzchar(orig)))
  body <- strsplit(top, "\n", fixed = TRUE)[[1]]
  body <- c(body[seq_len(max(which(nzchar(body))))], rep("", trail))
  wr(c(x[seq_len(t$a - 1L)], body, x[(t$b + 1L):length(x)]), fp)
  cat("split", length(splits[[nf]]), "sentences in", nf, "\n")
}
wrap_par <- function(lines) {
  lead <- if (startsWith(lines[1], "* ")) "* " else
    if (startsWith(lines[1], "  - ")) "  - " else "  "
  cont <- if (lead == "  - ") "    " else "  "
  body <- substring(lines[1], nchar(lead) + 1L)
  txt <- trimws(paste(c(body, trimws(lines[-1])), collapse = " "))
  w0 <- strsplit(txt, " +")[[1]]
  # a code span is one unit, so no wrap breaks it (`trunc(lb = )`), and
  # a line never starts with what markdown would read as a block quote,
  # a list item or a heading (">= k" ends a paragraph)
  words <- character(0)
  open <- FALSE
  for (w in w0) {
    if (open) words[length(words)] <- paste(words[length(words)], w)
    else words <- c(words, w)
    if (lengths(regmatches(w, gregexpr("`", w))) %% 2L == 1L) open <- !open
  }
  risky <- function(w) grepl("^([>#+*-]|[0-9]+[.)]( |$))", w)
  out <- character(0)
  cur <- paste0(lead, words[1])
  for (w in words[-1]) {
    cand <- paste(cur, w)
    if (nchar(cand) > 72 && !(risky(w) && nchar(cand) <= 80)) {
      out <- c(out, cur)
      cur <- paste0(cont, w)
    } else {
      cur <- cand
    }
  }
  c(out, cur)
}
rewrap_top <- function(x) {
  t <- top_of(x)
  # from the first "## " heading on: the lead above it is the
  # consolidation's text, not a bullet
  t$a <- t$a - 1L + grep("^## ", x[t$a:t$b])[1]
  sec <- x[t$a:t$b]
  out <- character(0)
  i <- 1L
  while (i <= length(sec)) {
    if (grepl("^(\\* |  - |  [^ -])", sec[i])) {
      j <- i
      while (j < length(sec) && nzchar(sec[j + 1L]) &&
             !grepl("^(\\* |  - |#)", sec[j + 1L])) j <- j + 1L
      out <- c(out, wrap_par(sec[i:j]))
      i <- j + 1L
    } else {
      out <- c(out, sec[i])
      i <- i + 1L
    }
  }
  c(x[seq_len(t$a - 1L)], out, x[(t$b + 1L):length(x)])
}
# no lane name, no "lane", no `dev/` path in a new section
strip_lanes <- lanes
mentions <- function(x) {
  pat <- paste0("\\blanes?\\b|dev/|\\b(",
                paste(strip_lanes, collapse = "|"), ")\\b")
  sum(vapply(gregexpr(pat, x, ignore.case = TRUE, perl = TRUE),
             function(m) sum(m > 0), 0))
}
for (nf in news_files) {
  fp <- file.path(rel, nf)
  x <- rewrap_top(rd(fp))
  t <- top_of(x)
  left <- mentions(x[t$a:t$b])
  if (left > 0L) stop(nf, ": ", left, " lane or dev/ mentions left")
  wr(x, fp)
  cat("rewrapped", nf, "| lane or dev/ mentions:", left,
      "| lines over 80:", sum(nchar(x[t$a:t$b]) > 80),
      "| em dashes:", sum(grepl("\u2014", x[t$a:t$b])), "\n")
}
m <- do.call(rbind, newsmap)
out <- c("| package | section (as the lane filed it) | lane | record | bullet |",
         "|---|---|---|---|---|",
         sprintf("| %s | %s | %s | %s | %s |", m[, 1], m[, 2], m[, 3],
                 m[, 4], gsub("|", "\\|", m[, 5], fixed = TRUE)))
wr(out, file.path(rel, "dev/rel073-news-map.md"))
cat("wrote dev/rel073-news-map.md:", nrow(m), "bullets\n")
cat("NEWS MERGE DONE\n")
