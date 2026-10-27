# Build the 0.76.0 NEWS sections from the lanes' development-version
# sections (read from each lane's worktree, which is not edited), the
# consolidation's leads (dev/rel076-news-<pkg>-lead.md) and the
# consolidation's own bullets (dev/rel076-news-<pkg>-extra.md, for a
# package a lane changed without a NEWS bullet), regrouped: breaking
# changes, new features, bug fixes, then any other section. Within a
# section the bullets keep the merge order of the lanes (mocode,
# acprior, loglik), then the consolidation's. The released sections
# below come from the HEAD blob (70476dc3), so the merge's conflict
# markers at the top of a file never reach the result, and the script
# is rerunnable. Adapted from dev/rel075-news-merge.R; only the lanes,
# the packages and their versions change.
#
# A bullet marked "BREAKING:" moves to "Breaking changes" and loses the
# prefix. A bullet under no heading goes to the package's `loose`
# section. A lane's "- " bullet is read as "* ", and a lane's
# "Extension interface" or "For extension authors" section joins "New
# features", as at 0.74.0. A lane's own lead paragraph is replaced by
# the consolidation's lead, which names the version (a lead of NULL
# keeps the lane's). Then the sentences over 25 words are split
# (dev/rel076-news-split.R), every bullet paragraph of the new sections
# is rewrapped to 72 columns, and the new sections name no lane and no
# `dev/` path (the merge stops if one is left): the trail is
# dev/rel076-news-map.md, which dev/round-20261027.md pastes.
#
#   Rscript dev/rel076-news-merge.R
root <- "C:/Users/adf44/source/r"
rel <- Sys.getenv("REL076_ROOT", file.path(root, "frmtmb-wt-rel076"))
git_dir <- Sys.getenv("REL076_GIT", rel)
setwd(rel)
lanes <- strsplit(Sys.getenv("REL076_LANES",
                             "mocode acprior loglik"), " ")[[1]]
records <- c(mocode = "`dev/mocode-findings.md`",
             acprior = "`dev/acprior-findings.md`",
             loglik = "`dev/loglik-findings.md`")

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
  end <- if (length(h) >= 2L) h[2] - 1L else length(x)
  out <- x[(h[1] + 1):end]
  # a top-level "- " bullet is a "* " bullet
  sub("^- ", "* ", out)
}
split_secs <- function(v, loose) {
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
      # a package name keeps its case ("brms's `bhaz()` ...")
      first <- if (grepl("^(brms|frmtmb)", r)) substr(r, 1, 1) else
        toupper(substr(r, 1, 1))
      b[1] <- paste0(sub("BREAKING: $", "", p), first,
                     substr(r, 2, nchar(r)))
      return(b)
    }
  }
  NULL
}
trim_blank <- function(x) {
  while (length(x) && x[1] == "") x <- x[-1]
  while (length(x) && x[length(x)] == "") x <- x[-length(x)]
  x
}
newsmap <- list()
# A section name a lane used that means the same as one of the release's
# headings
same_as <- c("Extension interface" = "New features",
             "For extension authors" = "New features")
build <- function(pkg_dir, pkg, version, prev_version, lead, loose,
                  order_first) {
  f <- file.path(pkg_dir, "NEWS.md")
  rel_f <- sub(paste0(rel, "/?"), "", f)
  sec <- list()
  count <- integer(0)
  lane_lead <- character(0)
  src <- c(lanes, "consolidation")
  for (l in src) {
    if (l == "consolidation") {
      ef <- file.path("dev", paste0("rel076-news-", sub("^frmtmb[.]?", "",
                                    pkg), "-extra.md"))
      if (pkg == "frmtmb") ef <- "dev/rel076-news-core-extra.md"
      d <- if (file.exists(ef)) rd(ef) else NULL
      rec <- "`dev/round-20261027.md`"
    } else {
      d <- dev_section(file.path(root, paste0("frmtmb-wt-", l), rel_f))
      rec <- records[[l]]
    }
    if (is.null(d)) next
    s <- split_secs(d, loose)
    if (!length(lane_lead)) lane_lead <- trim_blank(s$lead)
    n <- 0L
    for (k in names(s$secs)) {
      for (b in bullets(s$secs[[k]])) {
        kk <- if (k %in% names(same_as)) same_as[[k]] else k
        u <- unbreak(b)
        if (!is.null(u)) {
          b <- u
          kk <- "Breaking changes"
          cat("  moved to Breaking changes (", l, "): ", b[1], "\n", sep = "")
        }
        sec[[kk]] <- c(sec[[kk]], list(b))
        n <- n + 1L
        hd <- gsub("[*]", "", sub("^[*] ", "", b[1]))
        newsmap[[length(newsmap) + 1L]] <<- c(pkg, k, l, rec,
                                               substr(hd, 1, 60))
      }
    }
    count[l] <- n
  }
  if (is.null(lead)) lead <- lane_lead
  keys <- c(intersect(order_first, names(sec)),
            setdiff(names(sec), order_first))
  out <- c(paste0("# ", pkg, " ", version), "", lead, "")
  for (k in keys) {
    out <- c(out, paste("##", k), "")
    for (b in sec[[k]]) out <- c(out, b, "")
  }
  if (is.null(prev_version)) {
    res <- out
  } else {
    x <- sub("\r$", "", system2("git", c("-C", git_dir, "show",
                                         paste0("HEAD:", rel_f)),
                                stdout = TRUE))
    prev <- which(x == paste0("# ", pkg, " ", prev_version))
    stopifnot(length(prev) == 1L)
    res <- c(out, x[prev:length(x)])
  }
  res <- res[!(c(FALSE, res[-1] == "" & res[-length(res)] == ""))]
  while (res[length(res)] == "") res <- res[-length(res)]
  stopifnot(!any(grepl("^(<<<<<<<|=======|>>>>>>>)", res)),
            !any(grepl("^# .*development version", res)))
  wr(res, f)
  cat(pkg, version, ": sections", paste(keys, collapse = " | "), "\n")
  cat("  bullets per source:", paste(names(count), count, collapse = ", "),
      "| total:", sum(lengths(sec)), "\n")
  for (k in keys) cat("  ", k, ": ", length(sec[[k]]), "\n", sep = "")
}
order_first <- c("Breaking changes", "New features", "Bug fixes",
                 "Performance", "Documentation")
lead_of <- function(nm) rd(file.path(rel, "dev",
                                     paste0("rel076-news-", nm, "-lead.md")))
ext <- function(p) file.path(rel, "extensions", p)
build(rel, "frmtmb", "0.76.0", "0.75.0", lead_of("core"), "New features",
      order_first)
build(ext("frmtmb.sample"), "frmtmb.sample", "0.23.0", "0.22.0",
      lead_of("sample"), "Other changes", order_first)
build(ext("frmtmb.spatial"), "frmtmb.spatial", "0.3.0", "0.2.0",
      lead_of("spatial"), "New features", order_first)
news_files <- c("NEWS.md", file.path("extensions", paste0("frmtmb.",
  c("sample", "spatial")), "NEWS.md"))

top_of <- function(x) {
  h <- grep("^# ", x)
  list(a = h[1] + 1L, b = if (length(h) >= 2L) h[2] - 1L else length(x))
}
# the sentences over 25 words, split (dev/rel076-news-split.R); each
# edit is a fixed text matched once in the file's top section with any
# run of whitespace standing for a space, so it does not depend on
# where the lane wrapped its lines
ws_regex <- function(s) {
  parts <- strsplit(trimws(s), "\\s+")[[1]]
  esc <- gsub("([][{}()+*^$|\\\\?.])", "\\\\\\1", parts)
  paste(esc, collapse = "\\s+")
}
sp <- file.path(rel, "dev/rel076-news-split.R")
splits <- if (file.exists(sp)) source(sp)$value else list()
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
    top <- sub(re, gsub("\\", "\\\\", gsub("\\s+", " ", trimws(e[2])),
                        fixed = TRUE), top, perl = TRUE)
  }
  orig <- x[t$a:t$b]
  trail <- length(orig) - max(which(nzchar(orig)))
  body <- strsplit(top, "\n", fixed = TRUE)[[1]]
  body <- c(body[seq_len(max(which(nzchar(body))))], rep("", trail))
  rest <- if (t$b < length(x)) x[(t$b + 1L):length(x)] else character(0)
  wr(c(x[seq_len(t$a - 1L)], body, rest), fp)
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
  # from the first "## " heading on, or the first bullet: the lead above
  # it is the consolidation's text, not a bullet
  st <- grep("^(## |[*] )", x[t$a:t$b])[1]
  t$a <- t$a - 1L + st
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
  rest <- if (t$b < length(x)) x[(t$b + 1L):length(x)] else character(0)
  c(x[seq_len(t$a - 1L)], out, rest)
}
# no lane name, no "lane", no `dev/` path in a new section. Code spans
# are left out first, and a name followed by "(" is a function: lane
# bhaz is named after brms's bhaz() term, which the NEWS names often
mentions <- function(x) {
  x <- gsub("`[^`]*`", "CODE", paste(x, collapse = "\n"))
  pat <- paste0("\\blanes?\\b|dev/|(?<![.])\\b(",
                paste(lanes, collapse = "|"), ")\\b(?![(])")
  sum(gregexpr(pat, x, ignore.case = TRUE, perl = TRUE)[[1]] > 0)
}
for (nf in news_files) {
  fp <- file.path(rel, nf)
  x <- rewrap_top(rd(fp))
  t <- top_of(x)
  left <- mentions(x[t$a:t$b])
  if (left > 0L) stop(nf, ": ", left, " lane or dev/ mentions left")
  # a lane wrote "the frmtmb release that ..." where it could not write
  # the number; the splits put the number in
  vague <- grep(paste0("frmtmb release|next frmtmb|next release|",
                       "that frmtmb[.]|the frmtmb that"),
                x[t$a:t$b], value = TRUE)
  if (length(vague)) stop(nf, ": a release without its number: ", vague[1])
  wr(x, fp)
  cat("rewrapped", nf, "| lane or dev/ mentions:", left,
      "| lines over 80:", sum(nchar(x[t$a:t$b]) > 80),
      "| em dashes:", sum(grepl("\u2014", x[t$a:t$b])), "\n")
}
m <- do.call(rbind, newsmap)
out <- c("| package | section (as filed) | source | record | bullet |",
         "|---|---|---|---|---|",
         sprintf("| %s | %s | %s | %s | %s |", m[, 1], m[, 2], m[, 3],
                 m[, 4], gsub("|", "\\|", m[, 5], fixed = TRUE)))
wr(out, file.path(rel, "dev/rel076-news-map.md"))
cat("wrote dev/rel076-news-map.md:", nrow(m), "bullets\n")
cat("NEWS MERGE DONE\n")
