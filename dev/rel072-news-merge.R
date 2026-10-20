# Build the 0.72.0 NEWS sections from the lanes' development-version
# sections (read from each lane's worktree, which is not edited) and the
# consolidation's leads (dev/rel072-news-*-lead.md), regrouped:
# breaking changes, new features, bug fixes, then any other section.
# Within a section the bullets keep the merge order of the lanes
# (testiso, catmulti, prilabel, sampleprior, bin2port, ndaterm, illse2,
# adhess, drop1raw). The released sections below come from the HEAD
# blob (47cf659d, 0.71.0), so the merge's conflict markers at the top
# of the file never reach the result, and the script is rerunnable.
# Adapted from dev/rel071-news-merge.R.
#
# A bullet marked "BREAKING:" moves to "Breaking changes" and loses the
# prefix. A lane bullet outside any section is a new feature unless
# `place` below names its section. A lane's own lead paragraph is
# replaced by the consolidation's lead, which names the version. The
# fixed-string edits in `subst` replace a lane's "the next release"
# wording with the version, each matched exactly once.
#
#   Rscript dev/rel072-news-merge.R
root <- "C:/Users/adf44/source/r"
rel <- Sys.getenv("REL072_ROOT", file.path(root, "frmtmb-wt-release"))
git_dir <- Sys.getenv("REL072_GIT", rel)
setwd(rel)
lanes <- strsplit(Sys.getenv("REL072_LANES",
  "testiso catmulti prilabel sampleprior bin2port ndaterm illse2 adhess drop1raw"),
  " ")[[1]]

rd <- function(f) sub("\r$", "", readLines(f, warn = FALSE))
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
# a lane bullet outside any section, placed by hand (lane/package)
place <- c("drop1raw/frmtmb.latent" = "Bug fixes")
# lane wording that names a release to come; each must match once
subst <- list(
  # lane catmulti's bare multinomial() is in the release, and the bin-2
  # rows the parenthesis waited for hold (dev/rel072-bin2-moves.R)
  "NEWS.md" = list(
    c("`brmsfamily(\"dirichlet2\")` and `logistic_normal()`
  (`multinomial()` once it reads its category count from the data),
  shims", "`brmsfamily(\"dirichlet2\")`, `logistic_normal()` and
  `multinomial()`, shims")),
  "extensions/frmtmb.latent/NEWS.md" = list(
    c("This needs the frmtmb release with the new `update()`.",
      "This needs frmtmb 0.72.0.")))
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
    loose <- place[paste0(l, "/", pkg)]
    s <- split_secs(d, if (is.na(loose)) "New features" else loose)
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
        # the trail the NEWS no longer carries (dev/rel072-news-strip.R)
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
  s <- paste(res, collapse = "\n")
  for (e in subst[[rel_f]]) {
    k <- sum(gregexpr(e[1], s, fixed = TRUE)[[1]] > 0)
    if (k != 1L) stop(rel_f, ": expected 1 match, found ", k, ": ", e[1])
    s <- sub(e[1], e[2], s, fixed = TRUE)
  }
  res <- strsplit(s, "\n", fixed = TRUE)[[1]]
  stopifnot(!any(grepl("^(<<<<<<<|=======|>>>>>>>)", res)),
            !any(grepl("^# .*development version", res)),
            !any(grepl("frmtmb release|next frmtmb", res[seq_along(out)])))
  con <- file(f, "wb"); writeLines(res, con, sep = "\n"); close(con)
  cat(pkg, version, ": sections", paste(keys, collapse = " | "), "\n")
  cat("  bullets per lane:", paste(names(count), count, collapse = ", "),
      "| total:", sum(lengths(sec)), "\n")
  for (k in keys) cat("  ", k, ": ", length(sec[[k]]), "\n", sep = "")
}
order_first <- c("Breaking changes", "New features", "Bug fixes",
                 "Performance", "Documentation")
src <- function(nm) rd(file.path(rel, "dev", paste0("rel072-news-", nm, ".md")))

build(rel, "frmtmb", "0.72.0", "0.71.0", src("core-lead"), order_first)
build(file.path(rel, "extensions/frmtmb.sample"), "frmtmb.sample",
      "0.20.0", "0.19.0", src("sample-lead"), order_first)
build(file.path(rel, "extensions/frmtmb.eam"), "frmtmb.eam",
      "0.13.0", "0.12.0", src("eam-lead"), order_first)
build(file.path(rel, "extensions/frmtmb.latent"), "frmtmb.latent",
      "0.6.2", "0.6.1", src("latent-lead"), order_first)
build(file.path(rel, "extensions/frmtmb.coupling"), "frmtmb.coupling",
      "0.6.1", "0.6.0", src("coupling-lead"), order_first)
# the sentences over 25 words, split (dev/rel072-news-split.R)
splits <- source(file.path(rel, "dev/rel072-news-split.R"))$value
for (nf in names(splits)) {
  fp <- file.path(rel, nf)
  s <- paste(rd(fp), collapse = "\n")
  for (e in splits[[nf]]) {
    k <- sum(gregexpr(e[1], s, fixed = TRUE)[[1]] > 0)
    if (k != 1L) stop(nf, ": expected 1 match, found ", k, ": ",
                      substr(e[1], 1, 60))
    s <- sub(e[1], e[2], s, fixed = TRUE)
  }
  con <- file(fp, "wb")
  writeLines(strsplit(s, "\n", fixed = TRUE)[[1]], con, sep = "\n")
  close(con)
  cat("split", length(splits[[nf]]), "sentences in", nf, "\n")
}
# The release review (dev/reviews/2026-10-20-release.md): four bullets
# move to "Breaking changes" (dev/rel072-news-moves.R), "relevelled"
# is US "releveled", and every bullet paragraph of the new sections is
# rewrapped to 72 columns, since the splits above left lines over 80
# and wraps that break mid-sentence
moves <- source(file.path(rel, "dev/rel072-news-moves.R"))$value
spell <- list("NEWS.md" = "A relevelled or",
              "extensions/frmtmb.eam/NEWS.md" = "A relevelled factor")
top_of <- function(x) {
  h <- grep("^# ", x)
  list(a = h[1] + 1L, b = h[2] - 1L)
}
wrap_par <- function(lines) {
  lead <- if (startsWith(lines[1], "* ")) "* " else
    if (startsWith(lines[1], "  - ")) "  - " else "  "
  cont <- if (lead == "  - ") "    " else "  "
  body <- substring(lines[1], nchar(lead) + 1L)
  txt <- trimws(paste(c(body, trimws(lines[-1])), collapse = " "))
  words <- strsplit(txt, " +")[[1]]
  out <- character(0)
  cur <- paste0(lead, words[1])
  for (w in words[-1]) {
    cand <- paste(cur, w)
    if (nchar(cand) > 72) {
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
  # from the first "## " heading on: the lead above it (a lane list in
  # core) is the consolidation's text, not a bullet
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
for (nf in unique(c(names(moves), names(spell), names(splits)))) {
  fp <- file.path(rel, nf)
  s <- paste(rd(fp), collapse = "\n")
  for (m in moves[[nf]]) {
    k <- sum(gregexpr(m$remove, s, fixed = TRUE)[[1]] > 0)
    if (k != 1L) stop(nf, ": move, expected 1 match, found ", k, ": ",
                      substr(m$remove, 1, 60))
    s <- sub(m$remove, "", s, fixed = TRUE)
    s <- sub("\n## New features\n",
             paste0("\n", m$bullet, "\n\n## New features\n"), s,
             fixed = TRUE)
  }
  if (!is.null(spell[[nf]])) {
    stopifnot(sum(gregexpr(spell[[nf]], s, fixed = TRUE)[[1]] > 0) == 1L)
    s <- sub(spell[[nf]], sub("relevelled", "releveled", spell[[nf]]), s,
             fixed = TRUE)
  }
  x <- rewrap_top(strsplit(s, "\n", fixed = TRUE)[[1]])
  t <- top_of(x)
  long <- sum(nchar(x[t$a:t$b]) > 80)
  con <- file(fp, "wb")
  writeLines(x, con, sep = "\n")
  close(con)
  cat("moved", length(moves[[nf]]), "bullets and rewrapped", nf,
      "| lines over 80 columns in its new section:", long, "\n")
}
# The user's decision (option b): no lane name and no `dev/` path in
# the 0.72.0 sections (dev/rel072-news-strip.R); the trail goes to
# dev/rel072-news-map.md, which dev/round-20261020.md pastes
source(file.path(rel, "dev/rel072-news-strip.R"))
for (nf in c("NEWS.md", "extensions/frmtmb.sample/NEWS.md",
             "extensions/frmtmb.eam/NEWS.md",
             "extensions/frmtmb.latent/NEWS.md",
             "extensions/frmtmb.coupling/NEWS.md")) {
  fp <- file.path(rel, nf)
  x <- rd(fp)
  t <- top_of(x)
  top <- paste(x[t$a:t$b], collapse = "\n")
  for (e in strip_edits[[nf]]) {
    re <- ws_regex(e[1])
    k <- sum(gregexpr(re, top, perl = TRUE)[[1]] > 0)
    if (k != 1L) stop(nf, ": strip, expected 1 match, found ", k, ": ",
                      substr(e[1], 1, 60))
    top <- sub(re, e[2], top, perl = TRUE)
  }
  # strsplit() drops trailing empty strings: keep the blank lines that
  # end the section, before the next release's heading
  orig <- x[t$a:t$b]
  trail <- length(orig) - max(which(nzchar(orig)))
  body <- strsplit(top, "\n", fixed = TRUE)[[1]]
  body <- c(body[seq_len(max(which(nzchar(body))))], rep("", trail))
  x <- c(x[seq_len(t$a - 1L)], body, x[(t$b + 1L):length(x)])
  x <- rewrap_top(x)
  t <- top_of(x)
  left <- mentions(x[t$a:t$b])
  if (left > 0L) stop(nf, ": ", left, " lane or dev/ mentions left")
  con <- file(fp, "wb")
  writeLines(x, con, sep = "\n")
  close(con)
  cat("stripped", length(strip_edits[[nf]]), "provenance notes from", nf,
      "| lane or dev/ mentions left:", left, "| lines over 80:",
      sum(nchar(x[t$a:t$b]) > 80), "\n")
}
m <- do.call(rbind, newsmap)
out <- c("| package | section (as the lane filed it) | lane | record | bullet |",
         "|---|---|---|---|---|",
         sprintf("| %s | %s | %s | %s | %s |", m[, 1], m[, 2], m[, 3],
                 m[, 4], gsub("|", "\\|", m[, 5], fixed = TRUE)))
con <- file(file.path(rel, "dev/rel072-news-map.md"), "wb")
writeLines(out, con, sep = "\n")
close(con)
cat("wrote dev/rel072-news-map.md:", nrow(m), "bullets\n")
cat("NEWS MERGE DONE\n")
