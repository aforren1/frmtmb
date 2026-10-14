# Build the 0.71.0 NEWS sections from the lanes' development-version
# sections (read from each lane's worktree, which is not edited) and the
# consolidation's leads (dev/rel071-news-*-lead.md), regrouped: breaking
# changes, new features, bug fixes, performance, extension API. Within
# a section the bullets keep the merge order of the lanes (dblbar,
# famsimplex, famuni, postfix, optstart, illse), then the
# consolidation's. The released sections below come from the HEAD blob
# (253b3ef7, 0.70.1), so the merge's conflict markers at the top of the
# file never reach the result, and the script is rerunnable.
#
# A bullet marked "BREAKING:" moves to "Breaking changes" and loses the
# prefix, which the heading then says. A lane bullet outside any
# section is a new feature. A lane's own lead paragraph ("Needs the
# frmtmb ...") is replaced by the consolidation's lead, which names the
# version, and "the next frmtmb" in a bullet becomes "frmtmb 0.71.0".
#
#   Rscript dev/rel071-news-merge.R
root <- "C:/Users/adf44/source/r"
rel <- file.path(root, "frmtmb-wt-release")
setwd(rel)
lanes <- strsplit(Sys.getenv("REL071_LANES",
                             "dblbar famsimplex famuni postfix optstart remeasure"),
                  " ")[[1]]

rd <- function(f) sub("\r$", "", readLines(f, warn = FALSE))
dev_section <- function(f) {
  if (!file.exists(f)) return(NULL)
  x <- rd(f)
  h <- grep("^# ", x)
  if (!grepl("development version", x[h[1]])) return(NULL)
  x[(h[1] + 1):(h[2] - 1)]
}
split_secs <- function(v) {
  h <- grep("^## ", v)
  lead <- if (length(h)) v[seq_len(h[1] - 1)] else v
  secs <- list()
  # bullets before the first heading are new features
  s0 <- grep("^\\* ", lead)
  if (length(s0)) {
    secs[["New features"]] <- lead[s0[1]:length(lead)]
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
# "* **BREAKING: an alias" -> "* **An alias", "* BREAKING: `x`" ->
# "* `x`"; the line only gets shorter, so no rewrap is needed
unbreak <- function(b) {
  if (startsWith(b[1], "* **BREAKING: ")) {
    r <- sub("^[*] [*][*]BREAKING: ", "", b[1])
    b[1] <- paste0("* **", toupper(substr(r, 1, 1)), substr(r, 2, nchar(r)))
    return(b)
  }
  if (startsWith(b[1], "* BREAKING: ")) {
    r <- sub("^[*] BREAKING: ", "", b[1])
    b[1] <- paste0("* ", toupper(substr(r, 1, 1)), substr(r, 2, nchar(r)))
    return(b)
  }
  NULL
}
build <- function(pkg_dir, pkg, version, prev_version, lead, extra,
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
        # the phrase can break across lines
        b <- strsplit(gsub("by(\\s+)the next frmtmb", "by\\1frmtmb 0.71.0",
                           paste(b, collapse = "\n")), "\n",
                      fixed = TRUE)[[1]]
        sec[[kk]] <- c(sec[[kk]], list(b))
        n <- n + 1L
      }
    }
    count[l] <- n
  }
  ne <- 0L
  if (length(extra)) {
    e <- split_secs(extra)
    for (k in names(e$secs)) {
      for (b in bullets(e$secs[[k]])) {
        sec[[k]] <- c(sec[[k]], list(b))
        ne <- ne + 1L
      }
    }
  }
  keys <- c(intersect(order_first, names(sec)),
            setdiff(names(sec), c(order_first, "Extension API")),
            intersect("Extension API", names(sec)))
  out <- c(paste0("# ", pkg, " ", version), "", lead, "")
  for (k in keys) {
    out <- c(out, paste("##", k), "")
    for (b in sec[[k]]) out <- c(out, b, "")
  }
  x <- sub("\r$", "", system2("git", c("show", paste0("HEAD:", rel_f)),
                              stdout = TRUE))
  prev <- which(x == paste0("# ", pkg, " ", prev_version))
  stopifnot(length(prev) == 1L)
  res <- c(out, x[prev:length(x)])
  res <- res[!(c(FALSE, res[-1] == "" & res[-length(res)] == ""))]
  while (res[length(res)] == "") res <- res[-length(res)]
  stopifnot(!any(grepl("^(<<<<<<<|=======|>>>>>>>)", res)),
            !any(grepl("^# .*development version", res)),
            !any(grepl("next frmtmb|frmtmb of this round", res[seq_along(out)])))
  con <- file(f, "wb"); writeLines(res, con, sep = "\n"); close(con)
  cat(pkg, version, ": sections", paste(keys, collapse = " | "), "\n")
  cat("  bullets per lane:", paste(names(count), count, collapse = ", "),
      "| consolidation:", ne, "| total:", sum(lengths(sec)), "\n")
  for (k in keys) cat("  ", k, ": ", length(sec[[k]]), "\n", sep = "")
}
order_first <- c("Breaking changes", "New features", "Bug fixes",
                 "Performance")
src <- function(nm) {
  f <- file.path(rel, "dev", paste0("rel071-news-", nm, ".md"))
  if (file.exists(f)) rd(f) else character(0)
}
illse <- "illse" %in% lanes
sfx <- if (illse) "" else "-six"

build(rel, "frmtmb", "0.71.0", "0.70.1", src(paste0("core-lead", sfx)),
      src("core-extra"), order_first = order_first)
build(file.path(rel, "extensions/frmtmb.sample"), "frmtmb.sample",
      "0.19.0", "0.18.0", src("sample-lead"), character(0),
      order_first = order_first)
build(file.path(rel, "extensions/frmtmb.eam"), "frmtmb.eam",
      "0.12.0", "0.11.3", src("eam-lead"), character(0),
      order_first = order_first)
build(file.path(rel, "extensions/frmtmb.latent"), "frmtmb.latent",
      "0.6.1", "0.6.0", src("latent-lead"), character(0),
      order_first = order_first)
if (illse) {
  build(file.path(rel, "extensions/frmtmb.spline"), "frmtmb.spline",
        "0.10.2", "0.10.1", src("spline-lead"), character(0),
        order_first = order_first)
}
# the release review's sentence splits (section 10)
splits <- source(file.path(rel, "dev/rel071-news-split.R"))$value
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
cat("NEWS MERGE DONE\n")
