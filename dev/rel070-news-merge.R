# Build the 0.70.0 NEWS sections from the lanes' development-version
# sections (read from each lane's worktree, which is not edited) and the
# consolidation's own lead and bullets (dev/rel070-news-*.md), regrouped:
# breaking changes, new features, bug fixes, performance, extension
# API. Within a section the bullets keep the merge order of the lanes
# (moparity, optrap, seflags), then the consolidation's. The released
# sections below come from the HEAD blob (2884ccbd, 0.69.1), so the
# merge's conflict markers at the top of the file never reach the
# result, and the script is rerunnable.
#
# Lane seflags marks two bug-fix bullets "BREAKING:"; they move to
# "Breaking changes", and the prefix, which the heading now says, goes.
#
#   Rscript dev/rel070-news-merge.R
root <- "C:/Users/adf44/source/r"
rel <- file.path(root, "frmtmb-wt-release")
setwd(rel)
lanes <- c("moparity", "optrap", "seflags")

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
  for (i in seq_along(h)) {
    end <- if (i < length(h)) h[i + 1] - 1 else length(v)
    body <- if (end > h[i]) v[(h[i] + 1):end] else character(0)
    secs[[sub("^## ", "", v[h[i]])]] <- body
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
# "* **BREAKING: an alias ..." -> "* **An alias ...", rewrapped is not
# needed: the line only gets shorter
unbreak <- function(b) {
  if (!startsWith(b[1], "* **BREAKING: ")) return(NULL)
  r <- sub("^[*] [*][*]BREAKING: ", "", b[1])
  b[1] <- paste0("* **", toupper(substr(r, 1, 1)), substr(r, 2, nchar(r)))
  b
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
  stopifnot(!any(grepl("^(<<<<<<<|=======|>>>>>>>)", res)),
            !any(grepl("^# .*development version", res)))
  con <- file(f, "wb"); writeLines(res, con, sep = "\n"); close(con)
  cat(pkg, version, ": sections", paste(keys, collapse = " | "), "\n")
  cat("  bullets per lane:", paste(names(count), count, collapse = ", "),
      "| consolidation:", ne, "| total:", sum(lengths(sec)), "\n")
  for (k in keys) cat("  ", k, ": ", length(sec[[k]]), "\n", sep = "")
}
order_first <- c("Breaking changes", "New features", "Bug fixes",
                 "Performance")
src <- function(nm) rd(file.path(rel, "dev", paste0("rel070-news-", nm,
                                                     ".md")))

build(rel, "frmtmb", "0.70.0", "0.69.1", src("core-lead"),
      src("core-extra"), order_first = order_first)
build(file.path(rel, "extensions/frmtmb.sample"), "frmtmb.sample",
      "0.18.0", "0.17.0", src("sample-lead"), character(0),
      order_first = order_first)
cat("NEWS MERGE DONE\n")
