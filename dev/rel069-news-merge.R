# Build the 0.69.0 NEWS sections from the lanes' development-version
# sections (read from each lane's worktree, which is not edited) and the
# consolidation's own lead and bullets (dev/rel069-news-*.md), regrouped:
# breaking changes, new features, bug fixes, performance, extension
# API. Within a section the bullets keep the merge order of the lanes
# (ciharden, surface, setier, optima), then the consolidation's. Each
# package's top section (everything above its previous release heading)
# is replaced, so the script is rerunnable.
#
#   Rscript dev/rel069-news-merge.R
root <- "C:/Users/adf44/source/r"
rel <- file.path(root, "frmtmb-wt-release")
lanes <- c("ciharden", "surface", "setier", "optima")

dev_section <- function(f) {
  x <- sub("\r$", "", readLines(f, warn = FALSE))
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
starts <- function(b, s) startsWith(b[1], s)
build <- function(pkg_dir, pkg, version, lead, extra, moves = list(),
                  drop = character(0), order_first) {
  f <- file.path(pkg_dir, "NEWS.md")
  sec <- list()
  count <- integer(0)
  for (l in lanes) {
    lf <- file.path(root, paste0("frmtmb-wt-", l),
                    sub(paste0(rel, "/?"), "", pkg_dir), "NEWS.md")
    d <- dev_section(lf)
    if (is.null(d)) next
    s <- split_secs(d)
    n <- 0L
    for (k in names(s$secs)) {
      for (b in bullets(s$secs[[k]])) {
        if (any(vapply(drop, function(x) starts(b, x), NA))) {
          cat("  dropped (", l, "): ", b[1], "\n", sep = "")
          next
        }
        kk <- k
        for (m in moves) if (starts(b, m[1])) kk <- m[2]
        sec[[kk]] <- c(sec[[kk]], list(b))
        n <- n + 1L
      }
    }
    count[l] <- n
  }
  e <- split_secs(extra)
  ne <- 0L
  for (k in names(e$secs)) {
    for (b in bullets(e$secs[[k]])) {
      sec[[k]] <- c(sec[[k]], list(b))
      ne <- ne + 1L
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
  x <- sub("\r$", "", readLines(f, warn = FALSE))
  h <- grep("^# ", x)
  # the previous release heading: the first "# <pkg> <number>" line of
  # another version (the merge left conflict markers above it)
  prev <- h[grepl(paste0("^# ", pkg, " [0-9]"), x[h]) &
              x[h] != paste0("# ", pkg, " ", version)][1]
  res <- c(out, x[prev:length(x)])
  res <- res[!(c(FALSE, res[-1] == "" & res[-length(res)] == ""))]
  con <- file(f, "wb"); writeLines(res, con, sep = "\n"); close(con)
  cat(pkg, version, ": sections", paste(keys, collapse = " | "), "\n")
  cat("  bullets per lane:", paste(names(count), count, collapse = ", "),
      "| consolidation:", ne, "| total:",
      sum(lengths(sec)), "\n")
}
order_first <- c("Breaking changes", "New features", "Bug fixes",
                 "Performance")
src <- function(nm) {
  sub("\r$", "", readLines(file.path(rel, "dev",
                                     paste0("rel069-news-", nm, ".md")),
                           warn = FALSE))
}

build(rel, "frmtmb", "0.69.0", src("core-lead"), src("core-extra"),
      # lane ciharden's kriging bullet is a memory change
      moves = list(c("* The kriging covariance of a `gp()` term no longer",
                     "Performance")),
      order_first = order_first)
build(file.path(rel, "extensions/frmtmb.sample"), "frmtmb.sample",
      "0.17.0", src("sample-lead"), src("sample-extra"),
      # a defect of the development version only, which no release had
      drop = "* `check_laplace()` stopped with \"length(ml) == length(keep)",
      order_first = order_first)
cat("NEWS MERGE DONE\n")
