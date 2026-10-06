# Build the 0.68.0 NEWS sections from the lanes' development-version
# sections (read from each lane's worktree, which is not edited) and the
# consolidation's own bullets (below), regrouped as at 0.67.0: breaking
# changes first, then new features, then bug fixes, then the extension
# API. Each package's top section (everything above its previous
# release heading) is replaced.
#
#   Rscript dev/rel068-news-merge.R
root <- "C:/Users/adf44/source/r"
rel <- file.path(root, "frmtmb-wt-release")
# lane nanse, merged last, lists its breaking bullets first (the nanse
# review's merge recipe, steps 2 and 7)
lanes <- c("fixes", "gpby", "ordmix")

dev_section <- function(f) {
  x <- readLines(f, warn = FALSE)
  h <- grep("^# ", x)
  stopifnot(grepl("development version", x[h[1]]))
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
# a bullet moved to another section, keyed by its first line
move_to <- function(b, rules) {
  for (r in rules) if (startsWith(b[1], r$start)) return(r$to)
  NA_character_
}
edit_text <- function(lines, old, new) {
  s <- paste(lines, collapse = "\n")
  if (!grepl(old, s, fixed = TRUE)) stop("edit not found: ", substr(old, 1, 60))
  strsplit(sub(old, new, s, fixed = TRUE), "\n", fixed = TRUE)[[1]]
}
build <- function(pkg_dir, pkg, version, lead, extra, moves = list(),
                  edits = list(), order_first, lanes_here,
                  lead_to = NULL) {
  f <- file.path(pkg_dir, "NEWS.md")
  sec <- list()
  for (l in lanes_here) {
    lf <- file.path(root, paste0("frmtmb-wt-", l),
                    sub(paste0(rel, "/?"), "", pkg_dir), "NEWS.md")
    s <- split_secs(dev_section(lf))
    # a lane section with bullets and no "##" heading (lane nanse's
    # frmtmb.spline) files them under `lead_to`
    if (!is.null(lead_to)) {
      for (b in bullets(s$lead)) sec[[lead_to]] <- c(sec[[lead_to]], list(b))
    }
    for (k in names(s$secs)) {
      for (b in bullets(s$secs[[k]])) {
        to <- move_to(b, moves)
        kk <- if (is.na(to)) k else to
        sec[[kk]] <- c(sec[[kk]], list(b))
      }
    }
  }
  e <- split_secs(extra)
  for (k in names(e$secs)) {
    for (b in bullets(e$secs[[k]])) sec[[k]] <- c(sec[[k]], list(b))
  }
  keys <- c(intersect(order_first, names(sec)),
            setdiff(names(sec), c(order_first, "Extension API")),
            intersect("Extension API", names(sec)))
  out <- c(paste0("# ", pkg, " ", version), "", lead, "")
  for (k in keys) {
    out <- c(out, paste("##", k), "")
    for (b in sec[[k]]) out <- c(out, b, "")
  }
  for (ed in edits) out <- edit_text(out, ed[1], ed[2])
  x <- readLines(f, warn = FALSE)
  h <- grep("^# ", x)
  # the release's own file holds the merge's conflicted top, or this
  # script's own section from an earlier run; the previous release
  # heading is the first "# <pkg> <number>" line of another version
  prev <- h[grepl(paste0("^# ", pkg, " [0-9]"), x[h]) &
              x[h] != paste0("# ", pkg, " ", version)][1]
  res <- c(out, x[prev:length(x)])
  res <- res[!(c(FALSE, res[-1] == "" & res[-length(res)] == ""))]
  writeLines(res, f)
  cat(pkg, version, ": sections", paste(keys, collapse = " | "), "\n")
}
order_first <- c("Breaking changes", "New features", "Performance", "Bug fixes",
                 "Bug fixes and changes")
src <- function(nm) {
  readLines(file.path(rel, "dev", paste0("rel068-news-", nm, ".md")),
            warn = FALSE)
}

build(rel, "frmtmb", "0.68.0", src("core-lead"), src("core-extra"),
      moves = list(list(start = "* `ord_thres_linpred()` joins the sampling API",
                        to = "Extension API")),
      edits = list(
        c("`groups = `, a\n  `hurdle_cumulative()` component beside one without a hurdle, and\n  `cs()` on a `cumulative()` component are refused.",
          "`groups = ` and a\n  `hurdle_cumulative()` component beside one without a hurdle are\n  refused."),
        c("thresholds. The log density equals brms 2.23.0's compiled program\n  to at most 1.6 ulp",
          "thresholds, with brms's warning that the effects are experimental.\n  The log density equals brms 2.23.0's compiled program to at most\n  1.6 ulp"),
        c("  frmtmb.sample's `posterior_linpred(incl_thres = TRUE)` reads per\n  draw.",
          "  frmtmb.sample's `posterior_linpred(incl_thres = TRUE)` reads per\n  draw. An ordinal mixture is refused first, in brms's words:\n  \"'incl_thres' is not supported for mixture models.\""),
        # one form for a breaking item under "Breaking changes": the
        # section says it, so no bullet repeats "BREAKING:" (release
        # review, m8)
        c("* **BREAKING: `get_prior()`", "* **`get_prior()`"),
        c("* **BREAKING: `summary()$gp`", "* **`summary()$gp`"),
        # release review, m11
        c("  `order = \"mu\"`) lands where brms puts it. brms's equidistant mixture",
          "  `order = \"mu\"`) lands where brms puts it. A sum-to-zero component\n  lists no class `\"Intercept\"` rows, where brms lists one per\n  threshold, as for one sum-to-zero family. brms's equidistant mixture")),
      order_first = order_first, lanes_here = c("nanse", lanes))
build(file.path(rel, "extensions/frmtmb.sample"), "frmtmb.sample", "0.16.0",
      src("sample-lead"), src("sample-extra"),
      edits = list(
        c("`hurdle_cumulative()` is refused: brms returns\n  the hurdle probability beside the threshold predictors times\n  `1 - hu` there, which is the predictor of nothing.",
          "`hurdle_cumulative()` is refused, a divergence\n  the user decided on 2026-10-06: brms returns the hurdle probability\n  beside the threshold predictors times `1 - hu` there, which is the\n  predictor of nothing (`dev/upstream-bugs.md`, brms-21). An ordinal\n  mixture is refused in brms's words, \"'incl_thres' is not supported\n  for mixture models.\""),
        c("* Through frmtmb, a smooth's null-space draws (`bs_sx_1`) are on brms's\n  scale, and",
          "* **Through frmtmb, a smooth's null-space draws (`bs_sx_1`) are on\n  brms's scale**, and")),
      order_first = order_first, lanes_here = lanes)
build(file.path(rel, "extensions/frmtmb.spline"), "frmtmb.spline", "0.10.0",
      src("spline-lead"), character(0),
      edits = list(
        c("* **BREAKING: the `\"Sigma\"`", "* **The `\"Sigma\"`"),
        c("* **BREAKING: `frm_curve_deriv()`", "* **`frm_curve_deriv()`"),
        # this release's floor is the frmtmb that marks the rows
        c("Needs the frmtmb release that marks those rows\n  (`frm_lp_basis()$se_nonest`); with an older frmtmb nothing changes.",
          "frmtmb 0.68.0 marks those rows\n  (`frm_lp_basis()$se_nonest`).")),
      order_first = order_first, lanes_here = c("nanse", "gpby"),
      lead_to = "Bug fixes and changes")
cat("NEWS MERGE DONE\n")
