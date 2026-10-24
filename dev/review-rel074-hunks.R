# Release review of 0.74.0, item 2: is any lane's change lost in the
# merge? Independent of dev/rel074-merge.sh: a line-multiset test per
# file, not a re-run of git merge-file.
#   For each lane L and each tracked file F that L changed (not deleted):
#   - every line L added (counted with multiplicity) must occur in the
#     release file at least that many times more than in the base,
#     after subtracting what the other lanes removed;
#   - every line L removed must occur in the release file no more often
#     than base - removed_by_L + added_by_others.
#   Deleted files must be absent. Untracked lane files must equal the
#   release copy (CR stripped), or be listed.
#   NAMESPACE: every export a lane added is in the release NAMESPACE and
#   every export a lane removed is not.
#   Rscript dev/review-rel074-hunks.R
R <- "C:/Users/adf44/source/r"
rel <- file.path(R, "frmtmb-wt-release")
lanes <- c("splitsurv", "mvnlatent", "spatial")
base <- "19b162b9"
git <- function(dir, ...) {
  suppressWarnings(system2("git", c("-C", shQuote(dir), ...), stdout = TRUE,
                           stderr = FALSE))
}
rdl <- function(f) if (file.exists(f)) sub("\r$", "", readLines(f, warn = FALSE)) else character()
baseblob <- function(f) sub("\r$", "", git(rel, "show", paste0(base, ":", f)))
tab <- function(x) if (length(x)) table(x) else table(character())
cnt <- function(t, k) { v <- t[k]; v[is.na(v)] <- 0L; as.integer(v) }
st <- lapply(setNames(lanes, lanes), function(l) {
  x <- git(file.path(R, paste0("frmtmb-wt-", l)), "diff", "--name-status",
           "--no-renames", base)
  data.frame(st = sub("\t.*", "", x), f = sub("^[^\t]*\t", "", x))
})
gen <- "^(man/.*[.]Rd|extensions/[^/]*/man/.*[.]Rd|R/links-brms[.]R|extensions/frmtmb.sample/tests/testthat/helper-brms-suite[.]R)$"
# per lane, per file: added and removed multisets
delta <- list()
for (l in lanes) {
  for (i in seq_len(nrow(st[[l]]))) {
    f <- st[[l]]$f[i]
    if (st[[l]]$st[i] != "M" || grepl(gen, f)) next
    b <- baseblob(f)
    n <- rdl(file.path(R, paste0("frmtmb-wt-", l), f))
    tb <- tab(b); tn <- tab(n)
    keys <- union(names(tb), names(tn))
    d <- cnt(tn, keys) - cnt(tb, keys)
    names(d) <- keys
    delta[[f]][[l]] <- d[d != 0]
  }
}
nbad <- 0L
for (f in names(delta)) {
  b <- baseblob(f); r <- rdl(file.path(rel, f))
  tb <- tab(b); tr <- tab(r)
  for (l in names(delta[[f]])) {
    d <- delta[[f]][[l]]
    oth <- Reduce(function(a, o) {
      k <- union(names(a), names(o)); v <- cnt(a, k) + cnt(o, k)
      names(v) <- k; v
    }, delta[[f]][setdiff(names(delta[[f]]), l)], integer())
    add <- d[d > 0]; rem <- -d[d < 0]
    for (k in names(add)) {
      o <- if (length(oth)) cnt(oth, k) else 0L
      if (cnt(tr, k) < cnt(tb, k) + add[[k]] + min(o, 0L)) {
        nbad <- nbad + 1L
        cat("LOST ADD [", l, "] ", f, ": ", substr(k, 1, 110), "\n", sep = "")
      }
    }
    for (k in names(rem)) {
      o <- if (length(oth)) cnt(oth, k) else 0L
      if (cnt(tr, k) > cnt(tb, k) - rem[[k]] + max(o, 0L)) {
        nbad <- nbad + 1L
        cat("KEPT REMOVED [", l, "] ", f, ": ", substr(k, 1, 110), "\n",
            sep = "")
      }
    }
  }
  # lines in release that no lane added and the base lacks: resolve edits
  allk <- unique(unlist(lapply(delta[[f]], names)))
  addany <- unlist(lapply(delta[[f]], function(d) names(d[d > 0])))
  extra <- setdiff(names(tr)[cnt(tr, names(tr)) > cnt(tb, names(tr))], addany)
  if (length(extra)) {
    cat("NEW IN RELEASE (no lane added) ", f, ": ", length(extra), " lines\n",
        sep = "")
    cat(paste0("    ", substr(head(extra, 8), 1, 110)), sep = "\n")
  }
}
cat("tracked-file check: files", length(delta), "; problems", nbad, "\n")
# deletions
for (l in lanes) {
  dd <- st[[l]]$f[st[[l]]$st == "D"]
  pres <- dd[file.exists(file.path(rel, dd))]
  cat("lane", l, "deleted", length(dd), "; still present in release:",
      length(pres), paste(pres, collapse = " "), "\n")
}
# untracked
uskip <- "^(man/.*[.]Rd|extensions/[^/]*/man/.*[.]Rd|.*stan-cache.*|dev/ciharden-out/.*|dev/brms-suite/.*|dev/review-.*|dev/core-boundary[.]md)$"
for (l in lanes) {
  lt <- file.path(R, paste0("frmtmb-wt-", l))
  u <- git(lt, "ls-files", "--others", "--exclude-standard")
  sk <- grepl(uskip, u); u2 <- u[!sk]
  miss <- character(); dif <- character()
  for (f in u2) {
    rf <- file.path(rel, f)
    if (!file.exists(rf)) { miss <- c(miss, f); next }
    a <- readBin(file.path(lt, f), "raw", file.size(file.path(lt, f)))
    b <- readBin(rf, "raw", file.size(rf))
    a <- a[a != as.raw(13)]; b <- b[b != as.raw(13)]
    if (!identical(a, b)) dif <- c(dif, f)
  }
  cat("lane", l, "untracked", length(u), "skipped", sum(sk), "missing in release",
      length(miss), "differ", length(dif), "\n")
  if (length(miss)) cat(paste0("   MISSING ", miss), sep = "\n")
  if (length(dif)) cat(paste0("   DIFFER ", dif), sep = "\n")
  # skipped review scripts and other skipped files, for the record
  cat("   skipped:", paste(head(u[sk & !grepl("stan-cache|brms-suite", u)], 40), collapse = " "), "\n")
}
# NAMESPACE files
for (ns in c("NAMESPACE", file.path("extensions", c("frmtmb.latent", "frmtmb.sample", "frmtmb.learn", "frmtmb.survival", "frmtmb.spatial"), "NAMESPACE"))) {
  r <- rdl(file.path(rel, ns)); b <- baseblob(ns)
  for (l in lanes) {
    n <- rdl(file.path(R, paste0("frmtmb-wt-", l), ns))
    if (!length(n)) next
    add <- setdiff(n, b); rem <- setdiff(b, n)
    la <- setdiff(add, r); lr <- intersect(rem, r)
    if (length(la) || length(lr)) {
      cat("NAMESPACE", ns, "lane", l, ": added but absent:", la,
          "| removed but present:", lr, "\n")
    }
  }
}
cat("NAMESPACE check done\n")
