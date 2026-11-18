# Q3 of the codesize lane: duplication in core R/.
#   Rscript dev/codesize-dup.R [W]
# (a) Whole functions whose bodies are the same after normalization:
#     comments and whitespace dropped, every SYMBOL (variable) -> V,
#     every string -> S, every number -> N; call names kept ("strict"),
#     or also -> F ("loose").
# (b) Blocks: windows of W (default 60) normalized tokens, strict
#     normalization, that occur at two or more places (a rolling hash with
#     two moduli; a window may not cross a top-level expression). Adjacent
#     duplicated windows are merged into regions; a line is "duplicated"
#     when a duplicated window covers a token on it.
# (c) Counts of repeated patterns: condition construction and branches
#     on a family's name.
# Reads dev/codesize-out/parsed.rds and functions.tsv. Writes
# dev/codesize-out/dup-*.tsv and prints the blocks for the findings.
a <- commandArgs(trailingOnly = TRUE)
W <- if (length(a)) as.integer(a[1]) else 60L
out <- "C:/Users/adf44/source/r/frmtmb-wt-codesize/dev/codesize-out"
P <- readRDS(file.path(out, "parsed.rds"))
funs <- read.delim(file.path(out, "functions.tsv"), stringsAsFactors = FALSE)

norm <- function(pd, loose = FALSE) {
  t <- pd$token; x <- pd$text
  x[t == "SYMBOL"] <- "V"
  x[t == "STR_CONST"] <- "S"
  x[t == "NUM_CONST"] <- "N"
  x[t %in% c("SYMBOL_SUB", "SYMBOL_FORMALS")] <- "A"
  if (loose) x[t == "SYMBOL_FUNCTION_CALL"] <- "F"
  x
}

# ---- (a) whole functions ----
tok <- list()
for (f in names(P$core)) {
  x <- P$core[[f]]
  pd <- x$pd[x$pd$token != "COMMENT", ]
  t <- x$top[x$top$is_fun & !is.na(x$top$name), ]
  for (i in seq_len(nrow(t))) {
    k <- pd$line1 >= t$first[i] & pd$line1 <= t$last[i]
    q <- pd[k, ]
    # drop `name <-` so only the function itself is compared
    q <- q[-(1:2), ]
    tok[[paste0(f, ".R::", t$name[i])]] <- list(
      strict = paste(norm(q), collapse = " "),
      loose = paste(norm(q, TRUE), collapse = " "), n = nrow(q))
  }
}
fk <- data.frame(key = names(tok),
                 strict = vapply(tok, `[[`, "", "strict"),
                 loose = vapply(tok, `[[`, "", "loose"),
                 ntok = vapply(tok, `[[`, 1L, "n"))
fk$name <- sub("^.*::", "", fk$key); fk$file <- sub("::.*$", "", fk$key)
m <- match(paste(fk$file, fk$name), paste(funs$file, funs$name))
fk$code <- funs$code[m]
clusters <- function(col, min_code) {
  d <- fk[fk$code >= min_code, ]
  g <- split(d, d[[col]])
  g <- g[vapply(g, nrow, 1L) > 1]
  g[order(-vapply(g, function(z) sum(z$code) - max(z$code), 1))]
}
cat("BEGIN codesize-dup\n")
for (col in c("strict", "loose")) {
  g <- clusters(col, 5)
  sav <- sum(vapply(g, function(z) sum(z$code) - max(z$code), 1))
  cat(sprintf("\nWhole-function clones (%s, functions of 5+ code lines): %d clusters, %d functions, %d code lines beyond one copy each\n",
              col, length(g), sum(vapply(g, nrow, 1L)), sav))
  for (z in head(g, if (col == "strict") 25 else 15))
    cat(sprintf("  %d x %d lines: %s\n", nrow(z), max(z$code),
                paste0(z$name, " (", z$file, ")", collapse = ", ")))
}

# ---- (b) blocks ----
stream <- list()
for (f in names(P$core)) {
  x <- P$core[[f]]
  pd <- x$pd[x$pd$token != "COMMENT", ]
  # expression id per token, so no window crosses a top-level expression
  ex <- findInterval(pd$line1, x$top$first)
  stream[[f]] <- data.frame(file = f, line = pd$line1, ex = ex,
                            tok = norm(pd))
}
S <- do.call(rbind, stream); rownames(S) <- NULL
code <- match(S$tok, unique(S$tok))
n <- length(code); nw <- n - W + 1L
hash <- function(b, p) {
  h <- numeric(nw)
  for (k in seq_len(W)) h <- (h * b + code[k:(k + nw - 1L)]) %% p
  h
}
h1 <- hash(131, 1000000007); h2 <- hash(137, 998244353)
ok <- S$file[1:nw] == S$file[W:n] & S$ex[1:nw] == S$ex[W:n]
hk <- paste(h1, h2)
hk[!ok] <- NA
dupw <- !is.na(hk) & (duplicated(hk) | duplicated(hk, fromLast = TRUE))
cover <- integer(n + 1L)
st <- which(dupw)
cover[st] <- cover[st] + 1L
cover[st + W] <- cover[st + W] - 1L
dt <- cumsum(cover)[1:n] > 0
dl <- unique(S[dt, c("file", "line")])
codel <- unique(S[, c("file", "line")])
pf <- merge(aggregate(line ~ file, codel, length),
            aggregate(line ~ file, dl, length), by = "file", all.x = TRUE)
names(pf) <- c("file", "code_lines", "dup_lines")
pf$dup_lines[is.na(pf$dup_lines)] <- 0
pf$share <- round(100 * pf$dup_lines / pf$code_lines, 1)
pf <- pf[order(-pf$dup_lines), ]
write.table(pf, file.path(out, "dup-files.tsv"), sep = "\t", quote = FALSE,
            row.names = FALSE)
pf$subsystem <- P$sub_of[pf$file]
ps <- aggregate(cbind(code_lines, dup_lines) ~ subsystem, pf, sum)
ps$share <- round(100 * ps$dup_lines / ps$code_lines, 1)
ps <- ps[order(-ps$dup_lines), ]
cat(sprintf("\nBlock clones, W = %d tokens: %d of %d code lines (%.1f%%) lie in a window that occurs twice or more\n",
            W, nrow(dl), nrow(codel), 100 * nrow(dl) / nrow(codel)))
cat("Per subsystem:\n")
cat(sprintf("  %s: %d of %d (%.1f%%)\n", ps$subsystem, ps$dup_lines,
            ps$code_lines, ps$share), sep = "")
cat("Top 20 files:\n")
x <- head(pf, 20)
cat(sprintf("  %s.R: %d of %d (%.1f%%)\n", x$file, x$dup_lines, x$code_lines,
            x$share), sep = "")

# regions: maximal runs of duplicated window starts; a region's siblings
# are the other places of its first window
r <- rle(dupw)
e <- cumsum(r$lengths); s0 <- e - r$lengths + 1
reg <- data.frame(s = s0[r$values], e = e[r$values])
reg$last_tok <- reg$e + W - 1L
reg$file <- S$file[reg$s]
reg$l1 <- S$line[reg$s]; reg$l2 <- S$line[reg$last_tok]
reg$lines <- reg$l2 - reg$l1 + 1L
reg$h <- hk[reg$s]
grp <- split(seq_len(nrow(reg)), reg$h)
reg$copies <- vapply(reg$h, function(z) sum(hk == z, na.rm = TRUE), 1L)
reg$fun <- vapply(seq_len(nrow(reg)), function(i) {
  t <- funs[funs$file == paste0(reg$file[i], ".R") &
              funs$first <= reg$l1[i] & funs$last >= reg$l1[i], ]
  if (nrow(t)) t$name[1] else "(top level)"
}, "")
reg$sib <- vapply(seq_len(nrow(reg)), function(i) {
  j <- setdiff(which(hk == reg$h[i]), reg$s[i])
  paste(unique(paste0(S$file[j], ".R:", S$line[j])), collapse = " ")
}, "")
reg <- reg[order(-reg$lines), ]
write.table(reg[c("file", "l1", "l2", "lines", "copies", "fun", "sib")],
            file.path(out, "dup-regions.tsv"), sep = "\t", quote = FALSE,
            row.names = FALSE)
cat(sprintf("\nRegions (maximal duplicated runs): %d; 40+ lines: %d; 20-39 lines: %d\n",
            nrow(reg), sum(reg$lines >= 40), sum(reg$lines >= 20 & reg$lines < 40)))
cat("Largest 40 regions (file:lines, function, copies of the first window, where else):\n")
x <- head(reg, 40)
cat(sprintf("  %s.R:%d-%d (%d lines) in %s; %d copies; also %s\n", x$file,
            x$l1, x$l2, x$lines, x$fun, x$copies,
            substr(x$sib, 1, 160)), sep = "")

# ---- (c) patterns ----
allcode <- unlist(lapply(names(P$core), function(f) {
  x <- P$core[[f]]
  k <- unique(x$pd$line1[x$pd$token != "COMMENT"])
  x$lines[sort(k)]
}))
pat <- c(
  `stop(` = "\\bstop\\(",
  `warning(` = "\\bwarning\\(",
  `message(` = "\\bmessage\\(",
  `frm_stop(` = "\\bfrm_stop\\(",
  `frm_warning(` = "\\bfrm_warning\\(",
  `frm_message(` = "\\bfrm_message\\(",
  `stop_*/warn_* helper call` = "\\b(stop|warn)_[a-z_]+\\(",
  `structure(class = c(..., "condition"))` = "\"condition\"\\)",
  `call. = FALSE` = "call\\. = FALSE",
  `family name compared (== "<fam>" or %in% c("<fam>"...))` =
    "(family|fam|famname|fname)[\\]\\[\"a-z$]*\\s*(==|%in%)\\s*(c\\()?\"",
  `switch on a family name` = "switch\\((family|fam|famname|fname)",
  `inherits(x, "frmtmb_...")` = "inherits\\([^,]+, \"frmtmb_")
cat("\nPattern counts over code lines of core R/:\n")
for (p in names(pat)) cat(sprintf("  %s: %d lines\n", p,
                                  sum(grepl(pat[[p]], allcode, perl = TRUE))))

# lines spanned by condition calls (the message text usually spans lines)
cond_fun <- c("stop", "warning", "message", "frm_stop", "frm_warning",
              "frm_message", "warningCondition", "errorCondition",
              "simpleCondition", "signalCondition")
cl <- 0; cl_by <- list()
for (f in names(P$core)) {
  p <- parse(file.path(P$wt, "R", paste0(f, ".R")), keep.source = TRUE)
  pd <- utils::getParseData(p)
  helper <- grepl("^(stop|warn)_[a-z_]+$", pd$text)
  k <- which(pd$token == "SYMBOL_FUNCTION_CALL" &
               (pd$text %in% cond_fun | helper))
  e1 <- pd$parent[k]
  e2 <- pd$parent[match(e1, pd$id)]
  r <- match(e2, pd$id)
  ln <- unique(unlist(Map(seq, pd$line1[r], pd$line2[r])))
  cl_by[[f]] <- length(ln); cl <- cl + length(ln)
}
cb <- sort(unlist(cl_by), decreasing = TRUE)
cat(sprintf("\nLines inside a condition call (%s, stop_*/warn_* helpers): %d of %d code lines (%.1f%%)\n",
            paste(cond_fun[1:6], collapse = ", "), cl, length(allcode),
            100 * cl / length(allcode)))
cat("  by file, top 12:", paste0(names(cb)[1:12], ".R ", cb[1:12],
                                   collapse = ", "), "\n")
# the definitions of stop_*/warn_* helpers
hd <- funs[grepl("^(stop|warn)_[a-z_]+$", funs$name), ]
cat(sprintf("stop_*/warn_* helper definitions: %d, %d code lines\n",
            nrow(hd), sum(hd$code)))
cat("END codesize-dup\n")
