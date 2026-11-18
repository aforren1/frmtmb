# Q2 (static half) and Q4 of the codesize lane: a static call graph of
# core's top-level functions; functions no root reaches; single-use
# helpers; guard and check code against modeling code.
#   Rscript dev/codesize-deadcode.R
# Reads dev/codesize-out/parsed.rds and functions.tsv (codesize-parse.R,
# codesize-inventory.R). Writes dev/codesize-out/callgraph.tsv and prints
# the blocks pasted into dev/codesize-findings.md.
#
# An edge f -> g is any SYMBOL or SYMBOL_FUNCTION_CALL token `g` inside
# f's definition, so do.call(g, ...), match.fun(g), lapply(x, g) and a
# registry entry list(k = g) all count. Tokens after `$`, `@` or another
# package's `::` do not. A string equal to a function's name is a second
# kind of edge ("string"), reported separately because it may be a
# message or a column name and not a lookup.
#
# Roots: exports and registered S3 methods; methods g.cls of a generic g
# that the package or base/stats/utils/graphics define with UseMethod
# (dispatch); hooks .onLoad/.onAttach/.onUnload; every top-level
# expression that is not a function definition (registries filled at
# build time). External roots: any symbol or string in tests/, inst/ and
# extensions/ (their R/ and tests/).
out <- "C:/Users/adf44/source/r/frmtmb-wt-codesize/dev/codesize-out"
P <- readRDS(file.path(out, "parsed.rds"))
funs <- read.delim(file.path(out, "functions.tsv"), stringsAsFactors = FALSE)
# a function stored into a registry slot (x[["k"]]$f <- function) has no
# name; its body is a root, like any other build-time expression
funs <- funs[!is.na(funs$name), ]
fn <- unique(funs$name)
dupdef <- funs$name[duplicated(funs$name)]

ref_tokens <- function(pd) {
  if (is.null(pd)) return(pd)
  pd <- pd[pd$token != "COMMENT", ]
  prev <- c("", pd$token[-nrow(pd)])
  prev2 <- c("", "", pd$text[seq_len(max(0, nrow(pd) - 2))])
  ok <- pd$token %in% c("SYMBOL", "SYMBOL_FUNCTION_CALL", "SPECIAL") &
    !prev %in% c("'$'", "'@'") &
    !(prev %in% c("NS_GET", "NS_GET_INT") & prev2 != "frmtmb")
  sym <- pd[ok, ]
  st <- pd[pd$token == "STR_CONST", ]
  st$text <- substr(st$text, 2, nchar(st$text) - 1)
  list(sym = sym, str = st)
}

# per-function symbol and string references, and root expressions
edges <- list(); sedges <- list(); root_sym <- character(); root_str <- character()
for (f in names(P$core)) {
  x <- P$core[[f]]
  r <- ref_tokens(x$pd)
  for (i in seq_len(nrow(x$top))) {
    a <- x$top$first[i]; b <- x$top$last[i]
    s <- unique(r$sym$text[r$sym$line1 >= a & r$sym$line1 <= b])
    t <- unique(r$str$text[r$str$line1 >= a & r$str$line1 <= b])
    nm <- x$top$name[i]
    if (isTRUE(x$top$is_fun[i]) && !is.na(nm)) {
      key <- paste0(f, ".R::", nm)
      edges[[key]] <- setdiff(intersect(s, fn), nm)
      sedges[[key]] <- setdiff(intersect(t, fn), nm)
    } else {
      # a non-function top-level expression runs at build time
      root_sym <- c(root_sym, intersect(s, fn))
      root_str <- c(root_str, intersect(t, fn))
    }
  }
}
key_name <- sub("^.*::", "", names(edges))
# a name defined twice: an edge reaches every definition
expand <- function(nms) names(edges)[key_name %in% nms]

# S3: generics the package defines with UseMethod, plus base ones
pkg_generics <- unique(unlist(lapply(names(P$core), function(f) {
  x <- P$core[[f]]
  um <- x$pd$line1[x$pd$text == "UseMethod"]
  t <- x$top[x$top$is_fun, ]
  t$name[vapply(seq_len(nrow(t)), function(i)
    any(um >= t$first[i] & um <= t$last[i]), TRUE)]
})))
base_generics <- unique(unlist(lapply(
  c("base", "stats", "utils", "graphics", "methods"), function(p) {
    ns <- asNamespace(p)
    nm <- ls(ns, all.names = TRUE)
    nm[vapply(nm, function(n) {
      g <- get(n, ns)
      is.function(g) && !is.primitive(g) &&
        any(grepl("UseMethod", deparse(body(g)), fixed = TRUE))
    }, TRUE)]
  })))
base_generics <- c(base_generics, "print", "format", "summary", "plot",
                   "as.data.frame", "as.list", "as.matrix", "length",
                   "names", "[", "[[", "$", "c", "rep", "t", "dim", "+", "-",
                   "*", "/", "==", "!=", "$<-", "[<-", "[[<-", "Ops", "Math",
                   "Summary", "names<-", "dimnames", "unique", "all.equal")
gen_all <- unique(c(pkg_generics, base_generics))
s3_dispatch <- fn[vapply(fn, function(n) {
  parts <- strsplit(n, ".", fixed = TRUE)[[1]]
  if (length(parts) < 2) return(FALSE)
  any(vapply(seq_len(length(parts) - 1), function(k)
    paste(parts[1:k], collapse = ".") %in% gen_all, TRUE))
}, TRUE)]

hooks <- intersect(fn, c(".onLoad", ".onAttach", ".onUnload", ".onDetach"))
roots_pkg <- unique(c(intersect(fn, P$exports), intersect(fn, P$s3$method),
                      hooks, root_sym))
roots_disp <- setdiff(s3_dispatch, roots_pkg)

# external references
ext_sym <- character(); ext_str <- character(); ext_by <- list()
for (nm in names(P$ext)) {
  x <- P$ext[[nm]]
  if (is.null(x$pd)) {
    w <- unique(regmatches(x$lines, gregexpr("[A-Za-z.][A-Za-z0-9._]*",
                                             x$lines))[[1]])
    s <- intersect(unlist(w), fn); t <- character()
  } else {
    r <- ref_tokens(x$pd)
    s <- intersect(unique(r$sym$text), fn)
    t <- intersect(unique(r$str$text), fn)
  }
  grp <- if (grepl("^extensions/", nm)) "extensions" else
    if (grepl("^tests/", nm)) "tests" else "inst"
  ext_by[[grp]] <- unique(c(ext_by[[grp]], s, t))
  ext_sym <- c(ext_sym, s); ext_str <- c(ext_str, t)
}
ext_sym <- unique(ext_sym); ext_str <- unique(ext_str)

reach <- function(start, use_str = FALSE) {
  seen <- expand(start); frontier <- seen
  while (length(frontier)) {
    nxt <- unique(unlist(c(edges[frontier],
                           if (use_str) sedges[frontier])))
    nk <- setdiff(expand(nxt), seen)
    seen <- c(seen, nk); frontier <- nk
  }
  seen
}
r_pkg <- reach(roots_pkg)
r_disp <- reach(c(roots_pkg, roots_disp))
r_str <- reach(c(roots_pkg, roots_disp, root_str), use_str = TRUE)
r_ext <- reach(c(roots_pkg, roots_disp, root_str, ext_sym, ext_str),
               use_str = TRUE)

cg <- data.frame(key = names(edges), name = key_name,
                 file = sub("::.*$", "", names(edges)))
cg$status <- ifelse(cg$key %in% r_pkg, "live",
             ifelse(cg$key %in% r_disp, "live via S3 dispatch",
             ifelse(cg$key %in% r_str, "uncertain: string only",
             ifelse(cg$key %in% r_ext, "only tests/extensions/inst",
                    "dead"))))
m <- match(paste(cg$file, cg$name), paste(funs$file, funs$name))
cg$code <- funs$code[m]; cg$comment <- funs$comment[m]; cg$rox <- funs$rox[m]
cg$first <- funs$first[m]; cg$last <- funs$last[m]
cg$subsystem <- funs$subsystem[m]
cg$exported <- funs$exported[m]
# who references the not-live ones from outside
cg$ext_ref <- vapply(cg$name, function(n) paste(names(ext_by)[
  vapply(ext_by, function(v) n %in% v, TRUE)], collapse = ","), "")

# call sites in package code: token count of the name in other functions
site <- table(unlist(lapply(names(P$core), function(f) {
  x <- P$core[[f]]; r <- ref_tokens(x$pd)
  tf <- x$top[x$top$is_fun, ]
  # drop the defining token of each function
  own <- paste(tf$name, tf$first)
  s <- r$sym[!paste(r$sym$text, r$sym$line1) %in% own, ]
  s$text[s$text %in% fn]
})))
callers <- table(unlist(lapply(names(edges), function(k) edges[[k]])))
cg$sites <- as.vector(site[cg$name]); cg$sites[is.na(cg$sites)] <- 0
cg$callers <- as.vector(callers[cg$name]); cg$callers[is.na(cg$callers)] <- 0
cg$str_refs <- cg$name %in% unlist(sedges) | cg$name %in% root_str
write.table(cg, file.path(out, "callgraph.tsv"), sep = "\t", quote = FALSE,
            row.names = FALSE)

cat("BEGIN codesize-deadcode\n")
cat(sprintf("functions %d (names defined twice: %s); roots: %d exported or S3-registered, %d hooks, %d named in top-level non-function code, %d S3 methods by dispatch\n",
    nrow(cg), paste(unique(dupdef), collapse = " "),
    length(intersect(fn, c(P$exports, P$s3$method))), length(hooks),
    length(unique(root_sym)), length(roots_disp)))
st <- aggregate(cbind(n = 1, code = cg$code, comment = cg$comment,
                      rox = cg$rox) ~ status, cg, sum)
cat("| status | functions | code lines | other comments | roxygen |\n|---|---|---|---|---|\n")
for (i in seq_len(nrow(st))) with(st[i, ], cat(sprintf(
  "| %s | %d | %d | %d | %d |\n", status, n, code, comment, rox)))
for (s in c("dead", "uncertain: string only", "only tests/extensions/inst",
            "live via S3 dispatch")) {
  d <- cg[cg$status == s, ]
  d <- d[order(-d$code), ]
  cat(sprintf("\n%s (%d functions, %d code lines):\n", s, nrow(d), sum(d$code)))
  cat(sprintf("  %s %s:%d (%d code lines)%s\n", d$name, d$file, d$first,
              d$code, ifelse(nzchar(d$ext_ref), paste0(" refs: ", d$ext_ref),
                             "")), sep = "")
}
cat("\nNot live, per file:\n")
nl <- cg[cg$status != "live" & cg$status != "live via S3 dispatch", ]
pf <- aggregate(cbind(n = 1, code = nl$code) ~ file + status, nl, sum)
pf <- pf[order(pf$file), ]
cat(sprintf("  %s %s: %d functions, %d code lines\n", pf$file, pf$status,
            pf$n, pf$code), sep = "")

# Q4: single-use helpers
cat("\nSingle-use helpers (internal, one call site in package code):\n")
su <- cg[!cg$exported & cg$sites == 1 & cg$status == "live", ]
for (n in c(5, 10, 20)) {
  x <- su[su$code < n, ]
  cat(sprintf("  under %d code lines: %d functions, %d code lines, %d roxygen lines, %d comment lines\n",
              n, nrow(x), sum(x$code), sum(x$rox), sum(x$comment)))
}
cat(sprintf("  any size: %d functions, %d code lines\n", nrow(su), sum(su$code)))
lu <- cg[!cg$exported & cg$status == "live", ]
cat(sprintf("Internal live functions: %d; with 1 site %d, 2 sites %d, 3+ sites %d\n",
            nrow(lu), sum(lu$sites == 1), sum(lu$sites == 2),
            sum(lu$sites >= 3)))
fs <- aggregate(cbind(n = 1, code = su$code[su$code < 10]) ~ file,
                su[su$code < 10, ], sum)
fs <- fs[order(-fs$n), ][1:15, ]
cat("Files with the most single-use helpers under 10 code lines:\n")
cat(sprintf("  %s: %d (%d code lines)\n", fs$file, fs$n, fs$code), sep = "")
# median roxygen per small single-use helper: the doc cost of a split
x <- su[su$code < 10, ]
cat(sprintf("Small single-use helpers: median code %g, median roxygen %g, roxygen/code %.2f\n",
            median(x$code), median(x$rox), sum(x$rox) / sum(x$code)))
# A symbol edge also matches a local variable of the same name. So list
# the internal live functions that no file ever CALLS by name (no
# SYMBOL_FUNCTION_CALL token anywhere): each is either passed as a value
# (registry, do.call, lapply) or reached only through a namesake variable.
calls_anywhere <- unique(unlist(lapply(c(P$core, P$ext), function(x) {
  if (is.null(x$pd)) return(NULL)
  x$pd$text[x$pd$token %in% c("SYMBOL_FUNCTION_CALL", "SPECIAL")]
})))
# a value use: the bare name as an argument or a list entry
nv <- cg[cg$status == "live" & !cg$exported & !cg$name %in% calls_anywhere, ]
nv <- nv[order(-nv$code), ]
cat(sprintf("\nInternal live functions never called by name anywhere (value use or namesake variable): %d, %d code lines\n",
            nrow(nv), sum(nv$code)))
cat(sprintf("  %s %s:%d (%d code lines, %d sites)\n", nv$name, nv$file,
            nv$first, nv$code, nv$sites), sep = "")
write.table(nv, file.path(out, "callgraph-novalue.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)
cat("END codesize-deadcode\n")
