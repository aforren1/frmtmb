# Item 7: was the last change text only? The R CMD check source copy against the worktree, and
# the two installed objects. S2: frm_linpred(type = "response") on an ordinal body.
source("dev/review-nlfix-common.R")
ns <- asNamespace("frmtmb")
src <- "dev/nlfix-check/frmtmb/frmtmb.Rcheck/00_pkg_src/frmtmb"
cat("check copy built:", format(file.info(file.path(src, "DESCRIPTION"))$mtime), "\n")
fs <- list.files("R", "[.]R$")
for (f in fs) {
  a <- readLines(file.path("R", f), warn = FALSE); b <- if (file.exists(file.path(src, "R", f))) readLines(file.path(src, "R", f), warn = FALSE) else character(0)
  if (!identical(a, b)) cat("differs:", f, " lines changed:", sum(a != b[seq_along(a)], na.rm = TRUE), " sizes", length(a), length(b), "\n")
}
e <- new.env(parent = ns)
for (f in list.files("R", "[.]R$", full.names = TRUE)) sys.source(f, envir = e, keep.source = FALSE)
nm <- ls(e, all.names = TRUE); isf <- vapply(nm, function(x) is.function(e[[x]]), TRUE)
dif <- character(0)
for (x in nm[isf]) { g <- get0(x, envir = ns, inherits = FALSE); if (is.null(g) || !identical(paste(deparse(body(e[[x]])), collapse = ""), paste(deparse(body(g)), collapse = ""))) dif <- c(dif, x) }
cat("worktree against wt-nlfix-lib: functions", sum(isf), " differing", length(dif), if (length(dif)) paste(dif, collapse = ","), "\n")
e2 <- new.env(parent = ns)
for (f in list.files(file.path(src, "R"), "[.]R$", full.names = TRUE)) sys.source(f, envir = e2, keep.source = FALSE)
dif2 <- character(0); nm2 <- ls(e2, all.names = TRUE)
for (x in union(nm[isf], nm2)) { a <- e[[x]]; b <- e2[[x]]; if (is.function(a) != is.function(b) || (is.function(a) && !identical(paste(deparse(body(a)), collapse = ""), paste(deparse(body(b)), collapse = "")))) dif2 <- c(dif2, x) }
cat("worktree against the check copy, by function body:", length(dif2), if (length(dif2)) paste(dif2, collapse = ","), "\n")
if ("compat_hand_rules_tbl" %in% dif2) {
  t1 <- e$compat_hand_rules_tbl(); t2 <- e2$compat_hand_rules_tbl()
  t1 <- as.data.frame(t1); t2 <- as.data.frame(t2)
  cat("rows", nrow(t1), nrow(t2), " columns", names(t1), "\n")
  key <- function(d) paste(d[[1]], d[[2]], d[[3]], sep = "|")
  if (nrow(t1) == nrow(t2)) {
    cat("feature and status columns identical:", identical(t1[setdiff(names(t1), "note")], t2[setdiff(names(t2), "note")]), "\n")
    w <- which(t1$note != t2$note); cat("notes differing:", length(w), "\n"); for (i in w) cat(" ", unlist(t1[i, 1:3]), "\n  now:", substr(t1$note[i], 1, 220), "\n  was:", substr(t2$note[i], 1, 220), "\n")
  }
}
cat("\n--- logs of the suite, lib lines\n")
for (f in c("dev/nlfix-suite-core.log", "dev/nlfix-suite-ext.log", "dev/nlfix-suite-gated.log")) if (file.exists(f)) cat(f, format(file.info(f)$mtime), " | ", tail(readLines(f, warn = FALSE), 1), "\n")
fl <- list.files("dev/nlfix-suite-core", full.names = TRUE)
libs <- unique(unlist(lapply(head(fl, 400), function(p) grep("^lib", readLines(p, warn = FALSE, n = 3), value = TRUE))))
cat("core logs:", length(fl), " newest", format(max(file.info(fl)$mtime)), " lib lines:", head(libs, 4), "\n")
cat("\n--- S2: frm_linpred on an ordinal body, type = response\n")
set.seed(9202); n <- 300
d <- data.frame(x = round(rnorm(n), 3)); d$y <- ordered(cut(0.8 * d$x + rlogis(n), c(-Inf, -1, 0, 1, Inf), labels = 1:4))
fn <- frm(bf(y ~ b1 * x, b1 ~ 1, nl = TRUE), family = cumulative(), data = d); fl <- frm(y ~ 0 + x, family = cumulative(), data = d)
nd <- d[1:3, ]
cat("body link:", format(frm_linpred(fn, newdata = nd, type = "link"), digits = 5), "\n")
cat("body response:", format(frm_linpred(fn, newdata = nd, type = "response"), digits = 5), "\n")
cat("linear response:", format(frm_linpred(fl, newdata = nd, type = "response"), digits = 5), "\n")
x <- tt(fitted(fn, newdata = nd, dpar = "mu")); cat("fitted(dpar = mu) body:", if (iserr(x)) x$msg else format(x[, 1], digits = 5), "\n")
x <- tt(fitted(fl, newdata = nd, dpar = "mu")); cat("fitted(dpar = mu) linear:", if (iserr(x)) x$msg else format(x[, 1], digits = 5), "\n")
cat("brms posterior_epred(dpar = mu) would apply the family's link: plogis(body link)", format(plogis(frm_linpred(fn, newdata = nd, type = "link")), digits = 5), "\n")
