# False-alarm rates of the three formula refusals this lane adds,
# measured against brms 2.23.0 as the judge of what is a valid formula.
#
#   Rscript dev/priorform-falsealarm.R
#
# Two formula sets, neither hand-picked:
#  H  every formula literal in brms's own test suite, its purled
#     vignettes and its Rd examples, harvested by walking the parsed
#     AST (dev/brms-suite/brms, sha256 in dev/brms-suite-audit.md).
#  G  generated designs: every single, every unordered pair and 300
#     seeded triples from a pool of group-level terms; and every
#     two-atom arrangement of whole-term specials (cs(), smooths, gp(),
#     ar()) and plain variables over four operators and three wrappers.
#  T  twins and grouping-order designs: every pool term written twice,
#     reversed and slash-nested interactions, twin gr(cov = ) and mm(),
#     and twins in an nlpar, in sigma and in one response of an mvbf().
#
# A FALSE ALARM is a formula brms accepts that frmtmb refuses with the
# new message. A MISS is a formula brms refuses for the same reason that
# frmtmb does not. Formulas that fail in either package for any OTHER
# reason are counted separately and excluded from both rates.
#
# Writes dev/priorform-falsealarm.tsv (one row per formula and check)
# and prints a generated block for dev/priorform-findings.md.
.libPaths(c("C:/Users/adf44/source/r/priorform-lib",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
suppressMessages(requireNamespace("brms"))
cat("frmtmb from", find.package("frmtmb"), "\n")
cat("brms", as.character(packageVersion("brms")), "\n")
seed <- 20260916
set.seed(seed)

msg_of <- function(expr) {
  tryCatch({
    suppressWarnings(suppressMessages(force(expr)))
    ""
  }, error = function(e) conditionMessage(e))
}

# ---- harvest -----------------------------------------------------------
root <- "dev/brms-suite/brms"
src_files <- c(
  list.files(file.path(root, "tests", "testthat"), pattern = "[.]R$",
             full.names = TRUE),
  list.files(file.path(root, "inst", "doc"), pattern = "[.]R$",
             full.names = TRUE))
rd_files <- list.files(file.path(root, "man"), pattern = "[.]Rd$",
                       full.names = TRUE)
exprs <- list()
for (f in src_files) {
  ex <- tryCatch(parse(f, keep.source = FALSE), error = function(e) NULL)
  if (!is.null(ex)) exprs[[f]] <- ex
}
n_rd_parsed <- 0L
for (f in rd_files) {
  tf <- tempfile(fileext = ".R")
  ok <- tryCatch({
    tools::Rd2ex(f, tf, commentDontrun = FALSE, commentDonttest = FALSE)
    TRUE
  }, error = function(e) FALSE)
  if (!ok || !file.exists(tf)) next
  ex <- tryCatch(parse(tf, keep.source = FALSE), error = function(e) NULL)
  if (!is.null(ex)) {
    exprs[[f]] <- ex
    n_rd_parsed <- n_rd_parsed + 1L
  }
}
formulas <- list()
walk <- function(x) {
  if (is.call(x) && identical(x[[1L]], as.name("~"))) {
    formulas[[length(formulas) + 1L]] <<- x
    return(invisible())
  }
  if (!(is.call(x) || is.expression(x) || is.pairlist(x) || is.list(x))) {
    return(invisible())
  }
  # an empty argument (`x[, 1]`, `function(a)`) is the missing symbol,
  # which cannot even be bound to a name, so it is compared in place
  for (i in seq_along(x)) {
    if (identical(x[[i]], quote(expr = ))) next
    walk(x[[i]])
  }
}
for (ex in exprs) walk(ex)
keys <- vapply(formulas, deparse1, "")
formulas <- formulas[!duplicated(keys)]
keys <- keys[!duplicated(keys)]
formulas <- lapply(formulas, function(x) eval(x, baseenv()))
cat("source files parsed:", length(exprs), " (Rd:", n_rd_parsed, ")\n")
cat("distinct formula literals harvested:", length(formulas), "\n")

rows <- list()
record <- function(check, set, formula, brms_msg, frm_msg, brms_says,
                   frm_says) {
  rows[[length(rows) + 1L]] <<- data.frame(
    check = check, set = set, formula = formula, brms = brms_says,
    frmtmb = frm_says, brms_msg = gsub("[\t\n]", " ", brms_msg),
    frmtmb_msg = gsub("[\t\n]", " ", frm_msg), stringsAsFactors = FALSE)
}
verdict <- function(msg, pattern) {
  if (!nzchar(msg)) "accept" else if (grepl(pattern, msg)) "refuse" else
    "other"
}

# ---- check 1: '~~' -----------------------------------------------------
for (i in seq_along(formulas)) {
  f <- formulas[[i]]
  bm <- msg_of(brms::bf(f))
  fm <- msg_of(frmtmb:::refuse_nested_formula(f))
  # brms's own refusal is its as_formula() message; a formula brms::bf()
  # refuses for any other reason says nothing about this check
  record("nested", "H", keys[i], bm, fm,
         verdict(bm, "Nested formulas"), verdict(fm, "Nested formulas"))
}

# ---- check 2: a whole-term special inside a term ------------------------
# cs(), the smooths, gp(), the autocorrelation terms and mmc(): brms's
# own list, frmtmb:::brms_whole_term_specials, which the test suite
# compares with brms's regex_sp()
specials <- frmtmb:::brms_whole_term_specials
calls_special <- function(f) {
  any(vapply(specials, function(sp) frmtmb:::calls_function(f, sp), TRUE))
}
sp_msg_frm <- function(f) {
  rhs <- f[[length(f)]]
  for (tm in frmtmb:::split_plus(rhs)) {
    r <- frmtmb:::special_term_refusal(tm)
    if (!is.null(r)) return(r)
  }
  ""
}
sp_msg_brms <- function(f) {
  # cs() needs an ordinal family to reach the term check at all
  fam <- if (frmtmb:::calls_function(f, "cs")) brms::sratio() else
    stats::gaussian()
  msg_of(brms::brmsterms(brms::bf(f, family = fam)))
}
sp_pat <- "is invalid"
for (i in seq_along(formulas)) {
  f <- formulas[[i]]
  if (!calls_special(f)) next
  bm <- sp_msg_brms(f)
  fm <- sp_msg_frm(f)
  record("special", "H", keys[i], bm, fm, verdict(bm, sp_pat),
         verdict(fm, sp_pat))
}
atoms <- c("x", "z", "cs(x)", "cs(z)", "cs(x:z)", "cs(1)", "s(z)",
           "t2(x, z)", "s(z, by = x)", "gp(z)", "gp(x, by = g)",
           "ar(time)", "te(x, z)")
ops <- c("+", ":", "*", "/")
gen_sp <- character(0)
for (a in atoms) for (o in ops) for (b in atoms) {
  gen_sp <- c(gen_sp, paste("y ~", a, o, b))
}
for (a in atoms) for (w in c("I", "log", "exp")) {
  gen_sp <- c(gen_sp, paste0("y ~ x + ", w, "(", a, ")"))
}
gen_sp <- unique(gen_sp)
for (s in gen_sp) {
  f <- stats::as.formula(s, env = globalenv())
  bm <- sp_msg_brms(f)
  fm <- sp_msg_frm(f)
  record("special", "G", s, bm, fm, verdict(bm, sp_pat),
         verdict(fm, sp_pat))
}

# ---- check 3: duplicated group-level effects ---------------------------
n <- 120
dat <- data.frame(y = stats::rnorm(n), x = stats::rnorm(n),
                  z = stats::rnorm(n),
                  f = factor(sample(c("a", "b", "c"), n, TRUE)),
                  g = factor(rep(1:6, 20)), h = factor(rep(1:5, 24)))
dup_pat <- "Duplicated group-level effects"
dup_one <- function(s, set, dpar_form = NULL) {
  f <- stats::as.formula(s, env = globalenv())
  bb <- if (is.null(dpar_form)) brms::bf(f) else
    brms::bf(f, stats::as.formula(dpar_form, env = globalenv()))
  fb <- if (is.null(dpar_form)) frmtmb::bf(f) else
    frmtmb::bf(f, stats::as.formula(dpar_form, env = globalenv()))
  bm <- msg_of(brms::default_prior(bb, data = dat))
  fm <- msg_of(frmtmb::frm(fb, data = dat, dry_run = "frame"))
  lab <- if (is.null(dpar_form)) s else paste(s, "|", dpar_form)
  record("dup", set, lab, bm, fm, verdict(bm, dup_pat), verdict(fm, dup_pat))
}
pool <- c("(1 | g)", "(x | g)", "(0 + x | g)", "(1 + x | g)", "(x || g)",
          "(z | g)", "(0 + x + z | g)", "(0 + f | g)", "(f | g)",
          "(1 | h)", "(x | h)", "(1 | g:h)", "(1 | g/h)", "(1 | p | g)",
          "(x | p | g)", "(1 | mm(g, h))", "(x | mm(g, h))",
          "(0 + z || g)")
designs <- as.list(pool)
pairs <- utils::combn(pool, 2L, simplify = FALSE)
designs <- c(designs, pairs)
triples <- utils::combn(pool, 3L, simplify = FALSE)
designs <- c(designs, triples[sample.int(length(triples), 300L)])
for (d in designs) {
  dup_one(paste("y ~ x +", paste(d, collapse = " + ")), "G")
}
# the same factor in two DIFFERENT predictors is two parameters in both
# packages, so none of these may be refused
for (d in pool) {
  dup_one(paste("y ~ x +", d), "G", paste("sigma ~", d))
}
# exact twins, which brms's terms() collapses to one term, and grouping
# factors written in either order or through a slash. The generated
# pairs above combine DISTINCT pool terms, so they contain none of these
for (d in pool) {
  dup_one(paste("y ~ x +", d, "+", d), "T")
}
for (s in c("y ~ (1 | g:h) + (1 | h:g)", "y ~ (1 | g:h) + (1 | g/h)",
            "y ~ (1 | h:g) + (1 | g/h)", "y ~ (x | g/h) + (x | g/h)",
            "y ~ (1 | g/h) + (1 | g)", "y ~ (x || g) + (1 | g)")) {
  dup_one(s, "T")
}
for (d in c("(1 | g)", "(x | g)")) {
  dup_one("y ~ x", "T", paste("sigma ~", d, "+", d))
}
A <- diag(6)
dimnames(A) <- list(levels(dat$g), levels(dat$g))
for (s in c("y ~ (1 | gr(g, cov = A)) + (1 | gr(g, cov = A))",
            "y ~ (1 | gr(g, cov = A)) + (1 | g)")) {
  f <- stats::as.formula(s, env = globalenv())
  bm <- msg_of(brms::default_prior(brms::bf(f), data = dat,
                                   data2 = list(A = A)))
  fm <- msg_of(frmtmb::frm(frmtmb::bf(f), data = dat, data2 = list(A = A),
                           dry_run = "frame"))
  record("dup", "T", s, bm, fm, verdict(bm, dup_pat), verdict(fm, dup_pat))
}


bm <- msg_of(brms::default_prior(
  brms::bf(y ~ a + x, a ~ 1 + (1 | g) + (1 | g), nl = TRUE), data = dat))
fm <- msg_of(frmtmb::frm(
  frmtmb::bf(y ~ a + x, a ~ 1 + (1 | g) + (1 | g), nl = TRUE), data = dat,
  dry_run = "frame"))
record("dup", "T", "bf(y ~ a + x, a ~ 1 + (1 | g) + (1 | g), nl = TRUE)",
       bm, fm, verdict(bm, dup_pat), verdict(fm, dup_pat))
dat$y2 <- stats::rnorm(n)
bm <- msg_of(brms::default_prior(
  brms::bf(y ~ (1 | g) + (1 | g)) + brms::bf(y2 ~ (1 | g)) +
    brms::set_rescor(FALSE), data = dat))
fm <- msg_of(frmtmb::frm(
  frmtmb::mvbf(frmtmb::bf(y ~ (1 | g) + (1 | g)), frmtmb::bf(y2 ~ (1 | g))),
  data = dat, dry_run = "frame"))
record("dup", "T", "mvbf(bf(y ~ (1 | g) + (1 | g)), bf(y2 ~ (1 | g)))",
       bm, fm, verdict(bm, dup_pat), verdict(fm, dup_pat))
# harvested formulas with a bar, on data synthesized for their variables
syn_data <- function(f) {
  v <- all.vars(f)
  bars <- unlist(lapply(frmtmb:::split_plus(f[[length(f)]]), function(tm) {
    if (is.call(tm) && identical(tm[[1L]], as.name("("))) tm <- tm[[2L]]
    if (is.call(tm) && as.character(tm[[1L]])[1] %in% c("|", "||")) {
      all.vars(tm[[length(tm)]])
    }
  }))
  out <- data.frame(row.names = seq_len(n))
  for (nm in v) {
    out[[nm]] <- if (nm %in% bars) factor(rep(seq_len(6), length.out = n))
      else stats::rnorm(n)
  }
  out
}
n_h_bar <- 0L
for (i in seq_along(formulas)) {
  f <- formulas[[i]]
  if (length(f) != 3L || !("|" %in% all.names(f[[3L]]))) next
  n_h_bar <- n_h_bar + 1L
  dd <- syn_data(f)
  bm <- msg_of(brms::default_prior(brms::bf(f), data = dd))
  fm <- msg_of(frmtmb::frm(frmtmb::bf(f), data = dd, dry_run = "frame"))
  record("dup", "H", keys[i], bm, fm, verdict(bm, dup_pat),
         verdict(fm, dup_pat))
}

tab <- do.call(rbind, rows)
utils::write.table(tab, "dev/priorform-falsealarm.tsv", sep = "\t",
                   row.names = FALSE, quote = FALSE)

cat("\n<!-- priorform-falsealarm:begin -->\n")
cat(sprintf("seed %d; brms %s; %d source files (%d Rd); %d distinct",
            seed, as.character(packageVersion("brms")), length(exprs),
            n_rd_parsed, length(formulas)),
    "harvested formula literals\n\n")
cat("| check | set | formulas | brms accepts | false alarms |",
    "brms refuses (same reason) | misses | other error, excluded |\n")
cat("|---|---|---|---|---|---|---|---|\n")
for (ck in c("nested", "special", "dup")) for (st in c("H", "G", "T")) {
  t <- tab[tab$check == ck & tab$set == st, ]
  if (!nrow(t)) next
  acc <- t$brms == "accept"
  ref <- t$brms == "refuse"
  oth <- t$brms == "other" | (t$frmtmb == "other" & !ref)
  cat(sprintf("| %s | %s | %d | %d | %d | %d | %d | %d |\n", ck, st,
              nrow(t), sum(acc & !oth), sum(acc & !oth & t$frmtmb == "refuse"),
              sum(ref), sum(ref & t$frmtmb != "refuse"), sum(oth)))
}
cat("<!-- priorform-falsealarm:end -->\n")
fa <- tab[tab$brms == "accept" & tab$frmtmb == "refuse", ]
if (nrow(fa)) {
  cat("\nFALSE ALARMS:\n")
  print(fa[, c("check", "set", "formula", "frmtmb_msg")], row.names = FALSE)
}
miss <- tab[tab$brms == "refuse" & tab$frmtmb != "refuse", ]
if (nrow(miss)) {
  cat("\nMISSES:\n")
  print(miss[, c("check", "set", "formula", "frmtmb", "frmtmb_msg")],
        row.names = FALSE)
}
