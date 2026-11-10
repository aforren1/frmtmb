# Each core addition term against each feature only an ordinal family
# has (thres(), disc, the equidistant and sum_to_zero thresholds, cs()),
# spelled on an ordinal family, run through frm() and through brms
# 2.23.0's make_stancode(). weights() is the control: the ordinal
# families take it, so its row must fit on both sides, which shows the
# structure's spelling is not what refuses.
#   Rscript dev/compat-ordinal-measure.R [lib;lib...] [out.tsv]
# Seed 20261109.
args <- commandArgs(trailingOnly = TRUE)
lib <- if (length(args) >= 1) args[1] else
  "C:/Users/adf44/source/r/rellib-r19"
out <- if (length(args) >= 2) args[2] else
  "dev/compat-ordinal-out/measure.tsv"
.libPaths(c(strsplit(lib, ";", fixed = TRUE)[[1]],
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({
  library(frmtmb)
  requireNamespace("brms")
})
cat("lib:", find.package("frmtmb"), as.character(packageVersion("frmtmb")),
    "\n")
cat("brms", as.character(packageVersion("brms")), "\n")
dir.create(dirname(out), showWarnings = FALSE, recursive = TRUE)

set.seed(20261109)
n <- 80
d <- data.frame(x = rnorm(n), z = rnorm(n))
d$y <- as.integer(cut(0.8 * d$x + rlogis(n), c(-Inf, -1, 0, 1, Inf)))
d$nt <- 4L
d$cc <- 0L
d$s <- 1
d$k <- 1L
d$rr <- 0.5
d$tt <- 2
d$w <- 1
dm <- d
dm$y[c(3, 9)] <- NA

aterms <- list(
  "weights()" = "weights(w)", "trials()" = "trials(nt)",
  "cens()" = "cens(cc)", "trunc()" = "trunc(lb = 0)", "se()" = "se(s)",
  "mi()" = "mi()", "vint()" = "vint(k)", "vreal()" = "vreal(rr)",
  "rate()" = "rate(tt)")
# each structure: the formula text (AT is replaced by the term) and
# the family text, read in the package that evaluates it
structs <- list(
  "thres()" = list(f = "bf(y | thres(3) + AT ~ x)", fam = "cumulative()"),
  disc = list(f = "bf(y | AT ~ x, disc ~ 0 + z)", fam = "cumulative()"),
  equidistant = list(f = "bf(y | AT ~ x)",
                     fam = "cumulative(threshold = 'equidistant')"),
  sum_to_zero = list(f = "bf(y | AT ~ x)",
                     fam = "cumulative(threshold = 'sum_to_zero')"),
  "cs_pred()" = list(f = "bf(y | AT ~ x + cs(z))", fam = "sratio()"))

ev <- function(txt, pkg) eval(parse(text = txt)[[1]], envir = asNamespace(pkg))
first_err <- function(expr) {
  w <- character()
  v <- withCallingHandlers(
    tryCatch({ force(expr); NA_character_ },
             error = function(e) conditionMessage(e)),
    warning = function(cw) {
      w <<- c(w, conditionMessage(cw))
      invokeRestart("muffleWarning")
    })
  gsub("[\r\n\t]+", " ", v)
}

rows <- list()
for (sn in names(structs)) for (an in names(aterms)) {
  st <- structs[[sn]]
  ftxt <- sub("AT", aterms[[an]], st$f, fixed = TRUE)
  dat <- if (an == "mi()") dm else d
  e_frm <- first_err(suppressMessages(
    frm(ev(ftxt, "frmtmb"), family = ev(st$fam, "frmtmb"), data = dat)))
  e_brm <- first_err(suppressMessages(
    brms::make_stancode(ev(ftxt, "brms"), family = ev(st$fam, "brms"),
                        data = dat)))
  base <- sub("[(][)]$", "", an)
  rows[[length(rows) + 1L]] <- data.frame(
    aterm = an, structure = sn, formula = ftxt, family = st$fam,
    frm = if (is.na(e_frm)) "fits" else "refuses",
    frm_names_term = !is.na(e_frm) && grepl(base, e_frm, fixed = TRUE),
    brms = if (is.na(e_brm)) "builds" else "refuses",
    frm_error = substr(ifelse(is.na(e_frm), "", e_frm), 1, 300),
    brms_error = substr(ifelse(is.na(e_brm), "", e_brm), 1, 200), stringsAsFactors = FALSE)
  r <- rows[[length(rows)]]
  cat(sprintf("%-10s x %-12s frm %-8s names %-5s brms %-8s | %s\n", an, sn,
              r$frm, r$frm_names_term, r$brms, substr(r$frm_error, 1, 110)))
}
res <- do.call(rbind, rows)
utils::write.table(res, out, sep = "\t", row.names = FALSE, quote = TRUE,
                   qmethod = "double")
ctl <- res$aterm == "weights()"
cat("\ncontrol weights(): frm fits", sum(res$frm[ctl] == "fits"), "of",
    sum(ctl), "; brms builds", sum(res$brms[ctl] == "builds"), "of",
    sum(ctl), "\n")
cat("other terms: frm refuses", sum(res$frm[!ctl] == "refuses"), "of",
    sum(!ctl), "; names the term", sum(res$frm_names_term[!ctl]), "of",
    sum(!ctl), "; brms refuses", sum(res$brms[!ctl] == "refuses"), "of",
    sum(!ctl), "\n")
