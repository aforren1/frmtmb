# The extensions' addition terms against each feature only an ordinal
# family has, as dev/compat-ordinal-measure.R does for core's: frm() with
# every extension loaded, and brms 2.23.0's make_stancode() for the two
# terms brms has (bhaz() of its cox family, dec() of its wiener family).
# ndt_group(), reward(), payoff() and stage2() have no brms spelling.
#   Rscript dev/compat-ordinal-ext-measure.R <lib;lib> <out.tsv>
# Seed 20261110.
args <- commandArgs(trailingOnly = TRUE)
.libPaths(c(strsplit(args[1], ";", fixed = TRUE)[[1]],
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({
  library(frmtmb)
  for (p in c("frmtmb.survival", "frmtmb.eam", "frmtmb.learn")) {
    library(p, character.only = TRUE)
  }
  requireNamespace("brms")
})
cat("lib:", find.package("frmtmb"), as.character(packageVersion("frmtmb")),
    "\n")
dir.create(dirname(args[2]), showWarnings = FALSE, recursive = TRUE)
set.seed(20261110)
n <- 80
d <- data.frame(x = rnorm(n), z = rnorm(n), w = 1, two = 1L, g = "a",
                r1 = 1, r2 = 0, p1 = 1, p2 = 0, p3 = 0, p4 = 1, s1 = 1L,
                s2 = 1L)
d$y <- as.integer(cut(0.8 * d$x + rlogis(n), c(-Inf, -1, 0, 1, Inf)))
aterms <- c("weights()" = "weights(w)", "bhaz()" = "bhaz()",
            "dec()" = "dec(two)", "ndt_group()" = "ndt_group(g)",
            "reward()" = "reward(r1, r2)",
            "payoff()" = "payoff(p1, p2, p3, p4)",
            "stage2()" = "stage2(s1, s2)")
in_brms <- c("weights()", "bhaz()", "dec()")
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
  v <- withCallingHandlers(
    tryCatch({ force(expr); NA_character_ },
             error = function(e) conditionMessage(e)),
    warning = function(cw) invokeRestart("muffleWarning"))
  gsub("[\r\n\t]+", " ", v)
}
rows <- list()
for (sn in names(structs)) for (an in names(aterms)) {
  st <- structs[[sn]]
  ftxt <- sub("AT", aterms[[an]], st$f, fixed = TRUE)
  e_frm <- first_err(suppressMessages(
    frm(ev(ftxt, "frmtmb"), family = ev(st$fam, "frmtmb"), data = d)))
  e_brm <- if (an %in% in_brms) first_err(suppressMessages(
    brms::make_stancode(ev(ftxt, "brms"), family = ev(st$fam, "brms"),
                        data = d))) else "no brms spelling"
  base <- sub("[(][)]$", "", an)
  r <- data.frame(
    aterm = an, structure = sn,
    frm = if (is.na(e_frm)) "fits" else "refuses",
    frm_names_term = !is.na(e_frm) && grepl(base, e_frm, fixed = TRUE),
    brms = if (identical(e_brm, "no brms spelling")) "none" else
      if (is.na(e_brm)) "builds" else "refuses",
    frm_error = substr(ifelse(is.na(e_frm), "", e_frm), 1, 200),
    brms_error = substr(ifelse(is.na(e_brm), "", e_brm), 1, 120),
    stringsAsFactors = FALSE)
  rows[[length(rows) + 1L]] <- r
  cat(sprintf("%-12s x %-12s frm %-8s names %-5s brms %-8s | %s | %s\n", an,
              sn, r$frm, r$frm_names_term, r$brms, substr(r$frm_error, 1, 70),
              substr(r$brms_error, 1, 60)))
}
res <- do.call(rbind, rows)
utils::write.table(res, args[2], sep = "\t", row.names = FALSE, quote = TRUE,
                   qmethod = "double")
ctl <- res$aterm == "weights()"
cat("\ncontrol weights(): frm fits", sum(res$frm[ctl] == "fits"), "of",
    sum(ctl), "; brms builds", sum(res$brms[ctl] == "builds"), "\n")
cat("extension terms: frm refuses", sum(res$frm[!ctl] == "refuses"), "of",
    sum(!ctl), "; names the term", sum(res$frm_names_term[!ctl]), "of",
    sum(!ctl), "; brms refuses", sum(res$brms[!ctl] == "refuses"), "of",
    sum(res$brms[!ctl] != "none"), "with a brms spelling;",
    sum(res$brms[!ctl] == "none"), "have none\n")
