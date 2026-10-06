# Reviewer: brms's standata against frmtmb's frame for an ORDERED by
# factor under cmc = FALSE (contr.poly columns, nonzero on every row),
# a character by, and a logical by; Hilbert-space basis rows times Cgp.
.libPaths(c("C:/Users/adf44/source/r/wt-gpby-lib",
            "C:/Users/adf44/source/r/rellib-r5",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
wt <- "C:/Users/adf44/source/r/frmtmb-wt-gpby"
src <- parse(file.path(wt, "tests/testthat/test-gp-by.R"))
for (e in src) if (is.call(e) && identical(e[[1]], as.name("<-"))) eval(e)
d <- gpby_data()
d$fo <- factor(d$f, ordered = TRUE)
d$fch <- as.character(d$f)
d$flg <- d$f == "a"
for (fs in c("y ~ gp(x, by = fo, k = 8, cmc = FALSE)",
             "y ~ gp(x, by = fo, k = 8)", "y ~ gp(x, by = fch, k = 8)",
             "y ~ gp(x, by = flg, k = 8)")) {
  f <- stats::as.formula(fs)
  sdat <- suppressMessages(brms::standata(brms::bf(f), data = d,
                                          family = gaussian()))
  fr <- frm(bf(f), data = d, dry_run = "frame")
  gis <- fr$linpreds[["y.mu"]]$gps
  Z <- as.matrix(fr$linpreds[["y.mu"]]$Z)
  kgp <- as.integer(sdat$Kgp_1)
  cat(fs, "| Kgp", kgp, "frmtmb sub-GPs", length(gis), "\n")
  for (j in seq_len(kgp)) {
    sfx <- paste0("_1_", j)
    bk <- fr$re_blocks[[gis[[j]]$block_id]]
    ig <- as.integer(sdat[[paste0("Igp", sfx)]])
    xb <- sdat[[paste0("Xgp", sfx)]]
    jg <- sdat[[paste0("Jgp", sfx)]]
    if (!is.null(jg)) xb <- xb[jg, , drop = FALSE]
    cg <- sdat[[paste0("Cgp", sfx)]]
    if (!is.null(cg)) xb <- xb * as.numeric(cg)
    cat(sprintf("  sub %d: rows identical %s | max |Z - brms| %.2e | slambda max rel %.2e\n",
                j, identical(which(rowSums(abs(Z[, bk$c_idx, drop = FALSE])) > 0), ig),
                max(abs(Z[ig, bk$c_idx] - xb)),
                max(abs(as.vector(gis[[j]]$omega) /
                          as.vector(sdat[[paste0("slambda", sfx)]]) - 1))))
  }
}
