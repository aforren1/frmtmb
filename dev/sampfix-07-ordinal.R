# Lane sampfix, script 07: brms's names on ordinal and cs() draws.
#
#   Rscript dev/sampfix-07-ordinal.R ref    # writes 07-ref.rds
#   Rscript dev/sampfix-07-ordinal.R lane   # writes 07-lane.rds, compares
#
# Four fits (data seed 405, sampler seed 3, one chain of 150 draws):
# cumulative yo ~ x (the log-increment map), sratio yo ~ x + cs(fc)
# (identity map and cs() columns), cumulative with thres(gr = ) (the
# per-group map) and a two-response cumulative model (response-prefixed
# names). Per fit: variables(ds) against variables(fit); on the lane, the
# stored draws against the reference build's (same sampler, same seed)
# after the map; posterior_epred() and fixef() against the reference;
# and a draws object stored the old way (internal names, internal
# scale) read by the lane build, which must give the same answers.

arm <- commandArgs(trailingOnly = TRUE)[1L]
LANE <- "C:/Users/adf44/source/r/wt-sampfix-lib"
REF <- "C:/Users/adf44/source/r/rellib-r3"
USER <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
.libPaths(if (identical(arm, "lane")) c(LANE, REF, USER) else c(REF, USER))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressMessages(library(frmtmb))
suppressMessages(library(frmtmb.sample))
cat("ARM ", arm, "\n", sep = "")
LOG <- "C:/Users/adf44/source/r/frmtmb-wt-sampfix/dev/sampfix-log"

set.seed(405L)
n <- 300L
do <- data.frame(x = rnorm(n), z = rnorm(n),
                 fc = factor(sample(c("a", "b", "c"), n, TRUE)),
                 grp = factor(sample(c("p", "q"), n, TRUE)))
ct <- function(eta) {
  factor(cut(eta + rlogis(n), c(-Inf, -0.5, 0.7, Inf), labels = FALSE),
         ordered = TRUE)
}
do$yo <- ct(0.8 * do$x)
do$yo2 <- ct(-0.5 * do$z)
fits <- list(
  cumulative = frm(bf(yo ~ x), family = cumulative(), data = do),
  sratio_cs = frm(bf(yo ~ x + cs(fc)), family = sratio(), data = do),
  thres_gr = frm(bf(yo | thres(gr = grp) ~ x), family = cumulative(),
                 data = do),
  mv = frm(bf(yo ~ x) + bf(yo2 ~ z), family = cumulative(), data = do))

res <- list()
for (nm in names(fits)) {
  fit <- fits[[nm]]
  ds <- suppressWarnings(suppressMessages(
    frm_sample(fit, chains = 1, iter = 300, refresh = 0, seed = 3)))
  resp <- names(fit$spec$responses)[1L]
  res[[nm]] <- list(
    variables_fit = variables(fit), variables_ds = variables(ds),
    draws = ds$draws,
    epred = posterior_epred(ds, ndraws = 20, resp = resp),
    fixef = fixef(ds), ds = ds)
  cat("\n-- ", nm, "\n  variables(fit): ",
      paste(variables(fit), collapse = " "), "\n  variables(ds):  ",
      paste(variables(ds), collapse = " "), "\n", sep = "")
  vf <- variables(fit)
  cat("  every variables(fit) name that is not sd_/cor_ is a draws column: ",
      all(vf[!grepl("^(sd|cor)_", vf)] %in% variables(ds)), "\n", sep = "")
}
saveRDS(lapply(res, function(r) r[setdiff(names(r), "ds")]),
        file.path(LOG, paste0("07-", arm, ".rds")))

if (identical(arm, "lane") && file.exists(file.path(LOG, "07-ref.rds"))) {
  ref <- readRDS(file.path(LOG, "07-ref.rds"))
  cat("\n== lane against the reference build (same sampler, same seed)\n")
  for (nm in names(res)) {
    a <- res[[nm]]
    b <- ref[[nm]]
    ds <- a$ds
    # the reference's draws mapped the lane's way: same values?
    nc <- frmtmb.sample:::draws_natural_cols(ds$fit)
    bm <- b$draws
    for (o in nc$ordinal) {
      j <- match(o$internal, colnames(bm))
      bm[, j] <- frmtmb.sample:::draws_rowmap(bm[, j, drop = FALSE], o$map)
      colnames(bm)[j] <- o$names
    }
    same_names <- identical(colnames(bm), colnames(a$draws))
    dmax <- max(abs(bm - a$draws))
    emax <- max(abs(a$epred - b$epred))
    fmax <- max(abs(a$fixef - b$fixef))
    same_rows <- identical(rownames(a$fixef), rownames(b$fixef))
    # the old storage read by the new build
    old <- ds
    old$draws <- b$draws
    e_old <- posterior_epred(old, ndraws = 20,
                             resp = names(ds$fit$spec$responses)[1L])
    f_old <- fixef(old)
    cat(sprintf(paste0("  %-11s names %s  max|draws diff| %.3g  ",
                       "max|epred diff| %.3g  fixef rows same %s, ",
                       "max diff %.3g  old storage: epred %.3g, fixef %.3g\n"),
                nm, same_names, dmax, emax, same_rows, fmax,
                max(abs(e_old - b$epred)), max(abs(f_old - b$fixef))))
  }
}
cat("\nDONE\n")
