# Reviewer, priority 4: the lane's shapes against brms's own, on the
# SAME data. The brms side is dev/shapes-rev-brmsref.rds, written by a
# process that never loaded frmtmb.
#
#   Rscript dev/shapes-rev-brmscmp.R

.libPaths(c("C:/Users/adf44/source/r/shapes-lib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
TREE <- "C:/Users/adf44/source/r/frmtmb-wt-shapes"
source(file.path(TREE, "dev/shapes-rev-fixtures.R"))
B <- readRDS(file.path(TREE, "dev/shapes-rev-brmsref.rds"))

dd <- rev_data()
dd$y2 <- with(dd, rnorm(nrow(dd), 0.3 - 0.2 * dd$x, 1))

fits <- list(
  gaussian = quote(frm(bf(y ~ x + f) + gaussian(), data = dd)),
  binomial = quote(frm(bf(bin ~ x + z) + bernoulli(), data = dd)),
  mixed    = quote(frm(bf(ymix ~ x + (1 | g)) + gaussian(), data = dd)),
  ordinal  = quote(frm(bf(ord ~ x) + cumulative(), data = dd)),
  distreg  = quote(frm(bf(y ~ x + f, sigma ~ z) + gaussian(), data = dd)),
  multivar = quote(frm(bf(mvbind(y, y2) ~ x) + gaussian(), data = dd))
)

safe <- function(e) tryCatch(suppressWarnings(e), error = function(c)
  structure(list(msg = conditionMessage(c)), class = "revErr"))
fmt <- function(x) {
  if (inherits(x, "revErr")) return(paste0("REFUSED: ",
                                           substr(x$msg, 1, 90)))
  if (is.null(x)) return("NULL")
  if (is.list(x) && !is.data.frame(x))
    return(paste0("list(", paste(names(x), collapse = ","), ") = ",
                  paste(unlist(lapply(x, function(z)
                    paste(z, collapse = "/"))), collapse = " ")))
  if (is.data.frame(x)) return(paste0("df ", nrow(x), "x", ncol(x),
                                      " rows[", paste(rownames(x),
                                                      collapse = ","),
                                      "] cols[",
                                      paste(colnames(x), collapse = ","),
                                      "]"))
  if (is.array(x) || is.matrix(x)) return(paste0(
    "[", paste(dim(x), collapse = "x"), "] ",
    paste(vapply(dimnames(x), function(z)
      if (is.null(z)) "NULL" else paste0("{",
        paste(utils::head(z, 6), collapse = ","),
        if (length(z) > 6) ",..." else "", "}"), ""),
      collapse = " | ")))
  paste(utils::head(as.character(x), 12), collapse = ",")
}
row <- function(lab, bx, lx) {
  fb <- fmt(bx); fl <- fmt(lx)
  mark <- if (identical(fb, fl)) "  =  " else " DIFF"
  cat(sprintf("  %-16s %s\n", lab, mark))
  cat("      brms : ", substr(fb, 1, 170), "\n", sep = "")
  cat("      lane : ", substr(fl, 1, 170), "\n", sep = "")
}

for (nm in names(fits)) {
  if (is.null(B[[nm]])) next
  cat("\n################ ", nm, " ################\n", sep = "")
  f <- eval(fits[[nm]])
  b <- B[[nm]]
  row("fixef dimnames", b$fixef_dimnm, safe(dimnames(fixef(f))))
  row("fixef class", b$fixef_class, safe(class(fixef(f))))
  row("vcov dimnames", b$vcov_dimnm, safe(dimnames(vcov(f))))
  row("vcov corr dimnm", b$vcov_corr,
      safe(dimnames(vcov(f, correlation = TRUE))))
  row("summary names", b$summary_names, safe(names(summary(f))))
  row("summary class", b$summary_class, safe(class(summary(f))))
  row("$fixed rows", b$fixed_rows, safe(rownames(summary(f)$fixed)))
  row("$fixed cols", b$fixed_cols, safe(colnames(summary(f)$fixed)))
  row("$spec_pars", b$spec_pars, safe(summary(f)$spec_pars))
  row("$random", b$random, safe(summary(f)$random))
  row("$cor_pars", b$cor_pars, safe(summary(f)$cor_pars))
  row("ngrps", b$ngrps, safe(ngrps(f)))
  row("ngrps class", b$ngrps_class, safe(class(ngrps(f))))
  row("fitted dim", b$fitted_dim, safe(dim(fitted(f))))
  row("fitted dimnames", b$fitted_dimnm, safe(dimnames(fitted(f))))
  row("fitted class", b$fitted_class, safe(class(fitted(f))))
  row("resid dim", b$resid_dim, safe(dim(residuals(f))))
  row("resid dimnames", b$resid_dimnm, safe(dimnames(residuals(f))))
  row("predict dim", b$predict_dim, safe(dim(predict(f, ndraws = 200))))
  row("predict dimnames", b$predict_dimnm,
      safe(dimnames(predict(f, ndraws = 200))))
  row("predict draws dim", b$predict_nosum,
      safe(dim(predict(f, summary = FALSE, ndraws = 200))))
  row("predict ndraws=3", b$predict_nd3,
      safe(dim(predict(f, ndraws = 3))))
  row("variables", b$variables, safe(variables(f)))
}
