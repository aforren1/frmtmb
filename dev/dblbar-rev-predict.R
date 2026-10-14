# Reviewer, lane dblbar: prediction from a factor || fit against the same
# model written as separate blocks over numeric dummies (a different
# parse), at old levels, new levels (every sample_new_levels option),
# re_formula partial terms, newdata with dropped or character factor
# levels, and poly() at new data. Data: dev/dblbar-data.R, seed 20261013.
# Output: dev/dblbar-rev-log/predict.txt
.libPaths(c("C:/Users/adf44/source/r/wt-dblbar-lib",
            "C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
wt <- "C:/Users/adf44/source/r/frmtmb-wt-dblbar"
source(file.path(wt, "dev/dblbar-data.R"))
suppressPackageStartupMessages(library(frmtmb))
sink(file.path(wt, "dev/dblbar-rev-log/predict.txt"), split = TRUE)
cat("lib:", find.package("frmtmb"), "\n")
dd <- dblbar_data()
dd$fb <- as.numeric(dd$f == "b")
dd$fc <- as.numeric(dd$f == "c")
rel <- function(a, b) max(abs(a - b) / pmax(abs(b), 1e-8))
A <- frm(bf(y ~ f + (1 + f || g)), family = gaussian(), data = dd)
B <- frm(bf(y ~ f + (1 | g) + (0 + fb | g) + (0 + fc | g)),
         family = gaussian(), data = dd)
cat("blocks A", length(A$frame$re_blocks), "B", length(B$frame$re_blocks),
    "\n")
cat(sprintf("logLik A %.10f B %.10f rel %.2e\n", logLik(A), logLik(B),
            rel(as.numeric(logLik(A)), as.numeric(logLik(B)))))
nd_old <- dd[c(1, 5, 9, 40), ]
nd_new <- nd_old
nd_new$g <- c("NEW1", "NEW2", "NEW1", "NEW3")
cmp <- function(lab, a, b) {
  if (is.character(a) || is.character(b)) {
    cat(sprintf("%-34s A: %s\n%-34s B: %s\n", lab, substr(a, 1, 110), "",
                substr(b, 1, 110)))
    return(invisible())
  }
  cat(sprintf("%-34s Est rel %.2e  Err rel %.2e\n", lab,
              rel(a[, 1], b[, 1]), rel(a[, 2], b[, 2])))
}
tc <- function(expr) tryCatch(expr, error = function(e) conditionMessage(e))
cmp("fitted old levels", fitted(A, newdata = nd_old),
    fitted(B, newdata = nd_old))
cmp("fitted new gaussian",
    fitted(A, newdata = nd_new, allow_new_levels = TRUE,
           sample_new_levels = "gaussian"),
    fitted(B, newdata = nd_new, allow_new_levels = TRUE,
           sample_new_levels = "gaussian"))
set.seed(3)
oa <- fitted(A, newdata = nd_new, allow_new_levels = TRUE,
             sample_new_levels = "old_levels")
set.seed(3)
ob <- fitted(B, newdata = nd_new, allow_new_levels = TRUE,
             sample_new_levels = "old_levels")
cmp("fitted new old_levels (seed 3)", oa, ob)
# old_levels must equal the fit at the seen level it picked: the row
# matches one seen level's fitted value
fa_all <- fitted(A, newdata = transform(nd_new[rep(1, 24), ],
                                        g = levels(dd$g)))
cat("old_levels row 1 equals a seen level:",
    any(abs(fa_all[, 1] - oa[1, 1]) < 1e-12 * abs(oa[1, 1])), "\n")
cmp("fitted new uncertainty",
    tc(fitted(A, newdata = nd_new, allow_new_levels = TRUE,
              sample_new_levels = "uncertainty")),
    tc(fitted(B, newdata = nd_new, allow_new_levels = TRUE,
              sample_new_levels = "uncertainty")))
cmp("fitted new, allow FALSE", tc(fitted(A, newdata = nd_new)),
    tc(fitted(B, newdata = nd_new)))
set.seed(4)
pa <- predict(A, newdata = nd_new, allow_new_levels = TRUE,
              sample_new_levels = "gaussian", ndraws = 20000)
set.seed(4)
pb <- predict(B, newdata = nd_new, allow_new_levels = TRUE,
              sample_new_levels = "gaussian", ndraws = 20000)
cat("predict new gaussian A Est", sprintf("%.3f", pa[, 1]), "Err",
    sprintf("%.3f", pa[, 2]), "\n")
cat("predict new gaussian B Est", sprintf("%.3f", pb[, 1]), "Err",
    sprintf("%.3f", pb[, 2]), "\n")
set.seed(5)
qa <- predict(A, newdata = nd_new, allow_new_levels = TRUE,
              sample_new_levels = "old_levels", ndraws = 20000)
set.seed(5)
qb <- predict(B, newdata = nd_new, allow_new_levels = TRUE,
              sample_new_levels = "old_levels", ndraws = 20000)
cat("predict new old_levels A Est", sprintf("%.3f", qa[, 1]), "\n")
cat("predict new old_levels B Est", sprintf("%.3f", qb[, 1]), "\n")
cmp("re_formula ~ (1 | g)",
    fitted(A, newdata = nd_old, re_formula = ~ (1 | g)),
    fitted(B, newdata = nd_old, re_formula = ~ (1 | g)))
cmp("re_formula ~ (1 | g) + (0 + fb | g)",
    fitted(A, newdata = nd_old, re_formula = ~ (1 + f || g)),
    fitted(B, newdata = nd_old,
           re_formula = ~ (1 | g) + (0 + fb | g) + (0 + fc | g)))
cmp("re_formula NA", fitted(A, newdata = nd_old, re_formula = NA),
    fitted(B, newdata = nd_old, re_formula = NA))
cat("re_formula ~ (0 + f || g) on A:",
    substr(tc(fitted(A, newdata = nd_old, re_formula = ~ (0 + f || g))), 1,
           150), "\n")
cat("re_formula ~ (0 + f | g) on A:",
    substr(tc(fitted(A, newdata = nd_old, re_formula = ~ (0 + f | g))), 1,
           150), "\n")
cat("re_formula ~ (1 + f | g) on A same as full:",
    identical(fitted(A, newdata = nd_old, re_formula = ~ (1 + f | g)),
              fitted(A, newdata = nd_old)), "\n")
x <- tc(fitted(A, newdata = nd_old, re_formula = ~ (f | g)))
cat("re_formula ~ (f | g) on A same as full:",
    if (is.character(x)) x else identical(x, fitted(A, newdata = nd_old)),
    "\n")
# a partial lhs that names a column: (0 + fb | g)? brms reads 0 + f
x <- tc(fitted(A, newdata = nd_old, re_formula = ~ (1 + x || g)))
cat("re_formula naming a column not in the fit:", substr(
  if (is.character(x)) x else "no error", 1, 150), "\n")
# newdata whose factor has dropped levels, or is character
nd_c <- dd[dd$f == "c", ][1:4, ]
ref <- fitted(A, newdata = nd_c)
nd_c2 <- nd_c
nd_c2$f <- factor(as.character(nd_c2$f))
nd_c3 <- nd_c
nd_c3$f <- as.character(nd_c3$f)
cmp("newdata f dropped levels", tc(fitted(A, newdata = nd_c2)), ref)
cmp("newdata f character", tc(fitted(A, newdata = nd_c3)), ref)
# poly(): coefficients carried to new data
P <- frm(bf(y ~ poly(x, 2) + (1 + poly(x, 2) || g)), family = gaussian(),
         data = dd)
cat("poly blocks", vapply(P$frame$re_blocks, function(b) {
  paste0(b$covstruct, "[", paste(b$cnms, collapse = ","), "]")
}, ""), "\n")
cmp("poly fitted newdata subset", fitted(P, newdata = dd[1:20, ]),
    fitted(P)[1:20, , drop = FALSE])
cmp("poly re_formula NA subset",
    fitted(P, newdata = dd[1:20, ], re_formula = NA),
    fitted(P, re_formula = NA)[1:20, , drop = FALSE])
ndp <- dd[1:5, ]
ndp$g <- "NEWP"
cmp("poly new level vs re_formula NA",
    fitted(P, newdata = ndp, allow_new_levels = TRUE),
    fitted(P, newdata = ndp, re_formula = NA))
P2 <- frm(bf(y ~ poly(x, 2) + diag(1 + poly(x, 2) | g)),
          family = gaussian(), data = dd)
cat("poly || vs diag() identical:", identical(P$estimates, P2$estimates),
    "\n")
sink()
