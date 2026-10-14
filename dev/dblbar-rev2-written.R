# Reviewer round 2, lane dblbar: does written_text leak into, or break,
# update(), refit(), formula(), print(), save and reload, spec identity,
# the duplicate message, frmtmb.sample and as_tmbstan()?
# Data seed 20261013 (dev/dblbar-data.R); sampler seed 3.
# Usage: Rscript ... base|lane. Output: dev/dblbar-rev2-log/written-<arm>.txt
arm <- commandArgs(TRUE)[1]
lib <- if (identical(arm, "base")) character(0) else
  "C:/Users/adf44/source/r/wt-dblbar-lib"
.libPaths(c(lib, "C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
wt <- "C:/Users/adf44/source/r/frmtmb-wt-dblbar"
source(file.path(wt, "dev/dblbar-data.R"))
suppressPackageStartupMessages({library(frmtmb); library(frmtmb.sample)})
sink(file.path(wt, sprintf("dev/dblbar-rev2-log/written-%s.txt", arm)),
     split = TRUE)
cat("lib:", find.package("frmtmb"), "\n")
dd <- dblbar_data()
tc <- function(expr) tryCatch(expr, error = function(e) {
  paste("ERROR:", substr(gsub("\n", " ", conditionMessage(e)), 1, 200))
})
fit_of <- function(txt) frm(bf(stats::as.formula(txt)), family = gaussian(),
                            data = dd)
sp <- c(dbar = "y ~ f + (1 + f || g)",
        corF = "y ~ f + (1 + f | gr(g, cor = FALSE))",
        diag = "y ~ f + diag(1 + f | g)")
# literal calls, so update() can re-evaluate them
fits <- list(
  dbar = frm(bf(y ~ f + (1 + f || g)), family = gaussian(), data = dd),
  corF = frm(bf(y ~ f + (1 + f | gr(g, cor = FALSE))), family = gaussian(),
             data = dd),
  diag = frm(bf(y ~ f + diag(1 + f | g)), family = gaussian(), data = dd))
# spec identity: the parsed spec with and without the written text
specs <- lapply(sp, function(t) frm(bf(stats::as.formula(t)),
                                    family = gaussian(), data = dd,
                                    dry_run = "spec"))
strip <- function(x) {
  if (is.list(x)) {
    x[["written_text"]] <- NULL
    x[["written"]] <- NULL
    x[] <- lapply(x, strip)
  }
  if (inherits(x, "formula")) attributes(x) <- NULL
  x
}
re_of <- function(s) s$responses[[1]]$dpars$mu$re
cat("written_text per spelling:",
    vapply(specs, function(s) re_of(s)[[1]][["written_text"]] %||% "NULL",
           ""), "\n")
cat("re spec identical dbar/corF:", identical(re_of(specs$dbar),
                                              re_of(specs$corF)),
    " dbar/diag:", identical(re_of(specs$dbar), re_of(specs$diag)), "\n")
cat("re spec identical after dropping written*: dbar/corF",
    identical(strip(re_of(specs$dbar)), strip(re_of(specs$corF))),
    " dbar/diag", identical(strip(re_of(specs$dbar)),
                            strip(re_of(specs$diag))), "\n")
cat("estimates identical dbar/corF/diag:",
    identical(fits$dbar$estimates, fits$corF$estimates),
    identical(fits$dbar$estimates, fits$diag$estimates), "\n")
cat("logLik identical:", identical(logLik(fits$dbar), logLik(fits$corF)),
    identical(logLik(fits$dbar), logLik(fits$diag)), "\n")
A <- fits$dbar
cat("formula(A):", deparse1(formula(A)$formula %||% formula(A)), "\n")
cat("attr on stored formula:",
    paste(names(attributes(formula(A)$formula %||% formula(A))),
          collapse = ","), "\n")
pr <- utils::capture.output(print(A))
cat("print lines mentioning written/diag/||:",
    paste(grep("diag|[|][|]|written", pr, value = TRUE), collapse = " / "),
    "\n")
cat("summary random rows:", rownames(summary(A)$random$g), "\n")
u1 <- update(A, formula. = ~ . + x)
cat("update + x:", deparse1(formula(u1)$formula %||% formula(u1)), "|",
    grep("^sd_|^cor_", variables(u1), value = TRUE), "\n")
u2 <- update(A, newdata = dd[1:300, ])
cat("update newdata blocks:", vapply(u2$frame$re_blocks, function(b) {
  paste0(b$covstruct, "[", paste(b$cnms, collapse = ","), "]")
}, ""), "\n")
u3 <- tc(update(A, formula. = ~ . - (1 + f || g) + (1 + f | g)))
cat("update || to |:", if (is.character(u3)) u3 else
  paste(vapply(u3$frame$re_blocks, `[[`, "", "covstruct"),
        grep("^cor_", variables(u3), value = TRUE)), "\n")
u4 <- update(fits$corF, formula. = ~ . + x)
cat("update corF + x:", deparse1(formula(u4)$formula %||% formula(u4)),
    vapply(u4$frame$re_blocks, `[[`, "", "covstruct"), "\n")
set.seed(5)
ynew <- dd$y + rnorm(nrow(dd), 0, 0.3)
r1 <- refit(A, newresp = ynew)
d1 <- frm(bf(y ~ f + (1 + f || g)), family = gaussian(),
          data = transform(dd, y = ynew))
cat(sprintf("refit logLik %.10f direct %.10f\n", logLik(r1), logLik(d1)))
f <- file.path(tempdir(), "dblbar-rev2-fit.rds")
saveRDS(A, f)
B <- readRDS(f)
cat("reloaded: fitted identical", identical(fitted(B), fitted(A)),
    "predict newdata identical", identical(
      fitted(B, newdata = dd[1:5, ]), fitted(A, newdata = dd[1:5, ])),
    "\n")
ub <- tc(update(B, formula. = ~ . + x))
cat("reloaded update:", if (is.character(ub)) ub else
  identical(ub$estimates, u1$estimates), "\n")
# the duplicate message for each spelling, and mixed-spelling twins
for (t in c("y ~ f + (1 + f || g) + (1 | g)",
            "y ~ f + (1 + f | gr(g, cor = FALSE)) + (1 | g)",
            "y ~ f + diag(1 + f | g) + (1 | g)",
            "y ~ f + (1 + f || g) + (1 + f | gr(g, cor = FALSE))",
            "y ~ f + (1 + f || g) + diag(1 + f | g)")) {
  r <- tc(fit_of(t))
  cat(t, "=>", if (is.character(r)) r else
    paste("ok", length(r$frame$re_blocks), "block(s)"), "\n")
}
# frmtmb.sample and as_tmbstan
dr <- frm_sample(A, chains = 2, iter = 500, seed = 3, refresh = 0)
cat("draws VarCorr rows:", rownames(VarCorr(dr)$g$sd), "\n")
cat("draws summary ok:", !inherits(tc(summary(dr)), "character"), "\n")
dr2 <- frm_sample(fits$corF, chains = 2, iter = 500, seed = 3, refresh = 0)
cat("draws identical across spellings:",
    identical(as.matrix(dr), as.matrix(dr2)), "\n")
ts <- tc(as_tmbstan(A, chains = 1, iter = 300, seed = 1, refresh = 0))
cat("as_tmbstan:", if (is.character(ts)) ts else class(ts), "\n")
sink()
