# Reviewer: does bhaz() compose with update(), newdata, get_prior(),
# frm_sample(), other addition terms and save/reload?
#   Rscript dev/review-bhaz-compose.R > dev/review-bhaz-log/compose.txt 2>&1
.libPaths(c("C:/Users/adf44/source/r/wt-bhaz-lib",
            "C:/Users/adf44/source/r/rellib-r12",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true")
suppressPackageStartupMessages(library(frmtmb.survival))
cat("lib:", find.package("frmtmb"), find.package("frmtmb.survival"), "\n")
set.seed(2); n <- 300; x <- rnorm(n); z <- rnorm(n)
tt <- stats::rexp(n, exp(-0.5 + 0.7 * x)); ct <- stats::rexp(n, 0.2)
dd <- data.frame(time = pmin(tt, ct), cens = as.numeric(tt > ct), x = x,
                 z = z, w = runif(n, 0.5, 2))
dd$entry <- dd$time * runif(n, 0, 0.3)
dd$c2 <- ifelse(dd$cens == 1, 2, 0)
dd$upper <- ifelse(dd$cens == 1, dd$time * 1.5, dd$time)
cat("range time", format(range(dd$time), digits = 8), "\n")
q <- function(expr) {
  tryCatch(suppressWarnings(suppressMessages(expr)),
           error = function(e) structure(conditionMessage(e),
                                         class = "err"))
}
K <- function(f) if (inherits(f, "err")) NA else length(cox_baseline(f))
same <- function(lab, a, b) {
  if (inherits(a, "err") || inherits(b, "err")) {
    cat(sprintf("%-46s ERR a: %s | b: %s\n", lab,
                if (inherits(a, "err")) a else "ok",
                if (inherits(b, "err")) b else "ok"))
    return(invisible())
  }
  cat(sprintf("%-46s K %s/%s logLik %.10f/%.10f identical ll %s est %s\n",
              lab, K(a), K(b), as.numeric(logLik(a)), as.numeric(logLik(b)),
              identical(logLik(a), logLik(b)),
              identical(unname(a$estimates), unname(b$estimates))))
}
f6b <- q(frm(bf(time | cens(cens) + bhaz(df = 6) ~ x), family = cox(),
             data = dd))
f6c <- q(frm(bf(time | cens(cens) ~ x), family = cox(df = 6), data = dd))
f5 <- q(frm(bf(time | cens(cens) ~ x), family = cox(), data = dd))
same("bhaz(df=6) vs cox(df=6)", f6b, f6c)
cat("stored family df:", f6b$bform$family$df %||% "NULL",
    " written attr present:",
    !is.null(attr(f6b$bform, "frmtmb_written")), "\n")

cat("\n== update() ==\n")
same("update(+ z) vs cox(df=6) + z",
     q(update(f6b, . ~ . + z)),
     q(frm(bf(time | cens(cens) ~ x + z), family = cox(df = 6), data = dd)))
same("update(drop bhaz, full formula) vs cox()",
     q(update(f6b, time | cens(cens) ~ x)), f5)
same("update(drop bhaz, bf()) vs cox()",
     q(update(f6b, bf(time | cens(cens) ~ x))), f5)
same("update(newdata = 250 rows) vs cox(df=6)",
     q(update(f6b, newdata = dd[1:250, ])),
     q(frm(bf(time | cens(cens) ~ x), family = cox(df = 6),
           data = dd[1:250, ])))
same("update(family = cox()) vs bhaz(df=6)",
     q(update(f6b, family = cox())), f6c)
same("update(family = cox(df = 7))", q(update(f6b, family = cox(df = 7))),
     f6c)
same("update(cox(df=6) fit, + bhaz(df = 7))",
     q(update(f6c, time | cens(cens) + bhaz(df = 7) ~ x)),
     q(frm(bf(time | cens(cens) ~ x), family = cox(df = 7), data = dd)))
same("update(cox() fit, + bhaz(df = 7))",
     q(update(f5, time | cens(cens) + bhaz(df = 7) ~ x)),
     q(frm(bf(time | cens(cens) ~ x), family = cox(df = 7), data = dd)))
same("update twice (+z, then -z)",
     q(update(update(f6b, . ~ . + z), . ~ . - z)), f6c)

cat("\n== predict(newdata) ==\n")
nd <- dd[c(3, 50, 120), ]
pb <- q(predict(f6b, newdata = nd)); pc <- q(predict(f6c, newdata = nd))
cat("predict identical:", identical(pb, pc), "\n")
nd2 <- nd; nd2$time <- NULL; nd2$cens <- NULL
cat("predict without response columns:",
    if (inherits(q(predict(f6b, newdata = nd2)), "err"))
      q(predict(f6b, newdata = nd2)) else "ok", "\n")

cat("\n== get_prior() ==\n")
gb <- q(get_prior(bf(time | cens(cens) + bhaz(df = 6) ~ x), family = cox(),
                  data = dd))
gc <- q(get_prior(bf(time | cens(cens) ~ x), family = cox(df = 6),
                  data = dd))
cat("get_prior identical:", identical(gb, gc), "\n")
print(as.data.frame(gb)[as.data.frame(gb)$class == "sbhaz", ])

cat("\n== other addition terms ==\n")
same("weights(w) + bhaz(df=6)",
     q(frm(bf(time | cens(cens) + weights(w) + bhaz(df = 6) ~ x),
           family = cox(), data = dd)),
     q(frm(bf(time | cens(cens) + weights(w) ~ x), family = cox(df = 6),
           data = dd)))
same("bhaz first, then cens",
     q(frm(bf(time | bhaz(df = 6) + cens(cens) ~ x), family = cox(),
           data = dd)), f6c)
same("trunc(lb) + bhaz(knots, Boundary.knots)",
     q(frm(bf(time | cens(cens) + trunc(lb = entry) +
                bhaz(knots = c(0.5, 1.1), Boundary.knots = c(0, 30)) ~ x),
           family = cox(), data = dd)),
     q(frm(bf(time | cens(cens) + trunc(lb = entry) ~ x),
           family = cox(knots = c(0.5, 1.1), Boundary.knots = c(0, 30)),
           data = dd)))
same("interval cens + bhaz(df=6)",
     q(frm(bf(time | cens(c2, upper) + bhaz(df = 6) ~ x), family = cox(),
           data = dd)),
     q(frm(bf(time | cens(c2, upper) ~ x), family = cox(df = 6),
           data = dd)))
mx <- max(dd$time)
r <- q(frm(bf(time | cens(c2, upper) + bhaz(Boundary.knots = c(0, mx)) ~ x),
           family = cox(), data = dd))
cat("interval cens, upper bound past Boundary.knots[2] (",
    sum(dd$upper > mx), "rows ):",
    if (inherits(r, "err")) r else
      sprintf("fit, logLik %.6f", as.numeric(logLik(r))), "\n")
same("family = 'cox' string + bhaz(df=6)",
     q(frm(bf(time | cens(cens) + bhaz(df = 6) ~ x), family = "cox",
           data = dd)), f6c)
same("brmsfamily('cox') + bhaz(df=6)",
     q(frm(bf(time | cens(cens) + bhaz(df = 6) ~ x),
           family = brmsfamily("cox"), data = dd)), f6c)
same("bf(...) + cox() spelling",
     q(frm(bf(time | cens(cens) + bhaz(df = 6) ~ x) + cox(), data = dd)),
     f6c)
same("frailty (1|g) + bhaz",
     q(frm(bf(time | cens(cens) + bhaz(df = 6) ~ x + (1 | g)),
           family = cox(), data = transform(dd, g = factor(rep(1:15, 20))))),
     q(frm(bf(time | cens(cens) ~ x + (1 | g)), family = cox(df = 6),
           data = transform(dd, g = factor(rep(1:15, 20))))))

cat("\n== spellings brms takes ==\n")
same("bhaz(df=3, intercept=FALSE)",
     q(frm(bf(time | cens(cens) + bhaz(df = 3, intercept = FALSE) ~ x),
           family = cox(), data = dd)),
     q(frm(bf(time | cens(cens) ~ x),
           family = cox(df = 3, intercept = FALSE), data = dd)))
cat("cox(df = 3) + bhaz(intercept = FALSE):",
    q(frm(bf(time | cens(cens) + bhaz(intercept = FALSE) ~ x),
          family = cox(df = 3), data = dd))[1], "\n")
same("bhaz(d = 6) [brms: df = 6]",
     q(frm(bf(time | cens(cens) + bhaz(d = 6) ~ x), family = cox(),
           data = dd)), f6c)
same("bhaz(Boundary = c(0, 30)) [brms ignores it]",
     q(frm(bf(time | cens(cens) + bhaz(Boundary = c(0, 30)) ~ x),
           family = cox(), data = dd)), f5)
same("bhaz(deg = 2) [brms: degree = 2]",
     q(frm(bf(time | cens(cens) + bhaz(deg = 2) ~ x), family = cox(),
           data = dd)),
     q(frm(bf(time | cens(cens) ~ x), family = cox(degree = 2), data = dd)))
same("bhaz(gr = NA)",
     q(frm(bf(time | cens(cens) + bhaz(gr = NA) ~ x), family = cox(),
           data = dd)), f5)
cat("bhaz(6) [brms: one group '6', df = 5]:",
    q(frm(bf(time | cens(cens) + bhaz(6) ~ x), family = cox(),
          data = dd))[1], "\n")
same("bhaz(NA, 6) [brms: df = 6]",
     q(frm(bf(time | cens(cens) + bhaz(NA, 6) ~ x), family = cox(),
           data = dd)), f6c)
same("bhaz(df = 6L)",
     q(frm(bf(time | cens(cens) + bhaz(df = 6L) ~ x), family = cox(),
           data = dd)), f6c)
same("bhaz(df = 6) with cox(df = 6L)",
     q(frm(bf(time | cens(cens) + bhaz(df = 6) ~ x), family = cox(df = 6L),
           data = dd)), f6c)
same("bhaz(knots=) with cox(knots=) reordered",
     q(frm(bf(time | cens(cens) + bhaz(knots = c(1.1, 0.5)) ~ x),
           family = cox(knots = c(0.5, 1.1)), data = dd)),
     q(frm(bf(time | cens(cens) ~ x), family = cox(knots = c(0.5, 1.1)),
           data = dd)))
same("bhaz(knots = numeric(0))",
     q(frm(bf(time | cens(cens) + bhaz(knots = numeric(0)) ~ x),
           family = cox(), data = dd)),
     q(frm(bf(time | cens(cens) ~ x), family = cox(df = 4), data = dd)))
kn <- c(0.5, 1.1)
fk <- function() {
  kn <- c(0.2, 0.4, 0.8)
  frm(bf(time | cens(cens) + bhaz(knots = kn) ~ x), family = cox(),
      data = dd)
}
fkk <- q(fk())
cat("knots from a function's environment: K", K(fkk), "(3 knots -> 7)\n")
cat("update of that fit (+ z): K", K(q(update(fkk, . ~ . + z))), "\n")

cat("\n== save and reload ==\n")
rds <- "C:/Users/adf44/source/r/frmtmb-wt-bhaz/dev/review-bhaz-fit.rds"
saveRDS(f6b, rds)
child <- function(lines) {
  f <- tempfile(fileext = ".R")
  writeLines(c(sprintf(".libPaths(%s)", paste(deparse(.libPaths()),
                                               collapse = "")), lines), f)
  out <- system2(file.path(R.home("bin"), "Rscript"), c("--vanilla", f),
                 stdout = TRUE, stderr = TRUE)
  cat(paste0("   ", out), sep = "\n")
}
child(c("suppressMessages(library(frmtmb))",
        sprintf("f <- readRDS('%s')", rds),
        "r <- tryCatch(predict(f), error = function(e) conditionMessage(e))",
        "cat('frmtmb only, predict:', substr(paste(r, collapse=' '), 1, 300), '\\n')",
        "r <- tryCatch(update(f, . ~ . + z), error = function(e) conditionMessage(e))",
        "cat('frmtmb only, update:', substr(paste(r, collapse=' '), 1, 300), '\\n')"))
child(c("suppressMessages(library(frmtmb.survival))",
        sprintf("f <- readRDS('%s')", rds),
        "nd <- f$data[c(3, 50, 120), ]",
        "p <- predict(f, newdata = nd)",
        "cat('with survival, predict(newdata) rows', NROW(p), '\\n')",
        "u <- suppressWarnings(update(f, . ~ . + z))",
        "cat('with survival, update(+z) logLik', format(as.numeric(logLik(u)), digits = 12), ' K', length(cox_baseline(u)), '\\n')"))
unlink(rds)
