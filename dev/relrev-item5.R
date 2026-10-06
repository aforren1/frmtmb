# Reviewer, item 5: what reads a variable missing from newdata (or from
# data) out of the global environment, in frmtmb 0.68.0 and in brms
# 2.23.0. Seed 5.
#   Rscript dev/relrev-item5.R > dev/relrev-log/item5.txt 2>&1
.libPaths(c("C:/Users/adf44/source/r/rellib-r6",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(brms); library(frmtmb)})
cat("frmtmb", format(packageVersion("frmtmb")), find.package("frmtmb"), "\n")
try_ <- function(tag, expr) {
  r <- tryCatch(expr, error = function(e) paste("ERROR:", conditionMessage(e)))
  if (is.character(r) && length(r) == 1 && startsWith(r, "ERROR")) {
    cat(sprintf("%-55s %s\n", tag, substr(r, 1, 160)))
  } else {
    cat(sprintf("%-55s OK %s\n", tag,
                paste(format(head(as.numeric(r), 4), digits = 6), collapse = " ")))
  }
}
set.seed(5)
n <- 60
d <- data.frame(x = runif(n, 0, 10),
                f = factor(sample(c("a", "b"), n, TRUE)),
                g = factor(sample(letters[1:5], n, TRUE)),
                z = rnorm(n))
d$y <- sin(d$x) + (d$f == "b") + d$z + rnorm(n, 0, 0.3)
fit_gp <- frm(bf(y ~ gp(x, by = f, k = 8)), data = d)
fit_lin <- frm(bf(y ~ x + z), data = d)
fit_re <- frm(bf(y ~ x + (1 | g)), data = d)
fit_off <- frm(bf(y ~ x + offset(z)), data = d)
cat("\n== frmtmb newdata lacking a column, no global of that name\n")
try_("gp by: newdata lacks f", frm_linpred(fit_gp, newdata = data.frame(x = 1:3)))
try_("linear: newdata lacks z", frm_linpred(fit_lin, newdata = data.frame(x = 1:3)))
try_("re: newdata lacks g", frm_linpred(fit_re, newdata = data.frame(x = 1:3)))
try_("offset: newdata lacks z", frm_linpred(fit_off, newdata = data.frame(x = 1:3)))
cat("\n== the same, with globals f (a factor level), z and g defined\n")
f <- factor(c("b", "b", "b"), levels = c("a", "b"))
z <- c(100, 100, 100)
g <- factor(c("a", "a", "a"))
try_("gp by: newdata lacks f, global f = b", frm_linpred(fit_gp, newdata = data.frame(x = 1:3)))
try_("  reference: newdata f = b", frm_linpred(fit_gp, newdata = data.frame(x = 1:3, f = f)))
try_("linear: newdata lacks z, global z = 100", frm_linpred(fit_lin, newdata = data.frame(x = 1:3)))
try_("re: newdata lacks g, global g = a", frm_linpred(fit_re, newdata = data.frame(x = 1:3)))
try_("offset: newdata lacks z, global z = 100", frm_linpred(fit_off, newdata = data.frame(x = 1:3)))
f <- "C:/some/path/test-gp-by.R"
try_("gp by: global f = a path string (the runner's)", frm_linpred(fit_gp, newdata = data.frame(x = 1)))
rm(f, z, g)
cat("\n== frmtmb at fit time, data lacking the column, global present\n")
d2 <- d[, c("x", "y", "z")]
f <- d$f; g <- d$g
try_("frm gp(x, by = f), f global only", fixef(frm(bf(y ~ gp(x, by = f, k = 8)), data = d2)))
try_("frm (1 | g), g global only", fixef(frm(bf(y ~ x + (1 | g)), data = d2)))
w <- d$z
try_("frm y ~ x + w, w global only", fixef(frm(bf(y ~ x + w), data = d2)))
rm(f, g, w)

cat("\n== brms 2.23.0 (standata on newdata, no sampling)\n")
b_gp <- brm(brms::bf(y ~ gp(x, by = f, k = 8)), data = d, empty = TRUE)
b_lin <- brm(brms::bf(y ~ x + z), data = d, empty = TRUE)
b_re <- brm(brms::bf(y ~ x + (1 | g)), data = d, empty = TRUE)
sd_ <- function(fit, nd) {
  nd$y <- 0
  s <- standata(fit, newdata = nd)
  c(s$N, if (!is.null(s$X)) s$X[1, ] else NA)
}
try_("brms gp by: newdata lacks f", sd_(b_gp, data.frame(x = 1:3)))
try_("brms linear: newdata lacks z", sd_(b_lin, data.frame(x = 1:3)))
try_("brms re: newdata lacks g", sd_(b_re, data.frame(x = 1:3)))
f <- factor(c("b", "b", "b"), levels = c("a", "b"))
z <- c(100, 100, 100)
g <- factor(c("a", "a", "a"))
try_("brms gp by: newdata lacks f, global f", sd_(b_gp, data.frame(x = 1:3)))
try_("brms linear: newdata lacks z, global z = 100", sd_(b_lin, data.frame(x = 1:3)))
try_("brms re: newdata lacks g, global g", sd_(b_re, data.frame(x = 1:3)))
rm(f, z, g)
f <- d$f; g <- d$g; w <- d$z
try_("brms fit gp(x, by = f), f global only",
     standata(brm(brms::bf(y ~ gp(x, by = f, k = 8)), data = d2, empty = TRUE))$N)
try_("brms fit (1 | g), g global only",
     standata(brm(brms::bf(y ~ x + (1 | g)), data = d2, empty = TRUE))$N)
try_("brms fit y ~ x + w, w global only",
     standata(brm(brms::bf(y ~ x + w), data = d2, empty = TRUE))$N)
