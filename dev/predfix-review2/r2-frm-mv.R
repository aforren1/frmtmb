source("C:/Users/adf44/source/r/frmtmb-wt-predfix/dev/predfix-review2/r2-helpers.R")
# frmtmb side of r2-brms-mv.R, same data (seed 606, saved by that script).
d <- readRDS("C:/Users/adf44/source/r/frmtmb-wt-predfix/dev/predfix-review2/r2-brms-mv-data.rds")
br <- readRDS("C:/Users/adf44/source/r/frmtmb-wt-predfix/dev/predfix-review2/r2-brms-mv.rds")
# frmtmb refuses cumulative() in a multivariate fit ("Families with
# extra parameters"), so o1 is left out and brms's layers 5:10 compared
m <- frm(mvbf(bf(y1 ~ x, family = gaussian()),
              bf(c1 ~ x, family = categorical()),
              bf(c2 ~ x, family = categorical()), rescor = FALSE), data = d)
show <- function(lab, f) {
  if (inherits(f, "error")) { cat(lab, ": ERROR", conditionMessage(f), "\n"); return() }
  cat(lab, ": dim", dim(f), "\n   layers:",
      paste0("'", dimnames(f)[[3]], "'", collapse = " "), "\n")
}
tc <- function(e) tryCatch(e, error = function(x) x)
f <- tc(fitted(m)); show("fitted(m)", f)
show("fitted(m, resp = c('c2', 'y1'))", tc(fitted(m, resp = c("c2", "y1"))))
show("fitted(m, resp = c('c1', 'y1'))", tc(fitted(m, resp = c("c1", "y1"))))
show("fitted(m, resp = 'c1')", tc(fitted(m, resp = "c1")))
nd <- data.frame(x = c(-1, 0, 2))
fn <- tc(fitted(m, newdata = nd)); show("fitted(m, newdata = 3 rows)", fn)
cat("same layer names and order as brms:",
    identical(dimnames(f)[[3]], dimnames(br$f)[[3]][5:10]), "\n")
cat("Estimate at newdata, frmtmb minus brms (ML against posterior mean):\n")
print(round(fn[, "Estimate", ] - br$fn[, "Estimate", 5:10], 4))
cat("Est.Error at newdata, frmtmb / brms:\n")
print(round(fn[, "Est.Error", ] / br$fn[, "Est.Error", 5:10], 2))
# each layer must equal the single-response call
one <- fitted(m, newdata = nd, resp = "c2")
cat("c2 layers equal fitted(resp = 'c2'):",
    isTRUE(all.equal(unname(fn[, , 5:6]), unname(one))), "\n")
show("fitted(m, resp = c('c1', 'c2'))", tc(fitted(m, resp = c("c1", "c2"))))
show("fitted(m, resp = 'nope')", tc(fitted(m, resp = "nope")))
show("fitted(m, resp = c('y1', 'y1'))", tc(fitted(m, resp = c("y1", "y1"))))
