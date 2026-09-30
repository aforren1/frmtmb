# Reviewer, re-check of the punch-round fixes other than nobs():
# rate() on newdata (refused at <= 0 on every newdata path, NA allowed),
# the character subset refusal (and its absent case), cat(x = 4).
# Seeds: 606 (rate), 223 (subset), 505 (cat), sampler 3.
# Log: dev/aterms2-rev-log-23-punch.txt
.libPaths(c("C:/Users/adf44/source/r/wt-aterms2-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressPackageStartupMessages({
  library(frmtmb)
  library(frmtmb.sample)
})
options(mc.cores = 1)
q <- function(e) suppressWarnings(suppressMessages(e))
show <- function(label, expr) {
  r <- tryCatch(q(expr), error = function(e) paste("ERROR:",
                                                   conditionMessage(e)))
  cat(sprintf("%-44s ", label))
  if (is.character(r) && length(r) == 1L) cat(substr(r, 1, 120), "\n")
  else if (is.numeric(r) || is.matrix(r)) {
    cat("returns", paste(format(head(as.numeric(r), 4), digits = 5),
                         collapse = " "), "\n")
  } else cat("returns", class(r)[1], "\n")
  invisible(r)
}
set.seed(606)
n <- 200
d <- data.frame(x = rnorm(n), time = runif(n, 0.5, 4))
d$y <- rpois(n, exp(0.2 + 0.5 * d$x) * d$time)
fr <- q(frm(y | rate(time) ~ x, data = d, family = poisson()))
fn <- q(frm(y | rate(time) ~ x, data = d, family = negbinomial()))
nd <- d[1:3, ]
cat("== rate() on newdata\n")
for (tt in list(c(1, -2, 3), c(1, 0, 3), c(1, NA, 3), c(1, 1e-3, 3))) {
  nd$time <- tt
  lab <- paste0("time = ", paste(tt, collapse = ","))
  show(paste("fitted", lab), fitted(fr, newdata = nd)[, 1])
  show(paste("predict", lab), {set.seed(1); predict(fr, newdata = nd,
                                                    ndraws = 50)[, 1]})
  show(paste("frm_linpred response", lab),
       frm_linpred(fr, newdata = nd, type = "response"))
  show(paste("negbin fitted", lab), fitted(fn, newdata = nd)[, 1])
  show(paste("simulate(newdata)", lab),
       unlist(simulate(fr, newdata = nd, seed = 1)))
}
nd$time <- c(1, 1e-3, 3)
mu <- fitted(fr, newdata = nd, dpar = "mu")[, 1]
cat("  at 1e-3: fitted / (mu * time) =",
    format(fitted(fr, newdata = nd)[, 1] / (mu * nd$time), digits = 15),
    "\n")
show("conditional_effects(conditions time = -1)",
     conditional_effects(fr, conditions = data.frame(time = -1)))
ds <- q(frm_sample(fr, chains = 1, iter = 200, refresh = 0, seed = 3))
nd$time <- c(1, -2, 3)
show("sample posterior_epred time -2", posterior_epred(ds, newdata = nd))
show("sample posterior_predict time -2",
     posterior_predict(ds, newdata = nd))
show("sample posterior_linpred(transform) -2",
     posterior_linpred(ds, newdata = nd, transform = TRUE))
nd$time <- c(1, NA, 3)
show("sample posterior_epred time NA", posterior_epred(ds, newdata = nd))

cat("\n== character subset\n")
set.seed(223)
n <- 60
dd <- data.frame(x = rnorm(n), z = rnorm(n), s = rep(c(TRUE, FALSE), n / 2))
dd$y1 <- 1 + dd$x + rnorm(n)
dd$y2 <- dd$z + rnorm(n)
dd$sch <- ifelse(dd$s, "TRUE", "FALSE")
dd$sfac <- factor(dd$sch)
dd$s01 <- as.integer(dd$s)
ref <- q(frm(bf(y1 | subset(s) ~ x) + bf(y2 ~ z), data = dd,
             family = gaussian()))
show("mv character", frm(bf(y1 | subset(sch) ~ x) + bf(y2 ~ z), data = dd,
                         family = gaussian()))
show("univariate character", frm(y1 | subset(sch) ~ x, data = dd))
for (sv in c("sfac", "s01")) {
  f <- q(eval(bquote(frm(bf(y1 | subset(.(as.name(sv))) ~ x) + bf(y2 ~ z),
                         data = dd, family = gaussian()))))
  cat(sprintf("%-44s fn identical to logical: %s\n", paste("mv", sv),
              identical(f$obj$fn(ref$opt$par), ref$obj$fn(ref$opt$par))))
}
fu <- q(frm(y1 | subset(s) ~ x, data = dd))
nd2 <- dd[1:6, ]
nd2$s <- nd2$sch
show("newdata with a character subset column",
     fitted(fu, newdata = nd2)[, 1])

cat("\n== cat()\n")
set.seed(505)
n <- 150
dc <- data.frame(x = rnorm(n))
dc$yi <- as.integer(cut(dc$x + rnorm(n), c(-Inf, -1, 0, 1, Inf),
                        labels = FALSE))
a <- show("cat(x = 4)", frm(yi | cat(x = 4) ~ x, data = dc,
                            family = cumulative()))
b <- q(frm(yi | cat(4) ~ x, data = dc, family = cumulative()))
if (inherits(a, "frmtmb_fit")) {
  cat("  fn identical cat(x = 4) vs cat(4):",
      identical(a$obj$fn(b$opt$par), b$obj$fn(b$opt$par)), "\n")
}
show("cat(k = 4)", frm(yi | cat(k = 4) ~ x, data = dc,
                       family = cumulative()))
show("cat(4, 5)", frm(yi | cat(4, 5) ~ x, data = dc, family = cumulative()))
