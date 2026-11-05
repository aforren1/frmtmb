# A dirichlet response whose cbind() holds an expression: is the term dropped?
source("C:/Users/adf44/source/r/frmtmb-wt-latres/dev/review-latres-common.R")
set.seed(20261102)
d <- expand.grid(t = 1:5, g = factor(1:10))
d <- d[order(d$g, d$t), ]
n <- nrow(d)
d$x <- rnorm(n)
d$un <- rbeta(n, 2, 2)
r <- msg_of(frm(bf(cbind(un, 1 - un) ~ x + ar(t, g, cov = TRUE)), family = dirichlet(), data = d, dry_run = "frame"))
print(r$err); print(r$warn)
if (!is.null(r$val)) {
  fr <- r$val
  cat("autocor", length(fr$autocor), "latres", length(fi("frame_latres")(fr)), "\n")
  print(names(fr$linpreds))
  str(lapply(fr$linpreds, function(l) colnames(l$X)))
}
r2 <- msg_of(frm(bf(cbind(un, 1 - un) ~ x + ar(t, g, cov = TRUE)), family = dirichlet(), data = d))
print(r2$err)
if (!is.null(r2$val)) print(summary(r2$val))
if (rev_arm == "lane") {
  rb <- msg_of(frm(cbind(un, 1 - un) ~ x, family = dirichlet(), data = d))
  print(rb$err)
}
