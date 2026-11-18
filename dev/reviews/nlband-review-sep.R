# Review of lane nlband: follow-ups. (1) the separated bernoulli fit the
# hole start returned silently (BSEP seed 2): its standard errors, and a
# direct fit from inside the hole; (2) a pole 0.003 from the hole's edge
# row (nearer the wall than the median gap); (3) separated data, seeds 3
# to 10, counts. NLBAND_ARM=before|after Rscript nlband-review-sep.R
source("dev/nlband-common.R")
show <- function(tag, r) {
  if (iserr(r)) { cat(tag, "| ERROR", r$msg, "\n"); return(invisible()) }
  se <- tryCatch(sqrt(diag(suppressWarnings(vcov(r)))),
                 error = function(e) NA)
  cat(sprintf("%s | ll %.6g code %d (%s) beta %s se %s | warn %s msg %s\n",
              tag, llv(r), r$opt$convergence, r$opt$message,
              paste(signif(r$estimates$beta, 6), collapse = " "),
              paste(signif(se, 4), collapse = " "), wsum(r), msum(r)))
  for (w in attr(r, "w")) cat("   WARN", substr(w, 1, 300), "\n")
}
for (s in 1:10) {
  d <- band_data(s)
  d$y <- as.integer(d$x > 2)
  show(paste("SEP default seed", s),
       tt(frm(band_form, data = d, family = bernoulli())))
  if (s <= 2) {
    show(paste("SEP from (1,-2) nl_starts off seed", s),
         tt(frm(band_form, data = d, family = bernoulli(),
                start = list(beta = c(1, -2)),
                control = frmtmb_control(nl_starts = FALSE))))
    show(paste("SEP glm x > 2 check seed", s),
         tt(frm(y ~ I(x > 2), data = d, family = bernoulli())))
  }
}
# a pole 0.003 below the first row above the hole
for (s in 1:3) {
  d <- band_data(s)
  xe <- min(d$x[d$x > 2])
  p <- xe - 0.003
  set.seed(s + 700)
  d$y <- 1 / (d$x - p) + rnorm(nrow(d), 0, 0.2)
  r <- tt(frm(band_form, data = d))
  g <- tt(frm(band_form, data = d, start = list(beta = c(1, -p)),
              control = frmtmb_control(nl_starts = FALSE)))
  cat(sprintf(paste0("NEARWALL seed %d pole %.5f gap %.5f | default %.6f b %.5f ",
                     "msg %s warn %s | from pole %.6f b %.5f\n"), s, p,
              median(diff(sort(unique(d$x)))), llv(r),
              if (iserr(r)) NA else r$estimates$beta[2], msum(r), wsum(r),
              llv(g), if (iserr(g)) NA else g$estimates$beta[2]))
}
cat("DONE\n")
