# Review of lane quadesc: would a floor on the skew normal's tail argument
# remove the rule's NaN at the source? Runs the capture part of
# dev/reviews/quadesc-review-iso.R (base build, ALPHAX k = 9) to get the
# rule's anchor and x*, then tapes a replica rule for group 11 at that
# anchor with the skew normal's log density written out, once as RTMBdist
# writes it and once with pnorm's argument floored at -1e6 (where the
# density is exp(-5e11) = 0 already). RTMB refuses a comparison on its AD
# type, so the floor is the soft one -1e6 + log(1 + exp(w + 1e6)).
#   Rscript dev/reviews/quadesc-review-clamp.R
src <- readLines("dev/reviews/quadesc-review-iso.R")
stop_at <- grep("^tpl <- o$env$parameters", src, fixed = FALSE)
stop_at <- which(startsWith(src, "tpl <- o$env$parameters"))
commandArgs <- function(...) "9"
eval(parse(text = src[seq_len(stop_at)]))
cat("anchor b[11]", tpl$b[11], "\n")
gk <- ns$quad_gk_integrate()
mk <- function(dg, b0, floor) {
  nll <- function(p) {
    RTMB::getAll(p)
    mu <- beta[1] + beta[2] * dg$x + b[1]
    a <- betad[2] + betad[3] * dg$x
    s <- exp(betad[1])
    dl <- a / sqrt(1 + a^2)
    om <- s / sqrt(1 - 2 * dl^2 / pi)
    xi <- mu - om * dl * sqrt(2 / pi)
    z <- (dg$y - xi) / om
    w <- a * z
    if (floor) w <- -1e6 + RTMB::logspace_add(0 * w, w + 1e6)
    ld <- log(2) - log(om) + RTMB::dnorm(z, log = TRUE) +
      RTMB::pnorm(w, log.p = TRUE)
    -(RTMB::dnorm(b, 0, exp(theta), log = TRUE) + sum(ld))
  }
  RTMB::MakeADFun(nll, list(beta = unname(tpl$beta),
                            betad = unname(tpl$betad), b = b0,
                            theta = tpl$theta),
                  random = "b", integrate = gk, silent = TRUE)
}
dg <- d[d$g == "11", ]
for (fl in c(FALSE, TRUE)) {
  o2 <- tryCatch(mk(dg, tpl$b[11], fl), error = function(e) e)
  if (inherits(o2, "error")) {
    cat("floor", fl, "tape error:", conditionMessage(o2), "\n")
    next
  }
  for (lbl in c("xb", "x*", "anchor")) {
    x <- switch(lbl, xb = unname(xb), "x*" = unname(xs),
                anchor = c(unname(tpl$beta), unname(tpl$betad), tpl$theta))
    cat(sprintf("floor %-5s %-6s fn %.10f gr %s\n", fl, lbl, o2$fn(x),
                paste(format(o2$gr(x), digits = 6), collapse = " ")))
  }
}
cat("END\n")
