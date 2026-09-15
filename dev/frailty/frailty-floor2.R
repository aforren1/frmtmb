# Does rp_floored() read the CENTRE DEVIATIONS, or only the population
# gamma1? The floor never fired on its own over 75 fits (see
# frailty-floor.R and the gamma1 recovery row), so the guard is pinned
# by construction instead: push ONE centre's deviation past the
# population slope and see whether the rows that come back belong to
# that centre.
source("frailty-common.R")
suppressMessages({
  library(frmtmb); library(frmtmb.spline)
})

d <- g1_sim(20260910L, n = 800L, n_centre = 20L, sd_u = 0.25)
bk <- range(log(d$time[d$event == 1L]))
f <- suppressWarnings(frm(bf(time | cens(censored) ~ trt,
                            gamma1 ~ (1 | centre)),
                         family = royston_parmar(knots = numeric(0),
                                                 bknots = bk),
                         data = d, se = TRUE))
g1 <- unname(fixef(f)$gamma1[["(Intercept)"]])
u <- as.numeric(frmtmb::ranef(f)[["centre"]])
cat("population gamma1", g1, " min centre slope", min(g1 + u), "\n")
r0 <- rp_floored(f, action = "report")
cat("as fitted: n_nonmonotone", r0$n_nonmonotone, "\n")

bad <- f
j <- 3L
# [["b"]], not $b: names(fit$estimates) is beta, betad, b, theta, and
# three of the four start with "b". The exact match wins today and
# test-bracket-access.R exists because that is not a guarantee.
bad$estimates[["b"]][j] <- -g1 - 0.2   # centre 3's slope becomes -0.2
r <- rp_floored(bad, action = "report")
rows <- attr(r, "rows")[["nonmonotone"]]
cc <- sort(unique(as.character(d$centre[rows])))
cat("after pushing centre", j, ": n_nonmonotone", r$n_nonmonotone,
    " from centres ", paste(cc, collapse = " "), "\n")
cat("rows of centre", j, "that are events:",
    sum(d$centre == levels(d$centre)[j] & d$event == 1L), "\n")
cat("observed rows of centre", j, ":", sum(d$centre == levels(d$centre)[j]),
    "\n")
e <- try(rp_floored(bad), silent = TRUE)
cat("refuses:", inherits(e, "try-error"), "\n")
if (inherits(e, "try-error")) {
  cat("message:", substr(conditionMessage(attr(e, "condition")), 1, 90),
      "\n")
}
cat("\nb slot names:", paste(head(names(bad$estimates$b), 3),
                             collapse = " "), " length ",
    length(bad$estimates$b), "\n")
