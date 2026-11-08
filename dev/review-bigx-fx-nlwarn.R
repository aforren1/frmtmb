# A warning that belongs to the point the nl second descent keeps must
# reach the user: a body linear in a shifted column plus an exact alias
# inside the body (b2 and b3 enter only as b2 + b3), against the same body
# on the centered column.
source("C:/Users/adf44/source/r/frmtmb-wt-bigx/dev/review-bigx-common.R")
for (s in 1:3) {
  set.seed(s); n <- 150; g <- factor(rep(1:15, length.out = n)); z <- rnorm(n)
  t <- runif(n)
  d <- data.frame(z = z, w = z + 1e6, g = g, t = t,
                  yg = 1 + 0.4 * z + 0.5 * t + rnorm(15, 0, 0.7)[g] + rnorm(n, 0, 0.5))
  for (body in c("b0 + b1 * %s + (b2 + b3) * t", "b0 + b1 * %s + b2 * t")) {
    mk <- function(v) bf(as.formula(paste("yg ~", sprintf(body, v))),
                         b0 ~ 1 + (1 | g), b1 ~ 1, b2 ~ 1,
                         if (grepl("b3", body)) b3 ~ 1 else NULL, nl = TRUE)
    fw <- tryCatch(mk("w"), error = function(e) NULL)
    rw <- rfit(if (grepl("b3", body)) bf(as.formula(paste("yg ~", sprintf(body, "w"))), b0 ~ 1 + (1 | g), b1 ~ 1, b2 ~ 1, b3 ~ 1, nl = TRUE)
               else bf(as.formula(paste("yg ~", sprintf(body, "w"))), b0 ~ 1 + (1 | g), b1 ~ 1, b2 ~ 1, nl = TRUE),
               family = gaussian(), data = d)
    rz <- rfit(if (grepl("b3", body)) bf(as.formula(paste("yg ~", sprintf(body, "z"))), b0 ~ 1 + (1 | g), b1 ~ 1, b2 ~ 1, b3 ~ 1, nl = TRUE)
               else bf(as.formula(paste("yg ~", sprintf(body, "z"))), b0 ~ 1 + (1 | g), b1 ~ 1, b2 ~ 1, nl = TRUE),
               family = gaussian(), data = d)
    cat(sprintf("seed %d body %-30s dll %+.3e\n   w warns(%d): %s\n   z warns(%d): %s\n", s, body,
                if (rw$err || rz$err) NA else ll_of(rw$f) - ll_of(rz$f),
                length(rw$w), substr(paste(unique(paste0("[", rw$wc, "] ", rw$w)), collapse = " || "), 1, 300),
                length(rz$w), substr(paste(unique(paste0("[", rz$wc, "] ", rz$w)), collapse = " || "), 1, 300)))
  }
}
