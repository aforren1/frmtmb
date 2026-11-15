# Review of lane wcheck2, item 4993: how the structure checks name a
# term on multivariate and nonlinear models.
#   Rscript dev/reviews/wcheck2-review-label.R <base|lane>
a <- commandArgs(TRUE)
lib_arg <- a[1]
source("dev/wcheck2-common.R")
set.seed(1)
n <- 60
d <- data.frame(obs = factor(seq_len(n)), x = rnorm(n),
                g = factor(rep(1:12, each = 5)), one = factor(rep("u", n)))
d$y <- 1 + 0.5 * d$x + rnorm(n)
d$y2 <- d$y + rnorm(n)
d$y_a <- d$y2
d$y.b <- d$y2
d$yp <- exp(d$y / 3)
stop_on <- frmtmb_control(check_olre = "stop", check_nlev_1 = "stop")
say <- function(lab, expr) {
  m <- tryCatch({ q(expr); "no error" },
                error = function(e) gsub("\\s+", " ", conditionMessage(e)))
  cat(sprintf("%-44s %s\n", lab, substr(m, 1, 120)))
}
say("mv olre on y", frm(mvbf(bf(y ~ x + (1 | obs)), bf(y2 ~ x)), data = d,
                        control = stop_on))
say("mv olre on y_a", frm(mvbf(bf(y ~ x), bf(y_a ~ x + (1 | obs))), data = d,
                          control = stop_on))
say("mv olre on y.b", frm(mvbf(bf(y ~ x), bf(y.b ~ x + (1 | obs))), data = d,
                          control = stop_on))
say("mv shared ID (1 | p | obs) both", frm(mvbf(bf(y ~ x + (1 | p | obs)),
                                               bf(y2 ~ x + (1 | p | obs))),
                                          data = d, control = stop_on))
say("mv (x | obs) slope on y", frm(mvbf(bf(y ~ x + (1 + x | obs)),
                                        bf(y2 ~ x)), data = d,
                                   control = stop_on))
say("mv one level on y2", frm(mvbf(bf(y ~ x), bf(y2 ~ x + (1 | one))),
                              data = d, control = stop_on))
say("mv one level shared ID", frm(mvbf(bf(y ~ x + (1 | p | one)),
                                       bf(y2 ~ x + (1 | p | one))),
                                  data = d, control = stop_on))
say("mv one level, sigma of y", frm(mvbf(bf(y ~ x, sigma ~ (1 | one)),
                                         bf(y2 ~ x)), data = d,
                                    control = stop_on))
say("mv gr(obs) by", frm(mvbf(bf(y ~ x + (1 | gr(obs))), bf(y2 ~ x)),
                         data = d, control = stop_on))
say("nl univariate olre on a", frm(bf(yp ~ a * exp(b * x), a ~ 1 + (1 | obs),
                                      b ~ 1, nl = TRUE), data = d,
                                   family = gaussian(), control = stop_on))
say("nl mv olre on a of yp", frm(mvbf(bf(yp ~ a * exp(b * x),
                                         a ~ 1 + (1 | obs), b ~ 1, nl = TRUE),
                                      bf(y2 ~ x)), data = d,
                                 control = stop_on))
say("nl mv one level on a of yp", frm(mvbf(bf(yp ~ a * exp(b * x),
                                              a ~ 1 + (1 | one), b ~ 1,
                                              nl = TRUE), bf(y2 ~ x)),
                                      data = d, control = stop_on))
say("univariate olre", frm(y ~ x + (1 | obs), data = d, control = stop_on))
say("univariate one level", frm(y ~ x + (1 | one), data = d,
                                control = stop_on))
cat("DONE\n")
