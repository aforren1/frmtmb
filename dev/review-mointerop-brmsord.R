## Reviewer: brms 2.23.0 standata on ordered factors with unused levels
## and on one- and two-value mo() (no compile).
source("C:/Users/adf44/source/r/frmtmb-wt-mointerop/dev/review-mointerop-common.R")
suppressMessages(library(brms))
d4 <- data.frame(xo = factor(c(2, 3, 4, 2, 3, 4), levels = 1:4,
                             ordered = TRUE), y = c(1, 2, 3, 1.5, 2.5, 3.2))
sd4 <- brms::standata(y ~ mo(xo), data = d4)
show("brms empty bottom level: Jmo, Xmo", list(Jmo = sd4$Jmo,
                                               Xmo = sd4$Xmo_1))
f4 <- frm(y ~ mo(xo), data = d4)
show("frm empty bottom level: D, codes",
     f4$frame$linpreds[[1]]$mo[[1]][c("D", "codes")])
d5 <- d4
d5$xo <- factor(c(1, 2, 3, 1, 2, 3), levels = 1:4, ordered = TRUE)
sd5 <- brms::standata(y ~ mo(xo), data = d5)
show("brms empty top level: Jmo, Xmo", list(Jmo = sd5$Jmo, Xmo = sd5$Xmo_1))
f5 <- frm(y ~ mo(xo), data = d5)
show("frm empty top level: D, codes",
     f5$frame$linpreds[[1]]$mo[[1]][c("D", "codes")])
## a two-level ordered factor with one level used
d6 <- data.frame(xo = factor(c("a", "a", "a"), levels = c("a", "b"),
                             ordered = TRUE), y = c(1, 2, 3))
show("brms one used of two levels", brms::standata(y ~ mo(xo), data = d6)[
  c("Jmo", "Xmo_1")])
show("frm one used of two levels", {
  f <- frm(y ~ mo(xo), data = d6)
  f$frame$linpreds[[1]]$mo[[1]][c("D", "codes")]
})
