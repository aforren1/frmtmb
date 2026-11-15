# Short check of lane silent2's fix round, F2: anova()'s row key.
#   Rscript dev/reviews/silent2-short-f2.R base|lane
arm <- commandArgs(TRUE)[1]
libs <- c("C:/Users/adf44/source/r/rellib-r21",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
if (identical(arm, "lane")) libs <- c("C:/Users/adf44/source/r/wt-silent2-lib", libs)
.libPaths(libs)
suppressPackageStartupMessages(library(frmtmb))
cat("arm", arm, "frmtmb", find.package("frmtmb"), "\n")
q <- function(...) suppressWarnings(suppressMessages(frm(...)))
say <- function(lab, e) {
  r <- tryCatch({force(e); "ran"}, error = function(err) {
    paste("REFUSED:", substr(gsub("\n", " ", conditionMessage(err)), 1, 100))
  })
  cat(sprintf("%-58s %s\n", lab, r))
}
set.seed(11)
n <- 80
d <- data.frame(x = rnorm(n), z = rnorm(n))
d$y <- d$x + rnorm(n)
d$y2 <- d$y
d$yeps <- d$y + 1e-12
d$o <- cut(d$y, c(-Inf, -1, 0, 1, Inf), labels = FALSE)
d$of <- factor(d$o, ordered = TRUE)
d$b <- as.integer(d$y > 0)
d$bl <- d$y > 0
cat("## must run\n")
say("y ~ 1 vs I(y) ~ x", anova(q(y ~ 1, data = d), q(I(y) ~ x, data = d)))
say("y ~ 1 vs y2 ~ x (copy)", anova(q(y ~ 1, data = d), q(y2 ~ x, data = d)))
say("cumulative int vs ordered factor",
    anova(q(o ~ 1, family = cumulative(), data = d),
          q(of ~ x, family = cumulative(), data = d)))
say("bernoulli int vs logical",
    anova(q(b ~ 1, family = bernoulli(), data = d),
          q(bl ~ x, family = bernoulli(), data = d)))
cat("## must refuse\n")
say("rows 1..72 vs 9..80",
    anova(q(y ~ 1, data = d[1:72, ]), q(y ~ x, data = d[9:80, ])))
dn1 <- d
dn1$x[1:3] <- NA
dn1$z[4:6] <- NA
say("NA rows of different predictors",
    anova(q(y ~ x, data = dn1), q(y ~ z, data = dn1)))
dr <- d
dr$yr <- rev(d$y)
say("y vs rev(y)", anova(q(y ~ 1, data = dr), q(yr ~ x, data = dr)))
say("y vs y + 1e-12", anova(q(y ~ 1, data = d), q(yeps ~ x, data = d)))
cat("yeps - y range:", format(range(d$yeps - d$y)), "\n")
# factor responses: equal codes, different outcomes
k <- d$o
d$fa <- factor(c("a", "b", "c", "d")[k], levels = c("a", "b", "c", "d"))
d$fb <- factor(c("w", "x", "y", "z")[k], levels = c("w", "x", "y", "z"))
d$oa <- factor(c(1, 2, 3, 4)[k], levels = 1:4, ordered = TRUE)
d$ob <- factor(c(1, 2, 3, 5)[k], levels = c(1, 2, 3, 5), ordered = TRUE)
d$oc <- factor(c(2, 3, 4, 5)[k], levels = 2:5, ordered = TRUE)
cat("codes equal fa/fb:", identical(as.integer(d$fa), as.integer(d$fb)),
    " oa/ob:", identical(as.integer(d$oa), as.integer(d$ob)), "\n")
say("categorical fa (a..d) vs fb (w..z), equal codes",
    anova(q(fa ~ 1, family = categorical(), data = d),
          q(fb ~ x, family = categorical(), data = d)))
say("cumulative oa (1,2,3,4) vs ob (1,2,3,5), equal codes",
    anova(q(oa ~ 1, family = cumulative(), data = d),
          q(ob ~ x, family = cumulative(), data = d)))
say("cumulative oa (1..4) vs oc (2..5), equal codes",
    anova(q(oa ~ 1, family = cumulative(), data = d),
          q(oc ~ x, family = cumulative(), data = d)))
d$ba <- factor(c("no", "yes")[d$b + 1], levels = c("no", "yes"))
d$bb <- factor(c("lo", "hi")[d$b + 1], levels = c("lo", "hi"))
say("bernoulli factor no/yes vs lo/hi, equal codes",
    anova(q(ba ~ 1, family = bernoulli(), data = d),
          q(bb ~ x, family = bernoulli(), data = d)))
# the same outcomes, level order differs (codes differ)
d$fr <- factor(as.character(d$fa), levels = c("d", "c", "b", "a"))
say("categorical fa vs same labels, levels reversed",
    anova(q(fa ~ 1, family = categorical(), data = d),
          q(fr ~ x, family = categorical(), data = d)))
# integer vs factor of an integer that does not start at 1
d$oi <- c(2L, 3L, 4L, 5L)[k]
say("cumulative int 2..5 vs ordered factor oc 2..5",
    anova(q(oi ~ 1, family = cumulative(), data = d),
          q(oc ~ x, family = cumulative(), data = d)))
cat("DONE\n")
