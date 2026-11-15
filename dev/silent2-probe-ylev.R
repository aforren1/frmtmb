# Lane silent2 fix round 2: how a frame stores a factor response.
.libPaths(c("C:/Users/adf44/source/r/wt-silent2-lib",
            "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
set.seed(1)
d <- data.frame(x = rnorm(60))
d$o <- sample(1:4, 60, TRUE)
d$of <- factor(d$o + 1, ordered = TRUE)
d$fa <- factor(c("a", "b", "c", "d")[d$o])
d$bf <- factor(ifelse(d$x > 0, "yes", "no"))
d$bl <- d$x > 0
q <- function(...) suppressWarnings(suppressMessages(frm(...)))
for (f in list(q(o ~ 1, family = cumulative(), data = d),
               q(of ~ 1, family = cumulative(), data = d),
               q(fa ~ 1, family = categorical(), data = d),
               q(bf ~ 1, family = bernoulli(), data = d),
               q(bl ~ 1, family = bernoulli(), data = d))) {
  r <- names(f$frame$y)
  cat(r, ": class", class(f$frame$y[[r]]), " head",
      head(as.character(f$frame$y[[r]]), 4), " y_levels",
      paste(f$frame$y_levels[[r]], collapse = ","), "\n")
}
