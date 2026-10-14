# Reviewer data for lane dblbar: the lane's data (seed 20261013) plus
# columns for the cases the lane did not probe.
source("C:/Users/adf44/source/r/frmtmb-wt-dblbar/dev/dblbar-data.R")
dblbar_data_lane <- dblbar_data
dblbar_data <- function(n = 360, seed = 20261013) {
  d <- dblbar_data_lane(n, seed)
  set.seed(seed + 1)
  # a factor with an unused level
  d$fu <- factor(as.character(d$f), levels = c("a", "b", "c", "d"))
  # an ordered factor (contr.poly) and a character variable
  d$fo <- factor(as.character(d$f), levels = c("a", "b", "c"),
                 ordered = TRUE)
  d$fch <- as.character(d$f)
  # a factor with sum contrasts set on it
  d$fs <- d$f
  contrasts(d$fs) <- contr.sum(3)
  # character grouping
  d$gch <- as.character(d$g)
  # a third nesting level
  d$g3 <- factor(sample(sprintf("r%d", 1:3), n, TRUE))
  d$w <- round(rnorm(n), 3)
  d
}
