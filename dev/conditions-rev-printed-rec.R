# Reviewer, lane wt-conditions (round 2): tabulate record mode of
# dev/conditions-rev-printed.R, base against lane.
o <- "C:/Users/adf44/source/r/frmtmb-wt-conditions/dev/conditions-rev-log/printed"
b <- readRDS(file.path(o, "record-base.rds")); l <- readRDS(file.path(o, "record-lane.rds"))
nm <- union(names(b), names(l))
g <- function(x, n, f) if (is.null(x[[n]])) "-" else x[[n]][[f]]
tab <- data.frame(t = nm,
  base = vapply(nm, g, "", x = b, f = "class"),
  lane = vapply(nm, g, "", x = l, f = "class"),
  msg = vapply(nm, function(n) identical(g(b, n, "message"), g(l, n, "message")), NA),
  call = vapply(nm, function(n) identical(g(b, n, "call"), g(l, n, "call")), NA),
  text = substr(vapply(nm, g, "", x = l, f = "message"), 1, 45))
options(width = 250); print(tab, row.names = FALSE, right = FALSE)
cat("lane frmtmb-classed:", sum(grepl("frmtmb_", tab$lane)), "of", nrow(tab), "\n")
