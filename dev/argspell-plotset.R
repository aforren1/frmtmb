## The graphical-argument set for plot()'s S3 contract, generated.
pdf(NULL)
pars <- names(graphics::par(no.readonly = TRUE))
dev.off()
pd <- setdiff(names(formals(graphics::plot.default)), c("x", "..."))
# names base R's own plot methods pass that are neither
straggler <- c("add", "angle", "border", "ci", "density", "freq",
               "horiz", "labels", "leaflab", "verticals",
               "xy.labels", "xy.lines")
all <- sort(unique(c(pars, pd, straggler)))
cat("n =", length(all), "\n")
w <- ""
line <- "  "
for (a in all) {
  piece <- paste0('"', a, '", ')
  if (nchar(line) + nchar(piece) > 70) {
    w <- paste0(w, sub(" $", "", line), "\n")
    line <- "  "
  }
  line <- paste0(line, piece)
}
w <- paste0(w, sub(", $", "", line), "\n")
cat(w)
