# Reviewer checks (lane ceplot) of the two plot() methods:
#  1. plot = FALSE opens no device and writes no file;
#  2. the returned objects survive saveRDS()/readRDS() and draw the same
#     display list afterwards;
#  3. printing the returned list draws each element;
#  4. what frm_check_dots() does with base graphical parameters;
#  5. the devAskNewPage state and par() after a call;
#  6. hyp_limit_chars() against brms:::limit_chars() on random labels.
#   Rscript dev/ceplot-rev-plot.R > dev/ceplot-rev-log/plot.txt
# Data seed 5; fuzz seed 99.
.libPaths(c("C:/Users/adf44/source/r/wt-ceplot-lib",
            "C:/Users/adf44/source/r/rellib-r4",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("frmtmb from", find.package("frmtmb"), "\n")
say <- function(...) cat(sprintf(...), "\n", sep = "")
setwd(tempdir())
unlink("Rplots.pdf")
set.seed(5)
dd <- data.frame(x = rnorm(80), z = rnorm(80),
                 f = factor(rep(c("a", "b"), 40)))
dd$y <- rnorm(80, 1 + 0.5 * dd$x + (dd$f == "b") + 0.3 * dd$x * dd$z, 1)
fit <- frm(bf(y ~ x * f + x * z), family = gaussian(), data = dd)
ce <- conditional_effects(fit, effects = c("x", "f", "x:f"), resolution = 5)
cc <- conditional_effects(fit, "x", resolution = 5,
                          conditions = data.frame(f = c("a", "b")))
h <- hypothesis(fit, c("x > 0", "fb = 0", "x + z = 0"))

## 1
say("devices before: %s", paste(names(dev.list()), collapse = ",") )
p1 <- plot(ce, plot = FALSE)
p2 <- plot(cc, plot = FALSE, points = TRUE, rug = TRUE,
           line_args = list(colour = "red"))
p3 <- plot(h, plot = FALSE)
say("after plot = FALSE x3: dev.cur() = %d (%s); Rplots.pdf exists: %s",
    dev.cur(), names(dev.cur()), file.exists("Rplots.pdf"))

## 2
ops <- function(obj) {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  grDevices::dev.control("enable")
  print(obj)
  rp <- grDevices::recordPlot()
  lapply(rp[[1L]], function(e) e[[2L]])
}
f <- tempfile(fileext = ".rds")
saveRDS(list(p1 = p1, p2 = p2, p3 = p3), f)
back <- readRDS(f)
for (nm in c("p1", "p2", "p3")) {
  a <- get(nm)
  b <- back[[nm]]
  same <- vapply(seq_along(a), function(i) {
    identical(ops(a[[i]]), ops(b[[i]]))
  }, NA)
  say("%s: %d objects; identical() after readRDS: %s; display lists equal: %s",
      nm, length(a), identical(a, b), paste(same, collapse = ","))
}

## 3
n <- 0L
grDevices::pdf(NULL)
setHook("before.plot.new", function() {
  if (isTRUE(graphics::par("page"))) n <<- n + 1L
})
print(p1)
setHook("before.plot.new", NULL, "replace")
dev.off()
say("print() of the returned list of %d effects drew %d pages", length(p1), n)

## 4
grDevices::pdf(NULL)
for (a in list(list(col = "red"), list(main = "title"), list(xlab = "X"),
               list(lwd = 3), list(foo = 1))) {
  r <- tryCatch({
    w <- NULL
    withCallingHandlers(do.call(plot, c(list(ce, ask = FALSE), a)),
                        warning = function(x) {
                          w <<- conditionMessage(x)
                          invokeRestart("muffleWarning")
                        })
    if (is.null(w)) "accepted, silent" else paste("warned:", w)
  }, error = function(e) paste("ERROR:", conditionMessage(e)))
  say("plot(ce, %s = ...): %s", names(a), r)
}
dev.off()

## 5
grDevices::pdf(NULL)
op <- par(no.readonly = TRUE)
a0 <- grDevices::devAskNewPage()
invisible(plot(ce))
invisible(plot(cc, facet_args = list(nrow = 1)))
invisible(plot(h, nvariables = 2))
say("devAskNewPage before %s after %s; par() unchanged: %s", a0,
    grDevices::devAskNewPage(), identical(op, par(no.readonly = TRUE)))
dev.off()

## 6
if (requireNamespace("brms", quietly = TRUE)) {
  set.seed(99)
  ok <- 0L
  bad <- character(0)
  for (i in 1:2000) {
    len <- sample(1:60, 1)
    s <- paste(sample(c(letters, "+", "-", "(", ")", " "), len, TRUE),
               collapse = "")
    s <- paste0(s, sample(c(" > 0", " = 0", " < 0"), 1))
    ch <- sample(c(1:50), 1)
    a <- frmtmb:::hyp_limit_chars(s, ch)
    b <- brms:::limit_chars(s, ch)
    if (identical(a, b)) ok <- ok + 1L else bad <- c(bad, s)
  }
  say("hyp_limit_chars == brms:::limit_chars on %d of 2000 random labels", ok)
  if (length(bad)) print(head(bad))
}
say("done")
