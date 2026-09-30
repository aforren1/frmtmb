# The 0.67.0 consolidation: does lane ceplot's conditional_effects()
# (effect validity, the raw-variable frame, plot objects, new levels)
# work on lane formrobust's frame that holds an offset's variables?
# Each case is set against brms 2.23.0 where brms can answer without a
# fit (brms:::get_all_effects() is the set its effect check reads).
#
#   Rscript dev/rel067-ceoffset.R > dev/rel067-log/ceoffset.txt
.libPaths(c("C:/Users/adf44/source/r/rellib-r5",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
set.seed(20260930)
n <- 120
d <- data.frame(x = rnorm(n), z = rnorm(n), time = runif(n, 1, 5),
                f = factor(sample(c("a", "b", "c"), n, TRUE)),
                g = factor(sample(1:15, n, TRUE)),
                h = factor(sample(1:8, n, TRUE)))
d$yc <- rpois(n, d$time * exp(0.3 + 0.4 * d$x))
show <- function(lab, expr) {
  w <- character(0)
  r <- withCallingHandlers(
    tryCatch(expr, error = function(e) paste("ERROR:", conditionMessage(e))),
    warning = function(cw) {
      w <<- c(w, conditionMessage(cw)); invokeRestart("muffleWarning")
    })
  cat("==", lab, "\n")
  if (is.character(r)) cat(" ", r, "\n") else
    cat("  ok, class", class(r)[1], "effects",
        paste(names(r), collapse = " "), "\n")
  for (x in w) cat("  warning:", gsub("\n", " ", x), "\n")
  invisible(r)
}
fit <- frm(bf(yc ~ x + poly(z, 2) + f + offset(log(time))), data = d,
           family = poisson())
fb <- brms::bf(yc ~ x + poly(z, 2) + f + offset(log(time)))
cat("brms get_all_effects:",
    vapply(brms:::get_all_effects(brms::brmsterms(fb)), paste, "",
           collapse = ":"), "\n")
show("default display", conditional_effects(fit))
show("effects = 'time' (offset only)", conditional_effects(fit, "time"))
show("effects = c('x', 'time')", conditional_effects(fit, c("x", "time")))
show("effects = 'z' (inside poly)", conditional_effects(fit, "z"))
ce <- conditional_effects(fit, "x")
cat("  grid time held at", unique(ce$x$time), "; mean(time) =", mean(d$time),
    "\n")
pdf(NULL)
p <- show("plot(plot = FALSE)", plot(ce, plot = FALSE))
cat("  plot object classes:", vapply(p, function(o) class(o)[1], ""), "\n")
dev.off()
# a crossed new-level display (lane ceplot) on an offset model
fr <- frm(bf(yc ~ x + (1 | g) + (1 | h) + (1 | g:h) + offset(log(time))),
          data = d, family = poisson())
obs <- unique(paste(d$g, d$h))
gh <- expand.grid(g = levels(d$g), h = levels(d$h))
miss <- gh[!paste(gh$g, gh$h) %in% obs, ][1, ]
cat("unseen combination:", as.character(miss$g), as.character(miss$h), "\n")
show("crossed unseen, Wald, re_formula = NULL",
     conditional_effects(fr, "x", conditions = data.frame(g = miss$g,
                         h = miss$h), re_formula = NULL))
show("crossed unseen, boot", conditional_effects(fr, "x",
     conditions = data.frame(g = miss$g, h = miss$h), re_formula = NULL,
     band = "boot", boot = 20, seed = 1))
show("fitted old_levels on new data", fitted(fr, newdata = data.frame(
  x = 0, time = 2, g = "99", h = "1"), re_formula = NULL,
  allow_new_levels = TRUE, sample_new_levels = "old_levels"))
show("emmeans offset", emmeans::emmeans(fit, "f"))
cat("done\n")
