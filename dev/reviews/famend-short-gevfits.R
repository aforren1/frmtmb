# Re-check of lane famend after fix round 2: GEV fits whose optimizer
# tries points where some rows lie outside the support. Run once on the
# base library and once on the lane's; compare the two outputs with
# famend-short-gevfits-cmp.R. A trial point counts as "outside" when its
# objective is finite and above 1e6, which only the support's barrier
# (1e8 per scale unit outside, and exp(-h) up to exp(700)) gives here.
#   Rscript dev/reviews/famend-short-gevfits.R <base|lane> <out.tsv>
a <- commandArgs(trailingOnly = TRUE)
lanelib <- "C:/Users/adf44/source/r/wt-famend-lib"
.libPaths(c(if (a[1] == "lane") lanelib, "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib", dirname(find.package("frmtmb")), "\n")
rec <- new.env()
rec$v <- numeric()
.rec_obj <- function(v) rec$v <- c(rec$v, v)
assign(".rec_obj", .rec_obj, globalenv())
suppressMessages(trace("nlminb", where = asNamespace("stats"), print = FALSE,
  tracer = quote({
    .of <- objective
    objective <- function(par, ...) {
      v <- .of(par, ...)
      .rec_obj(as.numeric(v)[1])
      v
    }
  })))
rg <- function(n, mu, s, xi) mu + s * ((-log(stats::runif(n)))^(-xi) - 1) / xi
rows <- list()
k <- 0
for (xi0 in c(0.2, 0.5, 0.7, 0.9, -0.3)) {
  for (seed in 1:12) {
    for (form in c("1", "x", "g")) {
      set.seed(seed * 100 + round(xi0 * 10))
      n <- 200
      d <- data.frame(x = stats::rnorm(n), g = factor(rep(1:20, each = 10)))
      u <- stats::rnorm(20, 0, 0.5)
      mu <- switch(form, "1" = 0, x = 0.8 * d$x, g = u[d$g])
      d$y <- rg(n, mu, 1, xi0)
      fo <- switch(form, "1" = y ~ 1, x = y ~ x, g = y ~ 1 + (1 | g))
      rec$v <- numeric()
      w <- character()
      f <- tryCatch(withCallingHandlers(
        frm(fo, family = gen_extreme_value(), data = d),
        warning = function(e) {
          w <<- c(w, class(e)[1])
          invokeRestart("muffleWarning")
        },
        message = function(e) invokeRestart("muffleMessage")),
        error = function(e) e)
      k <- k + 1
      big <- sum(is.finite(rec$v) & rec$v > 1e6)
      if (inherits(f, "error")) {
        rows[[k]] <- data.frame(xi0 = xi0, seed = seed, form = form,
                                code = NA, objective = NA, est = "error",
                                n_eval = length(rec$v), n_outside = big,
                                warns = substr(conditionMessage(f), 1, 60))
      } else {
        rows[[k]] <- data.frame(
          xi0 = xi0, seed = seed, form = form, code = f$opt$convergence,
          objective = sprintf("%.10f", f$opt$objective),
          est = paste(sprintf("%.8g", f$opt$par), collapse = ","),
          n_eval = length(rec$v), n_outside = big,
          warns = paste(sort(unique(w)), collapse = ";"))
      }
    }
  }
}
out <- do.call(rbind, rows)
utils::write.table(out, a[2], sep = "\t", row.names = FALSE, quote = FALSE)
cat("fits", nrow(out), "with an outside trial point", sum(out$n_outside > 0),
    "\n")
