## Reviewer, claim 3: the design set the worker did NOT include.
## Nine shapes the field produces, 10 replicates each at two sample
## sizes, one seed per (design, n, replicate) and none shared.
## usage: Rscript gradcheck-rev-08-extend.R <core-lib> <out.rds>
a <- commandArgs(TRUE)
LIB <- a[1]; OUT <- a[2]
.libPaths(c(LIB, "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("CORE:", find.package("frmtmb"), "\n")
`%||%` <- function(x, y) if (is.null(x)) y else x
gw <- function(expr) {
  w <- character(0)
  val <- withCallingHandlers(tryCatch(expr, error = function(e) e),
                             warning = function(cond) {
    w <<- c(w, conditionMessage(cond)); invokeRestart("muffleWarning")
  })
  list(fit = val, w = w,
       grad = any(grepl("Large maximum absolute gradient", w, fixed = TRUE)))
}

designs <- list(
  ## complete separation in a logistic model
  separation = function(n, s) {
    set.seed(s); dd <- data.frame(x = rnorm(n), z = rnorm(n))
    dd$y <- as.integer(dd$x > 0)
    frm(bf(y ~ x + z), family = bernoulli(), data = dd)
  },
  ## near separation: one flipped point, so the MLE is finite but large
  near_sep = function(n, s) {
    set.seed(s); dd <- data.frame(x = rnorm(n))
    dd$y <- as.integer(dd$x > 0)
    dd$y[which.max(dd$x)] <- 0L
    frm(bf(y ~ x), family = bernoulli(), data = dd)
  },
  zero_infl = function(n, s) {
    set.seed(s); dd <- data.frame(x = rnorm(n))
    dd$y <- ifelse(rbinom(n, 1, 0.3) == 1, 0L,
                   rpois(n, exp(0.5 + 0.4 * dd$x)))
    frm(bf(y ~ x, zi ~ 1), family = zero_inflated_poisson(), data = dd)
  },
  ## cumulative with a rare top category (about 0.5 percent of rows)
  rare_top = function(n, s) {
    set.seed(s); dd <- data.frame(x = rnorm(n))
    dd$yo <- cut(0.8 * dd$x + rlogis(n), breaks = c(-Inf, 0, 5.3, Inf),
                 labels = FALSE)
    if (length(unique(dd$yo)) < 3L) return(simpleError("too few levels"))
    frm(bf(yo ~ x), family = cumulative(), data = dd)
  },
  ## student() on data with no heavy tails, so nu runs off its link
  student_nu = function(n, s) {
    set.seed(s); dd <- data.frame(x = rnorm(n))
    dd$y <- rnorm(n, 1 + 2 * dd$x, 1)
    frm(bf(y ~ x), family = student(), data = dd)
  },
  ## a two-component gaussian mixture on unimodal data: one component
  ## collapses onto the other
  mix_collapse = function(n, s) {
    set.seed(s); dd <- data.frame(x = rnorm(n))
    dd$y <- rnorm(n, 1 + 0.5 * dd$x, 1)
    frm(bf(y ~ x), family = mixture(gaussian(), gaussian()), data = dd)
  },
  gp = function(n, s) {
    set.seed(s); dd <- data.frame(t = sort(runif(min(n, 600), 0, 12)))
    dd$y <- rnorm(nrow(dd), sin(dd$t) + 0.2 * dd$t, 0.3)
    frm(bf(y ~ gp(t, k = 30)), family = gaussian(), data = dd)
  },
  reml = function(n, s) {
    set.seed(s); ng <- max(20L, n %/% 25L)
    dd <- data.frame(g = factor(rep(seq_len(ng), length.out = n)))
    dd$x <- rnorm(n)
    re <- rnorm(ng, 0, 0.8)
    dd$y <- rnorm(n, 1 + 0.6 * dd$x + re[dd$g], 1)
    frm(bf(y ~ x + (1 | g)), family = gaussian(), data = dd, REML = TRUE)
  },
  smooth = function(n, s) {
    set.seed(s); dd <- data.frame(x = sort(runif(n, 0, 8)), z = rnorm(n))
    dd$y <- rnorm(n, sin(dd$x) + 0.4 * dd$z, 0.5)
    frm(bf(y ~ s(x) + z), family = gaussian(), data = dd)
  })

sizes <- c(2000L, 20000L)
REPS <- 10L
rows <- list()
for (dn in names(designs)) {
  for (n in sizes) {
    for (k in seq_len(REPS)) {
      s <- 60000L + match(dn, names(designs)) * 1000L +
        match(n, sizes) * 100L + k
      r <- gw(designs[[dn]](n, s))
      if (inherits(r$fit, "condition")) {
        rows[[length(rows) + 1L]] <- data.frame(
          design = dn, n = n, rep = k, seed = s, err = TRUE,
          conv = NA_integer_, gmax = NA_real_, proj = NA_real_,
          head = NA_real_, warn = NA, logLik = NA_real_,
          stringsAsFactors = FALSE)
        next
      }
      f <- r$fit
      d <- tryCatch(diagnose(f, quiet = TRUE), error = function(e) NULL)
      rows[[length(rows) + 1L]] <- data.frame(
        design = dn, n = n, rep = k, seed = s, err = FALSE,
        conv = f$opt$convergence,
        gmax = if (is.null(d)) NA_real_ else d$max_grad,
        proj = if (is.null(d)) NA_real_ else (d$grad_proj %||% NA_real_),
        head = if (is.null(d)) NA_real_ else (d$grad_headroom %||% NA_real_),
        warn = r$grad, logLik = as.numeric(logLik(f)),
        stringsAsFactors = FALSE)
    }
  }
  cat("done design", dn, "\n"); flush.console()
}
df <- do.call(rbind, rows)
saveRDS(df, OUT)
cat("\nrows:", nrow(df), " distinct seeds:", length(unique(df$seed)),
    " fit errors:", sum(df$err), "\n\n")
agg <- aggregate(cbind(warn = as.integer(warn), conv1 = as.integer(conv != 0)) ~
                   design + n, data = df[!df$err, ], FUN = sum)
agg$reps <- aggregate(rep ~ design + n, data = df[!df$err, ],
                      FUN = length)$rep
agg$maxg <- aggregate(gmax ~ design + n, data = df[!df$err, ],
                      FUN = function(x) max(x, na.rm = TRUE))$gmax
print(agg[order(agg$design, agg$n), ], row.names = FALSE)
cat("\nTOTAL fits", sum(!df$err), " gradient warnings", sum(df$warn, na.rm = TRUE),
    " convergence != 0", sum(df$conv != 0, na.rm = TRUE), "\n")
cat("DONE extend\n")
