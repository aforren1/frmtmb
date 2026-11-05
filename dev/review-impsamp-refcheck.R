# Were the "fn far below exact" rows of the away grid my reference's failure?
# Recompute the exact objective there by piecewise integrate().
source("C:/Users/adf44/source/r/frmtmb-wt-impsamp/dev/review-impsamp-common.R")
set_libs("after")
suppressPackageStartupMessages(library(frmtmb))
source("C:/Users/adf44/source/r/frmtmb-wt-impsamp/dev/review-impsamp-designs.R")
chk <- function(nm, seed, pn) {
  de <- design_list[[nm]](seed)
  fit <- suppressWarnings(frm(de$formula, family = de$family, data = de$data,
                              importance = de$draws))
  aw <- read.delim(sprintf(
    "C:/Users/adf44/source/r/frmtmb-wt-impsamp/dev/review-impsamp-out/accuracy/%s-%d-away.tsv",
    nm, seed))
  r <- aw[aw$pt == pn, ]
  par0 <- unname(fit$opt$par)
  i <- as.integer(sub("c([0-9]+)_.*", "\\1", pn))
  k <- as.numeric(sub("c[0-9]+_", "", pn))
  # the review's se_exact was from the fit's exact Hessian; recover the point from
  # the log's value by solving for theta on the fit's objective
  th <- uniroot(function(s) fit$obj$fn(replace(par0, i, s)) - r$fn,
                lower = par0[i], upper = par0[i] + 60)$root
  p <- replace(par0, i, th)
  sd <- de$sdof(p)
  dl <- split(de$data, de$data$g)
  tot <- 0
  for (d in dl) {
    f <- function(u) exp(de$ll(p, d, u) + dnorm(u, 0, sd, log = TRUE))
    br <- seq(-10 * sd - 8, 10 * sd + 8, length.out = 2001L)
    v <- 0
    for (j in seq_len(length(br) - 1L)) {
      v <- v + integrate(f, br[j], br[j + 1L], rel.tol = 1e-8, abs.tol = 0,
                         stop.on.error = FALSE)$value
    }
    tot <- tot - log(v)
  }
  cat(nm, seed, pn, "theta", th, "sd", format(sd, digits = 4), "fn", r$fn, "my AGH", r$exact,
      "integrate", tot, "\n")
}
chk("smallpois", 5L, "c3_3")
chk("cloglog", 4L, "c3_8")
