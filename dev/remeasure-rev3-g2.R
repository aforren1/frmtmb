# Reviewer round 3: G2's inv_link_<ordinal>() blocks on arrays
# (tests.distributions.R:300, 8 assertion lines, and :479, 2) are filed
# as "brms-internal, decision". By the G1 test (does an exported frmtmb
# object answer the assertion through a thin read?), shim brms's
# inv_link_<fam>(x, link) over the lane's brms_adapt_dord(), which reads
# each ordinal family object's exported ord_lpdf_make(identity), and run
# brms's block text. Seed 300; brms's own reference helpers.
.libPaths(c("C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
`%||%` <- function(a, b) if (is.null(a)) b else a
source("dev/brmsport-bin2-adapters.R")
sd <- "dev/brms-suite/brms/tests/testthat/helpers"
source(file.path(sd, "inv_link_ordinal_ch.R"))
source(file.path(sd, "simopts_catlike.R"))
# brms's x is thres - eta for cumulative and sratio, eta - thres for
# cratio and acat (as in block 238's reference calls)
il <- function(fam, x, link) {
  d <- dim(x); ncat <- d[3] + 1L
  out <- array(NA_real_, c(d[1], d[2], ncat))
  sgn <- if (fam %in% c("cratio", "acat")) -1 else 1
  for (j in seq_len(d[2])) {
    th <- matrix(sgn * x[, j, ], nrow = d[1])
    out[, j, ] <- brms_adapt_dord(fam, seq_len(ncat), rep(0, d[1]), th,
                                  link = link)
  }
  out
}
tally <- list()
note <- function(id, ok) tally[[id]] <<- c(tally[[id]], isTRUE(ok))
eqv <- function(a, b) isTRUE(all.equal(c(a), c(b)))
set.seed(300)
suppressWarnings(
for (ndraws in ndraws_vec) for (nobsv in nobsv_vec) for (ncat in ncat_vec) {
  x_test <- array(rnorm(ndraws * nobsv * (ncat - 1)),
                  dim = c(ndraws, nobsv, ncat - 1))
  nx_test <- -x_test
  for (link in c("logit", "probit", "cauchit", "cloglog")) {
    for (f in c("cumulative", "sratio", "cratio", "acat")) {
      arg <- if (f %in% c("cratio", "acat")) nx_test else x_test
      a <- il(f, arg, link)
      ch <- get(paste0("inv_link_", f, "_ch"))(arg, link = link)
      note(paste("300", f, "equivalent"), eqv(a, ch))
      note(paste("300", f, "dim"), identical(dim(a),
                                             as.integer(c(ndraws, nobsv, ncat))))
    }
    s <- il("sratio", x_test, link); cr <- il("cratio", nx_test, link)
    note("479 symmetric", if (link != "cloglog") eqv(s, cr) else
      !isTRUE(all.equal(s, cr)))
  }
})
for (id in names(tally)) {
  cat(sprintf("%-28s holds in %d of %d iterations\n", id,
              sum(tally[[id]]), length(tally[[id]])))
}
