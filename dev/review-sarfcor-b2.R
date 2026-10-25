# Reviewer: part B2, the full_rank property.
# (1) every built-in's frmtmb_covstruct_info();
# (2) a student-t latent per row, gr(obs, dist = "student"), which the
#     built-in entry states full_rank = FALSE: does the weighted fit run
#     away (no maximum) while the check stays open?
# (3) a registered structure with a PD vcov at the start but a bounded
#     latent support: does the inferred full_rank refuse a model that has
#     a maximum?
# Rscript dev/review-sarfcor-b2.R
.libPaths(c("C:/Users/adf44/source/r/wt-sarfcor-lib",
            "C:/Users/adf44/source/r/rellib-r12",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), "\n")
frmtmb_local_registries(environment())
core <- names(frmtmb:::covstruct_core)
tab <- do.call(rbind, lapply(core, function(nm) {
  i <- frmtmb_covstruct_info(nm)
  data.frame(name = nm, field = i$field, term = i$term,
             coupled = i$levels_coupled, full_rank = i$full_rank)
}))
print(tab, row.names = FALSE)

set.seed(20261020)
ng <- 12; nt <- 6; n <- ng * nt
d <- data.frame(x = rnorm(n))
d$obs <- factor(seq_len(n))
eta <- 1 + 0.4 * d$x + rnorm(n, 0, 0.5)
d$y <- rgamma(n, shape = 3, rate = 3 / exp(eta))
d$w15 <- 1.5
G <- Gamma(link = "log")
show <- function(lab, f) {
  w <- character()
  r <- tryCatch(withCallingHandlers(f(), warning = function(x) {
    w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning")
  }), error = function(e) e)
  if (inherits(r, "error")) {
    cat(lab, "-> REFUSED:", substr(conditionMessage(r), 1, 150), "\n")
  } else {
    cat(lab, "-> FIT logLik", format(as.numeric(logLik(r)), digits = 10),
        "shape", format(exp(r$estimates$betad[1]), digits = 4),
        "warnings:", length(w), "\n")
    if (length(w)) cat("   first warning:", substr(w[1], 1, 160), "\n")
  }
  invisible(r)
}
show("weights 1.5, (1 | obs)",
     function() frm(y | weights(w15) ~ x + (1 | obs), family = G, data = d))
show("weights 1.5, (1 | gr(obs, dist = 'student'))",
     function() frm(y | weights(w15) ~ x + (1 | gr(obs, dist = "student")),
                    family = G, data = d))
show("weights 1.0, (1 | gr(obs, dist = 'student'))",
     function() frm(y ~ x + (1 | gr(obs, dist = "student")), family = G,
                    data = d))

# (3) a latent bounded to (-0.3, 0.3): b = 0.3 * tanh-like support through
# a density that is -Inf outside; vcov at the start is PD
frmtmb_register_covstruct(
  "boxre",
  lpdf = function(b, theta, blk) {
    # a scaled beta(2, 2) on (-c, c), c = 0.3: bounded support
    cc <- 0.3
    z <- (b + cc) / (2 * cc)
    sum(log(6 * z * (1 - z) / (2 * cc)))
  },
  vcov = function(theta, blk) {
    matrix(0.3^2 / 5, 1, 1, dimnames = list(blk[["cnms"]], blk[["cnms"]]))
  },
  npar = 0L)
cat("boxre info full_rank:", frmtmb_covstruct_info("boxre")$full_rank, "\n")
show("weights 1.5, boxre(1 | obs) [bounded latent]",
     function() frm(y | weights(w15) ~ x + boxre(1 | obs), family = G,
                    data = d))
frmtmb_register_covstruct(
  "boxre2",
  lpdf = function(b, theta, blk) {
    cc <- 0.3
    z <- (b + cc) / (2 * cc)
    sum(log(6 * z * (1 - z) / (2 * cc)))
  },
  vcov = function(theta, blk) {
    matrix(0.3^2 / 5, 1, 1, dimnames = list(blk[["cnms"]], blk[["cnms"]]))
  },
  npar = 0L, full_rank = FALSE)
show("weights 1.5, boxre2(1 | obs) [bounded latent, full_rank = FALSE]",
     function() frm(y | weights(w15) ~ x + boxre2(1 | obs), family = G,
                    data = d))
