# Lane prilabel: frmtmb's class "sd" rows and per-coefficient sd prior
# resolution against brms 2.23.0, case by case.
#   arg 1: "lane" (lane library first) or "base" (rellib-r9)
# For each case: (a) default_prior()'s sd rows (class, coef, group, resp,
# dpar, nlpar, in order) must be identical(); (b) for each per-coefficient
# prior, the set of standard deviations that get normal(0, 0.5) must be
# the same: brms's from its Stan code (sd_<id>[j] mapped through
# brms:::frame_re()), frmtmb's from resolve_priorlist() entries mapped
# through the block's theta_idx. A standard deviation is keyed
# "group|coef|resp|dpar|nlpar". Data: dev/dblbar-data.R, seed 20261013;
# x1, x2, w1, w2 seed 7.
arm <- commandArgs(TRUE)[1]
if (is.na(arm)) arm <- "lane"
.libPaths(c(if (arm == "lane") "C:/Users/adf44/source/r/wt-prilabel-lib",
            "C:/Users/adf44/source/r/rellib-r9",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
wt <- "C:/Users/adf44/source/r/frmtmb-wt-prilabel"
source(file.path(wt, "dev/dblbar-data.R"))
suppressPackageStartupMessages({
  library(brms)
  library(frmtmb)
})
dd <- dblbar_data()
set.seed(7)
dd$x1 <- rnorm(nrow(dd))
dd$x2 <- rnorm(nrow(dd))
dd$w1 <- runif(nrow(dd))
dd$w2 <- 1 - dd$w1
dd$g3 <- factor(sample(c("r1", "r2"), nrow(dd), TRUE))
dd$yz <- dd$y + rnorm(nrow(dd))
out <- file.path(wt, paste0("dev/prilabel-log/sdcmp-", arm, ".txt"))
sink(out, split = TRUE)
cat("frmtmb", format(packageVersion("frmtmb")), "from",
    dirname(find.package("frmtmb")), "; brms",
    format(packageVersion("brms")), "\n")

key_cols <- c("class", "coef", "group", "resp", "dpar", "nlpar")
brms_rows <- function(bf_b, pr = NULL) {
  p <- as.data.frame(brms::default_prior(bf_b, data = dd, prior = pr))
  p <- p[p$class == "sd", key_cols]
  rownames(p) <- NULL
  p
}
frm_rows <- function(bf_f, pr = NULL) {
  p <- as.data.frame(frmtmb::default_prior(bf_f, data = dd))
  p <- p[p$class == "sd", key_cols]
  rownames(p) <- NULL
  p
}
# brms: which sd keys a Stan program gives `dens`
brms_hits <- function(bf_b, prior, dens = "normal_lpdf") {
  sc <- brms::stancode(bf_b, data = dd, prior = prior)
  re <- as.data.frame(brms:::frame_re(brms:::brmsterms(bf_b), data = dd))
  re$j <- ave(seq_len(nrow(re)), re$id, FUN = seq_along)
  key <- paste(re$group, re$coef, re$resp, re$dpar, re$nlpar, sep = "|")
  ln <- grep(paste0("lprior.*", dens, "[(]sd_"), strsplit(sc, "\n")[[1]],
             value = TRUE)
  hit <- character(0)
  for (l in ln) {
    m <- regmatches(l, regexec("sd_([0-9]+)(\\[([0-9]+)\\])?", l))[[1]]
    id <- as.integer(m[2])
    j <- if (nzchar(m[4])) as.integer(m[4]) else re$j[re$id == id]
    hit <- c(hit, key[re$id == id & re$j %in% j])
  }
  sort(unique(hit))
}
# frmtmb: which sd keys the resolver gives `kind`
frm_hits <- function(bf_f, prior, kind = "normal") {
  des <- frmtmb:::prior_design(bf_f, dd, NULL, list(),
                               check_response = FALSE)
  ent <- frmtmb:::resolve_priorlist(des, prior)$entries
  idx <- unlist(lapply(ent, function(e) {
    if (identical(e$comp, "theta") && identical(e$dist$kind, kind)) e$idx
  }))
  hit <- character(0)
  for (bk in des$frame[["re_blocks"]]) {
    sdi <- frmtmb:::block_sd_idx(bk)
    sp <- frmtmb:::block_sd_prefix(des$spec, des$frame, bk)
    cf <- if (exists("block_sd_coefs", asNamespace("frmtmb"))) {
      frmtmb:::block_sd_coefs(des$frame, bk)
    } else {
      rep(NA_character_, length(sdi))
    }
    for (r in seq_len(nrow(sp))) {
      k <- sp$k[r]
      if (bk[["theta_idx"]][sdi[k]] %in% idx) {
        hit <- c(hit, paste(bk[["group_name"]], cf[k], sp$resp[r],
                            sp$dpar[r], sp$nlpar[r], sep = "|"))
      }
    }
  }
  sort(unique(hit))
}

cases <- list(
  list(f = y ~ x + (1 + x | g)),
  list(f = y ~ x + (1 + x || g)),
  list(f = y ~ x + (1 + x | gr(g, cor = FALSE))),
  list(f = y ~ x + (1 + x | g1/g2)),
  list(f = y ~ x + (1 + x | g1:g2)),
  list(f = y ~ x + (1 | g1) + (1 | g1:g2)),
  list(f = y ~ x + (1 | g1/g2/g3)),
  list(f = y ~ x + (1 + x || g1/g2)),
  list(f = y ~ x + (1 | g2:g1)),
  list(f = y ~ x + (1 | g1/g2) + (0 + x | g1:g2)),
  list(f = y ~ x + (1 + x | mm(m1, m2))),
  list(f = y ~ x + (1 | mm(m1, m2, weights = cbind(w1, w2)))),
  list(f = y ~ x + (1 + mmc(x1, x2) | mm(m1, m2))),
  list(f = y ~ x + (f | g)),
  list(f = y ~ x + (1 | g) + (0 + x | g)),
  list(f = y ~ x + (1 + x | gr(g, by = hg))),
  list(f = y ~ x + mo(xo) + (1 + mo(xo) | g)),
  list(f = y ~ x + (1 + x | g) + (1 | g1), s = sigma ~ (1 + x | g)),
  list(f = y ~ x + (1 | g), s = sigma ~ (1 + x | g)),
  list(f = y ~ x + (1 + x | p | g), s = sigma ~ (1 | p | g))
)
mkbf <- function(pkg, cs) {
  fn <- getExportedValue(pkg, "bf")
  if (is.null(cs$s)) fn(cs$f) else fn(cs$f, cs$s)
}
n_same <- 0L
n_case <- 0L
for (cs in cases) {
  n_case <- n_case + 1L
  lab <- paste(deparse1(cs$f), if (!is.null(cs$s)) deparse1(cs$s))
  b <- tryCatch(brms_rows(mkbf("brms", cs)),
                error = function(e) paste("ERROR", conditionMessage(e)))
  f <- tryCatch(frm_rows(mkbf("frmtmb", cs)),
                error = function(e) paste("ERROR", conditionMessage(e)))
  same <- identical(b, f)
  n_same <- n_same + same
  cat(sprintf("ROWS %-5s %s\n", if (same) "SAME" else "DIFF", lab))
  if (!same) {
    cat("  brms:\n"); print(b)
    cat("  frmtmb:\n"); print(f)
  }
}
cat(sprintf("ROWS SUMMARY %d of %d cases identical\n", n_same, n_case))

# per-coefficient priors: (case index, set_prior arguments)
pcs <- list(
  list(1, list(group = "g", coef = "x")),
  list(1, list(group = "g", coef = "Intercept")),
  list(2, list(group = "g", coef = "x")),
  list(3, list(group = "g", coef = "x")),
  list(4, list(group = "g1:g2", coef = "x")),
  list(4, list(group = "g1", coef = "Intercept")),
  list(4, list(group = "g1:g2")),
  list(7, list(group = "g1:g2:g3", coef = "Intercept")),
  list(8, list(group = "g1:g2", coef = "x")),
  list(10, list(group = "g1:g2", coef = "x")),
  list(11, list(group = "mmm1m2", coef = "x")),
  list(11, list(group = "mmm1m2")),
  list(12, list(group = "mmm1m2")),
  list(13, list(group = "mmm1m2", coef = "mmcx1x2")),
  list(14, list(group = "g", coef = "fb")),
  list(15, list(group = "g", coef = "x")),
  list(16, list(group = "g", coef = "x")),
  list(17, list(group = "g", coef = "moxo")),
  list(18, list(group = "g", coef = "x")),
  list(18, list(group = "g", coef = "x", dpar = "sigma")),
  list(19, list(group = "g", coef = "Intercept")),
  list(19, list(group = "g", coef = "x", dpar = "sigma")),
  list(20, list(group = "g", coef = "Intercept")),
  list(20, list(group = "g", coef = "Intercept", dpar = "sigma")),
  list(20, list(group = "g"))
)
n_same <- 0L
n_case <- 0L
for (pc in pcs) {
  n_case <- n_case + 1L
  cs <- cases[[pc[[1]]]]
  a <- pc[[2]]
  lab <- paste0(paste(deparse1(cs$f), if (!is.null(cs$s)) deparse1(cs$s)),
                " :: ", paste(names(a), unlist(a), sep = "=", collapse = ", "))
  pb <- do.call(brms::set_prior, c(list("normal(0, 0.5)", class = "sd"), a))
  b <- tryCatch(brms_hits(mkbf("brms", cs), pb),
                error = function(e) paste("ERROR", conditionMessage(e)))
  f <- tryCatch({
    pf <- do.call(frmtmb::set_prior,
                  c(list("normal(0, 0.5)", class = "sd"), a))
    frm_hits(mkbf("frmtmb", cs), pf)
  }, error = function(e) paste("ERROR", conditionMessage(e)))
  same <- identical(b, f) && !any(grepl("^ERROR", b)) && length(b) > 0L
  n_same <- n_same + same
  cat(sprintf("PRIOR %-5s %s\n", if (same) "SAME" else "DIFF", lab))
  if (same) cat("  both:  ", paste(b, collapse = " ; "), "\n")
  if (!same) {
    cat("  brms:  ", paste(b, collapse = " ; "), "\n")
    cat("  frmtmb:", paste(f, collapse = " ; "), "\n")
  }
}
cat(sprintf("PRIOR SUMMARY %d of %d specifications identical\n", n_same,
            n_case))
sink()
